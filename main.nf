#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended practices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concrete, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { ingress } from './lib/ingress'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "wf_common"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    """
}


process getParams {
    label "wf_aav"
    cpus 1
    output:
        path "params.json"
    script:
        String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process make_transgene_plasmid_with_ITR_ITR_orientations_reference {
    label "wf_aav"
    cpus 1
    input:
        input:
            path("aav_transgene_plasmid.fa")
            path("transgene_annotation.bed")
        output:
            path("transgene_plasmid_ITR-ITR_orientations.fasta")
    """
    workflow-glue make_ITR_orientation_reference \
        --fasta "aav_transgene_plasmid.fa" \
        --bed "transgene_annotation.bed" \
        --outfile "transgene_plasmid_ITR-ITR_orientations.fasta"
    """
}

process make_combined_reference {
    label "wf_aav"
    cpus 1
    input:
          path(ref_host)
          path("ref_helper.fa")
          path("ref_rep_cap.fa")
          path("transgene_plasmid_ITR-ITR_orientations.fa")
    output:
          path('combined_reference.fa.gz'),
          emit: combined_reference
    script:
    def compress_ref = (ref_host.getExtension() == 'gz') ? 'true':'false'
    """
    if [[ $compress_ref == 'true' ]]; then
        cp $ref_host combined_reference.fa.gz
    else
        cat $ref_host | gzip > combined_reference.fa.gz
    fi

    cat transgene_plasmid_ITR-ITR_orientations.fa \
        ref_helper.fa \
        ref_rep_cap.fa \
        | gzip >> combined_reference.fa.gz
    """
}

process map_to_combined_reference {
    label "wf_aav"
    cpus params.mapping_threads
    input:
        tuple val(meta),
              path(input)
        path("combined_reference.fa.gz")
    output:
        tuple val(meta),
              path("${meta['alias']}_align.bam"),
              path("${meta['alias']}_align.bam.bai"),
              emit: bam
        tuple val(meta),
              path("${meta['alias']}_bam_info.tsv"),
              emit: bam_info

    script:
        def usebam = false
        if (input.getExtension() in ['ubam', 'bam']){
            usebam = true
        }
    """
    if [[ "$usebam" == "true" ]]; then
        samtools fastq "$input" \
        | minimap2 -ax map-ont -secondary=no -t $task.cpus \
            combined_reference.fa.gz - \
        | samtools sort -o ${meta['alias']}_align.bam -
    else
        minimap2 -ax map-ont -secondary=no -t $task.cpus \
            combined_reference.fa.gz "$input" \
        | samtools sort -o ${meta['alias']}_align.bam -
    fi
    samtools index ${meta['alias']}_align.bam
    seqkit bam ${meta['alias']}_align.bam 2> ${meta['alias']}_bam_info.tsv
    """
}

process truncations {
    label "wf_aav"
    cpus 1
    input:
        tuple val(meta),
              path("bam_info.tsv")
        path("annotation.bed")

    output:
        tuple val(meta),
              path('truncations.tsv')
    """
    workflow-glue truncations \
        --bam_info bam_info.tsv \
        --annotation annotation.bed \
        --outfile truncations.tsv \
        --sample_id "$meta.alias"
    """
}


process contamination {
    label "wf_aav"
    cpus 1
    input:
        tuple val(meta),
              path("bam_info.tsv"),
              path('read_stats/')
        path('helper.fa')
        path('rep_cap.fa')

    output:
        path('contam_class_counts.tsv'), emit: contam_class_counts
        path('contam_read_lengths.tsv'), emit: contam_read_lengths

    """
    touch s
    # Get all the read ids - explain why
    if [[ -f read_stats/bamstats/per-read-stats.tsv ]]; then
        cut -f1 read_stats/bamstats/per-read-stats.tsv | tail -n +2 > read_ids.tsv
    else
        cut -f1 read_stats/fastcat_stats/per-read-stats.tsv | tail -n +2 > read_ids.tsv
    fi

    workflow-glue contamination \
        --bam_info bam_info.tsv \
        --sample_id "$meta.alias" \
        --helper_fasta helper.fa \
        --rep_cap_fasta rep_cap.fa \
        --read_ids read_ids.tsv \
        --contam_class_counts contam_class_counts.tsv \
        --contam_read_lengths contam_read_lengths.tsv
    """
}


process aav_structures {
    label "wf_aav"
    cpus params.mapping_threads
    input:
        tuple val(meta),
              path("bam_info.tsv")
        path("annotation.bed")
    output:
        tuple val(meta),
              path("*_aav_structure_counts.tsv"),
              emit: structure_counts
        tuple val(meta),
              path("*_aav_per_read_info.tsv"),
              emit: per_read_info

    """
    export POLARS_MAX_THREADS=${task.cpus}
    workflow-glue aav_structures \
        --bam_info bam_info.tsv \
        --annotation annotation.bed \
        --output_plot_data '${meta.alias}_aav_structure_counts.tsv' \
        --output_per_read '${meta.alias}_aav_per_read_info.tsv' \
        --sample_id "${meta.alias}"
    """
}


process itr_coverage {
    label "wf_aav"
    cpus 4
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai")
        path("annotation.bed")
    output:
        tuple val(meta),
              path('itr_coverage_trimmed.tsv')
    """
    echo -e "ref\tpos\tdepth\tstrand\tsample_id" > itr_coverage.tsv

    # mapping to forward strand
    samtools view -h -F 16 -h align.bam \
        | awk '\$1 ~ "^@" || \$3 == "flip_flip" || \$3 == "flip_flop" || \$3== "flop_flip" || \$3 == "flop_flop"' \
        | samtools depth -a -@ $task.cpus - | sed s'/\$/\tforward\t$meta.alias/' >>  itr_coverage.tsv

    # mapping to reverse strand
    samtools view -h -f 16 -h align.bam \
        | awk '\$1 ~ "^@" || \$3 == "flip_flip" || \$3 == "flip_flop" || \$3== "flop_flip" || \$3 == "flop_flop"' \
        | samtools depth -a -@ $task.cpus - | sed s'/\$/\treverse\t$meta.alias/' >> itr_coverage.tsv

    # Trim the depth TSV to ITR-ITR regions only (+/- 3bp)
    workflow-glue coverage \
        --annotation annotation.bed \
        --itr_coverage itr_coverage.tsv \
        --output itr_coverage_trimmed.tsv
    """
}

process lookup_clair3_model {
    label "wf_aav"
    cpus 1
    input:
        path("lookup_table")
        val basecall_model
    output:
        path("model/")
    shell:
    '''
    clair3_model=$(resolve_clair3_model.py lookup_table '!{basecall_model}')
    cp -r ${CLAIR_MODELS_PATH}/${clair3_model} model
    echo "Basecall model: !{basecall_model}"
    echo "Clair3 model  : ${clair3_model}"
    '''
}


process split_by_plasmid_orientation {
    // Get transgene plasmid alignemnts + plasmid references from specific orientations
    label "wf_aav"
    cpus 1
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai"),
              val(orientation)
        path "combined_reference.fa.gz"
    output:
        tuple val(meta),
              val(orientation),
              path('transgene_plasmid_orientation.fa'),
              path('itr_orientation.bam'),
              path('itr_orientation.bam.bai'),
              emit: ori_data
    """
    samtools view align.bam -bh $orientation > itr_orientation.bam
    samtools index itr_orientation.bam
    seqkit grep -r -p ^$orientation combined_reference.fa.gz -o transgene_plasmid_orientation.fa
    """
}



process call_variants {
    label "clair3"
    cpus params.mapping_threads
    input:
        tuple val(meta),
              val(orientation),
              path("transgene_plasmid_orientation.fa"),
              path("itr_orientation.bam"),
              path("itr_orientation.bam.bai"),
              path('clair3_model')
    output:
        tuple val(meta),
              val(orientation),
              path("merge_output.vcf.gz"),
              path("transgene_plasmid_orientation.fa"),
              emit: vcf
    """
    faidx transgene_plasmid_orientation.fa

    run_clair3.sh \
        --bam_fn itr_orientation.bam \
        --ref_fn transgene_plasmid_orientation.fa \
        --model_path "$clair3_model" \
        --var_pct_full=1 \
        --ref_pct_full=1 \
        --include_all_ctgs \
        --no_phasing_for_fa \
        --threads $task.cpus --chunk_size=1000 \
        --platform ont \
        -o .
    """
}

process make_consensus {
    label "wf_aav"
    cpus 4
    input:
        tuple val(meta),
              val(orientation),
              path('variants.vcf.gz'),
              path("transgene_plasmid_orientation.fa")
    output:
        tuple val(meta),
              path('transgene_plasmid_consensus_*.fa'),
              emit: consensus_fasta
    """
    bcftools index -f variants.vcf.gz
	bcftools consensus \
	    --mark-ins lc \
	    --mark-snv lc \
	    -o "transgene_plasmid_consensus_${orientation}.fa" \
	    -f transgene_plasmid_orientation.fa \
	    variants.vcf.gz
    """
}


process makeReport {
    label "wf_aav"
    input:
        val metadata
        path per_read_stats
        path 'truncations.tsv'
        path 'itr_coverage.tsv'
        path 'contam_class_counts.tsv'
        path 'contam_read_lengths.tsv'
        path 'structure_counts.tsv'
        path "versions/*"
        path "params.json"
    output:
        path "wf-aav-qc-*.html"
    script:
        String report_name = "wf-aav-qc-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        String stats_args = \
            (per_read_stats.name == OPTIONAL_FILE.name) ? "" : "--stats $per_read_stats"
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        $stats_args \
        --params params.json \
        --metadata metadata.json \
        --truncations truncations.tsv \
        --itr_coverage itr_coverage.tsv \
        --contam_class_counts contam_class_counts.tsv \
        --contam_read_lengths contam_read_lengths.tsv \
        --aav_structures structure_counts.tsv
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636. This is the only way to
// publish files from a workflow whilst decoupling the publish from the process steps.
// The process takes a tuple containing the filename and the name of a sub-directory to
// put the file into. If the latter is `null`, puts it into the top-level directory.
process output {
    // publish inputs to output directory
    label "wf_aav"
    publishDir (
        params.out_dir,
        mode: "copy",
        saveAs: { dirname ? "$dirname/$fname" : fname }
    )
    input:
        tuple path(fname), val(dirname)
    output:
        path fname
    """
    """
}

// workflow module
workflow pipeline {
    take:
        reads
        ref_host
        ref_helper
        ref_rep_cap
        ref_transgene_plamid
        transgene_plasmid_annotation
        clair3_model
    main:
        per_read_stats = reads.map {
            it[2] ? it[2].resolve('per-read-stats.tsv') : null
        }
            | collectFile ( keepHeader: true )
            | ifEmpty ( OPTIONAL_FILE )
        software_versions = getVersions()
        workflow_params = getParams()
        metadata = reads.map { it[0] }.toList()
//         reads_with_key = reads.map {[it[0].alias, *it]}
//         aav_meta_with_key = aav_meta.map {[it['barcode'], it ]} // Should we get user to use alias here?

        // Create a channel which gives tuples with the following contents:
        // 0: metadata map
        // 1: fastq
        // 2: fastcat stats
        // 3: aav_sample_metamap
//         sample_data = reads_with_key
//             | join(aav_meta_with_key)
//             // Now remove the barcode prefix merge key
//             | map {it -> it[1..-1]}
        make_transgene_plasmid_with_ITR_ITR_orientations_reference(
            ref_transgene_plamid,
            transgene_plasmid_annotation
        )

        make_combined_reference(
            ref_host,
            ref_helper,
            ref_rep_cap,
            make_transgene_plasmid_with_ITR_ITR_orientations_reference.out
        )

        map_to_combined_reference(
            reads.map {meta, reads, stats -> [meta, reads]},
            make_combined_reference.out
        )

        truncations(
            map_to_combined_reference.out.bam_info,
            transgene_plasmid_annotation
        )

        itr_coverage(
            map_to_combined_reference.out.bam,
            transgene_plasmid_annotation
        )

        contamination(
            map_to_combined_reference.out.bam_info
            | join(reads.map {meta, fastq, stats -> [meta, stats]}),
            ref_helper, ref_rep_cap
        )

        aav_structures(
            map_to_combined_reference.out.bam_info,
            transgene_plasmid_annotation
        )

        split_by_plasmid_orientation(
            map_to_combined_reference.out.bam
            | combine(channel.from('flip_flip', 'flip_flop', 'flop_flip', 'flop_flop')),
            make_combined_reference.out.combined_reference
        )

        call_variants(
            split_by_plasmid_orientation.out.ori_data
            | combine(clair3_model)
        )

        make_consensus(call_variants.out.vcf)

        report = makeReport(
            metadata, per_read_stats,
            truncations.out.collectFile(keepHeader: true),
            itr_coverage.out.collectFile(keepHeader: true),
            contamination.out.contam_class_counts.collectFile(keepHeader: true),
            contamination.out.contam_read_lengths.collectFile(keepHeader: true),
            aav_structures.out.structure_counts.collectFile(keepHeader: true),
            software_versions.collect(),
            workflow_params
        )


    emit:
        telemetry = workflow_params
        workflow_params
        report
        per_read_stats
        bam = map_to_combined_reference.out.bam
        bam_info = map_to_combined_reference.out.bam_info
        combined_reference = make_combined_reference.out.combined_reference
        consensus = make_consensus.out.consensus_fasta
        per_read_genome_types = aav_structures.out.per_read_info
}



// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }

         /// PREPARE INPUTS  ///
        // make sure that one of `--fastq` or `--ubam` was given
        def input_type = ['fastq', 'ubam'].findAll { params[it] }
        if (input_type.size() != 1) {
            error "Only provide one of '--fastq' or '--ubam'."
        }
        input_type = input_type[0]

    samples = ingress([
        "input":params[input_type],
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "fastcat_stats": true,
        "fastcat_extra_args": "",
        "input_type": input_type])

    ref_host = file(params.ref_host, checkIfExists: true)
    ref_helper = file(params.ref_helper, checkIfExists: true)
    ref_rep_cap = file(params.ref_rep_cap, checkIfExists: true)
    ref_transgene_plamid = file(params.ref_transgene_plasmid, checkIfExists: true)
    transgene_plasmid_annotation = file(params.transgene_annotation, checkIfExists: true)

    if(params.clair3_model_path) {
            log.warn "Overriding Clair3 model with ${params.clair3_model_path}."
            clair3_model = Channel.fromPath(params.clair3_model_path, type: "dir", checkIfExists: true)
    } else {
            // map basecalling model to clair3 model
            lookup_table = Channel.fromPath("${projectDir}/data/clair3_models.tsv", checkIfExists: true)
            // TODO basecaller_model_path
            clair3_model = lookup_clair3_model(lookup_table, params.basecaller_cfg)
    }

    pipeline(
        samples,
        ref_host,
        ref_helper,
        ref_rep_cap,
        ref_transgene_plamid,
        transgene_plasmid_annotation,
        clair3_model)

    pipeline.out.per_read_stats
    | map { [it, "fastq_ingress_results"] }
    | concat (
        pipeline.out.report.concat(pipeline.out.workflow_params)
            | map { [it, null] },
        pipeline.out.bam.flatMap {it ->
            bam_and_indxs = []
            for (i in it[1..2] ){
                bam_and_indxs.add([i, it[0].alias])
            }
            return bam_and_indxs
        },
        pipeline.out.bam_info.map {meta, item -> [item, meta.alias]},
        pipeline.out.consensus
            | concat(pipeline.out.per_read_genome_types)
            | map {[it[1], it[0].alias]},
        pipeline.out.combined_reference
            | map {[it, null]}

    )
    | output
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }
}
