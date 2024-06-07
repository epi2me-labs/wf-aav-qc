#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process medakaVersion {
    label "medaka"
    cpus 2
    memory "2 GB"
    output:
        path "medaka_version.txt"
    """
    medaka --version | sed 's/ /,/' >> "medaka_version.txt"
    """
}


process getVersions {
    label "wf_aav"
    cpus 1
    memory "2 GB"
    input:
        path "input_versions.txt"
    output:
        path "versions.txt"

    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bcftools -v | head -n 1 |head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    python -c "import polars; print(f'polars,{polars.__version__}')"
    """
}


process get_ref_names {
    label "wf_aav"
    cpus 1
    memory "2 GB"
    input:
        path("transgene_plasmid.fa")
    output:
        path("transgene_id.txt", emit: transgene_name)
    script:
    def transgene_name = ''
    """
    seqkit seq -ni transgene_plasmid.fa | tr -d '\n' > transgene_id.txt
    """
}


process getParams {
    label "wf_aav"
    cpus 1
    memory "2 GB"
    output:
        path "params.json"
    script:
        String paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process mask_transgene_reference {
    /*
    Mask the variable regions within the transgene cassette ITRs
    */
    label "wf_aav"
    cpus 2
    memory "2 GB"
    input:
        input:
            path("aav_transgene_plasmid.fa")
            val(transgene_plasmid_name)
        output:
            path("itr_masked_transgene_plasmid.fasta"),
            emit: masked_transgene_plasmid
    """
    workflow-glue mask_itrs \
        --transgene_plasmid_fasta "aav_transgene_plasmid.fa" \
        --itr_locations $params.itr1_start $params.itr1_end $params.itr2_start $params.itr2_end \
        --outfile "itr_masked_transgene_plasmid.fasta" \
        --transgene_plasmid_name "${transgene_plasmid_name}"
    """
}


process make_combined_reference {
    /*
    Make a reference file containing the sequences of the host reference and the plasmids used in the rAAV prep.
    */
    label "wf_aav"
    cpus 2
    memory "2 GB"
    input:
          path(ref_host)
          path("ref_helper.fa")
          path("ref_rep_cap.fa")
          path("masked_transgene_plasmid.fa")
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

    cat masked_transgene_plasmid.fa \
        ref_helper.fa \
        ref_rep_cap.fa \
        | gzip >> combined_reference.fa.gz
    """
}

process make_index {
    label "wf_aav"
    cpus Math.min(params.threads, 16)
    memory "16 GB"
    input:
        path "ref_genome.fasta"
    output:
        path "genome_index.mmi", emit: index
    """
     minimap2 -t ${task.cpus} -I 16G -d genome_index.mmi ref_genome.fasta
    """
}


process map_to_combined_reference {
    label "wf_aav"
    cpus Math.min(params.threads, 16)
    memory "16 GB"
    input:
        tuple val(meta),
              path("reads.fastq.gz")
        path("genome_index.mmi")
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
    if (params.bam){
        usebam = true
    }
    // Keep two threads for samtools
    def mm2_threads = Math.max(task.cpus - 2, 1)
    """
    if [[ "$usebam" == "true" ]]; then
        samtools fastq reads.fastq.gz \
        | minimap2 -ax map-ont --secondary=no -t ${mm2_threads} \
            genome_index.mmi - \
        | samtools sort -o ${meta['alias']}_align.bam -
    else
        minimap2 -ax map-ont --secondary=no -t ${mm2_threads} \
            genome_index.mmi reads.fastq.gz \
        | samtools sort -o ${meta['alias']}_align.bam -
    fi
    samtools index ${meta['alias']}_align.bam
    seqkit bam --bins $params.bins --img '${meta['alias']}_bam.png' ${meta['alias']}_align.bam 2> ${meta['alias']}_bam_info.tsv
    """
}


process truncations {
    /*
    Filter alignments for those that start and end within the ITR-ITR cassette.
    */
    label "wf_aav"
    cpus 2
    memory "2 GB"

    input:
        tuple val(meta),
              path("bam_info.tsv")
        path("transgene_plasmid.fa")
        val(transgene_plasmid_name)

    output:
        path('truncations.tsv'),
        emit: locations
    script:
    """
    workflow-glue truncations \
        --bam_info bam_info.tsv \
        --itr_range $params.itr1_start $params.itr2_end \
        --transgene_plasmid_name "${transgene_plasmid_name}" \
        --outfile truncations.tsv \
        --sample_id "$meta.alias"
    """
}


process contamination {
    /*
    Make plot data detailing the frequency of reads mapping to various references.
    */
    label "wf_aav"
    cpus 2
    memory '4 GB'
    input:
        tuple val(meta),
              path("bam_info.tsv"),
              path('read_stats/')
        path('transgene.fa')
        path('helper.fa')
        path('rep_cap.fa')
        path('host_cell_line.fa')

    output:
        path('contam_class_counts.tsv'), emit: contam_class_counts

    """
    # Get read IDs from either bamstats or fastcat stats file.
    zcat < read_stats/*/*stats.tsv.gz | cut -f1 | tail -n +2 > read_ids.tsv

    workflow-glue contamination \
        --bam_info bam_info.tsv \
        --sample_id "$meta.alias" \
        --transgene_fasta transgene.fa \
        --helper_fasta helper.fa \
        --rep_cap_fasta rep_cap.fa \
        --host_fasta  host_cell_line.fa \
        --read_ids read_ids.tsv \
        --contam_class_counts contam_class_counts.tsv
    """
}


process aav_structures {
    label "wf_aav"
    cpus params.threads
    memory '4 GB'
    input:
        tuple val(meta),
              path("bam_info.tsv")
        path("transgene_plasmid.fa")
        val(transgene_plasmid_name)
    output:
        path("aav_structure_counts.tsv"),
              emit: structure_counts
        tuple val(meta),
              path("*_aav_per_read_info.tsv"),
              emit: per_read_info
    """
    export POLARS_MAX_THREADS=${task.cpus}
    workflow-glue aav_structures \
        --bam_info bam_info.tsv \
        --itr_locations \
            $params.itr1_start $params.itr1_end $params.itr2_start $params.itr2_end \
        --output_plot_data 'aav_structure_counts.tsv' \
        --output_per_read '${meta.alias}_aav_per_read_info.tsv' \
        --sample_id "${meta.alias}" \
        --transgene_plasmid_name "${transgene_plasmid_name}" \
        --itr_fl_threshold ${params.itr_fl_threshold} \
        --itr_backbone_threshold ${params.itr_backbone_threshold} \
        --symmetry_threshold ${params.symmetry_threshold}
    """
}


process itr_coverage {
    /*
    Make data to plot coverage at each of the four ITR-ITR cassette orientation references.
    */
    label "wf_aav"
    memory "2 GB"
    cpus 3
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai")
        path("transgene_plasmid.fa")
        val(transgene_plasmid_name)
    output:
        path('itr_coverage_trimmed.tsv')
    """
    # mapping to forward strand
    samtools view -h -F 16 -h align.bam ${transgene_plasmid_name} \
        | samtools depth -a -@ 1 - | sed s'/\$/\tforward\t$meta.alias/' >  itr_coverage_forward.tsv

    # mapping to reverse strand
    samtools view -h -f 16 -h align.bam ${transgene_plasmid_name} \
        | samtools depth -a -@ 1 - | sed s'/\$/\treverse\t$meta.alias/' > itr_coverage_reverse.tsv

    # Trim the depth TSV to ITR-ITR regions only (+/- 10bp so that the expected dropoff either side is visible)
    echo -e "ref\tpos\tdepth\tstrand\tsample_id" > itr_coverage_trimmed.tsv
    cat itr_coverage_forward.tsv itr_coverage_reverse.tsv \
        | awk  '\$2 > ${params.itr1_start} - 10 && \$2 < ${params.itr2_end} + 10' >> itr_coverage_trimmed.tsv
    """
}


process lookup_medaka_variant_model {
    label "wf_aav"
    cpus 1
    memory "2 GB"
    input:
        path("lookup_table")
        val basecall_model
    output:
        stdout
    shell:
    '''
    medaka_model=$(workflow-glue resolve_medaka_model lookup_table '!{basecall_model}' "medaka_variant")
    echo $medaka_model
    '''
}


process medaka_consensus {
    /*
    Generate a consensus sequence and a VCF with the variant sites from alignments mapping to the transgene plasmid.
    */
    label "medaka"
    cpus 4
    memory '2 GB'
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai")
        val(medaka_model)
        path('transgene_plasmid.fa')
        val(transgene_plasmid_name)

    output:
        tuple val(meta),
              path("${meta.alias}.transgene_plasmsid_consensus.fasta.gz"),
              emit: consensus
       tuple val(meta),
              path("${meta.alias}.transgene_plasmsid_sorted.vcf.gz"),
              emit: variants
    script:
        def model = medaka_model
    """
    # Extract reads mapping to transgene plasmid
    samtools view align.bam -bh "${transgene_plasmid_name}" > transgene_reads.bam
    samtools index transgene_reads.bam

    echo ${model}
    echo ${medaka_model}

    medaka consensus transgene_reads.bam "consensus_probs.hdf" \
        --threads 2 --model ${model}

    medaka stitch \
        --threads 2 \
         consensus_probs.hdf \
         transgene_plasmid.fa \
         "${meta.alias}.transgene_plasmsid_consensus.fasta"
    bgzip "${meta.alias}.transgene_plasmsid_consensus.fasta"

    medaka variant \
         transgene_plasmid.fa \
         consensus_probs.hdf \
         "transgene_plasmid.vcf"

   bcftools sort transgene_plasmid.vcf > "${meta.alias}.transgene_plasmsid_sorted.vcf.gz"
    """
}


process combine_stats {
    label "wf_aav"
    cpus 2
    memory "2 GB"
    input:
        tuple val(meta),
              path('per-read-stats.tsv.gz')
        output:
            path('stats.tsv')
    """
    gunzip -c per-read-stats.tsv.gz |
        # Add sample_id column
        sed "s/\$/\t${meta.alias}/" > stats.tsv
    """
}


process makeReport {
    label "wf_aav"
    cpus 2
    memory '4 GB'
    input:
        val metadata
        path 'per_read_stats.tsv'
        path 'truncations.tsv'
        path 'itr_coverage.tsv'
        path 'contam_class_counts.tsv'
        path 'structure_counts.tsv'
        path "versions/*"
        path "params.json"
    output:
        path "wf-aav-qc-*.html"
    script:
        String report_name = "wf-aav-qc-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --versions versions \
        --stats per_read_stats.tsv \
        --params params.json \
        --metadata metadata.json \
        --truncations truncations.tsv \
        --itr_coverage itr_coverage.tsv \
        --contam_class_counts contam_class_counts.tsv \
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
    cpus 2
    memory "2 GB"
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
        samples
        ref_host
        ref_helper
        ref_rep_cap
        ref_transgene_plasmid
    main:
        per_read_stats = samples.flatMap {
            it[2] ? file(it[2].resolve('*read*.tsv.gz')) : null
        }

        get_ref_names(ref_transgene_plasmid)
        transgene_plasmid_name = get_ref_names.out.transgene_name.splitCsv().flatten().first()

        medaka_version = medakaVersion()
        software_versions = getVersions(medaka_version)

        workflow_params = getParams()
        metadata = samples.map { it[0] }.toList()

        mask_transgene_reference(
            ref_transgene_plasmid, transgene_plasmid_name
        )

        make_combined_reference(
            ref_host,
            ref_helper,
            ref_rep_cap,
            mask_transgene_reference.out.masked_transgene_plasmid
        )

        make_index(
            make_combined_reference.out.combined_reference
        )

        map_to_combined_reference(
            samples.map {meta, reads, bins, stats -> [meta, reads]},
            make_index.out.index
        )

        truncations(
            map_to_combined_reference.out.bam_info,
            ref_transgene_plasmid,
            transgene_plasmid_name
        )

        itr_coverage(
            map_to_combined_reference.out.bam,
            ref_transgene_plasmid,
            transgene_plasmid_name
        )

        contamination(
            map_to_combined_reference.out.bam_info
            | join(samples.map {meta, fastq, stats -> [meta, stats]}),
            ref_transgene_plasmid, ref_helper, ref_rep_cap, ref_host
        )

        aav_structures(
            map_to_combined_reference.out.bam_info,
            ref_transgene_plasmid,
            transgene_plasmid_name
        )

       if (params.medaka_model) {
            log.warn "Overriding Medaka model with ${params.medaka_model}."
            medaka_model = params.medaka_model
       } else {
            Path lookup_table = file(
                "${projectDir}/data/medaka_models.tsv", checkIfExists: true)
            medaka_model = lookup_medaka_variant_model(lookup_table, params.basecaller_cfg)
       }

        medaka_consensus(
            map_to_combined_reference.out.bam,
            medaka_model,
            ref_transgene_plasmid,
            transgene_plasmid_name
        )

        report = makeReport(
            metadata,
            samples.map { meta, reads, stats ->
            stats ? [meta, file(stats.resolve('*read*.tsv.gz'))] : [null, null]
        } | combine_stats | collectFile(keepHeader: true),

            truncations.out.locations.collectFile(keepHeader: true),
            itr_coverage.out.collectFile(keepHeader: true),
            contamination.out.contam_class_counts.collectFile(keepHeader: true),
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
        consensus = medaka_consensus.out.consensus
        variants = medaka_consensus.out.variants
        per_read_genome_types = aav_structures.out.per_read_info
}



// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)

    if (params.fastq) {
    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "stats": true,
        "fastcat_extra_args": ""])
    } else {
        // if we didn't get a `--fastq`, there must have been a `--bam` (as is codified
        // by the schema)
         samples = xam_ingress([
        "input":params.bam,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "analyse_unclassified":params.analyse_unclassified,
        "keep_unaligned": true,
        "stats": true
        ])
    }

    ref_host = file(params.ref_host, checkIfExists: true)
    ref_helper = file(params.ref_helper, checkIfExists: true)
    ref_rep_cap = file(params.ref_rep_cap, checkIfExists: true)
    ref_transgene_plasmid = file(params.ref_transgene_plasmid, checkIfExists: true)

    pipeline(
        samples,
        ref_host,
        ref_helper,
        ref_rep_cap,
        ref_transgene_plasmid)

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
            | concat(
                pipeline.out.variants,
                pipeline.out.per_read_genome_types
            )
            | map {[it[1], it[0].alias]},
        pipeline.out.combined_reference
            | map {[it, null]}
    )
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}