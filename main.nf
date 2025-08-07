#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { configure_igv } from './lib/common'

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
    cache false
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
            val(itr_locs)
        output:
            path("itr_masked_transgene_plasmid.fasta"),
            emit: masked_transgene_plasmid
    script:
    """
    workflow-glue mask_itrs \
        --transgene_plasmid_fasta "aav_transgene_plasmid.fa" \
        --itr_locations $itr_locs.itr1_start $itr_locs.itr1_end $itr_locs.itr2_start $itr_locs.itr2_end \
        --outfile "itr_masked_transgene_plasmid.fasta" \
        --transgene_plasmid_name "${transgene_plasmid_name}"
    """
}


process make_combined_reference {
    /*
    1) Make a reference file containing the sequences of the host reference and the plasmids used in the rAAV prep.
    2) Get a mapping of reference file names to sequence IDs.
    */
    label "wf_aav"
    cpus 2
    memory "2 GB"
    input:
          path("non_transgene_refs/*")
          path(masked_transgene_plasmid)
    output:
          tuple path('combined_reference.fa.gz'),
                path('combined_reference.fa.gz.fai'),
                path('combined_reference.fa.gz.gzi'),
                emit: combined_ref
          path "reference_ids.json",
                emit:ids_json
    script:
    def refs_in_dir = params.non_transgene_refs ? true:false
    """
    echo $refs_in_dir
    if [ "$refs_in_dir" != "true" ]; then
        # Ref files will be individual file paths
        refs=\$(find non_transgene_refs/*)
    else
        # Ref files will be in a subdirectory
        refdir=\$(find non_transgene_refs/*)
        refs=\$(find \$refdir/*)
    fi

    echo \$refs

    json="{\n"

    for ref in \$refs ${masked_transgene_plasmid}; do
        if [[ \$ref == *.gz ]]; then
            zcat \$ref | bgzip >> combined_reference.fa.gz
        else
            cat \$ref | bgzip >> combined_reference.fa.gz
        fi
        # Associate filename with sequence names
        contig_ids=\$(seqkit seq --name --only-id "\${ref}" | awk '{printf "\\"%s\\",", \$0}')

         # Remove the trailing comma from contig IDs
        contig_ids=\$(echo \$contig_ids | sed 's/,\$//')
        name=\$(basename \$ref)
        json+="  \\"\$name\\": [ \$contig_ids ],\n"
    done

    # Remove the trailing comma from the last entry and write to file.
    echo \$json | sed '\$ s/,\$//' > "reference_ids.json"
    echo "}" >> "reference_ids.json"

    samtools faidx combined_reference.fa.gz
    """
}

process make_mmi_index {
    label "wf_aav"
    cpus Math.min(params.threads, 16)
    memory "15 GB"
    input:
        path "ref_genome.fasta"
    output:
        path "genome_index.mmi"
    """
     minimap2 -t ${task.cpus} -I 16G -d genome_index.mmi ref_genome.fasta
    """
}


process map_to_combined_reference {
    label "wf_aav"
    cpus Math.min(params.threads, 16)
    memory "15 GB"
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
    // Keep two threads for samtools
    def mm2_threads = Math.max(task.cpus - 2, 1)
    """
    minimap2 -ax map-ont --secondary=no -t ${mm2_threads} \
        genome_index.mmi reads.fastq.gz \
        | samtools sort - \
            --write-index \
            -o ${meta['alias']}_align.bam##idx##${meta['alias']}_align.bam.bai

    seqkit bam ${meta['alias']}_align.bam 2> ${meta['alias']}_bam_info.tsv
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
        val(itr_locs)

    output:
        path('truncations.tsv'),
        emit: locations
    script:
    """
    workflow-glue truncations \
        --bam_info bam_info.tsv \
        --itr_range $itr_locs.itr1_start $itr_locs.itr2_end \
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
              path("read_stats/")
        path("transgene.fa")
        path("ref_ids.json")

    output:
        path('contam_class_counts.tsv'), emit: contam_class_counts
    script:
        def n_reads = meta['n_seqs']
    """
    workflow-glue contamination \
        --bam_info bam_info.tsv \
        --sample_id "$meta.alias" \
        --transgene_fasta transgene.fa \
        --n_reads $n_reads \
        --ref_ids ref_ids.json \
        --contam_class_counts contam_class_counts.tsv
    """
}


process aav_structures {
    label "wf_aav"
    cpus params.threads
    memory '4 GB'
    input:
        tuple val(meta),
              path("bam_info.tsv"),
              path("sorted.bam"),
              path("sorted.bam.bai")
        path("transgene_plasmid.fa")
        val(transgene_plasmid_name)
        val(itr_locs)
    output:
        path("aav_structure_counts.tsv"),
              emit: structure_counts
        tuple val(meta),
              path("*_aav_per_read_info.tsv"),
              emit: per_read_info
        tuple val(meta),
              path("tagged_bams/"),
              emit: tagged_bams
    """
    export POLARS_MAX_THREADS=${task.cpus}
    mkdir tagged_bams

    workflow-glue aav_structures \
        --bam_info bam_info.tsv \
        --bam_in sorted.bam \
        --bam_out tagged_bams/sorted.tagged.bam \
        --itr_locations \
            $itr_locs.itr1_start $itr_locs.itr1_end $itr_locs.itr2_start $itr_locs.itr2_end \
        --output_plot_data 'aav_structure_counts.tsv' \
        --output_per_read '${meta.alias}_aav_per_read_info.tsv' \
        --sample_id "${meta.alias}" \
        --transgene_plasmid_name "${transgene_plasmid_name}" \
        --itr_fl_threshold ${params.itr_fl_threshold} \
        --itr_backbone_threshold ${params.itr_backbone_threshold} \
        --symmetry_threshold ${params.symmetry_threshold} \
        --threads ${params.threads}

    cd tagged_bams

    # The BAMs output by aav_structures.py will contain a AV:Z tag that tags an alignment with its assigned genome type
    # If params.output_genometype_bams is true, then split these tagged BAMs by the AV tag
    if [ "${params.output_genometype_bams}" == "true" ]; then
        samtools split -d AV -@${task.cpus} sorted.tagged.bam
        rm sorted.tagged.bam
    fi
    samtools index -M *.bam
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
        val(itr_locs)
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
        | awk  '\$2 > ${itr_locs.itr1_start} - 10 && \$2 < ${itr_locs.itr2_end} + 10' >> itr_coverage_trimmed.tsv
    """
}


process run_medaka {
    /*
    Generate a consensus sequence and a VCF with the variant sites from alignments mapping to the transgene plasmid.
    */
    label "medaka"
    cpus 4
    memory '8 GB'
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai")
        path('transgene_plasmid.fa')
        val(transgene_plasmid_name)

    output:
        tuple val(meta),
              path("${meta.alias}.transgene_plasmid_sorted.vcf.gz"),
              emit: variants
        tuple val(meta),
              path("${meta.alias}.transgene_plasmid_consensus.fasta.gz"),
              emit: consensus
    script:
    // we use `params.override_basecaller_cfg` if present; otherwise use
    // `meta.basecall_models[0]` (there should only be one value in the list because
    // we're running ingress with `allow_multiple_basecall_models: false`; note that
    // `[0]` on an empty list returns `null`)
    String basecall_model = params.override_basecaller_cfg ?: meta.basecall_models[0]
    if (!basecall_model) {
        error "Found no basecall model information in the input data for " + \
            "sample '$meta.alias'. Please provide it with the " + \
            "`--override_basecaller_cfg` parameter."
    }
    """
    # Extract reads mapping to transgene plasmid
    samtools view align.bam -bh "${transgene_plasmid_name}" > transgene_reads.bam
    samtools index transgene_reads.bam

    medaka consensus transgene_reads.bam consensus_probs.hdf \
        --threads $task.cpus --model "${basecall_model}:variant"

    medaka variant \
        transgene_plasmid.fa \
        consensus_probs.hdf \
        "${meta.alias}.transgene_plasmid_sorted.vcf"

    bgzip "${meta.alias}.transgene_plasmid_sorted.vcf"
    bcftools index "${meta.alias}.transgene_plasmid_sorted.vcf.gz"
    bcftools consensus "${meta.alias}.transgene_plasmid_sorted.vcf.gz" \
        -f transgene_plasmid.fa \
        -o "${meta.alias}.transgene_plasmid_consensus.fasta.gz"
    """
}


process makeReport {
    label "wf_common"
    cpus 2
    memory '4 GB'
    input:
        val metadata
        path stats, stageAs: "stats_*"
        path 'truncations.tsv'
        path 'itr_coverage.tsv'
        path 'contam_class_counts.tsv'
        path 'structure_counts.tsv'
        path "versions/*"
        path "params.json"
        val wf_version
    output:
        path "wf-aav-qc-*.html"
    script:
        String report_name = "wf-aav-qc-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report $report_name \
        --wf_version $wf_version \
        --versions versions \
        --stats ${stats} \
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
        non_transgene_refs
        ref_transgene_plasmid
        itr_locs
    main:
        samples.multiMap{ meta, path, stats ->
            meta: meta
            stats: stats ?: OPTIONAL_FILE
        }.set { for_report }

        get_ref_names(ref_transgene_plasmid)
        transgene_plasmid_name = get_ref_names.out.transgene_name.splitCsv().flatten().first()

        medaka_version = medakaVersion()
        software_versions = getVersions(medaka_version)

        workflow_params = getParams()

        mask_transgene_reference(
            ref_transgene_plasmid, transgene_plasmid_name, itr_locs
        )

        make_combined_reference(
            non_transgene_refs,
            mask_transgene_reference.out.masked_transgene_plasmid
        )

        make_mmi_index(
            make_combined_reference.out.combined_ref.first()
            .map { ref, idx, idx_gz -> ref }
        )

        map_to_combined_reference(
            samples.map {meta, reads, stats -> [meta, reads]},
            make_mmi_index.out
        )

        truncations(
            map_to_combined_reference.out.bam_info,
            ref_transgene_plasmid,
            transgene_plasmid_name,
            itr_locs
        )

        itr_coverage(
            map_to_combined_reference.out.bam,
            ref_transgene_plasmid,
            transgene_plasmid_name,
            itr_locs
        )

        contamination(
            map_to_combined_reference.out.bam_info
                | join(samples.map {meta, fastq, stats -> [meta, stats]}),
               ref_transgene_plasmid,
               make_combined_reference.out.ids_json
        )

        aav_structures(
            map_to_combined_reference.out.bam_info
                .join(map_to_combined_reference.out.bam),
            ref_transgene_plasmid,
            transgene_plasmid_name,
            itr_locs
        )

        run_medaka(
            map_to_combined_reference.out.bam,
            ref_transgene_plasmid,
            transgene_plasmid_name
        )

        // For IGV, get the final paths of the published tagged BAMs.
        // We don't know the identities of the final tagged BAMs because we may or may not be splitting them
        // by AAV genome type. So get all of the files in the folder emitted by aav_structures
         if (params.igv) {
            bams_and_indexes = aav_structures.out.tagged_bams
            .flatMap { meta, dir ->
                    files = [];
                    dir.eachFile {files.add(meta['alias'] + ',' + meta['alias'] + '/tagged_bams/' + it.name)}
                    files }

            refs = make_combined_reference.out.combined_ref
                .flatMap {it.name}

            igv_paths = bams_and_indexes.mix(refs)
                .collectFile(name: 'igv_fofn.txt', newLine: true, sort: true)

            configure_igv(
                igv_paths,
                transgene_plasmid_name,
                [displayMode: "SQUISHED", colorBy: "strand"],
                 Channel.of(null), // Variant extra opts
                 Channel.of(false)) // Keep track order
         }

        metadata = for_report.meta.collect()
        stats = for_report.stats.collect()

        report = makeReport(
            metadata,
            stats,
            truncations.out.locations.collectFile(keepHeader: true),
            itr_coverage.out.collectFile(keepHeader: true),
            contamination.out.contam_class_counts.collectFile(keepHeader: true),
            aav_structures.out.structure_counts.collectFile(keepHeader: true),
            software_versions.collect(),
            workflow_params,
            workflow.manifest.version
        )
    emit:
        telemetry = workflow_params
        workflow_params
        report
        bam = aav_structures.out.tagged_bams
        bam_info = map_to_combined_reference.out.bam_info
        combined_reference = make_combined_reference.out.combined_ref.flatten()
        consensus = run_medaka.out.consensus
        variants = run_medaka.out.variants
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
        "fastcat_extra_args": "",
        "allow_multiple_basecall_models": false,
    ])
    } else {
        // if we didn't get a `--fastq`, there must have been a `--bam` (as is codified
        // by the schema)
        samples = xam_ingress([
            "input":params.bam,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "analyse_unclassified":params.analyse_unclassified,
            "keep_unaligned": true,
            "return_fastq": true,
            "stats": true,
            "allow_multiple_basecall_models": false,
        ])
    }

     // It's possible to supply any number of reference files in a folder, but we allow specifying individual
     // hist, rep-cap and helper plasmid for backwards compatibility.
    if (params.non_transgene_refs){
        non_transgene_refs = file(params.non_transgene_refs, checkIfExists: true)
    }
    else{
        ref_host = file(params.ref_host, checkIfExists: true)
        ref_helper = file(params.ref_helper, checkIfExists: true)
        ref_rep_cap = file(params.ref_rep_cap, checkIfExists: true)
        non_transgene_refs = Channel.of([ref_host, ref_helper, ref_rep_cap])
    }

    // Get the ITR positions
    def itr_locs = [:]
    if (params.transgene_bed){
        log.info('Getting ITR locations from supplied BED file.')
        valid_itrs = ['itr1', 'itr2']
        file(params.transgene_bed).text.splitEachLine('\t') { line ->
            def (ref, start, end, feature) = line

            try {
                start = start as Integer
                end = end as Integer
            } catch (Exception e) {
                throw new Exception("Start and End must be valid integers in the BED file. Line: '${line}'")
            }

            if (start >= end){
                throw new Exception("Start location must be lower than end location. Line: '${line}'")
            }

            String itr = feature.toLowerCase()
            if (valid_itrs.contains(itr)) {
                valid_itrs.remove(itr)
                itr_locs[itr +'_start'] = start
                itr_locs[itr +'_end'] = end
            }
        }
        if (itr_locs.size() != 4){
            throw new Exception("Transgene plasmid BED file must contain entries for ITR1 and ITR2")
        }


    } else {
        log.info('Using ITR locations supplied at the command line.')
        itr_locs['itr1_start'] = params.itr1_start
        itr_locs['itr1_end'] = params.itr1_end
        itr_locs['itr2_start'] = params.itr2_start
        itr_locs['itr2_end'] = params.itr2_end
    }
    ref_transgene_plasmid = file(params.ref_transgene_plasmid, checkIfExists: true)

    pipeline(
        samples,
        non_transgene_refs,
        ref_transgene_plasmid,
        itr_locs)

    pipeline.out.bam
    .mix(
        pipeline.out.bam_info,
        pipeline.out.consensus,
        pipeline.out.variants,
        pipeline.out.per_read_genome_types)
    .map {[it[1], it[0].alias]}
    .mix(
        pipeline.out.combined_reference
        .mix(pipeline.out.report)
            | map {[it, null]}) // Sample-aggregated results, so publish to root output
    | output
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}