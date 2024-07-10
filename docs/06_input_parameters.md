### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| analyse_unclassified | boolean | Analyse unclassified reads from input directory. By default the workflow will not process reads in the unclassified directory. | If selected and if the input is a multiplex directory the workflow will also process the unclassified directory. | False |
| itr_fl_threshold | integer | The maximum number of bases missing from an ITR in order for it to be classed as a full length ITR. | For ITR1, this many bases can be missing from the end of the ITR region. For ITR2, this many bases can be missing from the start of the ITR region. | 100 |
| itr_backbone_threshold | integer | The maximum number of bases and alignment is allowed to extended outside of the ITR-ITR region for an associated read to not be classed as `backbone`. | Reads mapping to the transgene plasmid sometimes extend beyond the ITRs. This parameter sets a maximum number or bases after which the read is classified as `backbone`. | 20 |
| itr1_start | integer | The start position of ITR1. |  |  |
| itr1_end | integer | The end position of ITR2. |  |  |
| itr2_start | integer | The start position of ITR2. |  |  |
| itr2_end | integer | The end position of ITR2. |  |  |
| symmetry_threshold | integer | The threshold to consider whether the start or end positions on opposite strands are classed as symmetrical or asymmetrical. | For certain categories of AAV genome type we want to test whether alignments on both strands are symmetrical or asymmetrical (i.e. whether the start and end positions are approximately the same or not) This parameter sets the threshold for this comparison. | 10 |
| ref_host | string | The reference FASTA file for the host organism (.fasta/fasta.gz). |  |  |
| ref_helper | string | The helper plasmid FASTA file. |  |  |
| ref_rep_cap | string | The rep/cap plasmid FASTA file. |  |  |
| ref_transgene_plasmid | string | The transgene plasmid FASTA file. |  |  |
| basecaller_cfg | string | Name of the basecaller model that processed the signal data; used to select an appropriate Medaka model. | The basecaller configuration is used to automatically select the appropriate Medaka model. The automatic selection can be overridden with the 'medaka_model' parameters. Available models are: 'dna_r10.4.1_e8.2_400bps_hac@v3.5.2', 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2', 'dna_r9.4.1_e8_fast@v3.4', 'dna_r9.4.1_e8_hac@v3.3', 'dna_r9.4.1_e8_sup@v3.3', 'dna_r10.4.1_e8.2_400bps_hac_prom', 'dna_r9.4.1_450bps_hac_prom', 'dna_r10.3_450bps_hac', 'dna_r10.3_450bps_hac_prom', 'dna_r10.4.1_e8.2_260bps_hac', 'dna_r10.4.1_e8.2_260bps_hac_prom', 'dna_r10.4.1_e8.2_400bps_hac', 'dna_r9.4.1_450bps_hac', 'dna_r9.4.1_e8.1_hac', 'dna_r9.4.1_e8.1_hac_prom'. | dna_r10.4.1_e8.2_400bps_sup@v3.5.2 |
| medaka_model | string | The name of the Medaka model to use. This will override the model automatically chosen based on the provided basecaller configuration. | The workflow will attempt to map the basecaller model (provided with 'basecaller_cfg') used to a suitable Medaka model. You can override this by providing a model with this option instead. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |
| output_genometype_bams | boolean | If true, output a BAM file per identified AAV genome structure type. Otherwise output a BAM file per sample. | Output individual BAM files by the assigned genome type. | False |


### Miscellaneous Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| threads | integer | Maximum number of CPU threads for a process to consume. Applies to the minimap2 mapping and the AAV structure determination stages. | A minimap2 and AAV structure determination process per sample will be will be run. This setting applies a maximum number of threads to be used for each of these. | 4 |


