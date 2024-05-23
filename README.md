# AAV QC workflow

Nextflow workflow for AAV vector quality control.



## Introduction

This workflow takes reads sequenced from adeno-associated virus (rAAV) vector preps and 
does some basic quality control checks. The main stages of the workflow are: 


+ Masking of transgene cassette variable ITR regions
+ Mapping reads to a combined reference sequence containing host cell reference genome, mask transgene plasmid and other AAV plasmids
+ Identification of reference which each reads maps to
+ Truncation hotspot identification
+ ITR transgene cassette coverage
+ Determination of AAV genome structure types
+ Calling transgene plasmid variants and creation of consensus 




## Compute requirements

Recommended requirements:

+ CPUs = 8
+ Memory = 32GB

Minimum requirements:

+ CPUs = 4
+ Memory = 16GB

Approximate run time: 15 minutes per sample - 150k reads and 8 cpus

ARM processor support: False




## Install and run

<!---Nextflow text remains the same across workflows, update example cmd and demo data sections.--->

These are instructions to install and run the workflow on command line. You can also access the workflow via the [EPI2ME application](https://labs.epi2me.io/downloads/).  

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and software resources, therefore nextflow will need to be installed before attempting to run the workflow. 

The workflow can currently be run using either [Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of 
the required software. Both methods are automated out-of-the-box provided 
either docker or singularity is installed. This is controlled by the [`-profile`](https://www.nextflow.io/docs/latest/config.html#config-profiles) parameter as exemplified below. 

It is not required to clone or download the git repository in order to run the workflow. 
More information on running EPI2ME workflows can be found on our [website](https://labs.epi2me.io/wfindex).

The following command can be used to obtain the workflow. This will pull the repository in to the assets folder of nextflow and provide a list of all parameters available for the workflow as well as an example command:

```
nextflow run epi2me-labs/wf-aav-qc –help 
```
A demo dataset is provided for testing of the workflow. It can be downloaded using: 
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-aav-qc/wf-aav-qc-demo.tar.gz
tar -xzvf wf-aav-qc-demo.tar.gz 
```
The workflow can be run with the demo data using: 
```
nextflow run epi2me-labs/wf-aav-qc \
    --fastq wf-aav-qc-demo/simulated_reads.fq \
    --ref_host wf-aav-qc-demo/cell_line.fasta.gz \
    --ref_helper wf-aav-qc-demo/helper.fasta \
    --ref_rep_cap  wf-aav-qc-demo/repcap.fasta \
    --ref_transgene_plasmid wf-aav-qc-demo/transgene.fasta \
    --itr1_start 11 \
    --itr1_end 156 \
    --itr2_start 2156 \
    --itr2_end 2286 \
    -profile standard
```
For further information about running a workflow on the command line interface, see https://labs.epi2me.io/wfquickstart/  



## Related protocols

This workflow is designed to take input sequences that have been produced from [Oxford Nanopore Technologies](https://nanoporetech.com/) devices using this protocol:

https://community.nanoporetech.com/docs/prepare/library_prep_protocols/ligation-sequencing-gdna-v14-adeno-associated-virus-sequencing/v/aav_9194_v114_reva_20sep2023




## Input example

<!---Example of input directory structure, delete and edit as appropriate per workflow.--->
This workflow accepts either FASTQ or BAM files as input.

The FASTQ or BAM input parameters for this workflow accept one of three cases: (i) the path to a single FASTQ or BAM file; (ii) the path to a top-level directory containing FASTQ or BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ or BAM files. In the first and second cases (i and ii), a sample name can be supplied with `--sample`. In the last case (iii), the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`.

```
(i)                     (ii)                 (iii)    
input_reads.fastq   ─── input_directory  ─── input_directory
                        ├── reads0.fastq     ├── barcode01
                        └── reads1.fastq     │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                              └── reads0.fastq
```



## Input parameters

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






## Outputs

Output files may be aggregated including information for all samples or provided per sample. Per-sample files will be prefixed with respective aliases and represented below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| Workflow report | ./wf-aav-qc-report.html | Report for all samples | aggregated |
| Combined reference sequence | ./combined_reference.fa.gz | Reference file containing all AAV plasmid and host genome sequences. | aggregated |
| Combined reference sequence index | ./combined_reference.fa.gz.fai | Index file for combined reference FASTA. | aggregated |
| Compressed combined reference sequence index | ./combined_reference.fa.gz.gzi | Extra index file for combined reference FASTA (required because combined reference is bgzip-compressed). | aggregated |
| Per read alignment info | ./{{ alias }}/{{ alias }}_bam_info.tsv | The result of `seqkit bam`. | per-sample |
| AAV structure assignment | ./{{ alias }}/{{ alias }}_aav_per_read_info.tsv | AAV per read genome subtypes. | per-sample |
| Transgene plasmid consensus | ./{{ alias }}/{{ alias }}_transgene_plasmsid_consensus.fasta.gz | The transgene plasmid consensus sequence generated by medaka. | per-sample |
| Transgene plasmid variants | ./{{ alias }}/{{ alias }}_transgene_plasmid_variants.vcf.gz | The transgene plasmid variants file generated by medaka. | per-sample |
| Alignment file | ./{{ alias }}/tagged_bams/sorted.tagged.bam | The resulting tagged BAM file from mapping reads to the combined reference. | per-sample |
| Alignment index file | ./{{ alias }}/tagged_bams/sorted.tagged.bam.bai | The index for the resulting tagged BAM file from mapping reads to the combined reference. | per-sample |
| backbone_contamination alignment | ./{{ alias }}/tagged_bams/backbone_contamination.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: backbone_contamination | per-sample |
| backbone_contamination alignment index | ./{{ alias }}/tagged_bams/backbone_contamination.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: backbone_contamination | per-sample |
| full_ssaav alignment | ./{{ alias }}/tagged_bams/full_ssaav.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: full_ssaav | per-sample |
| full_ssaav alignment index | ./{{ alias }}/tagged_bams/full_ssaav.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: full_ssaav | per-sample |
| partial_ssaav alignment | ./{{ alias }}/tagged_bams/partial_ssaav.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: partial_ssaav | per-sample |
| partial_ssaav alignment index | ./{{ alias }}/tagged_bams/partial_ssaav.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: partial_ssaav | per-sample |
| full_scaav alignment | ./{{ alias }}/tagged_bams/full_scaav.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: full_scaav | per-sample |
| full_scaav alignment index | ./{{ alias }}/tagged_bams/full_scaav.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: full_scaav | per-sample |
| itr_region_only alignment | ./{{ alias }}/tagged_bams/itr_region_only.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: itr_region_only | per-sample |
| itr_region_only alignment index | ./{{ alias }}/tagged_bams/itr_region_only.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: itr_region_only | per-sample |
| complex alignment | ./{{ alias }}/tagged_bams/complex.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: complex | per-sample |
| complex alignment index | ./{{ alias }}/tagged_bams/complex.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: complex | per-sample |
| partial_scaav alignment | ./{{ alias }}/tagged_bams/partial_scaav.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: partial_scaav | per-sample |
| partial_scaav alignment index | ./{{ alias }}/tagged_bams/partial_scaav.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: partial_scaav | per-sample |
| unknown alignment | ./{{ alias }}/tagged_bams/unknown.bam | The resulting tagged BAM file from mapping reads to the combined reference. This file contains alignmnets with the genotype assignment: unknown | per-sample |
| unknown alignment index | ./{{ alias }}/tagged_bams/unknown.bam.bai | The resulting tagged BAM index file from mapping reads to the combined reference. This indexes the file containing alignments with the genotype assignment: unknown | per-sample |
| IGV config JSON file | ./igv.json | JSON file with IGV config options to be used by the EPI2ME Desktop Application. | aggregated |




## Pipeline overview

<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->

The following (Fig.1) is a basic schematic of the workflow:
<figure>
<img src="docs/images/wf-aav-qc_overview_outline.png" alt="AAV QC overview"/>
<figcaption>Fig.1 wf-aav-qc workflow</figcaption>
</figure>

### 1: Make a combined reference sequence 
Reads can originate from the transgene cassette, but can also come from the other plasmids
used in the rAAV prep as well as host cell DNA. Therefore, a combined reference is created
that contains the following reference sequences:
* host reference genome
* Rep-Cap plasmid
* helper plasmid 
* transgene plasmid (variable ITR regions masked).

The transgene plasmid ITR cassette will naturally exist in four orientations. 
(termed flip-flip, flip-flop, flop-flop and flop-flip; see Fig.2)
This can lead to incorrect mapping of reads. To address this, the variable regions in the transgene cassette are masked.
This is done by taking the input transgene plasmid and locating the two ITR regions as defined in the `--transgene_annotation`
file. The `C'`, `C`, `B'` and `B` ITR regions are identified for each ITR. From these regions it can be determined which positions are constant between orientations and which are variable, and will be masked.

<figure>
<img src="docs/images/combined_reference.png" alt="Combined reference"/>
<figcaption>Fig.2 Making a combined reference</figcaption>
</figure>

### 2: Map to reference and get alignment summaries
The reads are mapped to the combined reference using minimap2 (secondary alignments are excluded).
`seqkit bam` is used to generate alignment summaries that are used in the rest of the workflow. 

### 3: Contamination
Reads that do not map to the transgene expression cassette are classified as contaminants. They can arise from
* The Rep-Cap or helper plasmids
* The host expression system
* None of the above reference sequences. The reads will be classified as `Unknown`. If there are a large proportion of reads 
in this category, it may warrant further investigation to identify the source.
<figure>
<img src="docs/images/contamination.png" alt="Contamination"/>
<figcaption>Fig.3 Contamination summary plot </figcaption>
</figure>


### 4: Per-base coverage of the transgene cassette
Depth of coverage is generated for the transgene cassette region using `samtools depth`. 
A plot of this data is shown which indicates whether sufficient coverage has been
achieved across the transgene cassette.
<figure>
<img src="docs/images/coverage.png" alt="Coverage" height="300"/>
<figcaption>Fig.4 ITR-ITR coverage</figcaption>
</figure>


### 5: Identification of transgene plasmid variants
Variants are called using [medaka](https://github.com/nanoporetech/medaka), and a consensus sequence is generated for the transgene plasmid sequence  
generating the following files: `output/{sample_id}.medaka_variants.vcf.gz`, `output/{sample_id}.medaka_consensus.fasta.gz`

A relevant user option for this part of the workflow is `--basecaller_cfg`. This is the name of the basecaller model that 
was used to process the sequencing signal data. This is used to select the correct medaka model.


### 6: Identification of truncated regions
The 'start' and 'end' positions of alignments that map within the transgene cassette are plotted to highlight potential
regions where sequences are becoming truncated.

<figure>
<img src="docs/images/truncations.png" alt="Truncations plot" height="300"/>
<figcaption>Fig.5 Plot of start and end positions of reads mapping to transgene cassette. </figcaption>
</figure>

### 7: rAAV structure determination
The rAAV transgene expression cassette will ideally exist as full length ITR-flanked regions. 
However, subgenomic particles will be present in any prep, and it can be useful to know the abundance of the various genome
types, which is the aim of this stage of the workflow. Genome types are assigned to each read by applying a series of heuristics that use the characteristics of each alignment from the read. 

There are two user-adjustable parameters relevant to this part of the workflow:
* `--itr_fl_threshold` (default 100). This parameter specifies the maximum number of bases missing from an ITR in order for it to be classed as a full length ITR. 
* `--itr_backbone_threshold` (default 20). Reads mapping to the transgene plasmid sometimes extend beyond the ITRs. This parameter sets a maximum number of bases after which the read is classified as `backbone`.

See the [AAV structures](#aav-structure-diagrams) section for some representative diagrams of AAV gene structures and how they are classified.

At this stage, the BAM alignment files are tagged with `AV:Z` which associates each alignment with an assigned genotype, in the format `AV:Z:full_ssaav`.

If --gtype_bams is set to `true`, these tagged BAMs are split on this tag into separate BAM files.



## Troubleshooting

<!---Any additional tips.--->
+ If the workflow fails please run it with the demo data set to ensure the workflow itself is working. This will help us determine if the issue is related to the environment, input parameters or a bug.
+ See how to interpret some common nextflow exit codes [here](https://labs.epi2me.io/trouble-shooting/).

+ Please ensure that the ITR-ITR cassette region spans and contiguous range in the transgene plasmid reference sequence.
(i.e. ITR1 sequence should precede the ITR2 sequence)




## FAQ's

<!---Frequently asked questions, pose any known limitations as FAQ's.--->

If your question is not answered here, please report any issues or suggestions on the [github issues](https://github.com/epi2me-labs/wf-aav-qc/issues) page or start a discussion on the [community](https://community.nanoporetech.com/). 



## Related blog posts

<!---Any other sections that are relevant specifically to this workflow and may be useful to users eg. ## Related blog posts. ## Learning center links.--->

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts


## AAV structure diagrams

The following diagrams illustrate some of the possible genome type configurations that can be 
found in an rAAV prep. Different types of subgenomic structures can be formed by numerous different combinations of
truncations and the examples are representative examples only.


The first of these is an annotated example describing the
components of the image.

<figure>
<img src="docs/images/diagram_notes.png" alt="Example" height="120"/>
<figcaption>Fig.6 Full and partial examples of a rAAV transgene expression cassette. The colour codes are preserved in the following figures </figcaption>
</figure>

#### Full ssAAV (single stranded AAV)
Contains a single alignment including both ITRs (up to `itr_fl_threshold` bases missing)
<figure>
<img src="docs/images/full_ss.png" alt="Example" height="50"/>

</figure>

### Genome deletion mutants (GDM)
A subgenomic type of ssAAV where part of the transgene expression cassette, internal to the ITRs, is deleted.
This class will have two alignments both on the same strand. 
<figure>
<img src="docs/images/gdm.png" alt="GDM" height="100"/>
</figure>

### Incomplete genome (ICG)
Another subgenomic type of ssAAV where one side contains a full ITR (up to `itr_fl_threshold` bases missing) and the ITR is partial or missing on the other side.
<figure>
<img src="docs/images/icg.png" alt="ICG" height="100"/>
</figure>


### Full scAAV (self complementary rAAV)
Contains a full or partial ITR (up to `itr_fl_threshold` bases missing) on both ends of the alignments
<figure>
<img src="docs/images/full_sc.png" alt="Full scAAV" height="80"/>
</figure>


### Snapback Genome (SBG)
An scAAV subtype where only the left or right ITR region is retained.  Reads of this category will have two alignments
on opposite strands. 
These can be symmetric or asymmetric based the relative starts and end positions at the non-ITR end of the 
transgene cassette.   
<figure>
<img src="docs/images/snapback.png" alt="SBG" height="140"/>
</figure>


### Unresolved SBG
A type of SBG genome (scAAV) in which ITRs present on one strand only or have a single 
ITR (no mid section) on second strand.
<figure>
<img src="docs/images/unresolved_sbg.png" alt="Unresolved SBG" height="130"/>
</figure>

### ITR-ITR concatemers
These scAAV subgenomic particles contain only ITR sequence and can be full or partial.
<figure>
<img src="docs/images/itr_concatemers.png" alt="ITR concatemers" height="100"/>
</figure>


### Backbone integration
Theis category contains regions from the plasmid backbone, where the start and/or end
positions of the alignment are found outside the ITR-ITR region.
<figure>
<img src="docs/images/bb_integration.png" alt="Backbone integration" height="200"/>
</figure>


### Complex
The complex category contains reads with n alignments >= 3



