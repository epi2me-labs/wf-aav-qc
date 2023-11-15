# wf-aav-qc

This document provides an overview of the wv-aav-qc analysis workflow. 
wf-aav-qc is a quality control workflow for recombinant adeno-associated virus (rAAV) preps.

Workflow parameters specified in this document are shown in the command line format:
`--param_name value`. Parameters can also be defined in a Nextflow config file like so: `param_name = value`. The former is used in this document. 
For more information see the [Configuration section of the Nextflow documentation](https://www.nextflow.io/docs/latest/config.html)


## Workflow overview
The following is a basic schematic of the workflow
<figure>
<img src="images/wf-aav-qc_overview_outline.png", alt="AAV QC overview">
<figcaption>Fig.1 wf-aav-qc workflow</figcaption>
</figure>


# Workflow Stages


## 1: Make a combined reference sequence 
Reads can originate from the transgene caseete, but can also come from the other plasmids
used in the rAAV prep as well as host cell DNA. Therefore a combined refernce is created
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
<img src="images/combined_reference.png", alt="Combined reference">
<figcaption>Fig.2 Making a combined reference</figcaption>
</figure>

## 2/3: Map to reference and get alignment summaries
The reads are mapped to the combined reference with minimap2 (secondary alignments are excluded).

`seqkit bam` is used to generate alignment summaries that are used in the rest of the workflow. 
* 
<figure>
<img src="images/map.png", alt="Map to reference">
<figcaption>Fig.3 Mapping reads to combined reference</figcaption>
</figure>


## 4: Contamination
Reads that do not map to the transgene expression cassette are classed as contaminants. They can arise from
* The Rep-Cap or helper plasmids
* The host expression system
* None of the above reference sequences. The reads will be classified as `Unknown`. If there are a large proportion of reads 
in this category, it may warrant further investigation to identify the source.
<figure>
<img src="images/contamination.png", alt="Contamination">
<figcaption>Fig.4 Contamination summary plot </figcaption>
</figure>


## 5: Per-base coverage of the transgene cassette
Depth of coverage is generated for the transgene cassette region using `samtools depth` 
A plot of this data is shown which can indicate whether sufficient coverage has been
achieved across the transgene cassette.
<figure>
<img src="images/coverage.png", alt="Coverage">
<figcaption>Fig.5 ITR_ITE coverage</figcaption>
</figure>


## 6: Identification of transgene plasmid variants
Using [medaka](https://github.com/nanoporetech/medaka), variants are called, and a consensus sequence is generated for the transgene plasmid sequence  
generating the following files: `output/{sample_id}.medaka_variants.vcf.gz`, `output/{sample_id}.medaka_consensus.fasta.gz`

A relevant user option for this part of the workflow is `--basecaller_cfg`. This is the name of the basecaller model that 
was used to process the sequencing signal data. This is used to select the correct medaka model.


## 7: Identification of truncated regions
The start and end positions of alignments that map within the transgene cassette are plotted to highlight potential
regions where sequences are becoming truncated.

<figure>
<img src="images/truncations.png", alt="Truncations plot">
<figcaption>Fig.6 Plot of start and end positions of reads mapping to transgene cassette. </figcaption>
</figure>

## 8: rAAV structure determination
The rAAV transgene expression cassette will ideally exist either as full length ssAAV or ssAAV
but can also fall into various subgenomic categories.

There are two user-adjustable parameters relevant to this part of the workflow:
* `--itr_fl_threshold` (default 100). This parameter specifies the maximum number of bases missing from an ITR in order for it to be classed as a full length ITR. 
* `--itr_backbone_threshold` (default 20). Reads mapping to the transgene plasmid sometimes extend beyond the ITRs. This parameter sets a maximum number or bases after which the read is classified as `backbone`.
Diagram here ->

The following diagrams illustrate some of the possible genome type configurations that can be 
found in an rAAV prep. Different types of subgenomic structures can be formed by numerous different combinations of
truncations and the examples are representative examples only.
 
The first of these is an annotated example describing the
components of the image.

<figure>
<img src="images/diagram_notes.png", alt="Example", height="120">
<figcaption>Fig.7 Full and partial examples of a rAAV transgene expression cassette. The colour codes are preserved in the following figures </figcaption>
</figure>

### Full ssAAV (single stranded AAV)
Contains a single alignment including both ITRs (up to `itr_fl_threshold` bases missing)
<figure>
<img src="images/full_ss.png", alt="Example", height="50">

</figure>

### Genome deletion mutants (GDM)
A subgenomic type of ssAAV where part of the transgene expression cassette, internal to the ITRs, is deleted.
This class will have two alignments both on the same strand. 
<figure>
<img src="images/gdm.png", alt="GDM", height="100">
</figure>

### Incomplete genome (ICG)
Another subgenomic type of ssAAV where one side contains a full ITR (up to `itr_fl_threshold` bases missing) and the ITR is partial or missing on the other side.
<figure>
<img src="images/icg.png", alt="ICG", height="100">
</figure>


### Full scAAV (self complementary rAAV)
Contains a full or partial ITR (up to `itr_fl_threshold` bases missing) on both ends of the alignments
<figure>
<img src="images/full_sc.png", alt="Full scAAV", height="80">
</figure>


### Snapback Genome (SBG)
An scAAV subtype where only the left or right ITR region is retained.  Reads of this category will have two alignments
on opposite strands. 
These can be symmetric or asymmetric based the relative starts and end positions at the non-ITR end of the 
transgene cassette.   
<figure>
<img src="images/snapback.png", alt="SBG", height="140">
</figure>


### Unresolved SBG
A type of SBG genome (scAAV) in which ITRs present on one strand only or have a single 
ITR (no mid section) on second strand.
<figure>
<img src="images/unresolved_sbg.png", alt="Unresolved SBG", height="130">
</figure>

### ITR-ITR concatemers
These scAAV subgenomic particles contain only ITR sequence and can be full or partial.
<figure>
<img src="images/itr_concatemers.png", alt="ITR concatemers", height="100">
</figure>


### Backbone integration
Theis category contains regions from the plasmid backbone, where the start and/or end
positions of the alignment are found outside the ITR-ITR region.
<figure>
<img src="images/bb_integration.png", alt="Backbone integration", height="200">
</figure>


### Complex
The complex category contains reads with n alignments >= 3



# Workflow outputs
For each sample an output folder is created in a subdirectory of `--out_dir`
and would include the following:
```
sample_1
├── sample_1_aav_per_read_info.tsv
├── sample_1_align.bam
├── sample_1_align.bam.bai
├── sample_1_bam_info.tsv
├── sample_1_transgene_plasmid_consensus.fasta.gz
└── sample_1_transgene_plasmid_sorted.vcf.gz
```

<b>sample_1_aav_per_read_info.tsv</b>: An example excerpt is shown below
<figure>
<img src="images/per_read_tsv.png", alt="ITR concatemers">
</figure>
For each read the file gives the assigned genome subtype and extra assigned info
and the sample id

<b>sample_1_align.bam</b>: 
This alignment file from mapping the reads to the combined reference.

<b>sample_1_bam_info.tsv</b>:
The raw alignment summaries from output `seqkit bam` 

<b>sample_1_transgene_plasmid_consensus.fasta.gz </b>: The consensus FASTA sequence
of the transgene plasmid.