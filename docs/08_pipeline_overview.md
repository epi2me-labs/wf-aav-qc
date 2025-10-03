<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->

The following (Fig.1) is a basic schematic of the workflow:
<figure>
<img src="docs/images/wf-aav-qc_overview_outline.png" alt="AAV QC overview"/>
<figcaption>Fig.1 wf-aav-qc workflow</figcaption>
</figure>

### 1: Make a combined reference sequence 
Reads can originate from the transgene cassette, but can also come from the other plasmids
used in the rAAV prep as well as host cell DNA. Therefore, a combined reference is created
that contains the following reference sequences (with relevant parameter):
* transgene plasmid with variable ITR regions masked; see below (`ref_transgene_plasmid`).
* host reference genome (`ref_host`)
* Rep-Cap plasmid (`ref_rep_cap`)
* helper plasmid (`ref_helper`)

An alternative to supplying the non-transgene plasmid sequences individually is to place them into a folder and reference
this folder with the parameter `non_transgene_refs`.



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
Transgene plasmid variants are called using [medaka](https://github.com/nanoporetech/medaka), producing a VCF file that is used to generate a consensus sequence using [bcftool concensus](https://samtools.github.io/bcftools/bcftools.html#consensus)
The workflow selects the appropriate Medaka models based on the basecaller configuration that was used to process the signal data.
By default, the workflow will attempt to determine the basecaller model from the input data.
When this fails (or when you wish to override the automatic selection), it can be provided with `--override_basecaller_cfg`.



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

At this stage, the BAM alignment files are tagged with `AV:Z` which associates each alignment with an assigned genotype, in the format `AV:Z:full_ssaav`. If the read does not map to the transgene plasmid, it will be assigned the tag `AV:Z:non_transgene`.

If --gtype_bams is set to `true`, these tagged BAMs are split on this tag into separate BAM files.