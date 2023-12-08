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