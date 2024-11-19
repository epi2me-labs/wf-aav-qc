<!---Any other sections that are relevant specifically to this workflow and may be useful to users eg. ## Related blog posts. ## Learning center links.--->

See the [EPI2ME website](https://labs.epi2me.io/) for lots of other resources and blog posts


## AAV structure definitions

The following diagrams illustrate representative examples of the main genome configurations identifiable by wf-aav-qc.

The definitions are derived from [this publication](https://doi.org/10.1016/j.omtn.2022.08.027) 

To begin with, here are some definitions of terms used in the following section:

**Complete ITR**:
> For a given genome structure, the ITR regions are classed as `complete` if the whole ITR sequence is
covered by an alignment. However, ITR regions are unstable and can often be missing parts of the terminal regions.
When the number of bases missing from an ITR is below a certain size, it may remain functional.
In order to not reject potentially functional partial ITRs, the workflow classes ITRs as `complete` if they are missing
less than `itr_fl_threshold` terminal bases (default 100).

**Partial ITR**:
> ITRs are classed as `partial` if there is alignment support, but the number of bases covered is less than or equal to `ITR length - itr_fl_threshold`.
For example, given an ITR of 145 nt and with the default `itr_fl_threshold` of 100, an ITR would be classed as `partial`
if less than 45 nt was covered.

<figure>
<img src="docs/images/genome_diagrams/complete_partial_itr.png" alt="Partial and complete ITRs"/>
</figure>


**Payload**:
> The genetic material of interest that is to be delivered into cells. In the diagrams below this is indicated by the 
bright and pale orange boxes that represent a promoter and gene of interest.

### 1. Single stranded AAV (ssAAV)
ssAAV is a single-stranded genome flanked by inverted terminal repeat sequences (ITRs).
The workflow classifies these genome types into either full or partial ssAAV.


#### 1.1 Full ssAAV
Contains a single alignment including both `complete` ITRs.

Reads are assigned as originating from a full ssAAV genome if it meets the following criteria:
* There exists a single alignment per read that maps to the transgene plasmid cassette 
* Both ITRs are present and `complete`


<figure>
<img src="docs/images/genome_diagrams/full_ss.png" alt="Example"/>
</figure>

#### 1.2 Partial ssAAV
There are several categories of partial ssAAV each defined by which parts of the transgene cassette are absent.

##### 1.2.1 Incomplete genome (ICG)
ICGs are identified as reads that result in a single transgene cassette alignment, 
but one or both ITRs are absent or not `complete`.
They are further categorized by which ITR is present (if any).

<figure>
<img src="docs/images/genome_diagrams/icg.png" alt="ICG"/>
</figure>


##### 1.2.2 Genome deletion mutants (GDM)
This genome type is characterized by the loss of the central part of the genome. 
A read will be assigned a GDM type if it has two alignments to the transgene cassette on the same strand with a gap in between. 
The ITRs may be `complete`, `partial` or absent.

<figure>
<img src="docs/images/genome_diagrams/gdm.png" alt="GDM"/>
</figure>


### 2 Self-complementary AAV (scAAV)

scAAV genomes are created using a modified ITR sequence that has a mutated terminal resolution site, which results in a 
self-complementary viral genome with an internal ITR and is flanked by terminal ITRs.  

The workflow identifies reads as originating from scAAV genomes if they meet the following criteria:
* Has two alignments to transgene cassette.
* Each alignment maps to opposite strands.


#### 2.1 Full scAAV
Contains a `complete` ITR on both ends of the alignments. Note that the internal ITR will be included in both 
alignments as the ITRs are partially palindromic.  
<figure>
<img src="docs/images/genome_diagrams/full_sc.png" alt="Full scAAV"/>
</figure>

Below is a screenshot from IGV showing an example of a read assigned the `full scAAV` type.
<figure>
<img src="docs/images/genome_diagrams/full_sc_igv.png" alt="scaav_igv"/>
</figure>


#### 2.2 Snapback Genome (SBG)
SBG genomes are a subtype of self-complementary AAV particle that have missing parts of the genome but contain
ITRs on both the 5' and 3' ends. 

They can be further classified as either symmetrical, asymmetrical or unresolved:
* Symmetrical: This is where the same amount of each complementary part of the genome has been deleted.
* Asymmetrical: This is where varying amounts of each complementary part of the genome has been deleted.
If the size of the non-overlapping regions of each alignment is greater than `symmetry_threshold` (default 10),
it is classed as asymmetrical.
<figure>
<img src="docs/images/genome_diagrams/sbg.png" alt="SBG"/>
</figure>


#### 2.3 Unresolved SBG
A type of self-complementary AAV characterised by pronounced genome loss of, and asymmetry between, the complementary parts of
the standard genome. Inferring the structure of these AAV genome subtypes can be difficult, so some representative
alignments (a-e) are shown rather than structure diagrams:
<figure>
<img src="docs/images/genome_diagrams/unresolved_sbg.png" alt="Unresolved SBG"/>
</figure>

### 3 Other 

#### 3.1 ITR region only
These small subgenomic particles contain no payload and can consist of either:
a. Single ITR sequence
b. Two concatenated ITRs. 

<figure>
<img src="docs/images/genome_diagrams/itr_only.png" alt="ITR concatemers"/>
</figure>


#### 3.2 Backbone integration
This genome category contains DNA sequence originating from the plasmid backbone.
The start and/or end positions of the alignment are found outside the ITR-ITR region.
An allowance can be made for some backbone plasmid contamination, and this can be set with the `itr_backbone_threshold` parameter (default 20)

<figure>
<img src="docs/images/genome_diagrams/bb_integration.png" alt="Backbone integration"/>
</figure>


#### 3.4 Complex
The complex category contains reads with 3 or more alignments.