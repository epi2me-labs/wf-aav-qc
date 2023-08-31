## Introduction

This workflow generates a reference geneome which contains the following
* Host reference genome (--ref_host) 
* Transgene plasmid (--transgene_plasmid)
* AAV helper plasmid (--ref_helper)
* AAV rep-cap (--ref_rep_cap)

The transgene plasmid ITR regions can exist in one of four orientations 
(FLIP/FLOP, FLOP/FLIP, FLIP/FLIP, or FLOP/FLOP) (refernce). Each of these four orientations are generated from the
given transgene plasmid sequence and added to the combined reference.

Reads are mapped to the combined reference using [minimap2](https://github.com/lh3/minimap2) and alignment 
statistics are then generated using [seqkit](https://bioinf.shenwei.me/seqkit/), 
which are used in the rest of the workflow.

Consensus sequence(s) are made using [clair3](https://github.com/HKU-BAL/Clair3) to generate VCF
variant files from which [bcftools](https://samtools.github.io/bcftools/bcftools.html) creates a consensus sequence. 
