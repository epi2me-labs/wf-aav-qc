# wf-aav-qc

The following provides an overview of the workflow. 
( This document is a work in progress) .


## Make reference
The transgene plasmid
#### Transgene plasmid with ITR-ITR orientations

#### Combined reference file with host, helper plasmid, rep cap plasmid and transgene plasmid



## Analyse
#### Map to reference

#### Get per-read stats (seqkit bam)

#### Per-base coverage (samtools depth)


## Consensus sequence
#### Create and map to orientation references

#### Call variants (clair3)

#### Call SV's (sniffles2)

#### Consensus generation (bcftools)

#### Map regions of original ref to consensus 

#### Accuracy (seqkit bam)

# Plots
## Truncations
The plot shows the frequncy of start and end possitions for reads that are fuly contained within the ITR-ITR region.

