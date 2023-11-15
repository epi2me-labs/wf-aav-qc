## Introduction

This is a quality control workflow for recombinant adeno-associated virus (rAAV) preps.
It provides information about the status of plasmid preparations including; contamination 
levels, integrity of individual AAV genomes, and unwanted sequence variation that may be present within the vectors.

The transgene plasmid ITR cassettes will naturally exist in one of four orientations 
To account for this, the variable regions of the supplied transgene plasmid are identified, masked and added to a 
combined reference containing the helper and Rep-Cap plasmid sequences and the host cell reference.

Reads are mapped to the combined reference using [minimap2](https://github.com/lh3/minimap2) and alignment 
statistics are then generated using [seqkit](https://bioinf.shenwei.me/seqkit/), 
which are used in the rest of the workflow to generate QC plots and tables.

Finally, transgene plasmid consensus sequences are generated using [medaka](https://github.com/nanoporetech/medaka) allowing genome integrity to be checked.
