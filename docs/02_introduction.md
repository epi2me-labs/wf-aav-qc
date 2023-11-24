This workflow takes reads sequenced from adeno-associated virus (rAAV) vector preps and 
does some basic quality control checks. The main stages of the workflow are: 


+ Masking of transgene cassette variable ITR regions
+ Mapping reads to a combined reference sequence containing host cell reference genome, mask transgene plasmid and other AAV plasmids
+ Identification of reference which each reads maps to
+ Truncation hotspot identification
+ ITR transgene cassette coverage
+ Determination of AAV genome structure types
+ Calling transgene plasmid variants and creation of consensus 
