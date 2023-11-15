## Quickstart

# Description
 A quality control workflow for recombinant adeno-associated virus (rAAV) preps.

# Summary 

In summary this workflow has the following steps:

+ Masking of transgene cassette variable ITR regions
+ Mapping reads to a combined reference sequence containing host cell reference genome and rAAV plasmids
+ Read mapping reference identity
+ Identification of truncation hotspots
+ ITR transgene cassette coverage
+ Determination of AAV genome structure types
+ Calling transgene plasmid variants and creation of consensus 

# Compute requirements

Compute requirements are based on running 6 samples with approximately 150k reads each

Minimum requirement

+ cpus = 4
+ memory = 16 GB

Recommended requirement

+ cpus 8
+ memory 32 GB
+ GPU for medaka variant calling / consensus generation

If you are running more samples, providing more memory and cpus will speed up the workflow.

# Install and run

## Nextflow requirements (not in app)

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.0/user-guide/index.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

### Run the workflow (not in app)   

To obtain the workflow, having installed `nextflow`, users can run:
```
nextflow run epi2me-labs/wf-aav-qc --help
```
Demo data:
```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-aav-qc/wf-aav-qc-demo.tar.gz
tar -xzvf wf-aav-qc-demo.tar.gz
```

    Note: The demo data is currently synthetic data containing randomly geenrated sequences, which currently 
    serves to check the workflow functions as expected. Expect the results gernated from using real data to look different.

Example cmd:
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
    --itr2_end 2286
```

# Related protocols
https://nanoporetech.atlassian.net/wiki/spaces/PMO/pages/264476502/New+AAV+sequencing+SQK-NBD114.24+end+to+end+workflow+release

# Inputs

* `fastq or bam` -  A fastq/BAM file or directory containing fastq/bam input files or directories of input files. Here is an example of the input formats that will be accepted.
```
reads0.fastq  ─── input_directory        ─── input_directory
                 ├── reads0.fastq            ├── barcode01
                 └── reads1.fastq            │   ├── reads0.fastq
                                             │   └── reads1.fastq
                                             ├── barcode02
                                             │   ├── reads0.fastq
                                             │   ├── reads1.fastq
                                             │   └── reads2.fastq
                                             └── barcode03
                                                 └── reads0.fastq
```


* `ref_host` - The reference genomic sequence of the host cell line used to create the rAAV prep (fastq/fast.gz).
* `ref_helper` - The helper plasmid sequence (fasta).
* `ref_rep_cap` - the Rep-Cap plasmid sequence (fasta).
* `ref_transgene_plasmid` - The transgene plasmid sequence (fasta).
* `basecaller_cfg` - Name of the basecaller model that processed the signal data; used to select an appropriate Medaka model.
* `itr1_start` - start position (0-based) of the ITR1.
* `itr1_end` - end position (0-based) of ITR1.
* `itr2_start` - start position (0-based) of the ITR2.
* `itr2_end` - end position (0-based) of ITR2.

# Outputs
The primary outputs of the workflow include:

* An HTML report document detailing the primary findings of the workflow.
* A simple text file providing a summary of sequencing reads.
* A TSV file detailing the genome type read assignments.
* A VCF file of transgene plasmid variants.
* Transgene plasmid consensus sequence.


# Pipeline overview
To see a detailed overview of the workflow see the [Workflow Overview](./docs/workflow_overview.md) at docs/workflow_overview.md


# Troubleshooting
If the ITR cassette locations are split across the reference, the wrokflow will currently
not be able to process the data correctly. For example:
```
itr1_start = 4000
itr1_end = 4156
itr2_start = 100
itr2_end = 156
```
In this case, the transgene reference sequence will have to be edited to make the
ITR cassette span an unbroken range of the reference sequence.


# FAQs and tips
