## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[Singularity](https://docs.sylabs.io/guides/latest/user-guide/) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either Docker or Singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit our website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-aav-qc --help
```

to see the options for the workflow.

**Download demonstration data**

A small test dataset is provided for the purposes of testing the workflow software,
it can be downloaded using: (Not available yet)

```
wget -O demo_data.tar.gz \
    https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-aav-qc/demo_data.tar.gz
tar -xzvf demo_data.tar.gz
```

The workflow can be run with the demo data as flollows: (to complete)
```
OUTPUT=output
nextflow run epi2me-labs/wf-aav-qc \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq demo_data/fastq \
    --ref_host\
    --ref_helper \
    --ref_rep_cap \
    --transgene_plasmid \
    --transgene_annotation \
    --basecaller_cfg 'dna_r10.4.1_e8.2_400bps_sup@v3.5.2'  \
    --out_dir ${OUTPUT}
```

**Workflow outputs**

The primary outputs of the workflow include:

* an HTML report document detailing the primary findings of the workflow.
* a simple text file providing a summary of sequencing reads.
* a TSV file detailing the genome type read assignemnts.
* A consensus sequence generated 



