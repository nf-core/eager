# ![nf-core/EAGER2](docs/images/EAGER2_logo.png)

[![Build Status](https://travis-ci.org/nf-core/EAGER2.svg?branch=master)](https://travis-ci.org/nf-core/EAGER2)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.27.0-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/EAGER2)

[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/eager2.svg)](https://hub.docker.com/r/nfcore/eager2/)
![Singularity Container available](
https://img.shields.io/badge/singularity-available-7E4C74.svg)

## Introduction
THIS PIPELINE IS A WORK IN PROGRESS. Thanks for checking it out! Hopefully it will be functional soon.

**nf-core/EAGER2** is a bioinformatics best-practice analysis pipeline for ancient DNA data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FastQ inputs, aligns the reads and performs extensive quality-control on the results. It comes with docker / singularity containers making installation trivial and results highly reproducible.

### Pipeline steps

* Make BWA reference genome index (optional)
* FastQC
* Clip & Merge / AdapterRemoval for read clipping and merging
* Align with BWA
* Samtools sort, index, stats & conversion to BAM/CRAM
* Samtools idxstats
* DeDup read deduplication or MarkDuplicates
* QualiMap
* Preseq estimation
* DamageProfiler damage profiling
* Schmutzi contamination assessment
* Genotyping using GATK or ANGSD or sequenceTools (optional)
* Consensus generation with VCF2Genome (optional)


### Documentation
The nf-core/EAGER2 pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

### Credits
This pipeline was written by Alexander Peltzer ([apeltzer](https://github.com/apeltzer)).
