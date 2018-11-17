# ![nf-core/eager](docs/images/eager_logo.png)

[![Build Status](https://travis-ci.org/nf-core/eager.svg?branch=master)](https://travis-ci.org/nf-core/eager)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/eager)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/eager.svg)](https://hub.docker.com/r/nfcore/eager/)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)
[![DOI](https://zenodo.org/badge/135918251.svg)](https://zenodo.org/badge/latestdoi/135918251)



## Introduction

**nf-core/eager** is a bioinformatics best-practice analysis pipeline for ancient DNA data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FASTQ inputs, aligns the reads and performs extensive quality-control on the results. It comes with docker / singularity containers making installation trivial and results highly reproducible.

## Pipeline steps

* Create reference genome indices (optional)
    * BWA 
    * Samtools Index
    * Sequence Dictionary
* QC with FastQC
* AdapterRemoval for read clipping and merging
* Read mapping with BWA, BWA Mem or CircularMapper
* Samtools sort, index, stats & conversion to BAM
* DeDup or MarkDuplicates read deduplication
* QualiMap BAM QC Checking
* Preseq Library Complexity Estimation
* DamageProfiler damage profiling
* BAM Clipping for UDG+/UDGhalf protocols
* PMDTools damage filtering / assessment

## Quick Start

1. Install [`nextflow`](docs/installation.md)
2. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)
3. Download the EAGER pipeline

```bash
nextflow pull nf-core/eager
```

4. Set up your job with default parameters

```bash
nextflow run nf-core -profile <docker/singularity/conda> --reads'*_R{1,2}.fastq.gz' --fasta '<REFERENCE.fasta'
```

5. See the overview of the run with under `<OUTPUT_DIR>/MultiQC/multiqc_report.html`

## Documentation
The nf-core/eager pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. Pipeline configuration
    * [Local installation](docs/configuration/local.md)
    * [Adding your own system](docs/configuration/adding_your_own.md)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [Troubleshooting](docs/troubleshooting.md)

## Credits
This pipeline was written by Alexander Peltzer ([apeltzer](https://github.com/apeltzer)), with major contributions from Stephen Clayton, ideas and documentation from James Fellows Yates, Raphael Eisenhofer and Judith Neukamm. If you want to contribute, please open an issue and ask to be added to the project - happy to do so and everyone is welcome to contribute here!