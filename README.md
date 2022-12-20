# ![nf-core/eager](docs/images/nf-core_eager_logo_outline_drop.png)

**A fully reproducible and state-of-the-art ancient DNA analysis pipeline.**

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/eager/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/eager)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23eager-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/eager)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/eager** is a scalable and reproducible bioinformatics best-practise processing pipeline for genomic NGS sequencing data, with a focus on ancient DNA (aDNA) data. It is ideal for the (palaeo)genomic analysis of humans, animals, plants, microbes and even microbiomes.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/eager/results).

## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

- (Optionally) create reference genome indices for mapping (`bwa`, `samtools`, and `picard`)
- Sequencing quality control (`FastQC`)
- Sequencing adapter removal, paired-end data merging (`AdapterRemoval`)
- Read mapping to reference using (`bwa aln`, `bwa mem`, `CircularMapper`, or `bowtie2`)
- Post-mapping processing, statistics and conversion to bam (`samtools`)
- Ancient DNA C-to-T damage pattern visualisation (`DamageProfiler`)
- PCR duplicate removal (`DeDup` or `MarkDuplicates`)
- Post-mapping statistics and BAM quality control (`Qualimap`)
- Library Complexity Estimation (`preseq`)
- Overall pipeline statistics summaries (`MultiQC`)

### Additional Steps

Additional functionality contained by the pipeline currently includes:

#### Input

- Automatic merging of complex sequencing setups (e.g. multiple lanes, sequencing configurations, library types)

#### Preprocessing

- Illumina two-coloured sequencer poly-G tail removal (`fastp`)
- Post-AdapterRemoval trimming of FASTQ files prior mapping (`fastp`)
- Automatic conversion of unmapped reads to FASTQ (`samtools`)
- Host DNA (mapped reads) stripping from input FASTQ files (for sensitive samples)

#### aDNA Damage manipulation

- Damage removal/clipping for UDG+/UDG-half treatment protocols (`BamUtil`)
- Damaged reads extraction and assessment (`PMDTools`)
- Nuclear DNA contamination estimation of human samples (`angsd`)

#### Genotyping

- Creation of VCF genotyping files (`GATK UnifiedGenotyper`, `GATK HaplotypeCaller` and `FreeBayes`)
- Creation of EIGENSTRAT genotyping files (`pileupCaller`)
- Creation of Genotype Likelihood files (`angsd`)
- Consensus sequence FASTA creation (`VCF2Genome`)
- SNP Table generation (`MultiVCFAnalyzer`)

#### Biological Information

- Mitochondrial to Nuclear read ratio calculation (`MtNucRatioCalculator`)
- Statistical sex determination of human individuals (`Sex.DetERRmine`)

#### Metagenomic Screening

- Low-sequenced complexity filtering (`BBduk`)
- Taxonomic binner with alignment (`MALT`)
- Taxonomic binner without alignment (`Kraken2`)
- aDNA characteristic screening of taxonomically binned data from MALT (`MaltExtract`)

#### Functionality Overview

A graphical overview of suggested routes through the pipeline depending on context can be seen below.

<p align="center">
    <img src="docs/images/eager2_metromap_complex.png" alt="nf-core/eager metro map" width="70%"
</p>

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/eager -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

      <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run nf-core/eager --input samplesheet.csv --outdir <OUTDIR> --reference '<your_reference>.fasta' -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

5. Once your run has completed successfully, clean up the intermediate files.

   ```bash
   nextflow clean -f -k
   ```

## Documentation

The nf-core/eager pipeline comes with documentation about the pipeline [usage](https://nf-co.re/eager/usage), [parameters](https://nf-co.re/eager/parameters) and [output](https://nf-co.re/eager/output).

## Credits

This pipeline was established by Alexander Peltzer ([apeltzer](https://github.com/apeltzer)) and [James A. Fellows Yates](https://github.com/jfy133). Version two had major contributions from [Stephen Clayton](https://github.com/sc13-bioinf), [Thiseas C. Lamnidis](https://github.com/TCLamnidis), [Maxime Borry](https://github.com/maxibor), [Zandra Fagernäs](https://github.com/ZandraFagernas), [Aida Andrades Valtueña](https://github.com/aidaanva) and [Maxime Garcia](https://github.com/MaxUlysse) and the nf-core community.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Alex Hübner](https://github.com/alexhbnr)
- [Alexandre Gilardet](https://github.com/alexandregilardet)
- Arielle Munters
- [Åshild Vågene](https://github.com/ashildv)
- [Charles Plessy](https://github.com/charles-plessy)
- [Elina Salmela](https://github.com/esalmela)
- [Fabian Lehmann](https://github.com/Lehmann-Fabian)
- [He Yu](https://github.com/paulayu)
- [Hester van Schalkwyk](https://github.com/hesterjvs)
- [Ido Bar](https://github.com/IdoBar)
- [Irina Velsko](https://github.com/ivelsko)
- [Işın Altınkaya](https://github.com/isinaltinkaya)
- [Johan Nylander](https://github.com/nylander)
- [Jonas Niemann](https://github.com/NiemannJ)
- [Katerine Eaton](https://github.com/ktmeaton)
- [Kathrin Nägele](https://github.com/KathrinNaegele)
- [Kevin Lord](https://github.com/lordkev)
- [Luc Venturini](https://github.com/lucventurini)
- [Mahesh Binzer-Panchal](https://github.com/mahesh-panchal)
- [Marcel Keller](https://github.com/marcel-keller)
- [Megan Michel](https://github.com/meganemichel)
- [Pierre Lindenbaum](https://github.com/lindenb)
- [Pontus Skoglund](https://github.com/pontussk)
- [Raphael Eisenhofer](https://github.com/EisenRa)
- [Roberta Davidson](https://github.com/roberta-davidson)
- [Rodrigo Barquera](https://github.com/RodrigoBarquera)
- [Selina Carlhoff](https://github.com/scarlhoff)
- [Torsten Günter](https://bitbucket.org/tguenther)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#eager` channel](https://nfcore.slack.com/channels/eager) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->

# <!-- If you use  nf-core/eager for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

If you use nf-core/eager for your analysis, please cite it using the following:

> Fellows Yates JA, Lamnidis TC, Borry M, Valtueña Andrades A, Fagernäs Z, Clayton S, Garcia MU, Neukamm J, Peltzer A. 2021. Reproducible, portable, and efficient ancient genome reconstruction with nf-core/eager. PeerJ 9:e10947. DOI: [10.7717/peerj.10947](https://doi.org/10.7717/peerj.10947).

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
