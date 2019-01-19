# nf-core/eager: Output

## Table of contents
* [Introduction](#introduction)
* [Directory Structure](#directory-structure)
* [Module Overview](#module-overview)

## Introduction
The output of EAGER2 consists of two main components: output files (e.g. BAM or FASTQ files), and summary statistics of the whole run presented in a [`MultiQC`](https://multiqc.info) report. (Some) intermediate files and module-specific statstics files are also retained.

## Directory Structure
The directory structure of EAGER2 is as follows

```
<RUN_OUTPUT_DIRECTORY>/
├── MultiQC/
├── <MODULE_1>/
├── <MODULE_2>/
├── <MODULE_3>/
├── pipeline_info/
├── reference_genome/
└── work/
```

* The parent directory `<RUN_OUTPUT_DIRECTORY` is the parent directory of the run, either the directory the pipeline was run from or as specified by the `--outdir` flag.

**Primary Output Directories**
These directories are the ones you will use on a day-to-day basis and are those which you should familiarise yourself with.

* The `MultiQC` directory is the most important directory and contains the main summary report of the run in HTML format, which can be viewed in a web-browser of your choice. The sub-diectory contains the MultiQC collected data used to build the HTML report. The Report allows you to get an overview of the sequencing and mapping quality as well as aDNA metrics. 
* A `<MODULE_1>` directory contains the (cleaned-up) output from a particular software module. This is the second most important set of directories. This contains output files such as FASTQ, BAM, statistics, and/or plot files of a specific module. The latter two are only needed when you need finer detail about that particular module.

**Secondary Output Directories**
These are less important directories which are used less often, normally in the context of bug-reporting.

* `pipeline_info` contains back-end reporting of the pipeline itself such as run times and computational statistics. You rarely need this information other than for curiosity or when bug-reporting.
* `reference_genome` contains either text files describing the location of specified reference genomes, and if not already supplied when running the pipeline, auxilary indexing files. This is often useful when re-running other samples using the same reference genome, but is otherwise often not otherwise important.
* The `work` directory contains all the `nextflow` processing directories. This is where `nextflow` actually does all the work, but in an efficient programatic procedure that is not intuitive to human-readers. Due to this, the directory is often not important to a user as all the useful output files are linked to the module directories (see above). Otherwise, this directory maybe useful when a bug-reporting.


## Module overview

In this section we will run through the output of each module as reported in a MultiQC output. This can be viewed by opening the HTML file in your `<RUN_OUTPUT_DIRECTORY>/MultiQC/` directory in a web browser. The section will also provide some basic tips on how to interpret the plots and values, although we highly recommend reading the READMEs or original papers of the tools used in the pipeline. A list of references can be seen on the [EAGER2 github repository](https://github.com/nf-core/eager/)

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## FastP
## AdapterRemoval
## BWA 
## Picard
## Samtools
## DeDup
## QualiMap
## DamageProfiler
## PMDTools
## Preseq
## BamUtil

