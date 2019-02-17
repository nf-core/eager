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

### FastQC
#### Background
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

### FastP
### AdapterRemoval
#### Background
AdapterRemoval a tool that does the post-sequencing clean up of your sequencing reads. It performs the following functions
  - 'Merges' (or 'collapses') forward and reverse reads of Paired End data
  - Removes remaining library indexing adapters
  - Trims low quality base tails from ends of reads
  - Removes too-short reads

In more detail merging is where the same read from the forward and reverse files of a single library (based on the flowcell coordinates), are compared to find a stretch of sequence that are the same. If this overlap reaches certain quality thresholds, the two reads are 'collapsed' into a single read, with the base quality scores are updated accordingly accounting for the increase quality call precision.

Adapter removal involves finding overlaps at the 5' and 3' end of reads for the artificial NGS library adapters (which connect the DNA molecule insert, and the index), and stretches that match each other are then removed from the read itself.  Note, by default AdapterRemoval does _not_ remove 'internal barcodes' (between insert and the adapter), so these statistics are not considered.

Quality trimming (or 'truncating') involves looking at ends of reads for low-confidence bases (i.e. where the FASTQ Phred score is below a certain threshold). These are then removed remove the read.

Length filtering involves removing any read that does not reach the number of bases specified by a particular value. 

#### Retained and Discarded Reads Plot

These stacked bars plots are unfortunately a little confusing, when displayed in MultiQC. However are relatively straight-forward once you understand each category. They can be displayed as counts of reads per AdapterRemoval read-category, or as percentages of the same values. Each forward(/reverse) file combination are displayed once.

The most important value is the **Retained Read Pairs** which gives you the final number of reads output into the file that goes into mapping. Note, however, this section of the stack bar _includes_ the other categories displayed (see below) in the calculation.

Other Categories:
  * If paired-end, the **Singleton [mate] R1(/R2)** cateogries represent reads which were unable to be collapsed, possibly due to the reads being too long to overlap.
  * If paired-end, **Full-length collapsed pairs** are reads which were collapsed and did not require low-quality bases at end of reads to be removed.
  * If paired-end, **Truncated collapsed pairs** are paired-end that were collapsed but did required the removal of low quality bases at the end of reads.
  * **Discarded [mate] R1/R2** represent reads which were a part of a pair, but one member of the pair did not reach other quality criteria and weas discarded. However the other member of the pair is still retained in the output file as it still reached other quality critea.
  
For ancient DNA, assuming a good quality run, you expect to see a the vast majority of your reads overlapping because we have such fragmented molecules. Large numbers of singletons suggest your molecules are too long and may not represent true ancient DNA. 

If you see high numbers of discarded or truncated reads, you should check your FastQC results for low sequencing quality of that particular run.

#### Length Distribution Plot

The length distribution plots show the number of reads at each read-length. You can change the plot to display different cateogories.

  * All represent the overall distribution of reads. In the case of paired-end sequencing You may see a peak at the turn around from forward to reverse cycles.
  * **Mate 1** and **Mate 2** represents the length of the forward and reverse read respectively prior collapsing
  * **Singleton** represent those reads that had a one member of a pair discarded
  * **Collapsed** and **Collapsed Truncated** represent reads that overlapped and able to merge into a single read, with the latter including base-quality trimming off ends of reads. These plots will start with a vertical rise representing where you are above the minimum-read threshold you set.
  * **Discarded** here represents the number of reads that did not each the read length filter. You will likely see a vertical drop at what your threshold was set to.

With paired-end ancient DNA sequencing runs You expect to see a slight increase in shorter fragments in the reverse (R2) read, as our fragments are so short we often don't reach the maximum number of cycles of that particular sequencing run. 

### Samtools
### DeDup
### QualiMap
### DamageProfiler
#### Background
DamageProfiler is a tool which calculates a variety of standard 'aDNA' metrics from a BAM file. The primary plots here are the misincorporation and length distribution plots. Ancient DNA undergoes depurination and hydrolysis, causing fragmentation of molecules into gradually shorter fragments, and cytosine to thymine deamination damage, that occur on the subsequent single-stranded overhangs at the ends of molecules.

Therefore, three main characteristics of ancient DNA are:
  * Short DNA fragments
  * Elevated G and As (purines) just before strand breaks
  * Increased C and Ts at ends of fragments
  
#### Misincorporation Plots
The MultiQC DamageProfiler module misincorporation plots shows the percent frequency of C to T mismatches at 5' read ends and complementary G to A mismatches at the 3' ends. The mismatches are when compared to the base of the reference genome at that position. 

When looking at the misincorporation plots, keep in mind the following:
  * As few-base single-stranded overhangs are more likely to occur than long overhangs, we expect to see a gradual decrease in the frequency of the modifications from position 1 to the inside of the reads. 
  * If your library has been **partially-UDG treated**, only the first one or two bases will display the the misincorporation frequency.
  * If your library has been **UDG treated** you will expect to see extremely-low to no misincorporations at read ends.
  * If your library is **single-stranded**, you will expect to see only C to T misincorporations at both 5' and 3' ends of the fragments.
  * We generally expect that the older the sample, or the less-ideal preservational environtment (hot/wet) the greater the frequency of C to T/G to A.
  * The curve will be not smooth then you have few reads informing the frequency calculation. Read counts of less than 500 are likely not reliable.

> **NB:** An important difference to note compared to the MapDamage tool, which DamageProfiler is an exact-reimplmentation of, is that the percent frequency on the Y axis is not fixed between 0 and 0.3, and will 'zoom' into small values the less damage there is

#### Length Distribution
The MultiQC DamageProfiler module length distribution plots show the frequency of read lengths across forward and reverse reads respectively.

When looking at the length distribution plots, keep in mind the following:
  * Your curves will likely not start at 0, and will start wherever your minimum read-length setting was when removing adapters.
  * You should typically see the bulk of the distribution falling between 40-120bp, which is normal for aDNA
  * You may see large peaks at paired-end turn arounds, due to very-long reads that could not overlap for merging being present, however this reads are normally from modern contamination.
  
### PMDTools
### Preseq
### BamUtil

