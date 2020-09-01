### Tutorial - How to set up nf-core/eager for pathogen screening

🛠️ Under Construction 🛠️

#### Tutorial Pathogens - Introduction

This tutorial will give a basic example on how to set up nf-core/eager to
perform reconstruct a bacterial genome from samples in the context of ancient
pathogenomic research.

> :warning: this tutorial does not describe how to install and setup
> nf-core/eager For this please see other documentation on the
> [nf-co.re](https://nf-co.re/usage/installation) website.

We will describe how to set up mapping ancient pathogen samples against the
reference of the targeted organism genome to allow sequencing and library
quality-control, calculation of depth and breath of coverage, check for damage
profiles, snp calling and producing an SNP alignment for its usage in
phylogenetic analysis. DATA!!!!

#### Tutorial Pathogens - Preparation

Prior setting up an nf-core/eager run for pathogen reconstruction, we will need:

1. Raw sequencing data in FASTQ format
2. Reference genome in FASTA format, with associated pre-made `bwa`, `samtools`
   and `picard SequenceDictionary` indices**

**Note: This files are not required per se. nf-core/eager will generate those
for the reference fasta if they are not present, however they will be removed
after the run if not specified otherwise. If your reference genome is large, is
best to generate this files once and keep them in a separate folder.

We should also ensure we have the very latest version of the nf-core/eager
pipeline so we have all latest bugfixes etc. In this case we will be using
nf-core/eager version 2.2.0. You should always check on the
[nf-core](https://nf-co.re/eager) website  whether a newer release has been made
(particularly point releases e.g. 2.2.1).

```bash
nextflow pull nf-core/eager -r 2.2.0
```

It is important to note that if you are planning on running multiple runs of
nf-core/eager for a given project, that the version should be **kept the same**
for all runs to ensure consistency in settings for all of your libraries.

#### Tutorial Pathogens - Inputs and Outputs

To start, lets make a directory where all your nf-core/eager related files for
this run will go, and change into it.

```bash
mkdir projectX_pestis20200720
cd projectX_pestis20200720
```

The first part of constructing any nf-core/eager run is specifying a few generic
parameters that will often be common across all runs. This will be which
pipeline, version and _profile_ we will use. We will also specify a unique name
of the run to help us keep track of all the nf-core/eager runs you may be
running.

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
<...>
```

For the `-profile` parameter, I have indicated that I wish to use singularity as
my software container environment, and I will use the MPI-SHH institutional
config as listed on
[nf-core/configs](https://github.com/nf-core/configs/blob/master/conf/shh.config),
and using the profile for the 'sdag' cluster. These profiles specify settings
optimised for the specific cluster/institution, such as maximum memory available
or which scheduler queues to submit to. More explanations about configs and
profiles can be seen in the [nf-core
website](https://nf-co.re/usage/configuration) and the [profile
tutorial](../tutorial_profiles.md).

Next we need to specify our input data. nf-core/eager can accept input FASTQs
files in two main ways, either with direct paths to files (with wildcards), or
with a Tab-Separate-Value (TSV) file which contains the paths and extra
metadata. In this example, we will use the TSV method, as to simulate a
realistic use-case. DESCRIPTION DATA

Our TSV file will look something like the following:

```bash
Sample_Name     Library_ID      Lane    Colour_Chemistry        SeqType Organism        Strandedness    UDG_Treatment   R1      R2      BAM
EGR001  EGR001.B0101.SG1        1       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L001_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L001_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        2       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L002_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L002_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        3       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L003_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L003_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        4       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L004_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.1/EGR001.B0101.SG1.1_S0_L004_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        5       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L001_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L001_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        6       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L002_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L002_R2_001.fastq.gz NA
EGR001  EGR001.B0101.SG1        7       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L003_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L003_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        8       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L004_R1_001.fastq.gz ../../02-raw_data/EGR001.B0101.SG1.2/EGR001.B0101.SG1.2_S0_L004_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        1       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L001_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L001_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        2       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L002_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L002_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        3       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L003_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L003_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        4       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L004_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.1/EGR002.B0201.SG1.1_S0_L004_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        5       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L001_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L001_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        6       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L002_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L002_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        7       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L003_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L003_R2_001.fastq.gz NA
EGR002  EGR002.B0201.SG1        8       2       PE      homo_sapiens    double  half    ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L004_R1_001.fastq.gz ../../02-raw_data/EGR002.B0201.SG1.2/EGR002.B0201.SG1.2_S0_L004_R2_001.fastq.gz NA
```

You can see that we have a single line for each pair of FASTQ files representing
each `Lane`, but the `Sample_Name` and `Library_ID` columns identify and group
them together accordingly. Secondly, as we have NextSeq data, we have specified
we have two `Colour_Chemistry`, which is important for downstream processing
(see below). The other columns are less important for this particular context of
metagenomic screening. See the nf-core/eager [usage](../usage.md) documentation
for more specifications on how to set up a TSV file (e.g. why despite NextSeqs
only having 4 lanes, we go up to 8 in the example above).

Alongside our input TSV file, we will also specify the paths to our reference
FASTA file and the corresponding indices.

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
--input 'pestis20200720.tsv' \
--fasta '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta' \
--bwa_index '../Reference/Y_pestis_CO92_NC_003143/' \
--fasta_index '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta.fai' \
--seq_dict '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.dict' \
<...>
```

We specify the paths to each reference genome and it's corresponding tool
specific index. Paths should always be encapsulated in quotes to ensure Nextflow
evaluates them, rather than your shell! Also note that as `bwa` generates
multiple index files, nf-core/eager takes a _directory_ that must contain these
indices instead.

> Note the difference between single and double `-` parameters. The former
> represent nextflow flags, while double are nf-core/eager specific flags.

Finally, we can also specify the output directory and the Nextflow `work/`
directory (which contains "intermediate" working files and directories).

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
--input 'pestis20200720.tsv' \
--fasta '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta' \
--bwa_index '../Reference/Y_pestis_CO92_NC_003143/' \
--fasta_index '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta.fai' \
--seq_dict '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.dict' \
--outdir './results/' \
- w './work/' \
<...>
```

#### Tutorial Pathogens - Pipeline Configuration

Now that we have specified the input data, we can start moving onto specifying
settings for each different module we will be running. As mentioned above, we
are pretending to run with NextSeq data, which is generated with a two-colour
imaging technique. What this means is when you have shorter molecules than the
number of cycles of the sequencing chemistry, the sequencer will repeatedly see
'G' calls (no colour) at the last few cycles, and you get long poly-G 'tails' on
your reads. We therefore will turn on the poly-G clipping functionality offered
by [`fastp`](https://github.com/OpenGene/fastp), and any pairs of files
indicated in the TSV file as having `2` in the `Colour_Chemistry` column will be
passed to `fastp`. We will not change the default minimum length of a poly-G
string to be clipped.

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
--input 'pestis20200720.tsv' \
--fasta '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta' \
--bwa_index '../Reference/Y_pestis_CO92_NC_003143/' \
--fasta_index '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta.fai' \
--seq_dict '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.dict' \
--outdir './results/' \
- w './work/' \
--complexity_filter_poly_g \
<...>
```

We will now specify the mapping parameters and we will filter all unmapped reads
in order to reduce the size of the files.

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
--input 'pestis20200720.tsv' \
--fasta '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta' \
--bwa_index '../Reference/Y_pestis_CO92_NC_003143/' \
--fasta_index '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta.fai' \
--seq_dict '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.dict' \
--outdir './results/' \
- w './work/' \

```

Finally, multiVCFAnalyzer

```bash
nextflow run nf-core/eager \
-r 2.2.0 \
-profile sdag,shh,singularity \
-name 'projectX_pestis20200720' \
--input 'pestis20200720.tsv' \
--fasta '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta' \
--bwa_index '../Reference/Y_pestis_CO92_NC_003143/' \
--fasta_index '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.fasta.fai' \
--seq_dict '../Reference/Y_pestis_CO92_NC_003143/Y_pestis_CO92_NC_003143.dict' \
--outdir './results/' \
- w './work/' \

```

With this, we are ready to submit! If running on a remote cluster/server, make
sure to run this in `screen` session or similar, so that if you get a `ssh`
signal drop or want to log off, Nextflow will not crash.

#### Tutorial Pathogens - Results

Assuming the run completed without any crashes (if problems do occur, check all
against [#usage](usage.md) all parameters are as expected, or check the
[FAQ](../faq.md)), we can now check our results in `results/`.

##### Tutorial Pathogens - MultiQC Report

In here there are many different directories containing different output files.
The first directory to check is the `MultiQC/` directory. In here you should
find a `multiqc_report.html` file. You will need to view this in a web browser,
so I recommend either mounting your server to your file browser, or downloading
it to your own local machine (PC/Laptop etc.).

Once you've opened this you can go through each section and evaluate all the
results. Now we can evaluate the quality of the sequencing data and how well
preserve the pathogen genome is.

For example, I normally look for things like:

General Stats Table:

- Do I see the expected number of raw sequencing reads (summed across each set
  of FASTQ files per library) that was requested for sequencing?
- Does the percentage of trimmed reads look normal for aDNA, and do lengths
  after trimming look short as expected of aDNA?
- Does ClusterFactor or 'Dups' look high (e.g. >2 or >10% respectively)
  suggesting over-amplified or badly preserved samples?

FastQC (pre-AdapterRemoval):

- Do I see any very early drop off of sequence quality scores suggesting
  problematic sequencing run?
- Do I see outlier GC content distributions?
- Do I see high sequence duplication levels?

AdapterRemoval:

- Do I see high numbers of singletons or discarded read pairs?

FastQC (post-AdapterRemoval):

- Do I see improved sequence quality scores along the length of reads?
- Do I see reduced adapter content levels?

General Genome Stats:

- Do you have a reasonable Mean Coverage (depth)?
  - To perform phylogenetic analysis, I would like to see a mean coverage >=3X
- How well do the reads cover the reference genome (breath of coverage)?
  - The percent of the reference coverage should be the closest possible to
    100%. However, we do not expect to cover the entire chromosome, since we
    remove reads that map equally well to more than 1 locations, thus excluding
    repetitive or duplicated regions.
  - We want to observe a gradually lost or no change in the percent of the
    reference covered. If a drastic drop in the percent of coverage is observed
    this may indicate cross-mapping from a related species.

Samtools Flagstat (pre/post Filter):

- Do I see outliers, e.g. with unusually high levels of pathogen DNA,
  (indicative of cross-mapping with a closer microorganism) that require
  downstream closer assessment?

DeDup/Picard MarkDuplicates):

- Do I see large numbers of duplicates being removed, possibly indicating
  over-amplified or badly preserved samples?

DamageProfiler:

- Do I see evidence of damage on pathogen DNA?
  - If high numbers of reads mapped to the pathogen reference but no damage may
    indicate missmapping from a related modern microorganism.

> Detailed documentation and descriptions for all MultiQC modules can be seen in
> the the 'Documentation' folder of the results directory or here in the [output
> documentation](output.md)

If you're happy everything looks good in terms of sequencing, we then look at
specific directories to find any files you might want to use for downstream
processing.

Note that when you get back to writing up your publication, all the versions of
the tools can be found under the 'nf-core/eager Software Versions' section of
the MultiQC report. Note that all tools in the container are listed, so you may
have to remove some of them that you didn't actually use in the set up.

For example, in this example, we have used: Nextflow, nf-core/eager, FastQC,
AdapterRemoval, fastP, BWA, Samtools, endorS.py, Picard Markduplicates,
Qualimap, PreSeq, DamageProfiler, MALT, Maltextract and MultiQC.

Citations to all used tools can be seen
[here](https://nf-co.re/eager#tool-references)

##### Tutorial Pathogens - Files for Downstream Analysis

IGV --> visualise the bam files

MultiVCFAnalyze results: -SNP alignment --> phylogeny -Genotype/blabla -->
structure -SNP table --> filter for unique snps of the samples analysed -SNP
table for snpEff --> to obtain the effects of the identified snps

#### Tutorial Pathogens - Clean up

Finally, I would recommend cleaning up your `work/` directory of any
intermediate files (if your `-profile` does not already do so). You can do this
by going to above your `results/` and `work/` directory, e.g.

```bash
cd /<path>/<to>/projectX_screening20200720
```

and running

```bash
nextflow clean -f -k
```

#### Tutorial Pathogens - Summary

In this this tutorial we have described an example on how to set up a pathogen
genome reconstruction run of ancient human samples where we have an indication
of a specific pathogen. We have covered how set up nf-core/eager to reconstruct
the pathogen genome, check for its ancient authenticity and set up an initial
multiVCFAnalyser run to contruct a phylogenetic tree and check for unique SNPs
and their effect. Finally we have also described what to look for in the MultiQC
run summary report and where to find output files that can be used for
downstream analysis.
