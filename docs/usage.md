# nf-core/eager: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/eager/usage](https://nf-co.re/eager/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with 3 columns, and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Supplying BAM input

It is possible to also supply BAM files as input to nf-core/eager. This can allow you to skip earlier steps of the pipeline (preprocessing and mapping) when desired - e.g. when re-processing public data. You can also convert input BAM files back to FASTQ files to re-undergo preprocessing and mapping. This may be desired when you want to standardise the mapping parameters between your own and previously published data.

You will still need to fill the `pairment` column in the input TSV sheet for the BAM files. If you do not convert the BAM files back to FASTQ, you must specify the column as `single`. If you do do the conversion, you must specify the type of reads the BAM file contains, i.e.:

- If the mapped reads in the BAM file are single end then specify `single`
- If the mapped reads in the BAM file are paired-end _but merged pairs_ (i.e. overlapping pairs collapsed to a single read), then you must also supply `single`
- If the mapped reads in the BAM file are paired-end and are _not_ merged (i.e., paired-end mapping was originally performed), then you must specify `paired`

Note that if you do not specify to merge BAM converted paired-end FASTQs (i.e., request paired-end mapping), only forward and reverse pairs will be used - singletons in the BAMs will be discarded!

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Reference input

nf-core/eager supports two methods of supplying reference FASTA files.
The first is a direct path to a single FASTA file with optional prebuilt indicies via `--fasta`, `--fasta_fai`, `--fasta_dict`, etc., and the second is via a reference sheet.

Providing a reference sheet to `--fasta_sheet` allows users to align their input reads to _multiple_ reference genomes in the same run. The reference sheet must be in the format of a comma- or tab-separated table, and the file extension must be `csv` or `tsv` respectively.

In addition to including the path to the FASTA, the reference sheet can also be used to specify paths to pre-built indices of each reference (namely, `.fai` from `samtools faidx`, `.dict` from `picard CreateSequenceDictionary`, and/or a directory pointing to a directory containing the indices for the given mapper - e.g. created with `bwa index`).

Note that passing a reference sheet to the pipeline with `--fasta_sheet` will _override_ any corresponding directly-supplied parameters specifying user-build indices (`--fasta_fai`, `--fasta_dict`).

An example of a reference sheet in `csv` format is as follows:

```txt
reference_name,fasta,fai,dict,mapper_index,circular_target, mitochondrion_header
Mammoth_MT_Krause,//<path>/<to>/data/Mammoth_MT_Krause.fasta,/<path>/<to>/data/Mammoth_MT_Krause.fasta.fai,/<path>/<to>/data/Mammoth_MT_Krause.dict,/<path>/<to>/data/,NC_007596.2,NC_007596.2
Human_MT,/<path>/<to>/data/Human_MT.fasta.gz,,,,,
```

Only the `reference_name`, and `fasta` columns are mandatory, whereas all other cells can be empty depending on your context.

| Header                      | Required | Description                                                                                                                |
| --------------------------- | -------- | -------------------------------------------------------------------------------------------------------------------------- |
| reference_name              | Yes      | Name of the reference to be used in file names and nextflow console log                                                    |
| fasta                       | Yes      | Path to FASTA of reference. Can be optionally gzipped                                                                      |
| fai                         | No       | Optional path to pre-build SAMtools `fai` index file corresponding to the FASTA                                            |
| dict                        | No       | Optional path to pre-build picard `dict` index file corresponding to the FASTA                                             |
| mapper_index                | No       | Optional path to _directory_ containing pre-build mapper index files corresponding to the FASTA                            |
| circular_target             | No       | Optional string, and only required for CircularMapper, with name of entry in FASTA to extend (up to first space in header) |
| mitochondrion_header        | No       | Optional string, and only required for MTNucRatio, with name of entry in FASTA of mitochondrion entry's header             |
| snpcapture_bed              | No       | Optional path to BED file with SNP capture positions, only required for QualiMap                                           |
| pileupcaller_bedfile        | No       | Optional path to BED file with SNP capture positions for genotyping with pileupCaller                                      |
| pileupcaller_snpfile        | No       | Optional path to EIGENSTRAT SNP panel file for genotyping with pileupCaller                                                |
| hapmap_file                 | No       | Optional path to HapMap files for contamination estimation with ANGSD                                                      |
| pmdtools_masked_fasta       | No       | Optional path to masked FASTA files for PMDtools                                                                           |
| pmdtools_bed_for_masking    | No       | Optional path to SNP capture BED file to mask the reference for PMDtools                                                   |
| sexdeterrmine_snp_bed       | No       | Optional path to SNP capture bed files for genetic sex estimation with SexDetERRmine                                       |
| bedtools_feature_file       | No       | Optional path to feature file for coverage calculation with bedtools                                                       |
| genotyping_reference_ploidy | No       | Optional integer to specify organism ploidy for genotyping with GATK or FreeBayes                                          |
| genotyping_gatk_dbsnp       | No       | Optional path to SNP annotation file for genotyping with GATK                                                              |

Files for `fai`, `dict`, `mapper_index` will be generated by the pipeline for you if not specified.

A real-world example could look as follows, where a user-supplied `.dict` file and `circular_target ` and `mitochondrion_header` are not specified:

```txt
reference_name,fasta,fai,dict,mapper_index,circular_target,mitochondrion
tannerella_forsythia,/home/james/data/reference/tannerella_forstythia.fasta,/home/james/data/reference/tannerella_forstythia.fasta.fai,,/home/james/data/reference/,,
```

> ⚠️ The names of the files in the `mapper_index` directory must have the same basename as the FASTA itself! e.g. if FASTA is Mammoth_MT_Krause.fasta, then each BWA index file should be Mammoth_MT_Krause.fasta.\<suffix\>

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/eager --input ./samplesheet.csv --outdir ./results --genome GRCh37 -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/eager -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/eager
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/eager releases page](https://github.com/nf-core/eager/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
