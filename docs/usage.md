# nf-core/eager: Usage

## Table of contents

<!-- Install Atom plugin markdown-toc-auto for this ToC to auto-update on save -->
<!-- TOC START min:2 max:3 link:true asterisk:true update:true -->
* [Table of contents](#table-of-contents)
* [Introduction](#general-nextflow-info)
* [Running the pipeline](#running-the-pipeline)
* [Updating the pipeline](#updating-the-pipeline)
* [Reproducibility](#reproducibility)
* [Main arguments](#mandatory-arguments)
* [Other command line parameters](#other-command-line-parameters)
* [Adjustable parameters for nf-core/eager](#adjustable-parameters-for-nf-coreeager)
* [Automatic resubmission](#automatic-resubmission)
* [Clean up](#clean-up)

## General Nextflow info
Nextflow handles job submissions on SLURM or other environments, and supervises running the jobs. Thus the Nextflow process must run until the pipeline is finished. We recommend that you put the process running in the background through `screen` / `tmux` or similar tool. Alternatively you can run nextflow within a cluster job submitted your job scheduler.

To create a screen session:

```bash
screen -R eager2
```
To disconnect, press `ctrl+a` then `d`.

To reconnect, type :

```bash
screen -r eager2
```
to end the screen session while in it type `exit`.


## Help Message
To access the nextflow help message run: `nextflow run -help`

## Running the pipeline

> Before you start you should change into the output directory you wish your results to go in. When you start the nextflow job, it will place all the 'working' folders in the current directory and NOT necessarily the directory the output files will be in.

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/eager --reads '*_R{1,2}.fastq.gz' --fasta 'some.fasta' -profile standard,docker
```
where the reads are from libraries of the same pairing.

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow.log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

To see the the EAGER pipeline help message run: `nextflow run nf-core/eager --help`

### Updating the pipeline
When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/eager
```
See [below](#other-command-line-parameters) for more details about EAGER2 versioning.

## Mandatory Arguments

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different computing environments (e.g. schedulers, software environments, memory limits etc). Note that multiple profiles can be loaded, for example: `-profile standard,docker` - the order of arguments is important! The first entry takes precendence over the others, e.g. if a setting is set by both the first and second profile, the first entry will be used and the second entry ignored. 

> *Important*: If running EAGER2 on a cluster - ask your system administrator what profile to use.

For more details on how to set up your own private profile, please see [installation](../configuration/adding_your_own.md).

**Basic profiles**
These are basic profiles which primarily define where you derive the pipeline's software packages from. These are typically the profiles you would use if you are running the pipeline on your **own PC** (vs. a HPC cluster - see below).

* `awsbatch`
  * A generic configuration profile to be used with AWS Batch.
* `conda`
  * A generic configuration profile to be used with [conda](https://conda.io/docs/)
  * Pulls most software from [Bioconda](https://bioconda.github.io/)
* `docker`
  * A generic configuration profile to be used with [Docker](http://docker.com/)
  * Pulls software from dockerhub: [`nfcore/eager`](http://hub.docker.com/r/nfcore/eager/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](http://singularity.lbl.gov/)
    * Pulls software from singularity-hub
* `conda`
    * A generic configuration profile to be used with [conda](https://conda.io/docs/)
    * Pulls most software from [Bioconda](https://bioconda.github.io/)
 * `awsbatch`
    * A generic configuration profile to be used with AWS Batch.
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `none`
    * No configuration at all. Useful if you want to build your own config from scratch and want to avoid loading in the default `base` config profile (not recommended).
 
**Institution Specific Profiles**
These are profiles specific to certain **HPC clusters**, and are centrally maintained at [nf-core/configs](https://github.com/nf-core/configs). Those listed below are regular users of EAGER2, if you don't see your own institution here check the [nf-core/configs](https://github.com/nf-core/configs) repository.

* `uzh`
    * A profile for the University of Zurich Research Cloud
    * Loads Singularity and defines appropriate resources for running the pipeline.
* `binac`
    * A profile for the BinAC cluster at the University of Tuebingen
    * Loads Singularity and defines appropriate resources for running the pipeline
* `shh`
   * A profiler for the SDAG cluster at the Department of Archaeogenetics of the Max-Planck-Institute for the Science of Human History
   * Loads Singularity and defines appropriate resources for running the pipeline

    

### `--reads`
Use this to specify the location of your input FastQ files. The files maybe either from a single, or multiple samples. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```
for a single sample, or

```bash
--reads 'path/to/data/*/sample_*_{1,2}.fastq'
```

for multiple samples, where each sample's FASTQs are in it's own directory (indicated by the first `*`).

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

**Note**: It is not possible to run a mixture of single-end and paired-end files in one run.

### `--singleEnd`
If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads 'path/to/data/*.fastq'
```
for a single sample, or

```bash
--singleEnd --reads 'path/to/data/*/*.fastq'
```

for multiple samples, where each sample's FASTQs are in it's own directory (indicated by the first `*`)

**Note**: It is not possible to run a mixture of single-end and paired-end files in one run.

### `--pairedEnd`
If you have paired-end data, you need to specify `--pairedEnd` on the command line when you launc hthe pipeline. 

A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--pairedEnd --reads '*.fastq'
```

### `--fasta`
You specify the full path to your reference genome here. The FASTA file can have any file suffix, such as `.fasta`, `.fna`, `.fa`, `.FastA` etc. You may also supply a gzipped reference files, which will be unzipped automatically for you. 

For example:

```bash
--fasta '/<path>/<to>/my_reference.fasta'
```
> If you don't specify appropriate `--bwa_index`, `--fasta_index` parameters (see [below](#optional-reference-options)), the pipeline will create these indices for you automatically. Note that you can save the indices created for you for later by giving the `--saveReference` flag.

### `--large_ref`

This parameter is required to be set for large reference genomes. If your reference genome is larger than 3.5GB, the `samtools index` calls in the pipeline need to generate `CSI` indices instead of `BAI` indices to accompensate for the size of the reference genome. This parameter is not required for smaller references (including a human `hg19` or `grch37`/`grch38` reference), but `>4GB` genomes have been shown to need `CSI` indices. 

### `--genome` (using iGenomes)

The pipeline config files come bundled with paths to the illumina iGenomes reference index files. If running with docker or AWS, the configuration is set up to use the [AWS-iGenomes](https://ewels.github.io/AWS-iGenomes/) resource.

There are 31 different species supported in the iGenomes references. To run the pipeline, you must specify which to use with the `--genome` flag.

You can find the keys to specify the genomes in the [iGenomes config file](../conf/igenomes.config). Common genomes that are supported are:

* Human
  * `--genome GRCh37`
  * `--genome GRCh38`
* Mouse
  * `--genome GRCm38`
* _Drosophila_
  * `--genome BDGP6`
* _S. cerevisiae_
  * `--genome 'R64-1-1'`

> There are numerous others - check the config file for more.

Note that you can use the same configuration setup to save sets of reference files for your own use, even if they are not part of the iGenomes resource. See the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for instructions on where to save such a file.

The syntax for this reference configuration is as follows:

```nextflow
params {
  genomes {
    'GRCh37' {
      fasta   = '<path to the genome fasta file>' // Used if no star index given
    }
    // Any number of additional genomes, key is used with --genome
  }
}
```

## Optional Reference Options

### Generating Fresh Indices

#### `--saveReference`

Use this if you do not have pre-made reference FASTA indices for `bwa`, `samtools` and `picard`. If you turn this on, the indices EAGER2 generates for you will be stored in the `<your_output_dir>/results/reference_genomes` for you. 

### Premade Indices

Supplying pre-made indices saves time in pipeline execution and is especially advised when running multiple times on the same cluster system for example. You can even add a resource specific profile that sets paths to pre-computed reference genomes, saving even time when specifying these.

#### `--bwa_index`

If you want to use pre-existing `bwa index` indices, please supply the path **and file** to the FASTA you also specified in `--fasta` (see above). EAGER2 will automagically detect the index files by searching for the FASTA filename with the corresponding `bwa` index file suffixes.

For example:

```
nextflow run nf-core/eager \
-profile test_fna,docker \
--pairedEnd \
--reads *{R1,R2}*.fq.gz
--fasta results/reference_genome/bwa_index/BWAIndex/Mammoth_MT_Krause.fasta \
--bwa_index results/reference_genome/bwa_index/BWAIndex/Mammoth_MT_Krause.fasta
```

> `bwa index` does not give you an option to supply alternative suffixes/names for these indices. Thus, the file names generated by this command _must not_ be changed, otherwise EAGER2 will not be able to find them.

#### `--seq_dict`

If you want to use a pre-existing `picard CreateSequenceDictionary` dictionary file, use this to specify the required `.dict` file for the selected reference genome.

For example:

```
--seq_dict Mammoth_MT_Krause.dict
```

#### `--fasta_index`

If you want to use a pre-existing `samtools faidx` index, Use this to specify the required FASTA index file for the selected reference genome. This should be generated by `samtools faidx` and has a file suffix of `.fai`

For example:

```
--fasta_index Mammoth_MT_Krause.fasta.fai
```


## Other command line parameters

### `-r`
By default, EAGER2 will use the latest version of the pipeline that is downloaded on your system. However, it's a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/eager releases page](https://github.com/nf-core/eager/releases) and find the latest version number - numeric only (eg. `2.0`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 2.0`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

Additionally, EAGER pipeline releases are named after Swabian German Cities. The first release V2.0 is named "Kaufbeuren". Future releases are named after cities named in the [Swabian league of Cities](https://en.wikipedia.org/wiki/Swabian_League_of_Cities).

### `--outdir`
The output directory where the results will be saved.

### `-w / -work-dir`
The output directory where _intermediate_ files will be saved. It is **highly recommended** that this is the same path as `--outdir`, otherwise you may 'lose' your intermediate files if you need to re-run a pipeline. By default, if this flag is not given, 
the intermediate files will be saved in a `work/` and `.nextflow/` directory from wherever you have run EAGER from.

### `--max_memory`
Use to set a top-limit for the default memory requirement for each process.
Should be a string in the format integer-unit. eg. `--max_memory '8.GB'`. If not specified, will be taken from the configuration in the `-profile` flag.

### `--max_time`
Use to set a top-limit for the default time requirement for each process.
Should be a string in the format integer-unit. eg. `--max_time '2.h'`. If not specified, will be taken from the configuration in the `-profile` flag.

### `--max_cpus`
Use to set a top-limit for the default CPU requirement for each **process**. This is not the maximum number of CPUs that can be used for the whole pipeline, but the maximum number of CPUs each program can use for each program submission (known as a process). Do not set this higher than what is available on your workstation or computing node can provide. If you're unsure, ask your local IT administrator for details on compute node capabilities! 
Should be a string in the format integer-unit. eg. `--max_cpus 1`. If not specified, will be taken from the configuration in the `-profile` flag.

### `--email`
Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits. If set in your user config file (`~/.nextflow/config`) then you don't need to specify this on the command line for every run.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

### `-c`
Specify the path to a specific nextflow config file (this is a core NextFlow command).

**NB:** Single hyphen (core Nextflow option)

Note - you can use this to override pipeline defaults.

### `--custom_config_version`
Provide git commit id for custom Institutional configs hosted at `nf-core/configs`. This was implemented for reproducibility purposes. Default is set to `master`.

```bash
## Download and use config file with following git commid id
--custom_config_version d52db660777c4bf36546ddb188ec530c3ada1b96
```
### `--plaintext_email`
Set to receive plain-text e-mails instead of HTML formatted.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file. MultiQC produces final pipeline reports.

# Adjustable parameters for nf-core/eager

This part of the documentation contains a list of user-adjustable parameters in nf-core/eager. You can specify any of these parameters on the command line when calling the pipeline by simply prefixing the respective parameter with a double dash `--`

## Step skipping parameters

Some of the steps in the pipeline can be executed optionally. If you specify specific steps to be skipped, there won't be any output related to these modules.

### `--skip_fastqc`

Turns off FastQC pre- and post-Adapter Removal, to speed up the pipeline. Use of this flag is most common when data has been previously pre-processed and the post-Adapter Removal mapped reads are being re-mapped to a new reference genome.

### `--skip_adapterremoval`

Turns off adaptor trimming and paired-end read merging. Equivalent to setting both `--skip_collapse` and `--skip_trim`.

### `--skip_preseq`

Turns off the computation of library complexity estimation.  

### `--skip_deduplication`

Turns off duplicate removal methods DeDup and MarkDuplicates respectively. No duplicates will be removed on any data in the pipeline.

### `--skip_damage_calculation`

Turns off the DamageProfiler module to compute DNA damage profiles. 

### `--skip_qualimap`

Turns off QualiMap and thus does not compute coverage and other mapping metrics.


## Complexity Filtering Options

### `--complexity_filter_poly_g`

Performs a poly-G tail removal step in the beginning of the pipeline, if turned on. This can be useful for trimming ploy-G tails from short-fragments sequenced on two-colour Illumina chemistry such as NextSeqs (where no-fluorescence is read as a G on two-colour chemistry), which can inflate reported GC content values.

### `--complexity_filter_poly_g_min`

This option can be used to define the minimum length of a poly-G tail to begin low complexity trimming. By default, this is set to a value of `10` unless the user has chosen something specifically using this option.

## Adapter Clipping and Merging Options

These options handle various parts of adapter clipping and read merging steps.

### `--clip_forward_adaptor` 

Defines the adapter sequence to be used for the forward read. By default, this is set to `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC`.

### `--clip_reverse_adaptor`

Defines the adapter sequence to be used for the reverse read in paired end sequencing projects. This is set to `AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA` by default.

### `--clip_readlength` 30

Defines the minimum read length that is required for reads after merging to be considered for downstream analysis after read merging. Default is `30`.

### `--clip_min_read_quality` 20
Defines the minimum read quality per base that is required for a base to be kept. Individual bases at the ends of reads falling below this threshold will be clipped off. Default is set to `20`. 

### `--clip_min_adap_overlap` 1
Sets the minimum overlap between two reads when read merging is performed. Default is set to `1` base overlap.

### `--skip_collapse`

Turns off the paired-end read merging.

For example
```bash
--pairedEnd --skip_collapse  --reads '*.fastq'
```

### `--skip_trim`

Turns off adaptor and quality trimming.

For example:
```bash
--pairedEnd --skip_trim  --reads '*.fastq'
```

## Read Mapping Parameters

## BWA (default)

These parameters configure mapping algorithm parameters. 

### `--bwaalnn`

Configures the `bwa aln -n` parameter, defining how many mismatches are allowed in a read. By default set to `0.04`, if you're uncertain what to set check out [this](https://apeltzer.shinyapps.io/bwa-mismatches/) Shiny App for more information on how to set this parameter efficiently.

### `--bwaalnk`

Configures the `bwa aln -k` parameter for the seeding phase in the mapping algorithm. Default is set to `2`. 

### `--bwaalnl`

Configures the length of the seed used in `bwa aln -l`. Default is set to BWA default of `32`.

## CircularMapper

### `--circularmapper`

This turns on the CircularMapper application, that enhances the mapping procedure with the BWA algorithm on circular references utilizing a extend-remap procedure (see Peltzer et al 2016, Genome Biology for details). 

### `--circularextension`

The number of bases to extend the reference genome with. By default this is set to `500` if not specified otherwise.

### `--circulartarget`

The chromosome in your FastA reference that you'd like to be treated as circular. By default this is set to `MT` but can be configured to match any other chromosome. 

### `--circularfilter`

If you want to filter out reads that don't map to a circular chromosome, turn this on. By default this option is turned off.

## BWA Mem

### `--bwamem`

Turn this on to utilize BWA Mem instead of `bwa aln` for alignment. Can be quite useful for modern DNA, but is rarely used in projects for ancient DNA.

## Mapped reads Stripping

These parameters are used for removing mapped reads from the original input FASTQ files, usually in the context of uploading the original FASTQ files to a public read archive (NCBI SRA/EBI ENA). 

These flags will produce FASTQ files almost identical to your input files, except that reads with the same read ID as one found in the mapped bam file, are either removed or 'masked' (every base replaced with Ns). 

This functionality allows you to provide other researchers who wish to re-use your data to apply their own adapter removal/read merging procedures, while maintaining anonyminity for sample donors - for example with microbiome research.

### `--strip_input_fastq`

Create pre-Adapter Removal FASTQ files without reads that mapped to reference (e.g. for public upload of privacy sensitive non-host data)

### `--strip_mode`

Read removal mode. Strip mapped reads completely (strip) or just replace mapped reads sequence by N (replace)

## Read Filtering and Conversion Parameters

Users can configure to keep/discard/extract certain groups of reads efficiently in the nf-core/eager pipeline. 

### `--bam_discard_unmapped`

Defines whether unmapped reads should be discarded and stored in FastQ and/or BAM format separately. The behaviour depends on the choice of the `--bam_unmapped_type`.

### `--bam_unmapped_type`

Defines how to proceed with unmapped reads: "discard" removes all unmapped reads, "bam" keeps unmapped reads as BAM file, "fastq" keeps unmapped reads as FastQ file, "both" keeps both BAM and FastQ files. Only effective when option `--bam_discard_unmapped` is turned on.

### `--bam_mapping_quality_threshold`

Specify a mapping quality threshold for mapped reads to be kept for downstream analysis. By default keeps all reads and is therefore set to `0` (basically doesn't filter anything).


## Read DeDuplication Parameters

### `--dedupper`
Sets the duplicate read removal tool. By default uses `dedup` an ancient DNA specific read deduplication tool. Users can also specify `markdup` and use Picard MarkDuplicates instead, which is advised when working with paired end data that is *not* merged beforehand. In all other cases, it is advised to use `dedup`. 

### `--dedup_all_merged`
Sets DeDup to treat all reads as merged reads. This is useful if reads are for example not prefixed with `M_` in all cases.

## Library Complexity Estimation Parameters

### `--preseq_step_size`

Can be used to configure the step size of Preseqs `c_curve` method. Can be useful when only few and thus shallow sequencing results are used for extrapolation.

## DNA Damage Assessment Parameters

### `--damageprofiler_length`

Specifies the length filter for DamageProfiler. By default set to `100`. 

### `--damageprofiler_threshold`

Specifies the length of the read start and end to be considered for profile generation in DamageProfiler. By default set to `15` bases. 

### `--run_pmdtools`

Specifies to run PMDTools for damage based read filtering and assessment of DNA damage in sequencing libraries. By default turned off. 

### `--udg` false

Defines whether Uracil-DNA glycosylase (UDG) treatment was used to repair DNA damage on the sequencing libraries. If set, the parameter is used by downstream tools such as PMDTools to estimate damage only on CpG sites that are left after such a treatment. 

### `--pmd_udg_type` \`half`

If you have UDGhalf treated data (Rohland et al 2016), specify `half` as option to this parameter to use a different model for DNA damage assessment in PMDTools. Specify the parameter with `full` and the DNA damage assesment will use CpG context only. If you don't specify the parameter at all, the library will be treated as non UDG treated.

### `--pmdtools_range`

Specifies the range in which to consider DNA damage from the ends of reads. By default set to `10`. 

### `--pmdtools_threshold `

Specifies the PMDScore threshold to use in the pipeline when filtering BAM files for DNA damage. Only reads which surpass this damage score are considered for downstream DNA analysis. By default set to `3` if not set specifically by the user. 

### `--pmdtools_reference_mask` ''

Can be used to set a reference genome mask for PMDTools. 

### `--pmdtools_max_reads`

The maximum number of reads used for damage assessment in PMDtools. Can be used to significantly reduce the amount of time required for damage assessment in PMDTools. Note that a too low value can also obtain incorrect results. 

## BAM Trimming Parameters

For some library preparation protocols, users might want to clip off damaged bases before applying genotyping methods. This can be done in nf-core/eager automatically by turning on the `--trim_bam` parameter.

### `--trim_bam`

Turns on the BAM trimming method. Trims off `[n]` bases from reads in the deduplicated BAM file. Damage assessment in PMDTools or DamageProfiler remains untouched, as data is routed through this independently.

### `--bamutils_clip_left` / `--bamutils_clip_right`

Default set to `1` and clipps off one base of the left or right side of reads. Note that reverse reads will automatically be clipped off at the reverse side with this (automatically reverses left and right for the reverse read).

### `--bamutils_softclip`

By default, nf-core/eager uses hard clipping and sets clipped bases to `N` with quality `!` in the BAM output. Turn this on to use soft-clipping instead, masking reads at the read ends respectively using the CIGAR string.

## Library-Type Parameters

These parameters are required in some cases, e.g. when performing in-solution SNP capture protocols (390K,1240K, ...) for population genetics for example. Make sure to specify the required parameters in such cases. 

### `--snpcapture` false

This is by default set to `false`, but can be turned on to calculate on target metrics automatically for you. Note, that this requires setting `--bedfile` with the target SNPs simultaneously. 

### `--bedfile` 

Can be used to set a path to a BED file (3/6 column format) to calculate capture target efficiency on the fly. Will not be used without `--bedfile` set as parameter.

## Genotyping Parameters

There are options for different genotypers to be used. We suggest you the documentation of each tool to find the ones that suit your needs.

### `--run_genotyping`

Turns on genotyping to run on all post-dedup and downstream BAMs. For example if `--run_pmdtools` and `--trim_bam` are both supplied, the genotyper will be run on all three deduplication, pmd and trim bam BAM files.

## `--genotyping_tool`

Specifies which genotyper to use. Current options are GATK (v3.5) UnifiedGenotyper or GATK (v4.xx) HaplotypeCaller, specified by 'ug' and 'hc' respectively. Note that while UnifiedGenotyper is more suitable for low-coverage ancient DNA (HaplotypeCaller does _de novo_ assembly around each variant site), it is officially deperecated by the Broad Institute and is only accessible by an archived version not properly avaliable on `conda`. Therefore specifying 'ug' will download the GATK 3.5 `-jar` for you.

## `--gatk_out_mode`

If selected a GATK genotyper, what type of VCF to create, i.e. produce calls for every site or just confidence sites. Options: EMIT_VARIANTS_ONLY, EMIT_ALL_CONFIDENT_SITES, EMIT_ALL_SITES. Default: EMIT_VARIANTS_ONLY.

## `--gatk_call_conf`

If selected a GATK genotyper phred-scaled confidence threshold of a given SNP/INDEL call. Default: 30

## `--gatk_ploidy`                 

If selected a GATK genotyper, what is the ploidy of your reference organism. E.g. do you want to allow heterozygous calls from >= diploid orgaisms. . Default: 2        

## `--gatk_dbsnp`                  

(Optional)Specify VCF file for output VCF SNP annotation e.g. if you want annotate your VCF file with 'rs' SNP IDs. Check GATK documentation for more information. Gzip not accepted.

## `--gatk_ug_genotype_model`

If selected GATK UnifiedGenotyper, which likelihood model to follow, i.e. whether to call use SNPs or INDELS etc. Options: SNP, INDEL, BOTH, GENERALPLOIDYSNP, GENERALPLOIDYINDEL. Default: SNP. 

## `--gatk_hc_emitrefconf`

If selected GATK HaplotypeCaller, mode for emitting reference confidence calls. Options: NONE, BP_RESOLUTION, GVCF. Default: GVCF

## Automatic Resubmission
By default, if a pipeline step fails, EAGER2 will resubmit the job with twice the amount of CPU and memory. This will occur two times before failing.

## Clean up

Once completed a run has completed, you will have _lots_ of (some very large) intermediate files in your output directory, within the directory named `work`. 

Once you have verified your run completed correctly and everything in the module output directories are present as you expect and need, you can perform a clean up.

> **Important**: Once clean up is completed, you will _not_ be able to re-rerun the pipline from an earlier step and you'll have to re-run from scratch.

While in your output directory, firstly verify you're only deleting files stored in `work/` with the dry run command:

```bash
nextflow clean -n
```
If you're ready, you can then remove the files with

```bash
nextflow clean -f
```
This will make your system administrator very happy as you will _halve_ the harddrive footprint of the run, so be sure to do this!

### `--monochrome_logs`
Set to disable colourful command line output and live life in monochrome.

### `--multiqc_config`
Specify a path to a custom MultiQC configuration file.
