# nf-core/eager: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [2.4.6] - 2022-11-14

### `Added`

- [#933](https://github.com/nf-core/eager/issues/933) Added support for customising --seq-length in mapDamage rescaling (♥ to @ashildv for requesting)

### `Fixed`

- Changed endors.py license from GPL to MIT (♥ to @aidaanva for fixing)
- Removed erroneous R2 in single-end example in input TSV of usage docs (♥ to @aidaanva for fixing)
- [#928](https://github.com/nf-core/eager/issues/928) Fixed read group incompatibility by re-adding picard AddOrReplaceReadGroups for MultiVCFAnalyzer (♥ to @aidaanva, @meganemichel for reporting)
- Fixed edge case of DamageProfiler occasionally requiring FASTA index (♥ to @asmaa-a-abdelwahab for reporting)
- [#834](https://github.com/nf-core/eager/issues/834) Increased significance values in general stats table for Qualimap mean/median coverages (♥ to @neija2611 for reporting)
- Fixed parameter documentation for `--snpcapture_bed` regarding on-target SNP stats to state these stats currently not displayed in MultiQC only in the Qualimap results (♥ to @meganemichel and @TCLamnidis for reporting)
- [#934](https://github.com/nf-core/eager/issues/934) Fixed broken parameter setting in mapDamage2 rescale length (♥ to @ashildv for reporting)

### `Dependencies`

- Updated MultiQC to official 1.13 version (rather than alpha)
- Added pinned MALT dependency to ensure working version in future versions of eager

### `Deprecated`

## [2.4.5] - 2022-08-02

### `Added`

### `Fixed`

- [#882](https://github.com/nf-core/eager/pull/882) Define DSL1 execution explicitly, as new versions Nextflow made DSL2 default (♥ to & fix from @Lehmann-Fabian)
- [#879](https://github.com/nf-core/eager/issues/879) Add missing threads parameter for pre-clipping FastQC for single end data that caused insufficient memory in some cases (♥ to @marcel-keller for reporting)
- [#880](https://github.com/nf-core/eager/issues/880) Fix failure of endorSpy to be cached or reexecuted on resume (♥ to @KathrinNaegele, @TCLamnidis, & @mahesh-panchal for reporting and debugging)
- [#885](https://github.com/nf-core/eager/issues/885) Specify task memory for all tools in get_software_versions to account for incompatibilty of java with some SGE clusters causing hanging of the process (♥ to @maxibor for reporting)
- [#887](https://github.com/nf-core/eager/issues/887) Clarify what is considered 'ultra-short' reads in the help text of clip_readlength, for when you may wish to turn of length filtering during AdapterRemoval (♥ to @TCLamnidis for reporting)
- [#889](https://github.com/nf-core/eager/issues/889) Remove/update parameters from benchmarking test profiles (♥ to @TCLamnidis for reporting)
- [#895](https://github.com/nf-core/eager/issues/895) Output documentation typo fix and added location of output docs in pipeline summary (♥ to @RodrigoBarquera for reporting)
- [#897](https://github.com/nf-core/eager/issues/897) Fix pipeline crash if no Kraken2 results generated (♥ to @alexandregilardet for reporting)
- [#899](https://github.com/nf-core/eager/issues/897) Fix pipeline crash for circulargenerator if reference file does not end in .fasta (♥ to @scarlhoff for reporting)
- Fixed some missing default values in the nextflow parameter schema JSON
- [#789](https://github.com/nf-core/eager/issues/789) Substantial speed and memory optimisation of the `extract_map_reads.py` script (♥ to @ivelsko for reporting, @maxibor for optimisation)
- Fix staging of input bams for genotyping_pileupcaller process. Downstream changes from changes introduced when fixing endorspy caching.
- Made slight correction on metro map diagram regarding input data to SexDeterrmine (only BAM trimming output files)

### `Dependencies`

- Updated MultiQC to latest stable alpha version on bioconda, correcting the previously nonsensical AdapterRemoval plots (♥ to @NiemannJ for fixing in MultiQC)

### `Deprecated`

## [2.4.4] - 2022-04-08

### `Added`

### `Fixed`

- Fixed some auxiliary files (adapater list, snpcapture/pileupcaller/sexdeterrmine BED files, and pileupCaller SNP file, PMD reference mask) in some cases only be used against one sample (❤ to @meganemichel for reporting, fix by @jfy133)

### `Dependencies`

### `Deprecated`

## [2.4.3] - 2022-03-24

### `Added`

### `Fixed`

- [#828](https://github.com/nf-core/eager/issues/828) Improved error message if required metagenomic screening parameters not set correctly
- [#836](https://github.com/nf-core/eager/issues/836) Remove deprecated parameters from test profiles
- [#838](https://github.com/nf-core/eager/issues/838) Fix --snpcapture_bed files not being picked up by Nextflow (❤ to @meganemichel for reporting)
- [#843](https://github.com/nf-core/eager/issues/843) Re-add direct piping of AdapterRemovalFixPrefix to pigz
- [#844](https://github.com/nf-core/eager/issues/844) Fixed reference masking prior to pmdtools
- [#845](https://github.com/nf-core/eager/issues/845) Updates parameter documention to specify `-s` preseq parameter also applies to lc_extrap
- [#851](https://github.com/nf-core/eager/issues/851) Fixes a file-name clash during additional_library_merge, post-BAM trimming of different UDG treated libraries of a sample (❤ to @alexandregilardet for reporting)
- Renamed a range of MultiQC general stats table headers to improve clarity, documentation has been updated accordingly
- [#857](https://github.com/nf-core/eager/issues/857) Corrected samtools fastq flag to _retain_ read-pair information when converting off-target BAM files to fastq in paired-end mapping (❤ to @alexhbnr for reporting)
- [#866](https://github.com/nf-core/eager/issues/866) Fixed a typo in the indexing step of BWA mem when not-collapsing (❤ to @alexhbnr for reporting)
- Corrected tutorials to reflect updated BAM trimming flags (❤ to @marcel-keller for reporting and correcting)

### `Dependencies`

- [#829](https://github.com/nf-core/eager/issues/829) Bumped sequencetools: 1.4.0.5 -> 1.5.2
- Bumped MultiQC: 1.11 -> 1.12 (for run-time optimisation and tool citation information)

### `Deprecated`

## [2.4.2] - 2022-01-24

### `Added`

### `Fixed`

- [#824](https://github.com/nf-core/eager/issues/824) Fixes large memory footprint of bedtools coverage calculation.
- [#822](https://github.com/nf-core/eager/issues/822) Fixed post-adapterremoval trimmed files not being lane-merged and included in downstream analyses
- Fixed a couple of software version reporting commands

### `Dependencies`

### `Deprecated`

## [2.4.1] - 2021-11-30

### `Added`

- [#805](https://github.com/nf-core/eager/issues/805) Changes to bam_trim options to allow flexible trimming by library strandedness (in addition to UDG treatment). (@TCLamnidis)
- [#808](https://github.com/nf-core/eager/issues/808) Retain read group information across bam merges. Sample set to sample name (rather than library name) in bwa output 'RG' readgroup tag. (@TCLamnidis)
- Map and base quality filters prior to genotyping with pileupcaller can now be specified. (@TCLamnidis)
- [#774](https://github.com/nf-core/eager/issues/774) Added support for multi-threaded Bowtie2 build reference genome indexing (@jfy133)
- [#804](https://github.com/nf-core/eager/issues/804) Improved output documentation description to add how 'cluster factor' is calculated (thanks to @meganemichel)

### `Fixed`

- [#803](https://github.com/nf-core/eager/issues/803) Fixed mistake in metro-map diagram (`samtools index` is now correctly `samtools faidx`) (@jfy133)

### `Dependencies`

### `Deprecated`

## [2.4.0] - Wangen - 2021-09-14

### `Added`

- [#317](https://github.com/nf-core/eager/issues/317) Added bcftools stats for general genotyping statistics of VCF files
- [#651](https://github.com/nf-core/eager/issues/651) - Adds removal of adapters specified in an AdapterRemoval adapter list file
- [#642](https://github.com/nf-core/eager/issues/642) and [#431](https://github.com/nf-core/eager/issues/431) adds post-adapter removal barcode/fastq trimming
- [#769](https://github.com/nf-core/eager/issues/769) - Adds lc_extrap mode to preseq (suggested by @roberta-davidson)

### `Fixed`

- Fixed some missing or incorrectly reported software versions
- [#771](https://github.com/nf-core/eager/issues/771) Remove legacy code
- Improved output documentation for MultiQC general stats table (thanks to @KathrinNaegele and @esalmela)
- Improved output documentation for BowTie2 (thanks to @isinaltinkaya)
- [#612](https://github.com/nf-core/eager/issues/612) Updated BAM trimming defaults to 0 to ensure no unwanted trimming when mixing half-UDG with no-UDG (thanks to @scarlhoff)
- [#722](https://github.com/nf-core/eager/issues/722) Updated BWA mapping mapping parameters to latest recommendations - primarily alnn back to 0.01 and alno to 2 as per Oliva et al. 2021 (10.1093/bib/bbab076)
- Updated workflow diagrams to reflect latest functionality
- [#787](https://github.com/nf-core/eager/issues/787) Adds memory specification flags for the GATK UnifiedGenotyper and HaplotyperCaller steps (thanks to @nylander)
- Fixed issue where MultiVCFAnalyzer would not pick up newly generated VCF files, when specifying additional VCF files.
- [#790](https://github.com/nf-core/eager/issues/790) Fixed kraken2 report file-name collision when sample names have `.` in them
- [#792](https://github.com/nf-core/eager/issues/792) Fixed java error messages for AdapterRemovalFixPrefix being hidden in output
- [#794](https://github.com/nf-core/eager/issues/794) Aligned default test profile with nf-core standards (`test_tsv` is now `test`)

### `Dependencies`

- Bumped python: 3.7.3 -> 3.9.4
- Bumped markdown: 3.2.2 -> 3.3.4
- Bumped pymdown-extensions: 7.1 -> 8.2
- Bumped pyments: 2.6.1 -> 2.9.0
- Bumped adapterremoval: 2.3.1 -> 2.3.2
- Bumped picard: 2.22.9 -> 2.26.0
- Bumped samtools 1.9 -> 1.12
- Bumped angsd: 0.933 -> 0.935
- Bumped gatk4: 4.1.7.0 -> 4.2.0.0
- Bumped multiqc: 1.10.1 -> 1.11
- Bumped bedtools 2.29.2 -> 2.30.0
- Bumped libiconv: 1.15 -> 1.16
- Bumped preseq: 2.0.3 -> 3.1.2
- Bumped bamutil: 1.0.14 -> 1.0.15
- Bumped pysam: 0.15.4 -> 0.16.0
- Bumped kraken2: 2.1.1 -> 2.1.2
- Bumped pandas: 1.0.4 -> 1.2.4
- Bumped freebayes: 1.3.2 -> 1.3.5
- Bumped biopython: 1.76 -> 1.79
- Bumped xopen: 0.9.0 -> 1.1.0
- Bumped bowtie2: 2.4.2 -> 2.4.4
- Bumped mapdamage2: 2.2.0 -> 2.2.1
- Bumped bbmap: 38.87 -> 38.92
- Added bcftools: 1.12

### `Deprecated`

## [2.3.5] - 2021-06-03

### `Added`

- [#722](https://github.com/nf-core/eager/issues/722) - Adds bwa `-o` flag for more flexibility in bwa parameters
- [#736](https://github.com/nf-core/eager/issues/736) - Add printing of multiqc run report location on successful completion
- New logo that is more visible when a user is using darkmode on GitHub or nf-core website!

### `Fixed`

- [#723](https://github.com/nf-core/eager/issues/723) - Fixes empty fields in TSV resulting in uninformative error
- Updated template to nf-core/tools 1.14
- [#688](https://github.com/nf-core/eager/issues/688) - Clarified the pipeline is not just for humans and microbes, but also plants and animals, and also for modern DNA
- [#751](https://github.com/nf-core/eager/pull/751) - Added missing label to mtnucratio
- General code cleanup and standardisation of parameters with no default setting
- [#750](https://github.com/nf-core/eager/issues/750) - Fixed piped commands requesting the same number of CPUs at each command step
- [#757](https://github.com/nf-core/eager/issues/757) - Removed confusing 'Data Type' variable from MultiQC workflow summary (not consistent with TSV input)
- [#759](https://github.com/nf-core/eager/pull/759) - Fixed malformed software scraping regex that resulted in N/A in MultiQC report
- [#761](https://github.com/nf-core/eager/pull/759) - Fixed issues related to instability of samtools filtering related CI tests

### `Dependencies`

### `Deprecated`

## [2.3.4] - 2021-05-05

### `Added`

- [#729](https://github.com/nf-core/eager/issues/729) - Added Bowtie2 flag `--maxins` for PE mapping modern DNA mapping contexts

### `Fixed`

- Corrected explanation of the "--min_adap_overlap" parameter for AdapterRemoval in the docs
- [#725](https://github.com/nf-core/eager/pull/725) - `bwa_index` doc update
- Re-adds gzip piping to AdapterRemovalFixPrefix to speed up process after reports of being very slow
- Updated DamageProfiler citation from bioRxiv to publication

### `Dependencies`

- Removed pinning of `tbb` (upstream bug in bioconda fixed)
- Bumped `pigz` to 2.6 to fix rare stall bug when compressing data after AdapterRemoval
- Bumped Bowtie2 to 2.4.2 to fix issues with `tbb` version

### `Deprecated`

## [2.3.3] - 2021-04-08

### `Added`

- [#349](https://github.com/nf-core/eager/issues/349) - Added option enabling platypus formatted output of pmdtools misincorporation frequencies.

### `Fixed`

- [#719](https://github.com/nf-core/eager/pull/719) - Fix filename for bam output of `mapdamage_rescaling`
- [#707](https://github.com/nf-core/eager/pull/707) - Fix typo in UnifiedGenotyper IndelRealigner command
- Fixed some Java tools not following process memory specifications
- Updated template to nf-core/tools 1.13.2
- [#711](https://github.com/nf-core/eager/pull/711) - Fix conditional execution preventing multivcfanalyze to run
- [#714](https://github.com/nf-core/eager/issues/714) - Fixes bug in nuc contamination by upgrading to latest MultiQC v1.10.1 bugfix release

### `Dependencies`

### `Deprecated`

## [2.3.2] - 2021-03-16

### `Added`

- [#687](https://github.com/nf-core/eager/pull/687) - Adds Kraken2 unique kmer counting report
- [#676](https://github.com/nf-core/eager/issues/676) - Refactor help message / summary message formatting to automatic versions using nf-core library
- [#682](https://github.com/nf-core/eager/issues/682) - Add AdapterRemoval `--qualitymax` flag to allow FASTQ Phred score range max more than 41

### `Fixed`

- [#666](https://github.com/nf-core/eager/issues/666) - Fixed input file staging for `print_nuclear_contamination`
- [#631](https://github.com/nf-core/eager/issues/631) - Update minimum Nextflow version to 20.07.1, due to unfortunate bug in Nextflow 20.04.1 causing eager to crash if patch pulled
- Made MultiQC crash behaviour stricter when dealing with large datasets, as reported by @ashildv
- [#652](https://github.com/nf-core/eager/issues/652) - Added note to documentation that when using `--skip_collapse` this will use _paired-end_ alignment mode with mappers when using PE data
- [#626](https://github.com/nf-core/eager/issues/626) - Add additional checks to ensure pipeline will give useful error if cells of a TSV column are empty
- Added note to documentation that when using `--skip_collapse` this will use _paired-end_ alignment mode with mappers when using PE data
- [#673](https://github.com/nf-core/eager/pull/673) - Fix Kraken database loading when loading from directory instead of compressed file
- [#688](https://github.com/nf-core/eager/issues/668) - Allow pipeline to complete, even if Qualimap crashes due to an empty or corrupt BAM file for one sample/library
- [#683](https://github.com/nf-core/eager/pull/683) - Sets `--igenomes_ignore` to true by default, as rarely used by users currently and makes resolving configs less complex
- Added exit code `140` to re-tryable exit code list to account for certain scheduler wall-time limit fails
- [#672](https://github.com/nf-core/eager/issues/672) - Removed java parameter from picard tools which could cause memory issues
- [#679](https://github.com/nf-core/eager/issues/679) - Refactor within-process bash conditions to groovy/nextflow, due to incompatibility with some servers environments
- [#690](https://github.com/nf-core/eager/pull/690) - Fixed ANGSD output mode for beagle by setting `-doMajorMinor 1` as default in that case
- [#693](https://github.com/nf-core/eager/issues/693) - Fixed broken TSV input validation for the Colour Chemistry column
- [#695](https://github.com/nf-core/eager/issues/695) - Fixed incorrect `-profile` order in tutorials (originally written reversed due to [nextflow bug](https://github.com/nextflow-io/nextflow/issues/1792))
- [#653](https://github.com/nf-core/eager/issues/653) - Fixed file collision errors with sexdeterrmine for two same-named libraries with different strandedness

### `Dependencies`

- Bumped MultiQC to 1.10 for improved functionality
- Bumped HOPS to 0.35 for MultiQC 1.10 compatibility

### `Deprecated`

## [2.3.1] - 2021-01-14

### `Added`

### `Fixed`

- [#654](https://github.com/nf-core/eager/issues/654) - Fixed some values in JSON schema (used in launch GUI) not passing validation checks during run
- [#655](https://github.com/nf-core/eager/issues/655) - Updated read groups for all mappers to allow proper GATK validation
- Fixed issue with Docker container not being pullable by Nextflow due to version-number inconsistencies

### `Dependencies`

### `Deprecated`

## [2.3.0] - Aalen - 2021-01-11

### `Added`

- [#640](https://github.com/nf-core/eager/issues/640) - Added a pre-metagenomic screening filtering of low-sequence complexity reads with `bbduk`
- [#583](https://github.com/nf-core/eager/issues/583) - Added `mapDamage2` rescaling of BAM files to remove damage
- Updated usage (merging files) and workflow images reflecting new functionality.

### `Fixed`

- Removed leftover old DockerHub push CI commands.
- [#627](https://github.com/nf-core/eager/issues/627) - Added de Barros Damgaard citation to README
- [#630](https://github.com/nf-core/eager/pull/630) - Better handling of Qualimap memory requirements and error strategy.
- Fixed some incomplete schema options to ensure users supply valid input values
- [#638](https://github.com/nf-core/eager/issues/638#issuecomment-748877567) Fixed inverted circularfilter filtering (previously filtering would happen by default, not when requested by user as originally recorded in documentation)
- [DeDup:](https://github.com/apeltzer/DeDup/commit/07d47868f10a6830da8c9161caa3755d9da155bf) Fixed Null Pointer Bug in DeDup by updating to 0.12.8 version
- [#650](https://github.com/nf-core/eager/pull/650) - Increased memory given to FastQC for larger files by making it multithreaded

### `Dependencies`

- Update: DeDup v0.12.7 to v0.12.8

### `Deprecated`

## [2.2.2] - 2020-12-09

### `Added`

- Added large scale 'stress-test' profile for AWS (using de Barros Damgaard et al. 2018's 137 ancient human genomes).
  - This will now be run automatically for every release. All processed data will be available on the nf-core website: <https://nf-co.re/eager/results>
    - You can run this yourself using `-profile test_full`

### `Fixed`

- Fixed AWS full test profile.
- [#587](https://github.com/nf-core/eager/issues/587) - Re-implemented AdapterRemovalFixPrefix for DeDup compatibility of including singletons
- [#602](https://github.com/nf-core/eager/issues/602) - Added the newly available GATK 3.5 conda package.
- [#610](https://github.com/nf-core/eager/issues/610) - Create bwa_index channel when specifying circularmapper as mapper
- Updated template to nf-core/tools 1.12.1
- General documentation improvements

### `Deprecated`

- Flag `--gatk_ug_jar` has now been removed as GATK 3.5 is now avaliable within the nf-core/eager software environment.

## [2.2.1] - 2020-10-20

### `Fixed`

- [#591](https://github.com/nf-core/eager/issues/591) - Fixed offset underlines in lane merging diagram in docs
- [#592](https://github.com/nf-core/eager/issues/592) - Fixed issue where supplying Bowtie2 index reported missing bwamem_index error
- [#590](https://github.com/nf-core/eager/issues/592) - Removed redundant dockstore.yml from root
- [#596](https://github.com/nf-core/eager/issues/596) - Add workaround for issue regarding gzipped FASTAs and pre-built indices
- [#589](https://github.com/nf-core/eager/issues/582) - Updated template to nf-core/tools 1.11
- [#582](https://github.com/nf-core/eager/issues/582) - Clarify memory limit issue on FAQ

## [2.2.0] - Ulm - 2020-10-20

### `Added`

- **Major** Automated cloud tests with large-scale data on [AWS](https://aws.amazon.com/)
- **Major** Re-wrote input logic to accept a TSV 'map' file in addition to direct paths to FASTQ files
- **Major** Added JSON Schema, enabling web GUI for configuration of pipeline available [here](https://nf-co.re/launch?pipeline=eager&release=2.2.0)
- **Major** Lane and library merging implemented
  - When using TSV input, one library with the multiple _lanes_ will be merged together, before mapping
  - Strip FASTQ will also produce a lane merged 'raw' but 'stripped' FASTQ file
  - When using TSV input, one sample with multiple (same treatment) libraries will be merged together
  - Important: direct FASTQ paths will not have this functionality. TSV is required.
- [#40](https://github.com/nf-core/eager/issues/40) - Added the pileupCaller genotyper from [sequenceTools](https://github.com/stschiff/sequenceTools)
- Added validation check and clearer error message when `--fasta_index` is provided and filepath does not end in `.fai`.
- Improved error messages
- Added ability for automated emails using `mailutils` to also send MultiQC reports
- General documentation additions, cleaning, and updated figures with CC-BY license
- Added large 'full size' dataset test-profiles for ancient fish and human contexts human
- [#257](https://github.com/nf-core/eager/issues/257) - Added the bowtie2 aligner as option for mapping, following Poullet and Orlando 2020 doi: [10.3389/fevo.2020.00105](https://doi.org/10.3389/fevo.2020.00105)
- [#451](https://github.com/nf-core/eager/issues/451) - Adds ANGSD genotype likelihood calculations as an alternative to typical 'genotypers'
- [#566](https://github.com/nf-core/eager/issues/466) - Add tutorials on how to set up nf-core/eager for different contexts
- Nuclear contamination results are now shown in the MultiQC report
- Tutorial on how to use profiles for reproducible science (i.e. parameter sharing between different groups)
- [#522](https://github.com/nf-core/eager/issues/522) - Added post-mapping length filter to assist in more realistic endogenous DNA calculations
- [#512](https://github.com/nf-core/eager/issues/512) - Added flexible trimming of BAMs by library type. 'half' and 'none' UDG libraries can now be trimmed differentially within a single eager run.
- Added a `.dockstore.yml` config file for automatic workflow registration with [dockstore.org](https://dockstore.org/)
- Updated template to nf-core/tools 1.10.2
- [#544](https://github.com/nf-core/eager/pull/544) - Add script to perform bam filtering on fragment length
- [#456](https://github.com/nf-core/eager/pull/546) - Bumps the base (default) runtime of all processes to 4 hours, and set shorter time limits for test profiles (1 hour)
- [#552](https://github.com/nf-core/eager/issues/552) - Adds optional creation of MALT SAM files alongside RMA6 files
- Added eigenstrat snp coverage statistics to MultiQC report. Process results are published in `genotyping/*_eigenstrat_coverage.txt`.

### `Fixed`

- [#368](https://github.com/nf-core/eager/issues/368) - Fixed the profile `test` to contain a parameter for `--paired_end`
- Mini bugfix for typo in line 1260+1261
- [#374](https://github.com/nf-core/eager/issues/374) - Fixed output documentation rendering not containing images
- [#379](https://github.com/nf-core/eager/issues/378) - Fixed insufficient memory requirements for FASTQC edge case
- [#390](https://github.com/nf-core/eager/issues/390) - Renamed clipped/merged output directory to be more descriptive
- [#398](https://github.com/nf-core/eager/issues/498) - Stopped incompatible FASTA indexes being accepted
- [#400](https://github.com/nf-core/eager/issues/400) - Set correct recommended bwa mapping parameters from [Schubert et al. 2012](https://doi.org/10.1186/1471-2164-13-178)
- [#410](https://github.com/nf-core/eager/issues/410) - Fixed nf-core/configs not being loaded properly
- [#473](https://github.com/nf-core/eager/issues/473) - Fixed bug in sexdet_process on AWS
- [#444](https://github.com/nf-core/eager/issues/444) - Provide option for preserving realigned bam + index
- Fixed deduplication output logic. Will now pass along only the post-rmdup bams if duplicate removal is not skipped, instead of both the post-rmdup and pre-rmdup bams
- [#497](https://github.com/nf-core/eager/issues/497) - Simplifies number of parameters required to run bam filtering
- [#501](https://github.com/nf-core/eager/issues/501) - Adds additional validation checks for MALT/MaltExtract database input files
- [#508](https://github.com/nf-core/eager/issues/508) - Made Markduplicates default dedupper due to narrower context specificity of dedup
- [#516](https://github.com/nf-core/eager/issues/516) - Made bedtools not report out of memory exit code when warning of inconsistent FASTA/Bed entry names
- [#504](https://github.com/nf-core/eager/issues/504) - Removed uninformative sexdeterrmine-snps plot from MultiQC report.
- Nuclear contamination is now reported with the correct library names.
- [#531](https://github.com/nf-core/eager/pull/531) - Renamed 'FASTQ stripping' to 'host removal'
- Merged all tutorials and FAQs into `usage.md` for display on [nf-co.re](https://www.nf-co.re)
- Corrected header of nuclear contamination table (`nuclear_contamination.txt`).
- Fixed a bug with `nSNPs` definition in `print_x_contamination.py`. Number of SNPs now correctly reported
- `print_x_contamination.py` now correctly converts all NA values to "N/A"
- Increased amount of memory MultiQC by default uses, to account for very large nf-core/eager runs (e.g. >1000 samples)

### `Dependencies`

- Added sequenceTools (1.4.0.6) that adds the ability to do genotyping with the 'pileupCaller'
- Latest version of DeDup (0.12.6) which now reports mapped reads after deduplication
- [#560](https://github.com/nf-core/eager/issues/560) Latest version of Dedup (0.12.7), which now correctly reports deduplication statistics based on calculations of mapped reads only (prior denominator was total reads of BAM file)
- Latest version of ANGSD (0.933) which doesn't seg fault when running contamination on BAMs with insufficient reads
- Latest version of MultiQC (1.9) with support for lots of extra tools in the pipeline (MALT, SexDetERRmine, DamageProfiler, MultiVCFAnalyzer)
- Latest versions of Pygments (7.1), Pymdown-Extensions (2.6.1) and Markdown (3.2.2) for documentation output
- Latest version of Picard (2.22.9)
- Latest version of GATK4 (4.1.7.0)
- Latest version of sequenceTools (1.4.0.6)
- Latest version of fastP (0.20.1)
- Latest version of Kraken2 (2.0.9beta)
- Latest version of FreeBayes (1.3.2)
- Latest version of xopen (0.9.0)
- Added Bowtie 2 (2.4.1)
- Latest version of Sex.DetERRmine (1.1.2)
- Latest version of endorS.py (0.4)

## [2.1.0] - Ravensburg - 2020-03-05

### `Added`

- Added Support for automated tests using [GitHub Actions](https://github.com/features/actions), replacing travis
- [#40](https://github.com/nf-core/eager/issues/40), [#231](https://github.com/nf-core/eager/issues/231) - Added genotyping capability through GATK UnifiedGenotyper (v3.5), GATK HaplotypeCaller (v4.1) and FreeBayes
- Added MultiVCFAnalyzer module
- [#240](https://github.com/nf-core/eager/issues/240) - Added human sex determination module
- [#226](https://github.com/nf-core/eager/issues/226) - Added `--preserve5p` function for AdapterRemoval
- [#212](https://github.com/nf-core/eager/issues/212) - Added ability to use only merged reads downstream from AdapterRemoval
- [#265](https://github.com/nf-core/eager/issues/265) - Adjusted full markdown linting in Travis CI
- [#247](https://github.com/nf-core/eager/issues/247) - Added nuclear contamination with angsd
- [#258](https://github.com/nf-core/eager/issues/258) - Added ability to report bedtools stats to features (e.g. depth/breadth of annotated genes)
- [#249](https://github.com/nf-core/eager/issues/249) - Added metagenomic classification of unmapped reads with MALT and aDNA authentication with MaltExtract
- [#302](https://github.com/nf-core/eager/issues/302) - Added mitochondrial to nuclear ratio calculation
- [#302](https://github.com/nf-core/eager/issues/302) - Added VCF2Genome for consensus sequence generation
- Fancy new logo from [ZandraFagernas](https://github.com/ZandraFagernas)
- [#286](https://github.com/nf-core/eager/issues/286) - Adds pipeline-specific profiles (loaded from nf-core configs)
- [#310](https://github.com/nf-core/eager/issues/310) - Generalises base.config
- [#326](https://github.com/nf-core/eager/pull/326) - Add Biopython and [xopen](https://github.com/marcelm/xopen/) dependencies
- [#336](https://github.com/nf-core/eager/issues/336) - Change default Y-axis maximum value of DamageProfiler to 30% to match popular (but slower) mapDamage, and allow user to set their own value.
- [#352](https://github.com/nf-core/eager/pull/352) - Add social preview image
- [#355](https://github.com/nf-core/eager/pull/355) - Add Kraken2 metagenomics classifier
- [#90](https://github.com/nf-core/eager/issues/90) - Added endogenous DNA calculator (original repository: [https://github.com/aidaanva/endorS.py/](https://github.com/aidaanva/endorS.py/))

### `Fixed`

- [#227](https://github.com/nf-core/eager/issues/227) - Large re-write of input/output process logic to allow maximum flexibility. Originally to address [#227](https://github.com/nf-core/eager/issues/227), but further expanded
- Fixed Travis-Ci.org to Travis-Ci.com migration issues
- [#266](https://github.com/nf-core/eager/issues/266) - Added sanity checks for input filetypes (i.e. only BAM files can be supplied if `--bam`)
- [#237](https://github.com/nf-core/eager/issues/237) - Fixed and Updated script scrape_software_versions
- [#322](https://github.com/nf-core/eager/pull/322) - Move extract map reads fastq compression to pigz
- [#327](https://github.com/nf-core/eager/pull/327) - Speed up strip_input_fastq process and make it more robust
- [#342](https://github.com/nf-core/eager/pull/342) - Updated to match nf-core tools 1.8 linting guidelines
- [#339](https://github.com/nf-core/eager/issues/339) - Converted unnecessary zcat + gzip to just cat for a performance boost
- [#344](https://github.com/nf-core/eager/issues/344) - Fixed pipeline still trying to run when using old nextflow version

### `Dependencies`

- adapterremoval=2.2.2 upgraded to 2.3.1
- adapterremovalfixprefix=0.0.4 upgraded to 0.0.5
- damageprofiler=0.4.3 upgraded to 0.4.9
- angsd=0.923 upgraded to 0.931
- gatk4=4.1.2.0 upgraded to 4.1.4.1
- mtnucratio=0.5 upgraded to 0.6
- conda-forge::markdown=3.1.1 upgraded to 3.2.1
- bioconda::fastqc=0.11.8 upgraded to 0.11.9
- bioconda::picard=2.21.4 upgraded to 2.22.0
- bioconda::bedtools=2.29.0 upgraded to 2.29.2
- pysam=0.15.3 upgraded to 0.15.4
- conda-forge::pandas=1.0.0 upgraded to 1.0.1
- bioconda::freebayes=1.3.1 upgraded to 1.3.2
- conda-forge::biopython=1.75 upgraded to 1.76

## [2.0.7] - 2019-06-10

### `Added`

- [#189](https://github.com/nf-core/eager/pull/189) - Outputting unmapped reads in a fastq files with the --strip_input_fastq flag
- [#186](https://github.com/nf-core/eager/pull/186) - Make FastQC skipping [possible](https://github.com/nf-core/eager/issues/182)
- Merged in [nf-core/tools](https://github.com/nf-core/tools) release V1.6 template changes
- A lot more automated tests using Travis CI
- Don't ignore DamageProfiler errors any more
- [#220](https://github.com/nf-core/eager/pull/220) - Added post-mapping filtering statistics module and corresponding MultiQC statistics [#217](https://github.com/nf-core/eager/issues/217)

### `Fixed`

- [#152](https://github.com/nf-core/eager/pull/152) - DamageProfiler errors [won't crash entire pipeline any more](https://github.com/nf-core/eager/issues/171)
- [#176](https://github.com/nf-core/eager/pull/176) - Increase runtime for DamageProfiler on [large reference genomes](https://github.com/nf-core/eager/issues/173)
- [#172](https://github.com/nf-core/eager/pull/152) - DamageProfiler errors [won't crash entire pipeline any more](https://github.com/nf-core/eager/issues/171)
- [#174](https://github.com/nf-core/eager/pull/190) - Publish DeDup files [properly](https://github.com/nf-core/eager/issues/183)
- [#196](https://github.com/nf-core/eager/pull/196) - Fix reference [issues](https://github.com/nf-core/eager/issues/150)
- [#196](https://github.com/nf-core/eager/pull/196) - Fix issues with PE data being mapped incompletely
- [#200](https://github.com/nf-core/eager/pull/200) - Fix minor issue with some [typos](https://github.com/nf-core/eager/pull/196)
- [#210](https://github.com/nf-core/eager/pull/210) - Fix PMDTools [encoding issue](https://github.com/pontussk/PMDtools/issues/6) from `samtools calmd` generated files by running through `sa]mtools view` first
- [#221](https://github.com/nf-core/eager/pull/221) - Fix BWA Index [not being reused by multiple samples](https://github.com/nf-core/eager/issues/219)

### `Dependencies`

- Added DeDup v0.12.5 (json support)
- Added mtnucratio v0.5 (json support)
- Updated Picard 2.18.27 -> 2.20.2
- Updated GATK 4.1.0.0 -> 4.1.2.0
- Updated damageprofiler 0.4.4 -> 0.4.5
- Updated r-rmarkdown 1.11 -> 1.12
- Updated fastp 0.19.7 -> 0.20.0
- Updated qualimap 2.2.2b -> 2.2.2c

## [2.0.6] - 2019-03-05

### `Added`

- [#152](https://github.com/nf-core/eager/pull/152) - Clarified `--complexity_filter` flag to be specifically for poly G trimming.
- [#155](https://github.com/nf-core/eager/pull/155) - Added [Dedup log to output folders](https://github.com/nf-core/eager/issues/154)
- [#159](https://github.com/nf-core/eager/pull/159) - Added Possibility to skip AdapterRemoval, skip merging, skip trimming fixing [#64](https://github.com/nf-core/eager/issues/64),[#137](https://github.com/nf-core/eager/issues/137) - thanks to @maxibor, @jfy133

### `Fixed`

- [#151](https://github.com/nf-core/eager/pull/151) - Fixed [post-deduplication step errors](https://github.com/nf-core/eager/issues/128)
- [#147](https://github.com/nf-core/eager/pull/147) - Fix Samtools Index for [large references](https://github.com/nf-core/eager/issues/146)
- [#145](https://github.com/nf-core/eager/pull/145) - Added Picard Memory Handling [fix](https://github.com/nf-core/eager/issues/144)

### `Dependencies`

- Picard Tools 2.18.23 -> 2.18.27
- GATK 4.0.12.0 -> 4.1.0.0
- FastP 0.19.6 -> 0.19.7

## [2.0.5] - 2019-01-28

### `Added`

- [#127](https://github.com/nf-core/eager/pull/127) - Added a second test case for testing the pipeline properly
- [#129](https://github.com/nf-core/eager/pull/129) - Support BAM files as [input format](https://github.com/nf-core/eager/issues/41)
- [#131](https://github.com/nf-core/eager/pull/131) - Support different [reference genome file extensions](https://github.com/nf-core/eager/issues/130)

### `Fixed`

- [#128](https://github.com/nf-core/eager/issues/128) - Fixed reference genome handling errors

### `Dependencies`

- Picard Tools 2.18.21 -> 2.18.23
- R-Markdown 1.10 -> 1.11
- FastP 0.19.5 -> 0.19.6

## [2.0.4] - 2019-01-09

### `Added`

- [#111](https://github.com/nf-core/eager/pull/110) - Allow [Zipped FastA reference input](https://github.com/nf-core/eager/issues/91)
- [#113](https://github.com/nf-core/eager/pull/113) - All files are now staged via channels, which is considered best practice by Nextflow
- [#114](https://github.com/nf-core/eager/pull/113) - Add proper runtime defaults for multiple processes
- [#118](https://github.com/nf-core/eager/pull/118) - Add [centralized configs handling](https://github.com/nf-core/configs)
- [#115](https://github.com/nf-core/eager/pull/115) - Add DamageProfiler MultiQC support
- [#122](https://github.com/nf-core/eager/pull/122) - Add pulling from Dockerhub again

### `Fixed`

- [#110](https://github.com/nf-core/eager/pull/110) - Fix for [MultiQC Missing Second FastQC report](https://github.com/nf-core/eager/issues/107)
- [#112](https://github.com/nf-core/eager/pull/112) - Remove [redundant UDG options](https://github.com/nf-core/eager/issues/89)

## [2.0.3] - 2018-12-12

### `Added`

- [#80](https://github.com/nf-core/eager/pull/80) - BWA Index file handling
- [#77](https://github.com/nf-core/eager/pull/77) - Lots of documentation updates by [@jfy133](https://github.com/jfy133)
- [#81](https://github.com/nf-core/eager/pull/81) - Renaming of certain BAM options
- [#92](https://github.com/nf-core/eager/issues/92) - Complete restructure of BAM options

### `Fixed`

- [#84](https://github.com/nf-core/eager/pull/85) - Fix for [Samtools index issues](https://github.com/nf-core/eager/issues/84)
- [#96](https://github.com/nf-core/eager/issues/96) - Fix for [MarkDuplicates issues](https://github.com/nf-core/eager/issues/96) found by [@nilesh-tawari](https://github.com/nilesh-tawari)

### Other

- Added Slack button to repository readme

## [2.0.2] - 2018-11-03

### `Changed`

- [#70](https://github.com/nf-core/eager/issues/70) - Uninitialized `readPaths` warning removed

### `Added`

- [#73](https://github.com/nf-core/eager/pull/73) - Travis CI Testing of Conda Environment added

### `Fixed`

- [#72](https://github.com/nf-core/eager/issues/72) - iconv Issue with R in conda environment

## [2.0.1] - 2018-11-02

### `Fixed`

- [#69](https://github.com/nf-core/eager/issues/67) - FastQC issues with conda environments

## [2.0.0] "Kaufbeuren" - 2018-10-17

Initial release of nf-core/eager:

### `Added`

- FastQC read quality control
- (Optional) Read complexity filtering with FastP
- Read merging and clipping using AdapterRemoval v2
- Mapping using BWA / BWA Mem or CircularMapper
- Library Complexity Estimation with Preseq
- Conversion and Filtering of BAM files using Samtools
- Damage assessment via DamageProfiler, additional filtering using PMDTools
- Duplication removal via DeDup
- BAM Clipping with BamUtil for UDGhalf protocols
- QualiMap BAM quality control analysis

Furthermore, this already creates an interactive report using MultiQC, which will be upgraded in V2.1 "Ulm" to contain more aDNA specific metrics.
