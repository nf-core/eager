# nf-core/eager: Changelog

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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

## [2.1.0] - 2020-03-05 - "Ravensburg"

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
