# nf-core/eager: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [dev]

### `Added`

* Added Support for automated tests using [GitHub Actions](https://github.com/features/actions)
* Added genotyping capability through GATK UnifiedGenotyper (v3.5), GATK HaplotypeCaller (v4.1) and FreeBayes
* Added MultiVCFAnalyzer module

### `Fixed`

* [#227](https://github.com/nf-core/eager/issues/227) - Large re-write of input/output process logic to allow maximum flexibility. Originally to address [#227](https://github.com/nf-core/eager/issues/227), but further expanded

### `Dependencies`

* adapterremoval=2.2.2->2.3.1
* adapterremovalfixprefix=0.0.4->0.0.5
* picard=2.20.2->2.20.7
* angsd=0.923->0.931
* gatk4=4.1.2.0->4.1.3.0
* conda-forge::r-rmarkdown=1.12->1.15
* pysam=0.15.2->0.15.3
* python=3.6.3->3.7.3

## [2.0.7] - 2019-06-10

### `Added`

* [#189](https://github.com/nf-core/eager/pull/189) - Outputing unmapped reads in a fastq files with the --strip_input_fastq flag
* [#186](https://github.com/nf-core/eager/pull/186) - Make FastQC skipping [possible](https://github.com/nf-core/eager/issues/182)
* Merged in [nf-core/tools](https://github.com/nf-core/tools) release V1.6 template changes  
* A lot more automated tests using Travis CI
* Don't ignore DamageProfiler errors anymore
* [#220](https://github.com/nf-core/eager/pull/220) - Added post-mapping filtering statistics module and corresponding MultiQC statistics [#217](https://github.com/nf-core/eager/issues/217)

### `Fixed`

* [#152](https://github.com/nf-core/eager/pull/152) - DamageProfiler errors [won't crash entire pipeline anymore](https://github.com/nf-core/eager/issues/171)
* [#176](https://github.com/nf-core/eager/pull/176) - Increase runtime for DamageProfiler on [large reference genomes](https://github.com/nf-core/eager/issues/173)
* [#172](https://github.com/nf-core/eager/pull/152) - DamageProfiler errors [won't crash entire pipeline anymore](https://github.com/nf-core/eager/issues/171)
* [#174](https://github.com/nf-core/eager/pull/190) - Publish DeDup files [properly](https://github.com/nf-core/eager/issues/183)
* [#196](https://github.com/nf-core/eager/pull/196) - Fix reference [issues](https://github.com/nf-core/eager/issues/150)
* [#196](https://github.com/nf-core/eager/pull/196) - Fix issues with PE data being mapped incompletely
* [#200](https://github.com/nf-core/eager/pull/200) - Fix minor issue with some [typos](https://github.com/nf-core/eager/pull/196)
* [#210](https://github.com/nf-core/eager/pull/210) - Fix PMDTools [encoding issue](https://github.com/pontussk/PMDtools/issues/6) from `samtools calmd` generated files by running through `sa]mtools view` first
* [#221](https://github.com/nf-core/eager/pull/221) - Fix BWA Index [not being reused by multiple samples](https://github.com/nf-core/eager/issues/219)

### `Dependencies`

* Added DeDup v0.12.5 (json support)
* Added mtnucratio v0.5 (json support)
* Updated Picard 2.18.27 -> 2.20.2
* Updated GATK 4.1.0.0 -> 4.1.2.0
* Updated damageprofiler 0.4.4 -> 0.4.5
* Updated r-rmarkdown 1.11 -> 1.12
* Updated fastp 0.19.7 -> 0.20.0
* Updated qualimap 2.2.2b -> 2.2.2c

## [2.0.6] - 2019-03-05

### `Added`

* [#152](https://github.com/nf-core/eager/pull/152) - Clarified `--complexity_filter` flag to be specifically for poly G trimming.
* [#155](https://github.com/nf-core/eager/pull/155) - Added [Dedup log to output folders](https://github.com/nf-core/eager/issues/154)
* [#159](https://github.com/nf-core/eager/pull/159) - Added Possibility to skip AdapterRemoval, skip merging, skip trimming fixing [#64](https://github.com/nf-core/eager/issues/64),[#137](https://github.com/nf-core/eager/issues/137) - thanks to @maxibor, @jfy133

### `Fixed`

* [#151](https://github.com/nf-core/eager/pull/151) - Fixed [post-deduplication step errors](https://github.com/nf-core/eager/issues/128)
* [#147](https://github.com/nf-core/eager/pull/147) - Fix Samtools Index for [large references](https://github.com/nf-core/eager/issues/146)
* [#145](https://github.com/nf-core/eager/pull/145) - Added Picard Memory Handling [fix](https://github.com/nf-core/eager/issues/144)

### `Dependencies`

* Picard Tools 2.18.23 -> 2.18.27
* GATK 4.0.12.0 -> 4.1.0.0
* FastP 0.19.6 -> 0.19.7

## [2.0.5] - 2019-01-28

### `Added`

* [#127](https://github.com/nf-core/eager/pull/127) - Added a second testcase for testing the pipeline properly
* [#129](https://github.com/nf-core/eager/pull/129) - Support BAM files as [input format](https://github.com/nf-core/eager/issues/41)
* [#131](https://github.com/nf-core/eager/pull/131) - Support different [reference genome file extensions](https://github.com/nf-core/eager/issues/130)

### `Fixed`

* [#128](https://github.com/nf-core/eager/issues/128) - Fixed reference genome handling errors

### `Dependencies`

* Picard Tools 2.18.21 -> 2.18.23
* R-Markdown 1.10 -> 1.11
* FastP 0.19.5 -> 0.19.6

## [2.0.4] - 2019-01-09

### `Added`

* [#111](https://github.com/nf-core/eager/pull/110) - Allow [Zipped FastA reference input](https://github.com/nf-core/eager/issues/91)
* [#113](https://github.com/nf-core/eager/pull/113) - All files are now staged via channels, which is considered best practice by Nextflow
* [#114](https://github.com/nf-core/eager/pull/113) - Add proper runtime defaults for multiple processes
* [#118](https://github.com/nf-core/eager/pull/118) - Add [centralized configs handling](https://github.com/nf-core/configs)
* [#115](https://github.com/nf-core/eager/pull/115) - Add DamageProfiler MultiQC support
* [#122](https://github.com/nf-core/eager/pull/122) - Add pulling from Dockerhub again

### `Fixed`
* [#110](https://github.com/nf-core/eager/pull/110) - Fix for [MultiQC Missing Second FastQC report](https://github.com/nf-core/eager/issues/107)
* [#112](https://github.com/nf-core/eager/pull/112) - Remove [redundant UDG options](https://github.com/nf-core/eager/issues/89)

## [2.0.3] - 2018-12-12

### `Added`

* [#80](https://github.com/nf-core/eager/pull/80) - BWA Index file handling
* [#77](https://github.com/nf-core/eager/pull/77) - Lots of documentation updates by [@jfy133](https://github.com/jfy133)
* [#81](https://github.com/nf-core/eager/pull/81) - Renaming of certain BAM options
* [#92](https://github.com/nf-core/eager/issues/92) - Complete restructure of BAM options

### `Fixed`

* [#84](https://github.com/nf-core/eager/pull/85) - Fix for [Samtools index issues](https://github.com/nf-core/eager/issues/84)
* [#96](https://github.com/nf-core/eager/issues/96) - Fix for [MarkDuplicates issues](https://github.com/nf-core/eager/issues/96) found by [@nilesh-tawari](https://github.com/nilesh-tawari)

### Other

* Added Slack button to repository readme

## [2.0.2] - 2018-11-03

### `Changed`

* [#70](https://github.com/nf-core/eager/issues/70) - Uninitialized `readPaths` warning removed

### `Added`

* [#73](https://github.com/nf-core/eager/pull/73) - Travis CI Testing of Conda Environment added

### `Fixed`

* [#72](https://github.com/nf-core/eager/issues/72) - iconv Issue with R in conda environment

## [2.0.1] - 2018-11-02

### `Fixed`

* [#69](https://github.com/nf-core/eager/issues/67) - FastQC issues with conda environments

## [2.0.0] "Kaufbeuren" - 2018-10-17

Initial release of nf-core/eager:

### `Added`

* FastQC read quality control
* (Optional) Read complexity filtering with FastP
* Read merging and clipping using AdapterRemoval v2
* Mapping using BWA / BWA Mem or CircularMapper
* Library Complexity Estimation with Preseq
* Conversion and Filtering of BAM files using Samtools
* Damage assessment via DamageProfiler, additional filtering using PMDTools
* Duplication removal via DeDup
* BAM Clipping with BamUtil for UDGhalf protocols
* QualiMap BAM quality control analysis

Furthermore, this already creates an interactive report using MultiQC, which will be upgraded in V2.1 "Ulm" to contain more aDNA specific metrics.
