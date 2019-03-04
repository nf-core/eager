# nf-core/eager: Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## [Unpublished / Dev Branch]

### `Added`

* [#152](https://github.com/nf-core/eager/pull/152) - Clarified `--complexity_filter` flag to be specifically for poly G trimming.
* [#155](https://github.com/nf-core/eager/pull/155) - Added [Dedup log to output folders](https://github.com/nf-core/eager/issues/154)

### `Fixed`

* [#151](https://github.com/nf-core/eager/pull/151) - Fixed [post-deduplication step errors](https://github.com/nf-core/eager/issues/128
* [#147](https://github.com/nf-core/eager/pull/147) - Fix Samtools Index for [large references](https://github.com/nf-core/eager/issues/146)
* [#145](https://github.com/nf-core/eager/pull/145) - Added Picard Memory Handling [fix](https://github.com/nf-core/eager/issues/144)

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
* [#118](https://github.com/nf-core/eager/pull/118) - Add centralized configs handling by https://github.com/nf-core/configs
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
