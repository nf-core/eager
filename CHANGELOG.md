# nf-core/eager: Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

## unpublished

### `Added`
* [#80](https://github.com/nf-core/eager/pull/80) - BWA Index file handling 
* [#77](https://github.com/nf-core/eager/pull/77) - Lots of documentation updates by [@jfy133](https://github.com/jfy133)

### `Changed`
* [#81](https://github.com/nf-core/eager/pull/81) - Renaming of certain BAM options

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