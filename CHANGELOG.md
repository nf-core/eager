# nf-core/eager: Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](http://semver.org/spec/v2.0.0.html).

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