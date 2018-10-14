# nf-core/eager: Changelog

## 2.0 "Kaufbeuren" - 2018-10-14

Initial release of nf-core/eager featuring:

* FastQC read quality control
* (Optional) Read complexity filtering with FastP
* Read merging and clipping using AdapterRemoval v2
* Mapping using BWA
* Library Complexity Estimation with Preseq
* Conversion and Filtering of BAM files using Samtools
* Damage assessment via DamageProfiler, additional filtering using PMDTools
* Duplication removal via DeDup 
* BAM Clipping with BamUtil for UDGhalf protocols
* QualiMap BAM quality control analysis

Furthermore, this already creates an interactive report using MultiQC, which will be upgraded in V2.1 "Ulm" to contain more aDNA specific metrics.