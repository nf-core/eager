# ![nf-core/eager](docs/images/eager_logo.png)

**A fully reproducible ancient and modern DNA pipeline in Nextflow and with cloud support.**.

[![Build Status](https://travis-ci.com/nf-core/eager.svg?branch=master)](https://travis-ci.com/nf-core/eager)
[![GitHub Actions CI Status](https://github.com/nf-core/eager/workflows/nf-core%20CI/badge.svg)](https://github.com/nf-core/eager/actions)
[![GitHub Actions Linting Status](https://github.com/nf-core/eager/workflows/nf-core%20linting/badge.svg)](https://github.com/nf-core/eager/actions)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![Slack Status](https://nf-core-invite.herokuapp.com/badge.svg)](https://nf-core-invite.herokuapp.com)[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg)](http://bioconda.github.io/)
[![Docker Container available](https://img.shields.io/docker/automated/nfcore/eager.svg)](https://hub.docker.com/r/nfcore/eager/)
![Singularity Container available](https://img.shields.io/badge/singularity-available-7E4C74.svg)
[![DOI](https://zenodo.org/badge/135918251.svg)](https://zenodo.org/badge/latestdoi/135918251)

## Introduction

**nf-core/eager** is a bioinformatics best-practice analysis pipeline for NGS sequencing based ancient DNA (aDNA) data analysis.

The pipeline uses [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool. It pre-processes raw data from FASTQ inputs, or preprocessed BAM inputs, and can align reads and performs extensive general NGS and aDNA specific quality-control on the results. It comes with docker, singularity or conda containers making installation trivial and results highly reproducible.

## Pipeline steps

By default the pipeline currently performs the following:

* Create reference genome indices for mapping (`bwa`, `samtools`, and `picard`)
* Sequencing quality control (`FastQC`)
* Sequencing adapter removal and for paired end data merging (`AdapterRemoval`)
* Read mapping to reference using (`bwa aln`, `bwa mem` or `CircularMapper`)
* Post-mapping processing, statistics and conversion to bam (`samtools`)
* Ancient DNA C-to-T damage pattern visualisation (`DamageProfiler`)
* PCR duplicate removal (`DeDup` or `MarkDuplicates`)
* Post-mapping statistics and BAM quality control (`Qualimap`)
* Library Complexity Estimation (`preseq`)
* Overall pipeline statistics summaries (`MultiQC`)

Additional functionality contained by the pipeline currently includes:

* Illumina two-coloured sequencer poly-G tail removal (`fastp`)
* Automatic conversion of unmapped reads to FASTQ (`samtools`)
* Damage removal/clipping for UDG+/UDG-half treatment protocols (`BamUtil`)
* Damage reads extraction and assessment (`PMDTools`)

## Quick Start

1. Install [`nextflow`](https://nf-co.re/usage/installation)

2. Install one of [`docker`](https://docs.docker.com/engine/installation/), [`singularity`](https://www.sylabs.io/guides/3.0/user-guide/) or [`conda`](https://conda.io/miniconda.html)

3. Download the EAGER pipeline

        nextflow pull nf-core/eager

4. Test the pipeline using the provided test data

        nextflow run nf-core/eager -profile <docker/singularity/conda>,test --pairedEnd

5. Start running your own ancient DNA analysis!

        nextflow run nf-core/eager -profile <docker/singularity/conda> --reads'*_R{1,2}.fastq.gz' --fasta '<REFERENCE>.fasta'

6. Once your run has completed successfully, clean up the intermediate files.

        nextflow clean -k


NB. You can see an overview of the run in the MultiQC report located at `<OUTPUT_DIR>/MultiQC/multiqc_report.html`

Modifications to the default pipeline are easily made using various options
as described in the documentation.

## Documentation

The nf-core/eager pipeline comes with documentation about the pipeline, found in the `docs/` directory or on the main homepage of the nf-core project:

1. [Nextflow Installation](https://nf-co.re/usage/installation)
2. Pipeline configuration
    * [Pipeline installation](https://nf-co.re/usage/local_installation)
    * [Adding your own system config](https://nf-co.re/usage/adding_own_config)
    * [Reference genomes](https://nf-co.re/usage/reference_genomes)
3. [Running the pipeline](docs/usage.md)
4. [Output and how to interpret the results](docs/output.md)
5. [EAGER2 Code Contribution Guidelines](code_contribution.md)
6. [nf-core/nextflow Troubleshooting](https://nf-co.re/usage/troubleshooting)
7. [EAGER Troubleshooting](docs/troubleshooting.md)

## Credits

This pipeline was written by Alexander Peltzer ([apeltzer](https://github.com/apeltzer)), with major contributions from Stephen Clayton, ideas and documentation from James A. Fellows Yates, Raphael Eisenhofer, Maxime Borry and Judith Neukamm. If you want to contribute, please open an issue and ask to be added to the project - happy to do so and everyone is welcome to contribute here!

## Contributors

* [James A. Fellows-Yates](https://github.com/jfy133)
* [Stephen Clayton](https://github.com/sc13-bioinf)
* [Maxime Borry](https://github.com/maxibor)
* [Judith Neukamm](https://github.com/JudithNeukamm)
* [Raphael Eisenhofer](https://github.com/EisenRa)
* [Maxime Garcia](https://github.com/MaxUlysse)
* [Luc Venturini](https://github.com/lucventurini)
* [Hester van Schalkwyk](https://github.com/hesterjvs)

If you've contributed and you're missing in here, please let me know and I'll add you in.

## Tool References

* **EAGER v1**, CircularMapper, DeDup* Peltzer, A., Jäger, G., Herbig, A., Seitz, A., Kniep, C., Krause, J., & Nieselt, K. (2016). EAGER: efficient ancient genome reconstruction. Genome Biology, 17(1), 1–14. [https://doi.org/10.1186/s13059-016-0918-z](https://doi.org/10.1186/s13059-016-0918-z)  Download: [https://github.com/apeltzer/EAGER-GUI](https://github.com/apeltzer/EAGER-GUI) and [https://github.com/apeltzer/EAGER-CLI](https://github.com/apeltzer/EAGER-CLI)
* **FastQC** download: [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* **AdapterRemoval v2** Schubert, M., Lindgreen, S., & Orlando, L. (2016). AdapterRemoval v2: rapid adapter trimming, identification, and read merging. BMC Research Notes, 9, 88. [https://doi.org/10.1186/s13104-016-1900-2](https://doi.org/10.1186/s13104-016-1900-2) Download: [https://github.com/MikkelSchubert/adapterremoval](https://github.com/MikkelSchubert/adapterremoval)
* **bwa** Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. Bioinformatics , 25(14), 1754–1760. [https://doi.org/10.1093/bioinformatics/btp324](https://doi.org/10.1093/bioinformatics/btp324) Download: [http://bio-bwa.sourceforge.net/bwa.shtml](http://bio-bwa.sourceforge.net/bwa.shtml)
* **SAMtools** Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., … 1000 Genome Project Data Processing Subgroup. (2009). The Sequence Alignment/Map format and SAMtools. Bioinformatics , 25(16), 2078–2079. [https://doi.org/10.1093/bioinformatics/btp352](https://doi.org/10.1093/bioinformatics/btp352) Download: [http://www.htslib.org/](http://www.htslib.org/)
* **DamageProfiler** Judith Neukamm (Unpublished)
* **QualiMap** Okonechnikov, K., Conesa, A., & García-Alcalde, F. (2016). Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics , 32(2), 292–294. [https://doi.org/10.1093/bioinformatics/btv566](https://doi.org/10.1093/bioinformatics/btv566) Download: [http://qualimap.bioinfo.cipf.es/](http://qualimap.bioinfo.cipf.es/)
* **preseq** Daley, T., & Smith, A. D. (2013). Predicting the molecular complexity of sequencing libraries. Nature Methods, 10(4), 325–327. [https://doi.org/10.1038/nmeth.2375](https://doi.org/10.1038/nmeth.2375). Download: [http://smithlabresearch.org/software/preseq/](http://smithlabresearch.org/software/preseq/)
* **PMDTools** Skoglund, P., Northoff, B. H., Shunkov, M. V., Derevianko, A. P., Pääbo, S., Krause, J., & Jakobsson, M. (2014). Separating endogenous ancient DNA from modern day contamination in a Siberian Neandertal. Proceedings of the National Academy of Sciences of the United States of America, 111(6), 2229–2234. [https://doi.org/10.1073/pnas.1318934111](https://doi.org/10.1073/pnas.1318934111) Download: [https://github.com/pontussk/PMDtools](https://github.com/pontussk/PMDtools)
* **MultiQC** Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354) Download: [https://multiqc.info/](https://multiqc.info/)
* **BamUtils** Jun, G., Wing, M. K., Abecasis, G. R., & Kang, H. M. (2015). An efficient and scalable analysis framework for variant extraction and refinement from population-scale DNA sequence data. Genome Research, 25(6), 918–925. [https://doi.org/10.1101/gr.176552.114](https://doi.org/10.1101/gr.176552.114) Download: [https://genome.sph.umich.edu/wiki/BamUtil](https://genome.sph.umich.edu/wiki/BamUtil)
* **FastP** Chen, S., Zhou, Y., Chen, Y., & Gu, J. (2018). fastp: an ultra-fast all-in-one FASTQ preprocessor. Bioinformatics , 34(17), i884–i890. [https://doi.org/10.1093/bioinformatics/bty560](https://doi.org/10.1093/bioinformatics/bty560) Download: [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp)
* **GATK 3.8** DePristo, M. A., Banks, E., Poplin, R., Garimella, K. V., Maguire, J. R., Hartl, C., … Daly, M. J. (2011). A framework for variation discovery and genotyping using next-generation DNA sequencing data. Nature Genetics, 43(5), 491–498. [https://doi.org/10.1038/ng.806](https://doi.org/10.1038/ng.806.) [Download](https://software.broadinstitute.org/gatk/download/)
* **GATK 4.X** - no citation available yet
