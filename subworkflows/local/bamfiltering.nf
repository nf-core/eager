//
// Filter BAMs for mapping quality, length, unmapped etc.
//

include { SAMTOOLS_VIEW as SAMTOOLS_MAPPED    } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_UNMAPPED  } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ                      } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_STATS                      } from '../../modules/nf-core/samtools/stats/main'

workflow BAM_FILTER {

    take:
    bam // [ [meta], [bai] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    /* possible conditions


    filtering  mapped fastq  (a.k.a. metagenomic mapped)
    filtering  mapped bam
    filtering  mapped both

    filtering  unmapped fastq (a.k.a. metagenomic unmapped)
    filtering  unmapped bam
    filtering  unmapped both

    filtering  mapped/unmapped fastq (a.k.a. metagenomic unmapped)
    filtering  mapped/unmapped bam
    filtering  mapped/unmapped both

    withfiltering length filter
    nofiltering length filter

    */

    // Mapping quality filtering and separation of metagenomics
    if (  ( params.run_metagenomicscreening && metagenomicscreening_input == 'mapped' ) ) {
        SAMTOOLS_MAPPED ( bam )
        SAMTOOLS_FASTQ  ( SAMTOOLS_MAPPED.out.bam )

        ch_versions = ch_versions.mix(SAMTOOLS_MAPPED.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        // TODO account for `--reference_largeref`
        ch_bam_for_genomics     = SAMTOOLS_MAPPED.out.bam.join(SAMTOOLS_MAPPED.out.bam)
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ.out.fastq

    } else if ( params.run_metagenomicscreening && metagenomicscreening_input == 'unmapped' ) {
        SAMTOOLS_MAPPED   ( bam )
        SAMTOOLS_UNMAPPED ( bam )
        SAMTOOLS_FASTQ    ( SAMTOOLS_UNMAPPED.out.bam )
        ch_versions = ch_versions.mix(SAMTOOLS_MAPPED.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_UNMAPPED.out.versions)
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ.out.versions)

        // TODO account for `--reference_largeref`
        ch_bam_for_genomics     = SAMTOOLS_MAPPED.out.bam.join(SAMTOOLS_MAPPED.out.bam)
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ.out.fastq

    } else if ( !params.run_metagenomicscreening ) {
        SAMTOOLS_MAPPED   ( bam )
        ch_versions = ch_versions.mix(SAMTOOLS_MAPPED.out.versions)

        // TODO account for `--reference_largeref`
        ch_bam_for_genomics     = SAMTOOLS_MAPPED.out.bam.join(SAMTOOLS_MAPPED.out.bam)
        ch_bam_for_metagenomics = Channel.empty
    }

    // Length filtering
    if ( params.bamfiltering_readlength != 0  ) {}



    emit:
    genomics_bam     = ch_bam_for_genomics
    genomics_bai     = ch_bam_for_genomics
    metagenomics     = ch_bam_for_metagenomics

}
