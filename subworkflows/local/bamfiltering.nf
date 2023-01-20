//
// Filter BAMs for mapping quality, length, unmapped etc.
//

include { SAMTOOLS_VIEW as SAMTOOLS_MAPPED          } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_VIEW as SAMTOOLS_UNMAPPED        } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED   } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_STATS                            } from '../../modules/nf-core/samtools/stats/main'

workflow BAM_FILTER {

    take:
    bam // [ [meta], [bai] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    /* possible conditions


    filtering  mapped bam
    filtering  mapped fastq (a.k.a. metagenomic mapped)
    filtering  mapped both (Confiugred in publishDir !)

    filtering  filtering && unmapped bam
    filtering  filtering && unmapped fastq (a.k.a. metagenomic unmapped)
    filtering  filtering && unmapped both

    filtering  mapped/unmapped fastq (a.k.a. metagenomic unmapped)
    filtering  mapped/unmapped bam
    filtering  mapped/unmapped both

    withfiltering length filter
    nofiltering length filter

    */


    // TODO: add a check in eager.nf: if ( params.run_metagenomicscreening && !params.bamfiltering_generatefastq ) { exit 1, "If you wish to run metagneomic screening yo umust include --bamfiltering_generatefastq" }
    // TODO: add check in eager.nf: if ( params.run_metagenomicscreening && --metagenomicscreening_input 'unmapped' && !params.bamfiltering_generateunmapped ) { exit 1, "If you want unmapped input, you need to also generate this with `--bamfiltering_generatedunmmaped`  }
    // TODO: write length filter module!
    // TODO: add filter flags (-q, -f, -F) to modules.conf
    // TODO: add all params to nextflow.config & schema

    // Run map quality filtering
    SAMTOOLS_MAPPED ( bam )
    SAMTOOLS_INDEX ( SAMTOOLS_MAPPED.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_MAPPED.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    // (Filtered) Mapped BAM Generation
    if ( params.bamfiltering_readlength != 0  ) {
        ch_bam_for_mappedfqconvert = LENGTH_FILTER_MAPPED ( SAMTOOLS_MAPPED.out.bam )
        ch_versions = ch_versions.mix( LENGTH_FILTER_MAPPED.out.versions )
    } else {
        ch_bam_for_mappedfqconvert = SAMTOOLS_MAPPED.out.bam
    }

    ch_bai_for_genomics = SAMTOOLS_INDEX ( ch_bam_for_mappedfqconvert )

    // (Filtered) Mapped FASTQ Generation
    if ( params.bamfiltering_generatefastq ) {
        ch_bam_for_genomics = SAMTOOLS_FASTQ_MAPPED ( ch_bam_for_mappedfqconvert, false )
        ch_versions = ch_versions.mix(SAMTOOLS_FASTQ_MAPPED.out.versions)
    } else {
        ch_bam_for_genomics = ch_bam_for_mappedfqconvert
    }

    // (Filtered) Unmapped BAM Generation
    if ( params.bamfiltering_generateunmapped  && params.bamfiltering_readlength != 0 ) {
        LENGTH_FILTER_UNMAPPED ( bam )
        ch_unmapped_bam = SAMTOOLS_UNMAPPED ( LENGTH_FILTER_UNMAPPED.out.bam )
        ch_versions = ch_versions.mix( LENGTH_FILTER_UNMAPPED.out.versions )
        ch_versions = ch_versions.mix( SAMTOOLS_UNMAPPED.out.versions )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_UNMAPPED.out.versions )
    } else ( params.bamfiltering_generateunmapped  && params.bamfiltering_readlength == 0 ) {
        ch_unmapped_bam = SAMTOOLS_UNMAPPED ( bam )
        ch_versions = ch_versions.mix( SAMTOOLS_UNMAPPED.out.versions )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_UNMAPPED.out.versions )
    }

    // (Filtered) Unmapped FASTQ Generation
    if ( params.bamfiltering_generatefastq ) {
        ch_unmapped_fastq = SAMTOOLS_FASTQ_UNMAPPED ( ch_unmapped_bam, false )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_UNMAPPED.out.versions )
    } else {
        // Metagenomic screening requires FASTQ files
        ch_bam_for_metagenomics = Channel.empty()
    }

    // Routing for metagenomic screening -> TODO CHECK!
    if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) {
        ch_bam_for_metagenomics = ch_bam_for_metagenomics
    } else if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'mapped') {
        ch_bam_for_metagenomics = ch_bam_for_genomics
    } else if ( !params.run_metagenomicscreening ) {
        ch_bamfiltered_for_metagenomics = Channel.empty()
    }

    emit:
    genomics_bam     = ch_bam_for_genomics
    genomics_bai     = ch_bai_for_genomics // TODO Join here or later?
    metagenomics     = ch_bam_for_metagenomics

}
