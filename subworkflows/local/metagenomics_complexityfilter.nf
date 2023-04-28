include { BBMAP_BBDUK                 } from '../../modules/nf-core/bbmap/bbduk/main'
include { PRINSEQPLUSPLUS             } from '../../modules/nf-core/prinseqplusplus/main'


workflow METAGENOMICS_COMPLEXITYFILTER {
    take: ch_bamfiltered_for_metagenomics // [meta, fastq]

    main:

    //
    // Take the selected complexity filter tool
    //

    if (params.metagenomics_complexity_tool == 'bbduk') {
        BBMAP_BBDUK( ch_bamfiltered_for_metagenomics, [] )
        ch_versions = BBMAP_BBDUK.out.versions.first()
        ch_reads_for_metagenomics = BBMAP_BBDUK.out.reads
    }
    else if ( params.metagenomics_complexity_tool == 'prinseq' ) {
        // check if e.g. dustscore is set but entropy enabled
        PRINSEQPLUSPLUS ( ch_bamfiltered_for_metagenomics )
        ch_versions = PRINSEQPLUSPLUS.out.versions.first()
        ch_reads_for_metagenomics = PRINSEQPLUSPLUS.out.good_reads
    }

    emit:
        fastq    = ch_reads_for_metagenomics
        versions = ch_versions
}