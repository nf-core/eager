//
//  Generate genomic region intervals from samtools idxstats output
//

include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/samtools/idxstats/main'
include { BUILD_INTERVALS   } from '../../modules/local/build_intervals'

workflow BUILD_INTERVALS_FOR_BAM_SPLIT {
    take:
    ch_bam_bai // [ [ meta ], bam , bai ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    SAMTOOLS_IDXSTATS(ch_bam_bai)
    ch_versions   = ch_versions.mix ( SAMTOOLS_IDXSTATS.out.versions )

    BUILD_INTERVALS(SAMTOOLS_IDXSTATS.out.idxstats)
    ch_versions   = ch_versions.mix ( BUILD_INTERVALS.out.versions )

    emit:
    bed = BUILD_INTERVALS.out.bed  // [ [ meta ], bed ]
    versions = ch_versions
}
