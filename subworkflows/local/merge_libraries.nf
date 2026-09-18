//
// Merge libraries of the same sample, then sort, index, and flagstat the merged bam
//

include { addNewMetaFromAttributes                                } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

include { SAMTOOLS_MERGE    as SAMTOOLS_MERGE_LIBRARIES           } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT     as SAMTOOLS_SORT_MERGED_LIBRARIES     } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_MERGED_LIBRARIES    } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_MERGED_LIBRARIES } from '../../modules/nf-core/samtools/flagstat/main'

workflow MERGE_LIBRARIES {
    take:
    ch_bam_bai  // [ [ meta ], bam , bai ]

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    ch_library_merge_input = ch_bam_bai
        .map { addNewMetaFromAttributes( it, ["id", "sample_id", "strandedness", "reference"], ["id", "sample_id", "strandedness", "reference"], false ) }
        .groupTuple(by: 0)
        // Discrad library-level metas, and bais. Add single_end: true to all metas (no SE/PE distinction from here on)
        .map {
            meta, lib_metas, bam, bai ->
            [ meta + [ 'single_end': true ], bam ]
        }

    SAMTOOLS_MERGE_LIBRARIES ( ch_library_merge_input, [[], []], [[], []] )
    ch_versions = ch_versions.mix( SAMTOOLS_MERGE_LIBRARIES.out.versions.first() )

    SAMTOOLS_SORT_MERGED_LIBRARIES ( SAMTOOLS_MERGE_LIBRARIES.out.bam )
    ch_versions = ch_versions.mix( SAMTOOLS_SORT_MERGED_LIBRARIES.out.versions.first() )

    SAMTOOLS_INDEX_MERGED_LIBRARIES ( SAMTOOLS_SORT_MERGED_LIBRARIES.out.bam )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_MERGED_LIBRARIES.out.versions.first() )

    // Join merged sample-level bams and their bais for genotyping
    ch_merged_bams = SAMTOOLS_SORT_MERGED_LIBRARIES.out.bam
        .join( SAMTOOLS_INDEX_MERGED_LIBRARIES.out.bai )

    // Not sure if FLAGSTAT is really needed, but added here for completeness
    SAMTOOLS_FLAGSTAT_MERGED_LIBRARIES ( ch_merged_bams )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTAT_MERGED_LIBRARIES.out.versions.first() )

    ch_merged_flagstat = SAMTOOLS_FLAGSTAT_MERGED_LIBRARIES.out.flagstat
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_MERGED_LIBRARIES.out.flagstat )

    emit:
    bam_bai  = ch_merged_bams       // [ [ meta ], bam , bai ]
    flagstat = ch_merged_flagstat   // [ [ meta ], flagstat ]
    versions = ch_versions
    mqc      = ch_multiqc_files // Same as flagstat
}
