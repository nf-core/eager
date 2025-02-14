//
// Prepare reference indexing for downstream
//

include { SAMTOOLS_MERGE    as SAMTOOLS_MERGE_LANES_BAMINPUT        } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT     as SAMTOOLS_SORT_MERGED_LANES_BAMINPUT  } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_MAPPED_BAMINPUT    } from '../../modules/nf-core/samtools/flagstat/main'

workflow MERGE_LANES_INPUTBAM {
    take:
    bams // [ [meta], [bams] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    ch_input_for_lane_merge = bams
                                .map { meta, bam -> [ meta.clone().findAll{ it.key !in ['lane', 'colour_chemistry', 'shard_number'] }, bam ] }
                                .groupTuple()
                                .branch {
                                    meta, bam ->
                                        merge: bam.size() > 1
                                        skip: true
                                }

    SAMTOOLS_MERGE_LANES_BAMINPUT ( ch_input_for_lane_merge.merge, [[], []], [[], []] )
    ch_versions.mix( SAMTOOLS_MERGE_LANES_BAMINPUT.out.versions )

    // Then mix back merged and single lane libraries for everything downstream
    ch_mergedlanes_for_sorting = ch_input_for_lane_merge.skip
                                    .mix( SAMTOOLS_MERGE_LANES_BAMINPUT.out.bam )

    SAMTOOLS_SORT_MERGED_LANES_BAMINPUT ( ch_mergedlanes_for_sorting )
    ch_mapped_bam = SAMTOOLS_SORT_MERGED_LANES_BAMINPUT.out.bam
    ch_versions.mix( SAMTOOLS_SORT_MERGED_LANES_BAMINPUT.out.versions )

    SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT( ch_mapped_bam )
    ch_mapped_bai =  params.fasta_largeref ? SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT.out.csi : SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT.out.bai
    ch_versions.mix( SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT.out.versions )

    ch_input_for_flagstat = SAMTOOLS_SORT_MERGED_LANES_BAMINPUT.out.bam.join( SAMTOOLS_INDEX_MERGED_LANES_BAMINPUT.out.bai, failOnMismatch: true )

    SAMTOOLS_FLAGSTAT_MAPPED_BAMINPUT ( ch_input_for_flagstat )
    ch_versions.mix( SAMTOOLS_FLAGSTAT_MAPPED_BAMINPUT.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_MAPPED_BAMINPUT.out.flagstat )

    emit:
    bam        = ch_mapped_bam                            // [ [ meta ], bam ]
    bai        = ch_mapped_bai                            // [ [ meta ], bai/csi ]
    flagstat   = SAMTOOLS_FLAGSTAT_MAPPED_BAMINPUT.out.flagstat    // [ [ meta ], stats ]
    mqc        = ch_multiqc_files
    versions   = ch_versions

}
