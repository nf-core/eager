//
// Calculate different percent on target as well as complexity from samtools flagstats files
//

include { ENDORSPY as ENDORSPY_COMPLETE              } from '../../modules/nf-core/endorspy/main'
include { ENDORSPY as ENDORSPY_QUALITY_FILTER        } from '../../modules/nf-core/endorspy/main'
include { ENDORSPY as ENDORSPY_RAW                   } from '../../modules/nf-core/endorspy/main'
include { ENDORSPY as ENDORSPY_NOQF_DEDUP            } from '../../modules/nf-core/endorspy/main'

workflow ENDORSPY {
    take:
    flagstats_raw // [ [meta], [flagstats_raw] ]
    flagstats_filtered // [ [meta], [flagstats_filtered] ]
    flagstats_duplicated // [ [meta], [flagstats_duplicated] ]
//outputs for the flagstats will look like:     tuple val(meta), path("*.flagstat"), emit: flagstat

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.run_bamfiltering && !params.skip_deduplication ) {
        ch_for_endorspy_complete = flagstats_raw
                                                                 .join (flagstats_filtered)
                                                                 .join (flagstats_duplicated)
        ENDORSPY_COMPLETE ( ch_for_endorspy_complete )
    } else if ( params.run_bamfiltering && params.skip_deduplication ) {
        ch_for_endorspy_quality_filtered = flagstats_raw.join (flagstats_filtered)
                                                        .map {
                                                            meta, flags_raw, flags_filtered ->
                                                            [ meta, flags_raw, flags_filtered, [] ]
                                                        }
        ENDORSPY_QUALITY_FILTER ( ch_for_endorspy_quality_filtered)
    } else if ( !params.run_bamfiltering && !params.skip_deduplication) {
        ch_for_endorspy_noqf_dedup = flagstats_raw.join (flagstats_duplicated)
                                                        . map{
                                                            meta, flags_raw, flags_dedup ->
                                                            [ meta, flags_raw, [], flags_dedup ]
                                                            }
        ENDORSPY_QUALITY_FILTER ( ch_for_endorspy_noqf_dedup )
    } else {
        ch_for_endorspy_raw = flagstats_raw.map {
                                                    meta, flags_raw ->
                                                    [ meta, flags_raw, [], [] ]
        }
        ENDORSPY_RAW ( ch_for_endorspy_raw )
    }

}
