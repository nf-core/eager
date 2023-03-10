//
// Calculate different percent on target as well as complexity from samtools flagstats files
//

include { ENDORSPY_COMPLETE              } from '../../modules/nf-core/endorspy'
include { ENDORSPY_QUALITY_FILTER        } from '../../modules/nf-core/endorspy'
include { ENDORSPY_RAW                   } from '../../modules/nf-core/endorspy'

workflow ENDORSPY {
//outputs for the flagstats will look like:     tuple val(meta), path("*.flagstat"), emit: flagstat

    raw=MAP.out.flagstat
    dedup=DEDUPLICATE.out.flagstat
    filtered=?

    params.whateverfilteringiscalled
    params.skip_deduplication

}
