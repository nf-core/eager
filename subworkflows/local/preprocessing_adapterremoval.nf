//
// Process short raw reads with AdapterRemoval
//

include { ADAPTERREMOVAL as ADAPTERREMOVAL_SINGLE       } from '../../modules/nf-core/adapterremoval/main'
include { ADAPTERREMOVAL as ADAPTERREMOVAL_PAIRED       } from '../../modules/nf-core/adapterremoval/main'
include { CAT_FASTQ as CAT_FASTQ_ADAPTERREMOVAL         } from '../../modules/nf-core/cat/fastq/main'

workflow PREPROCESSING_ADAPTERREMOVAL {

    take:
    reads // [[meta], [reads]]
    adapterlist

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files      = Channel.empty()

    ch_input_for_adapterremoval = reads
                                    .branch{
                                        single: it[0].single_end
                                        paired: !it[0].single_end
                                    }

    ADAPTERREMOVAL_SINGLE ( ch_input_for_adapterremoval.single, adapterlist )
    ch_versions = ch_versions.mix( ADAPTERREMOVAL_SINGLE.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_SINGLE.out.settings)

    ADAPTERREMOVAL_PAIRED ( ch_input_for_adapterremoval.paired, adapterlist )
    ch_versions = ch_versions.mix( ADAPTERREMOVAL_PAIRED.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix(ADAPTERREMOVAL_PAIRED.out.settings)

    /*
     * Due to the ugly output implementation of the current AdapterRemoval2 version, each file
     * has to be exported in a separate channel and we must manually recombine when necessary.
     */

    // Merge and keep singletons
    if ( !params.preprocessing_skippairmerging && !params.preprocessing_excludeunmerged ) {

        ch_concat_fastq = Channel.empty()
            .mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated,
                ADAPTERREMOVAL_PAIRED.out.singles_truncated,
                ADAPTERREMOVAL_PAIRED.out.paired_truncated
            )
            .map { meta, reads ->
                def meta_new = meta.clone()
                meta_new.single_end = true
                [meta_new, reads]
            }
            .groupTuple()
            // Paired-end reads cause a nested tuple during grouping.
            // We want to present a flat list of files to `CAT_FASTQ`.
            .map { meta, fastq -> [meta, fastq.flatten()] }


        CAT_FASTQ_ADAPTERREMOVAL(ch_concat_fastq)

        ch_adapterremoval_reads_prepped = CAT_FASTQ_ADAPTERREMOVAL.out.reads
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    // Merge and exclude singletons
    } else if ( !params.preprocessing_skippairmerging && params.preprocessing_excludeunmerged ) {

        ch_concat_fastq = Channel.empty()
            .mix(
                ADAPTERREMOVAL_PAIRED.out.collapsed,
                ADAPTERREMOVAL_PAIRED.out.collapsed_truncated
            )
            .map { meta, reads ->
                def meta_new = meta.clone()
                meta_new.single_end = true
                [meta_new, reads]
            }
            .groupTuple()
            .map { meta, fastq -> [meta, fastq.flatten()] }


        CAT_FASTQ_ADAPTERREMOVAL(ch_concat_fastq)

        ch_adapterremoval_reads_prepped = CAT_FASTQ_ADAPTERREMOVAL.out.reads
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    // Don't merge (no singletons will exist)
    } else {

        ch_adapterremoval_reads_prepped = ADAPTERREMOVAL_PAIRED.out.paired_truncated
            .mix(ADAPTERREMOVAL_SINGLE.out.singles_truncated)

    }

    emit:
    reads    = ch_adapterremoval_reads_prepped  // channel: [ val(meta), [ reads ] ]
    versions = ch_versions  // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

