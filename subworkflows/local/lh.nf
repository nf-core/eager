include { LEEHOM    } from '../../modules/nf-core/modules/leehom/main'
include { CAT_FASTQ as CAT_LEEHOM } from '../../modules/nf-core/modules/cat/fastq/main'

workflow CLIPMERGE_LH {

    take:
    fastq // [ meta, fastqs ]

    main:
    ch_versions = Channel.empty()
    ch_logs_for_mqc     = Channel.empty()

    LEEHOM ( fastq )
    ch_versions = ch_versions.mix(LEEHOM.out.versions)
    ch_logs_for_mqc = ch_logs_for_mqc.mix(LEEHOM.out.log)


    if ( params.clipmerge_mergedonly ) {
        ch_reads_out = LEEHOM.out.fq_pass
    } else {
        // TODO breaks on single-end data as it some reason has a input1 element in `reads` when it gets to cat?
        // Probably due to this https://github.com/nf-core/modules/blob/251015c8bac16ecb55d738362503f17a84c45d18/modules/cat/fastq/main.nf#L10-L11
        // SOLUTION: MAKE OWN CAT LIKE AR2
        ch_reads_for_cat = LEEHOM.out.fq_pass.mix(LEEHOM.out.unmerged_r1_fq_pass, LEEHOM.out.unmerged_r2_fq_pass)
            .groupTuple()
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new['single_end'] = true

                [ meta_new, reads]
            }
            .dump(tag: "in_cat_leehom_metaupdate")


        CAT_LEEHOM (ch_reads_for_cat)

        ch_reads_out = CAT_LEEHOM.out.reads
        ch_versions = ch_versions.mix(CAT_LEEHOM.out.versions)
    }

    emit:
    reads    = ch_reads_out
    versions = ch_versions
    mqc      = ch_logs_for_mqc

}
