include { LEEHOM     } from '../../modules/nf-core/modules/leehom/main'
include { CAT_LEEHOM } from '../../modules/local/cat_leehom'

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
        // TODO ask Gabriel/Kay if this makes sense?
        // Only combine singletons into final fastq for PE data -
        // SE shouldn't have singletons
        ch_reads_for_cat = LEEHOM.out.fq_pass.mix(LEEHOM.out.unmerged_r1_fq_pass, LEEHOM.out.unmerged_r2_fq_pass)
            .groupTuple()
            .branch{
                se: it[0]['single_end']
                pe: !it[0]['single_end']
            }

        CAT_LEEHOM (ch_reads_for_cat.pe)

        // if we are not keep reads separate, should set as single so we
        // know downstream to do single-end mapping, so set as SE
        ch_reads_out = CAT_LEEHOM.out.reads
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new['single_end'] = true

                [ meta_new, reads]
            }
            .mix(ch_reads_for_cat.se)
            .dump(tag: "out_mix_pe_se_lh_withsingles")

        ch_reads_out = ch_reads_out
        ch_versions  = ch_versions.mix(CAT_LEEHOM.out.versions)
    }

    emit:
    reads    = ch_reads_out
    versions = ch_versions
    mqc      = ch_logs_for_mqc

}
