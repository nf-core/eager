include { ADAPTERREMOVAL     } from '../../modules/nf-core/modules/adapterremoval/main'
include { CAT_ADAPTERREMOVAL } from '../../modules/local/cat_adapterremoval'

workflow CLIPMERGE_AR {

    take:
    fastq // [ meta, fastqs ]

    main:
    ch_versions = Channel.empty()
    ch_logs_for_mqc     = Channel.empty()

    ADAPTERREMOVAL ( fastq )
    ch_versions = ch_versions.mix(ADAPTERREMOVAL.out.versions)
    ch_logs_for_mqc = ch_logs_for_mqc.mix(ADAPTERREMOVAL.out.log)

    ch_adapterremoval_into_cat = ADAPTERREMOVAL.out.singles_truncated.mix(
        ADAPTERREMOVAL.out.pair1_truncated,
        ADAPTERREMOVAL.out.pair2_truncated,
        ADAPTERREMOVAL.out.collapsed,
        ADAPTERREMOVAL.out.collapsed_truncated
    )
    .groupTuple()
    .dump(tag: "in_cat_ar")

    CAT_ADAPTERREMOVAL ( ch_adapterremoval_into_cat )
    ch_reads_out = CAT_ADAPTERREMOVAL.out.reads
    ch_versions = ch_versions.mix(CAT_ADAPTERREMOVAL.out.versions)


    emit:
    reads    = ch_reads_out
    versions = ch_versions
    mqc      = ch_logs_for_mqc

}
