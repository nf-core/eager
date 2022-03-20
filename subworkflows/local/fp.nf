include { FASTP } from '../../modules/nf-core/modules/fastp/main'

workflow CLIPMERGE_FP {

    take:
    fastq // [ meta, fastqs ]

    main:
    ch_versions      = Channel.empty()
    ch_logs_for_mqc  = Channel.empty()

    FASTP ( fastq, true, true )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    ch_logs_for_mqc = ch_logs_for_mqc.mix(FASTP.out.json)

    if ( params.clipmerge_skipcollapse ) {
        // PE mapping mode: SE trimmed as normal, PE unpaired reads (from reads),
        // so typical usage but with no meta update
        // TODO switch to filters from outputs then mix
        ch_reads_out =  FASTP.out.reads
                            .dump(tag: "out_mix_pe_se_fp_skipcollapse")
    } else {
        // Standard: SE trimmed, PE merged
        // (see modules.conf for when params.clipmerge_mergedonly)
        println("Bar")
        ch_reads_fp_se = FASTP.out.reads.filter { meta, reads -> meta['single_end'] }

        ch_reads_fp_pe = FASTP.out.reads_merged
            .filter { meta, reads -> !meta['single_end'] }
            .map {
                meta, reads ->
                    def meta_new = meta.clone()
                    meta_new['single_end'] = true

                [ meta_new, reads]
            }

        ch_reads_out = ch_reads_fp_se
                            .mix( ch_reads_fp_pe )
                            .dump(tag: "out_mix_pe_se_fp_standard")

    }

    emit:
    reads    = ch_reads_out
    versions = ch_versions
    mqc      = ch_logs_for_mqc

}

/*
        ch_reads_out = CAT_FASTP.out.reads

            .mix(ch_reads_for_cat.se)
            */
