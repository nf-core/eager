// Prepare various reference FASTA index files for downstream steps

include { FASTP as FASTP_POLYG_TRIM         } from '../../modules/nf-core/modules/fastp/main'
include { FASTP as FASTP_POLYX_TRIM         } from '../../modules/nf-core/modules/fastp/main'

include { CLIPMERGE_LH                      } from './lh.nf'
include { CLIPMERGE_AR                      } from './ar.nf'

// TODO add fastp as an actual trim/merger?
include { FASTP as FASTP_ENDTRIM            } from '../../modules/nf-core/modules/fastp/main'
include { FASTQC as FASTQC_AFTER_PROCESSING } from '../../modules/nf-core/modules/fastqc/main'


workflow FASTQ_PROCESSING {

    take:
    fastq

    main:
    // Set versions channel
    ch_versions = Channel.empty()
    ch_logs_for_mqc = Channel.empty()

    // polyG trimming
    if ( params.fastp_trim_polyg_run ) {
        ch_fastq_for_polygtrim = fastq
                                    .branch{
                                        twocol: it[0]['colour_chemistry'] == '2'  // Nextseq/Novaseq data with possible sequencing artefact
                                        fourcol: it[0]['colour_chemistry'] == '4' // HiSeq/MiSeq data where polyGs would be true
                                    }

    } else {
        ch_fastq_for_polygtrim = fastq
                                    .branch{
                                        twocol: it[0]['colour_chemistry'] == 'dummy'                                   // No clipping requested, so no need to send to fastp
                                        fourcol: it[0]['colour_chemistry'] == '4' || it[0]['colour_chemistry'] == '2'  // HiSeq/MiSeq data where polyGs would be true
                                    }

    }

    // Trim polyG tails only when 2 colour data present and requested
    FASTP_POLYG_TRIM ( ch_fastq_for_polygtrim.twocol, true, false )  // keeping fails as presumably shouldn't exist?
    ch_versions = ch_versions.mix( FASTP_POLYG_TRIM.out.versions )
    ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTP_POLYG_TRIM.out.json )

    ch_fastq_for_polygtrim.fourcol
        .mix( FASTP_POLYG_TRIM.out.reads )
        .dump(tag: "out_mix_post_fastp_pg")
        .set{ ch_polyg_for_polyx }

    // Trim poly X tails when any data requested
    if ( params.fastp_trim_polyx_run ) {
        FASTP_POLYX_TRIM ( ch_polyg_for_polyx, true, false )
        ch_polyx_out = FASTP_POLYX_TRIM.out.reads
        ch_versions = ch_versions.mix( FASTP_POLYX_TRIM.out.versions )
        ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTP_POLYX_TRIM.out.json )
    } else {
        ch_polyx_out = ch_polyg_for_polyx
    }

    // Adapter removing and merging
    // TODO Adapterremoval
    // TODO LeeHom
    // TODO FastP
    // Switch to when? As secondary processes (cat/fixprefix) only exectued if files go in anyway?
    //ADAPTERREMOVAL( ch_polyx_out )
        //ADAPTERREMOVAL_COMBINE () // to write
        //if ( params.deduplication_tool = 'dedup' ) // missingmodule

    if ( params.clipmerge_tool == 'adapterremoval' ) {
        CLIPMERGE_AR ( ch_polyx_out )
        ch_clipmerge_out = CLIPMERGE_AR.out.reads

        ch_versions = ch_versions.mix( CLIPMERGE_AR.out.versions )
        ch_logs_for_mqc = ch_logs_for_mqc.mix( CLIPMERGE_AR.out.mqc )
    } else if ( params.clipmerge_tool == 'leehom' ) {
        CLIPMERGE_LH ( ch_polyx_out )
        ch_clipmerge_out = CLIPMERGE_LH.out.reads

        ch_versions = ch_versions.mix( CLIPMERGE_LH.out.versions )
        ch_logs_for_mqc = ch_logs_for_mqc.mix( CLIPMERGE_LH.out.mqc )
    }


    // Final stats
    ch_final_fastq = ch_clipmerge_out
    ch_final_fastq.dump(tag: "out_clipmerge")

    if ( !params.skip_fastqc ) {
        FASTQC_AFTER_PROCESSING ( ch_final_fastq )
        ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTQC_AFTER_PROCESSING.out.zip )
        ch_versions = ch_versions.mix( FASTQC_AFTER_PROCESSING.out.versions )
    }

    emit:
    fastq       = ch_final_fastq
    mqc         = ch_logs_for_mqc
    versions    = ch_versions

}
