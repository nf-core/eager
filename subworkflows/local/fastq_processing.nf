// Prepare various reference FASTA index files for downstream steps

include { FASTP as FASTP_POLYG_TRIM         } from '../../modules/nf-core/modules/fastp/main'
include { FASTP as FASTP_POLYX_TRIM         } from '../../modules/nf-core/modules/fastp/main'
include { ADAPTERREMOVAL                    } from '../../modules/nf-core/modules/adapterremoval/main'
include { LEEHOM                            } from '../../modules/nf-core/modules/leehom/main'
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
        .dump(tag: "post_fastp_pg_mix")
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

    ch_final_fastq = ch_polyx_out

    if ( !params.skip_fastqc ) {
        FASTQC_AFTER_PROCESSING ( ch_final_fastq )
        ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTQC_AFTER_PROCESSING.out.zip )
        ch_versions = ch_versions.mix( FASTQC_AFTER_PROCESSING.out.versions ).dump(tag: "fastq_processing_versions")
    }

    emit:
    fastq       = ch_final_fastq
    mqc         = ch_logs_for_mqc
    versions    = ch_versions

}
