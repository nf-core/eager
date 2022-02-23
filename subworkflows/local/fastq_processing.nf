// Prepare various reference FASTA index files for downstream steps

include { FASTP as FASTP_POLYG_TRIM         } from '../../modules/nf-core/modules/fastp/main'
//include { FASTP as FASTP_POLYX_TRIM         } from '../../modules/nf-core/modules/fastp/main'
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
                                        twocol: it[0]['colour_chemistry'] == 'dummy'               // No clipping requested, so no need to send to fastp
                                        fourcol: it[0]['colour_chemistry'] == '4' || it[3] == '2'  // HiSeq/MiSeq data where polyGs would be true
                                    }

    }

    // TODO add new params to schema
    // TODO use only single FASTP with parameter condition to run polyX and/or polyG in `modules.config`?
    FASTP_POLYG_TRIM ( ch_fastq_for_polygtrim.twocol, true, false )  // keeping fails as presumably shouldn't exist?
    ch_versions = ch_versions.mix( FASTP_POLYG_TRIM.out.versions )
    ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTP_POLYG_TRIM.out.json )

    ch_fastq_for_polygtrim.fourcol
        .dump(tag: "post_fastp_merge")
        .mix( FASTP_POLYG_TRIM.out.reads )
        .set{ ch_final_fastq }

    // TODO Make optional
    FASTQC_AFTER_PROCESSING ( ch_final_fastq )//.dump(tag: "fastqc_processing_versions")
    ch_logs_for_mqc = ch_logs_for_mqc.mix( FASTQC_AFTER_PROCESSING.out.zip )
    ch_versions = ch_versions.mix( FASTQC_AFTER_PROCESSING.out.versions ).dump(tag: "fastq_processing_versions")

    emit:
    fastq       = ch_final_fastq
    mqc         = ch_logs_for_mqc
    versions    = ch_versions

}
