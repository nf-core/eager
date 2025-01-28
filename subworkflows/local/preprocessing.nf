//
// Perform read trimming and merging
//


include { PREPROCESSING_FASTP          } from './preprocessing_fastp'
include { PREPROCESSING_ADAPTERREMOVAL } from './preprocessing_adapterremoval'
include { FASTQC as FASTQC_PROCESSED   } from '../../modules/nf-core/fastqc/main'
include { FALCO as FALCO_PROCESSED     } from '../../modules/nf-core/falco/main'

workflow PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]
    adapterlist

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.preprocessing_tool == "fastp" ) {
        ch_processed_reads = PREPROCESSING_FASTP ( reads, adapterlist ).reads
        ch_versions        =  ch_versions.mix( PREPROCESSING_FASTP.out.versions )
        ch_multiqc_files   =  ch_multiqc_files.mix( PREPROCESSING_FASTP.out.mqc )
    } else if ( params.preprocessing_tool == "adapterremoval" ) {
        ch_processed_reads = PREPROCESSING_ADAPTERREMOVAL ( reads, adapterlist ).reads
        ch_versions        = ch_versions.mix( PREPROCESSING_ADAPTERREMOVAL.out.versions )
        ch_multiqc_files   = ch_multiqc_files.mix( PREPROCESSING_ADAPTERREMOVAL.out.mqc )
    } else {
        ch_processed_reads = reads
    }

    if ( params.sequencing_qc_tool == "falco" ) {
        FALCO_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FALCO_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FALCO_PROCESSED.out.txt )
    } else {
        FASTQC_PROCESSED ( ch_processed_reads )
        ch_versions = ch_versions.mix( FASTQC_PROCESSED.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )
    }

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc  = ch_multiqc_files
}

