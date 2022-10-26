//
// Perform read trimming and merging
//


include { PREPROCESSING_FASTP          } from './preprocessing_fastp'
include { PREPROCESSING_ADAPTERREMOVAL } from './preprocessing_adapterremoval'
//include { LEEHOM                     } from './preprocessing_leehom'
include { FASTQC as FASTQC_PROCESSED } from '../../modules/nf-core/fastqc/main'

workflow PREPROCESSING {
    take:
    reads //  [ [ meta ], [ reads ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.preprocessing_tool == "fastp" ) {
        ch_processed_reads = PREPROCESSING_FASTP ( reads ).reads
        ch_versions        =  ch_versions.mix( PREPROCESSING_FASTP.out.versions )
        ch_multiqc_files   =  ch_multiqc_files.mix( PREPROCESSING_FASTP.out.mqc )
    } else if ( params.preprocessing_tool == "adapterremoval" ) {
        ch_processed_reads = PREPROCESSING_ADAPTERREMOVAL ( reads ).reads
        ch_versions        = ch_versions.mix( PREPROCESSING_ADAPTERREMOVAL.out.versions )
        ch_multiqc_files   = ch_multiqc_files.mix( PREPROCESSING_ADAPTERREMOVAL.out.mqc )
    } else {
        ch_processed_reads = reads
    }

    FASTQC_PROCESSED ( ch_processed_reads )
    ch_versions = ch_versions.mix( FASTQC_PROCESSED.out.versions )
    ch_multiqc_files = ch_multiqc_files.mix( FASTQC_PROCESSED.out.zip )

    emit:
    reads    = ch_processed_reads   // channel: [ val(meta), [ reads ] ]
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}

