//
// Prepare reference indexing for downstream
//

include { FASTQ_ALIGN_BWAALN } from '../../subworkflows/nf-core/fastq_align_bwaaln/main'
include { SAMTOOLS_FLAGSTAT  } from '../../modules/nf-core/samtools/flagstat/main'

workflow MAP {
    take:
    reads // [ [meta], [read1, reads2] ] or [ [meta], [read1] ]
    index // [ [meta], [ index ] ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc  = Channel.empty()

    if ( params.mapping_tool == 'bwaaaln' ) {
        ch_mapped_bam = FASTQ_ALIGN_BWAALN ( reads, index )
    }

    SAMTOOLS_FLAGSTAT ( ch_mapped_bam )

    emit:
    bam   = ch_mapped_bam                         // [ [ meta ], bam ]
    bai   = ch_mapped_bai                         // [ [ meta ], bai ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat     // [ [ meta ], stats ]
    versions  = ch_versions

}
