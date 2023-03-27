//
// Alignment with BWA-MEM, sort and index
//

include { BWA_MEM            } from '../../modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX     } from '../../modules/nf-core/samtools/index/main'

workflow FASTQ_ALIGN_BWAMEM {

    take:
    ch_reads     // channel (mandatory): [ val(meta), [ reads ] ]
    ch_index     // channel (mandatory): [ val(meta), [ index ] ]
    val_sort_bam // boolean (mandatory): set to true (samtools sort) to be able to index afterwards

    main:
    ch_versions = Channel.empty()

    BWA_MEM ( ch_reads, ch_index, val_sort_bam )
    ch_versions = ch_versions.mix( BWA_MEM.out.versions.first() )

    ch_bam_for_index = BWA_MEM.out.bam

    SAMTOOLS_INDEX ( ch_bam_for_index )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam      = ch_bam_for_index           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai     // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi     // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}
