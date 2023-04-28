//
// Remove mapping reads to the host from the fastq
//

include { HOST_REMOVAL } from '../../modules/local/host_removal'

workflow HOST_REMOVAL {

    take:
    bam      // [ [ meta ], bam ]
    ch_fastq // [ [ meta ], [fastq_r1, fastq_r2 ] ] from INPUT_CHECK.out.fastqs?


    main:
    ch_versions           = Channel.empty()

    ch_for_hostremoval = bam.dump(tag: "bam").join( ch_fastq.dump(tag: "fastq") )

    HOST_REMOVAL ( ch_for_hostremoval )

    ch_versions = ch_versions.mix( HOST_REMOVAL.out.versions )

    emit:
    versions = ch_versions

}
