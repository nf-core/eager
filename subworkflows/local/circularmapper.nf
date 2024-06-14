//
// Run circularmapper
//

include { CIRCULARMAPPER_CIRCULARGENERATOR      } from '../../modules/nf-core/circularmapper/circulargenerator/main'
include { CIRCULARMAPPER_REALIGNSAMFILE         } from '../../modules/nf-core/circularmapper/realignsamfile/main'
include { BWA_ALN as BWA_ALN_CIRCULARMAPPER     } from '../../modules/nf-core/bwa/aln/main'
include { BWA_INDEX as BWA_INDEX_CIRCULARMAPPER } from '../../modules/nf-core/bwa/index/main'
include { BWA_SAMSE as BWA_SAMSE_CIRCULARMAPPER } from '../../../modules/nf-core/bwa/samse/main'

workflow CIRCULARMAPPER {

    take:
    fasta_reference // channel (mandatory): [ val(meta), path(reference) ]
    eval // channel (mandatory): val(elongation value)
    fastq_reads // channel (mandatory): [ val(meta), path(reads) ]

    main:

    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.run_circularmapper ) {

        ch_reference = fasta_reference
        ch_eval =  eval

        CIRCULARMAPPER_CIRCULARGENERATOR(ch_reference, ch_eval)
        ch_versions = ch_versions.mix( CIRCULARMAPPER_CIRCULARGENERATOR.out.versions.first() )

        BWA_INDEX_CIRCULARMAPPER(CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta)
        ch_versions = ch_versions.mix( BWA_INDEX_CIRCULARMAPPER.out.versions.first() )

        ch_input_bwa_aln = fastq_reads
            .map {
                // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .groupTuple(by:0)
            .combine( BWA_INDEX_CIRCULARMAPPER.out.index, by: 0 ) // [ [meta], fastq, bai ]
            .multiMap {
                combo_meta, metas, fastq, ref_bai, bai ->
                def ids = metas.collect { meta -> meta.id }
                fastqs: [ combo_meta + [id: ids], fastq ]
                bai:  [ ref_bai, bai ]
            }

        BWA_ALN_CIRCULARMAPPER(ch_input_bwa_aln)
        ch_versions = ch_versions.mix( BWA_ALN_CIRCULARMAPPER.out.versions.first() )

        ch_input_bwa_samse = ch_input_bwa_aln
        .combine(  BWA_ALN_CIRCULARMAPPER.out.sai, by: 0 ) // [ [meta], fastq, bai, sai ]
        .multiMap {
                metas, fastq, ref_bai, bai, ref_sai, sai ->
                fastqs: [ metas, fastq, sai ]
                bai:  [ ref_bai, bai ]
            }

        BWA_SAMSE_CIRCULARMAPPER(ch_input_bwa_samse)
        ch_versions = ch_versions.mix( BWA_SAMSE_CIRCULARMAPPER.out.versions.first() )

        ch_input_realignsamfile = BWA_SAMSE_CIRCULARMAPPER.out.bam
        .combine(CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta, by: 0)
        .combine(ch_eval)
        .multiMap {
            ref_bam, bam, ref_fasta, fasta, ch_eval ->
            bam:   [ ref_bam, bam ]
            fasta: [ ref_fasta, fasta ]
            eval:  [ ch_eval ]
        }

        CIRCULARMAPPER_REALIGNSAMFILE(ch_input_realignsamfile)
        ch_versions = ch_versions.mix( CIRCULARMAPPER_REALIGNSAMFILE.out.versions.first() )

        emit:

        bam      =   CIRCULARMAPPER_REALIGNSAMFILE.out.bam  // channel: [ val(meta), path(bam) ]
        versions = ch_versions         // channel: [ path(versions.yml) ]

    }
}
