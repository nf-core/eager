//
// Elongate a reference genome by circularising the target sequence by a given elongation factor.
//

include { CIRCULARMAPPER_CIRCULARGENERATOR      } from '../../modules/nf-core/circularmapper/circulargenerator/main'
include { BWA_INDEX as BWA_INDEX_CIRCULARISED } from '../../modules/nf-core/bwa/index/main'

workflow ELONGATE_REFERENCE {
    take:
    ch_reference            // [ meta, fasta, fai ]
    ch_elongated_reference  // [ meta, elongated_fasta, elongated_fai ]
    elongation_factor        // [ int ]
    // TODO CIRCULARMAPPER_CIRCULARGENERATOR module needs updating. `-s` option is the circular target and not the output file >.<

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    /*
        Check what fasta files we have:
        There are four options:
            1. Elongated reference with index (ignore circular target)
            2. Elongated reference without index (ignore circular target)
            3. No elongated reference, but circular target
            4. None of the above -> Throw error (should go in parameter validation)
    */
    ch_circulargenerator_input = ch_elongated_reference
                            .branch{
                                meta, elongated_fasta_index, elongated_fasta, circular_target ->
                                        ready:            elongated_fasta != "" && elongated_fasta_index != ""
                                        needs_index:      elongated_fasta != "" && elongated_fasta_index == ""
                                        needs_elongation: elongated_fasta == "" && circular_target != ""
                            }

    // Elongate references that need it
    // Join the original references to the branch of needs_elongation, to get the original fasta files, and elongate them.
    ch_references_to_elongate = ch_circulargenerator_input.needs_elongation
                            .join( ch_reference )
                            .map {
                                meta, elongated_fasta_index, elongated_fasta, circular_target, meta2, index, fasta ->
                                    [ meta, fasta ]
                            }

    CIRCULARMAPPER_CIRCULARGENERATOR(ch_circulargenerator_input.needs_elongation, elongation_value)
    ch_versions = ch_versions.mix( CIRCULARMAPPER_CIRCULARGENERATOR.out.versions.first() )

    // Collect newly generated circular references and provided ones without an index, and index them.
    ch_input_for_circular_indexing = ch_circulargenerator_input.needs_index
                            .map {
                                meta, elongated_fasta_index, elongated_fasta, circular_target ->
                                    [ meta, elongated_fasta ]
                            }
                            .mix( CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta )

    BWA_INDEX_CIRCULARISED(ch_input_for_circular_indexing)
    ch_versions = ch_versions.mix( BWA_INDEX_CIRCULARISED.out.versions.first() )

    ch_indexed_references = ch_input_for_circular_indexing
                            .join( BWA_INDEX_CIRCULARISED.out.index )

    // Then put all the indexed elongated references together and emit them
    ch_circular_reference = ch_circulargenerator_input.ready
                            .map {
                                meta, elongated_fasta_index, elongated_fasta, circular_target ->
                                    [ meta, elongated_fasta, elongated_fasta_index ]
                            }
                            .mix( ch_indexed_references )

    emit:
    circular_reference = ch_circular_reference // [ meta, fasta, fai ]
    versions           = ch_versions
    mqc                = ch_multiqc_files

}
