//
// Elongate a reference genome by circularising the target sequence by a given elongation factor.
//

include { GUNZIP as GUNZIP_ELONGATED_FASTA    } from '../../modules/nf-core/gunzip/main'
include { CIRCULARMAPPER_CIRCULARGENERATOR    } from '../../modules/nf-core/circularmapper/circulargenerator/main'
include { BWA_INDEX as BWA_INDEX_CIRCULARISED } from '../../modules/nf-core/bwa/index/main'
include { addNewMetaFromAttributes            } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

workflow ELONGATE_REFERENCE {
    take:
    ch_reference            // [ meta, fasta, fai, dict, mapindex ]
    ch_elongated_reference  // [ meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ]

    main:
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()
    ch_circular_reference = Channel.empty()
    ch_elongated_unzipped = Channel.empty()
    ch_elongated_chr      = Channel.empty()

    // Check if the provided elongated reference is gzipped, and if so, unzip it.
    ch_elongated_branches = ch_elongated_reference
                            .branch {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                    for_gunzip: circularmapper_elongated_fasta != '' && circularmapper_elongated_fasta.extension == "gz"
                                    skip_gunzip: true
                            }

    ch_elongated_for_gunzip = ch_elongated_branches.for_gunzip
                            .map {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                    [ meta, circularmapper_elongated_fasta ]
                            }

    GUNZIP_ELONGATED_FASTA( ch_elongated_for_gunzip )
    ch_versions = ch_versions.mix( GUNZIP_ELONGATED_FASTA.out.versions.first() )

    ch_elongated_unzipped = ch_elongated_reference
                            .join( GUNZIP_ELONGATED_FASTA.out.gunzip )
                            .map {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index, unzipped_fasta ->
                                    def final_fasta = unzipped_fasta ?: circularmapper_elongated_fasta
                                    [ meta, circular_target, unzipped_fasta, circularmapper_elongated_index ]
                            }
                            .mix( ch_elongated_branches.skip_gunzip )

    /*
        Check what fasta files we have.
        There are four options:
            1. Elongated reference with index (ignore circular target) -> Pass through
            2. Elongated reference without index (ignore circular target) -> Index and emit
            3. No elongated reference, but circular target -> Elongate, index and emit.
            4. None of the above -> Throw error and stop execution during parameter validation
    */

    ch_circulargenerator_input = ch_elongated_unzipped
                            .branch{
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                        ready:            circularmapper_elongated_fasta != "" && circularmapper_elongated_index != ""
                                        needs_index:      circularmapper_elongated_fasta != "" && circularmapper_elongated_index == ""
                                        needs_elongation: circularmapper_elongated_fasta == "" && circular_target != ""
                            }

    /* References that are already elongated, need ch_elongated_chr to be created from the circular target information
        1) Get the reference information ready for joinin with the new channel.
        2) Take all subchannels from the multiMap that do not go through CircularGenerator (.ready,.needs_index) and infer the name of the elongated_chr_list expected by RealignSAMFile
        3) Put the circular target in a file of that name FOR EACH REFERENCE. The resulting channel has no meta, so we need to add it.
        4) Add meta, and use to merge back to the reference channel. This way we can take the original reference's meta.

        This is a bit convoluted, but it should work. Would be simpler if I could create the meta within collectFile, but I did not find a way to do that.
    */
    ch_ref_for_chr_list = ch_reference
                            .map {
                                addNewMetaFromAttributes( it, "id", "id", false )
                            }

    ch_chr_list_for_already_elongated_ref = ch_circulargenerator_input.ready
                            .mix( ch_circulargenerator_input.needs_index )
                            .join( ch_reference )
                            .map {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index, fasta, fai, dict, mapindex ->
                                    [ meta, fasta, circular_target ]
                            }
                            .collectFile {
                                meta, fasta, circular_target ->
                                    [ "${fasta.name}_500_elongated", circular_target + '\n' ]
                            }
                            .map {
                                file ->
                                    def id = file.getSimpleName()
                                    [ [id: id ], file ]
                            }
                            .join(ch_ref_for_chr_list)
                            .map {
                                ignore_me, chr_list, meta, fasta, fai, dict, mapindex ->
                                    [ meta, chr_list ]
                            }
                            .dump(tag: "collected_files", pretty:true)

    // Elongate references that need it
    // Join the original references to the branch of needs_elongation, to get the original fasta files, and elongate them.
    ch_references_to_elongate = ch_circulargenerator_input.needs_elongation
                            .join( ch_reference )
                            .multiMap {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index, fasta, fai, dict, mapindex ->

                                    def elongation_factor = params.fasta_circularmapper_elongationfactor

                                    fasta:              [ meta, fasta ]
                                    elongation_factor : [ meta, elongation_factor ]
                                    target:             [ meta, circular_target ]
                            }

    CIRCULARMAPPER_CIRCULARGENERATOR(
        ch_references_to_elongate.fasta,
        ch_references_to_elongate.elongation_factor,
        ch_references_to_elongate.target
    )
    ch_elongated_chr = ch_chr_list_for_already_elongated_ref.mix(CIRCULARMAPPER_CIRCULARGENERATOR.out.elongated)
    ch_versions      = ch_versions.mix( CIRCULARMAPPER_CIRCULARGENERATOR.out.versions.first() )

    // Collect newly generated circular references and provided ones without an index, and index them.
    ch_input_for_circular_indexing = ch_circulargenerator_input.needs_index
                            .map {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                    [ meta, circularmapper_elongated_fasta ]
                            }
                            .mix( CIRCULARMAPPER_CIRCULARGENERATOR.out.fasta )

    BWA_INDEX_CIRCULARISED(ch_input_for_circular_indexing)
    ch_versions = ch_versions.mix( BWA_INDEX_CIRCULARISED.out.versions.first() )

    ch_indexed_references = ch_input_for_circular_indexing
                            .join( BWA_INDEX_CIRCULARISED.out.index )

    // Then put all the indexed elongated references together, replace any zipped ones with the unzipped version, and emit them
    ch_circular_reference = ch_circulargenerator_input.ready
                            .map {
                                meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                    [ meta, circularmapper_elongated_fasta, circularmapper_elongated_index ]
                            }
                            .mix( ch_indexed_references )

    emit:
    circular_reference = ch_circular_reference // [ meta, fasta, index ]
    elongated_chr_list = ch_elongated_chr      // [ meta, elongated_chr_list ]
    versions           = ch_versions
    mqc                = ch_multiqc_files
}
