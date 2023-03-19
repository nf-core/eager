//
// Index input reference as required
//

include { GUNZIP                          } from "../../modules/nf-core/gunzip/main"
include { BWA_INDEX                       } from "../../modules/nf-core/bwa/index/main"
include { BOWTIE2_BUILD                   } from "../../modules/nf-core/bowtie2/build/main"
include { SAMTOOLS_FAIDX                  } from "../../modules/nf-core/samtools/faidx/main"
include { PICARD_CREATESEQUENCEDICTIONARY } from "../../modules/nf-core/picard/createsequencedictionary/main"
// TODO missing: circulargeneraotr?

workflow REFERENCE_INDEXING_MULTI {

    take:
    referencesheet // file: /path/to/name.{csv,tsv}

    main:
    ch_versions = Channel.empty()

    // TODO geneal: try with fai/dict files also in referneces.csv
    // TODO versions!
    // TODO FINISH!


    // Parse CSV and detect files to load
    if ( referencesheet.extension == "tsv" ){
        ch_splitreferencesheet_for_branch = Channel.fromPath(referencesheet)
                                                .splitCsv ( header:true, sep:"\t" )
                                                .map {
                                                    row ->
                                                        def meta            = [:]
                                                        meta.id             = row["reference_name"]
                                                        def fasta           = file(row["fasta"], checkIfExists: true) // mandatory parameter!
                                                        def fai             = row["fai"] != "" ? file(row["fai"], checkIfExists: true) : ""
                                                        def dict            = row["dict"] != "" ? file(row["dict"], checkIfExists: true) : ""
                                                        def mapper_index    = row["mapper_index"] != "" ? file(row["mapper_index"], checkIfExists: true) : ""
                                                        def circular_target = row["circular_target"] // TODO Additional tests?
                                                        def mitochondrion   = row["mitochondrion"] // TODO Additional tests?
                                                        [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                                                }
        // TODO check same output as CSV (csv main testing)
    } else {
        ch_splitreferencesheet_for_branch = Channel.fromPath(referencesheet)
                                                .splitCsv ( header:true )
                                                .map {
                                                    row ->
                                                        def meta            = [:]
                                                        meta.id             = row["reference_name"]
                                                        def fasta           = file(row["fasta"], checkIfExists: true) // mandatory parameter!
                                                        def fai             = row["fai"] != "" ? file(row["fai"], checkIfExists: true) : ""
                                                        def dict            = row["dict"] != "" ? file(row["dict"], checkIfExists: true) : ""
                                                        def mapper_index    = row["mapper_index"] != "" ? file(row["mapper_index"], checkIfExists: true) : ""
                                                        def circular_target = row["circular_target"]
                                                        def mitochondrion   = row["mitochondrion"]
                                                        [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                                                }
    }

    // GENERAL DESCRIPTION FOR NEXT SECTIONS
    // This will be the same scheme for all other generation steps, i.e.
    // for those that need to be processed send them to a  multiMap (others skip via branch) 
    // take only those files needed for the generation step in one sub-channel, and all other
    // channel elements go into 'remainder'. After generation, join these back, and then mix back the
    // ones that don't skip

    //
    // DECOMPRESSION
    //

    // Detect if fasta is gzipped or not
    ch_fasta_for_gunzip = ch_splitreferencesheet_for_branch
                            .branch {
                                meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                                    forgunzip: fasta.extension == "gz"
                                    skip: true
                            }

    // Pull out name/file to match cardinality for GUNZIP module
    ch_gunzip_input = ch_fasta_for_gunzip
        .forgunzip
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                gunzip:    [ meta, fasta ]
                remainder: [ meta, fai, dict, mapper_index, circular_target, mitochondrion ]
        }


    GUNZIP ( ch_gunzip_input.gunzip )

    // Mix back gunzipped fasta with remaining files, and then mix back with pre-gunzipped references
    ch_gunzippedfasta_formix = GUNZIP.out.gunzip.join( ch_gunzip_input.remainder )
    ch_fasta_for_faiindexing = ch_fasta_for_gunzip.skip.mix(ch_gunzippedfasta_formix)

    //
    // INDEXING: fai
    //

    // Separate out non-faidxed references
    ch_fasta_for_faidx = ch_fasta_for_faiindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                forfaidx: fai == ""
                skip: true
        }

    // Split channel to ensure cardindality matching
    ch_faidx_input = ch_fasta_for_faidx
        .forfaidx
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                faidx:      [ meta, fasta ]
                remainder:  [ meta, fasta, dict, mapper_index, circular_target, mitochondrion ] // we drop fai here as we are going to make it
        }

    SAMTOOLS_FAIDX ( ch_faidx_input.faidx )

    // Rejoin output channel with main reference indicies channel elements
    // TODO this FAIDX was producing a nested ID for some reason, should work out why:  [['id':['id':'mammoth']], so we can drop the first map
    ch_faidxed_formix =  SAMTOOLS_FAIDX.out.fai
                            .join( ch_faidx_input.remainder )
                            .map {
                                meta, fai, fasta, dict, mapper_index, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }

    // Mix back newly faidx'd references with the pre-indexed ones
    ch_fasta_for_dictindexing = ch_fasta_for_faidx.skip.mix(ch_faidxed_formix)

    //
    // INDEXING: dict
    //

    ch_fasta_for_dict = ch_fasta_for_dictindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                fordict: dict == ""
                skip: true
        }

    ch_dict_input = ch_fasta_for_dict
        .fordict
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                dict:      [ meta, fasta ]
                remainder:  [ meta, fasta, fai, mapper_index, circular_target, mitochondrion ]
        }

    PICARD_CREATESEQUENCEDICTIONARY ( ch_dict_input.dict )

    ch_dicted_formix =  PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
                            .join( ch_dict_input.remainder )
                            .map {
                                meta, dict, fasta, fai, mapper_index, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }

    ch_dict_formapperindexing = ch_fasta_for_dict.skip.mix(ch_dicted_formix)

    //
    // INDEXING: Mapping indicies
    //

    // Generate mapper indicies if not supplied, and if supplied generate meta
    if ( params.mapping_tool == "bwaaln" ){

        ch_fasta_for_bwaindex = ch_dict_formapperindexing
            .branch {
                meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                    forindex: mapper_index == ""
                    skip: true
            }

        ch_mapindex_input = ch_fasta_for_bwaindex
            .forindex
            .multiMap {
                meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                    index:      [ meta, fasta ]
                    remainder:  [ meta, fasta, fai, dict, circular_target, mitochondrion ]
            }

        BWA_INDEX ( ch_mapindex_input.index.dump(tag: "looksOK?") )

        ch_indexed_formix = BWA_INDEX.out.index
                                .join( ch_mapindex_input.remainder )
                                .dump(tag: "post_join")
                                .map {
                                    meta, mapper_index, fasta, fai, dict, circular_target, mitochondrion ->

                                    [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                                }
                                .dump(tag: "post_remap")
        ch_indexmapper_for_reference = ch_fasta_for_bwaindex.skip.mix(ch_indexed_formix)

    } else if ( params.mapping_tool == "bowtie2" ) {

        println("Not yet implemented")
    }
    // TODO: document that the "base name" of all indicies must be the same, i.e. correspond to the FASTA

    // Join all together into a single map. Include failOnMismatch as a check if
    // a user supplies indicies with different "base" names.
    // ch_reference_for_mapping = ch_ungz_ref.join(ch_fasta_fai, failOnMismatch: true).join(ch_fasta_dict, failOnMismatch: true).join(ch_fasta_mapperindexdir, failOnMismatch: true)

    emit:
    reference = ch_indexmapper_for_reference //ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex, <etc.> ]
    versions = ch_versions
}

// def create_index_channel(LinkedHashMap row) {
//     def meta = [:]
//     // TODO create spanning main metadata
//     meta.reference_id = row.reference_name
//     meta.file_type    =

//     // MOVE FILE TYPE TO META AND TRANSPOSE TO SINGLE FILES?


//     def array = []
//     if (!file(row.bam).exists()) {
//         exit 1, "[nf-core/eager] error: Please check input referencesheet. BAM file does not exist!\nFile: ${row.bam}"
//     } else {
//         array = [ meta, file(row.bam) ]
//     }
//     return array
// }
