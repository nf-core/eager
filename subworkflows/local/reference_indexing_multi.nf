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
    samplesheet // file: /path/to/name.{csv,tsv}

    main:
    ch_versions = Channel.empty()

    // TODO geneal: try with fai/dict files also in referneces.csv
    // TODO versions!

    // Parse CSV and detect files to load
    if ( samplesheet.extension == "tsv" ){
        ch_splitsamplesheet_for_branch = samplesheet
                                            .splitCsv ( header:true, sep:"\t" )
        // TODO copy over from csv
    } else {
        ch_splitsamplesheet_for_branch = Channel.fromPath(samplesheet)
                                            .splitCsv ( header:true )
                                            .dump(tag: "splitcsv")
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

    //
    // DECOMPRESSION
    //

    // Detect if fasta is gzipped or not
    ch_fasta_for_gunzip = ch_splitsamplesheet_for_branch
                            .dump(tag: "ch_splitsamplesheet_for_branch")
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
    ch_gunzippedfasta_formix = GUNZIP.out.gunzip.join( ch_gunzip_input.remainder ).dump(tag: "ch_gunzippedfasta_formix")
    ch_fasta_for_faiindexing = ch_fasta_for_gunzip.skip.mix(ch_gunzippedfasta_formix).dump(tag: "prepped")

    //
    // INDEXING: fai
    //

    ch_fasta_for_faidx = ch_fasta_for_faiindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                forfaidx: fai == ""
                skip: true
        }

    ch_faidx_input = ch_fasta_for_faidx
        .forfaidx
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                faidx:      [ meta, fasta ]
                remainder:  [ meta, fasta, dict, mapper_index, circular_target, mitochondrion ] // we drop fai here as we are going to make it
        }

    SAMTOOLS_FAIDX ( ch_faidx_input.faidx )

    // TODO this FAIDX was producing a nested ID for some reason, should work out why:  [['id':['id':'mammoth']], so we can drop the first map
    ch_faidxed_formix =  SAMTOOLS_FAIDX.out.fai
                            .map {
                                meta, fai ->
                                [ meta['id'], fai ]
                            }
                            .dump(tag: "ch_faidx_input_remainder")
                            .join( ch_faidx_input.remainder)
                            .dump(tag: "ch_faidx_input_remainder_join")
                            .map {
                                meta, fai, fasta, dict, mapper_index, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }
                            .dump(tag: "ch_faidx_input_remainder_map")


    ch_faidxed_formix.dump(tag: "ch_faidxed_formix")
    ch_fasta_for_dictindexing = ch_fasta_for_faidx.skip.mix(ch_faidxed_formix).dump(tag: "ch_fasta_for_faidx_skip_mix")

    //
    // INDEXING: dict
    //

    ch_fasta_for_dict = ch_fasta_for_dictindexing
        .dump(tag: "ch_fasta_for_dictindexing")
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                fordict: dict == ""
                skip: true
        }

    ch_dict_input = ch_fasta_for_dict
        .fordict
        .dump(tag: "ch_fasta_for_dict_fordict")
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                dict:      [ meta, fasta ]
                remainder:  [ meta, fasta, fai, mapper_index, circular_target, mitochondrion ]
        }

    PICARD_CREATESEQUENCEDICTIONARY ( ch_dict_input.dict )

    ch_dicted_formix =  PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
                            .dump(tag: "picard_createsequencedictionary_out_reference_dict")
                            .join( ch_dict_input.remainder.dump(tag: "ch_dict_input_remainder") )
                            .dump(tag: "ch_dict_input_remainder_join")
                            .map {
                                meta, dict, fasta, fai, mapper_index, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }
                            .dump(tag: "ch_dict_input_remainder_map")

    ch_dict_formapperindexing = ch_fasta_for_dict.skip.dump(tag: "ch_fasta_for_dict_skip").mix(ch_dicted_formix).dump(tag: "ch_fasta_for_dict_skip_mix")


    // // Generate mapper indicies if not supplied, and if supplied generate meta
    // if ( params.mapping_tool == "bwaaln" ){

    //     if ( !fasta_mapperindexdir ) {
    //         ch_fasta_mapperindexdir = BWA_INDEX ( ch_ungz_ref ).index
    //     } else {
    //         ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
    //     }

    // } else if ( params.mapping_tool == "bowtie2" ) {

    //     if ( !fasta_mapperindexdir ) {
    //         ch_fasta_mapperindexdir = BOWTIE2_BUILD ( ch_ungz_ref ).index
    //     } else {
    //         ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
    //     }

    // }


    // TODO: document that the "base name" of all indicies must be the same, i.e. correspond to the FASTA

    // Join all together into a single map. Include failOnMismatch as a check if
    // a user supplies indicies with different "base" names.
    // ch_reference_for_mapping = ch_ungz_ref.join(ch_fasta_fai, failOnMismatch: true).join(ch_fasta_dict, failOnMismatch: true).join(ch_fasta_mapperindexdir, failOnMismatch: true)

    emit:
    reference = Channel.empty() //ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex ]
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
//         exit 1, "[nf-core/eager] error: Please check input samplesheet. BAM file does not exist!\nFile: ${row.bam}"
//     } else {
//         array = [ meta, file(row.bam) ]
//     }
//     return array
// }
