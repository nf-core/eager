//
// Index input reference as required
//

include { GUNZIP as GUNZIP_FASTA          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
// TODO missing: circulargeneraotr?

workflow REFERENCE_INDEXING_MULTI {

    take:
    referencesheet // file: /path/to/name.{csv,tsv}

    main:
    ch_versions = Channel.empty()

    // Parse CSV and detect files to load
    if ( referencesheet.extension == "tsv" ){
        ch_splitreferencesheet_for_map = Channel.fromPath(referencesheet)
                                                .splitCsv ( header:true, sep:"\t" )
    } else {
        ch_splitreferencesheet_for_map = Channel.fromPath(referencesheet)
                                                .splitCsv ( header:true )
    }

    ch_splitreferencesheet_for_branch = ch_splitreferencesheet_for_map
                                            .map {
                                                row ->
                                                    def meta            = [:]
                                                    meta.id             = row["reference_name"]
                                                    def fasta           = file(row["fasta"], checkIfExists: true) // mandatory parameter!
                                                    def fai             = row["fai"] != "" ? file(row["fai"], checkIfExists: true) : ""
                                                    def dict            = row["dict"] != "" ? file(row["dict"], checkIfExists: true) : ""
                                                    def mapper_index    = row["mapper_index"] != "" ? file(row["mapper_index"], checkIfExists: true) : ""
                                                    def circular_target = row["circular_target"]
                                                    def mitochondrion   = row["mitochondrion_header"]
                                                    [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
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
    ch_gunzip_input = ch_fasta_for_gunzip.forgunzip
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                gunzip:    [ meta, fasta ]
                remainder: [ meta, fai, dict, mapper_index, circular_target, mitochondrion ]
        }


    GUNZIP_FASTA ( ch_gunzip_input.gunzip )
    ch_version = ch_versions.mix( GUNZIP_FASTA.out.versions )

    // Mix back gunzipped fasta with remaining files, and then mix back with pre-gunzipped references
    ch_gunzippedfasta_formix = GUNZIP_FASTA.out.gunzip.join( ch_gunzip_input.remainder, failOnMismatch: true )
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
    ch_version = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    // Rejoin output channel with main reference indicies channel elements
    ch_faidxed_formix =  SAMTOOLS_FAIDX.out.fai
                            .join( ch_faidx_input.remainder, failOnMismatch: true )
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
    ch_version = ch_versions.mix( PICARD_CREATESEQUENCEDICTIONARY.out.versions )

    ch_dicted_formix =  PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
                            .join( ch_dict_input.remainder, failOnMismatch: true )
                            .map {
                                meta, dict, fasta, fai, mapper_index, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }

    ch_dict_formapperindexing = ch_fasta_for_dict.skip.mix(ch_dicted_formix)

    //
    // INDEXING: Mapping indicies
    //

    // Generate mapper indicies if not supplied, and if supplied generate meta

    ch_fasta_for_mapperindex = ch_dict_formapperindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                forindex: mapper_index == ""
                skip: true
        }

    ch_mapindex_input = ch_fasta_for_mapperindex
        .forindex
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ->
                index:      [ meta, fasta ]
                remainder:  [ meta, fasta, fai, dict, circular_target, mitochondrion ]
        }

    if ( params.mapping_tool == "bwaaln" || params.mapping_tool == "bwamem" ) {
        BWA_INDEX ( ch_mapindex_input.index )
        ch_version = ch_versions.mix( BWA_INDEX.out.versions )
        ch_indexed_forremap = BWA_INDEX.out.index
    } else if ( params.mapping_tool == "bowtie2" ) {
        BOWTIE2_BUILD ( ch_mapindex_input.index )
        ch_version = ch_versions.mix( BOWTIE2_BUILD.out.versions )
        ch_indexed_forremap = BOWTIE2_BUILD.out.index
    } else if ( params.mapping_tool == "circularmapper" ) {
        println("CircularMapper Indexing Not Yet Implemented")
    }

    ch_indexed_formix = ch_indexed_forremap
                            .join( ch_mapindex_input.remainder, failOnMismatch: true )
                            .map {
                                meta, mapper_index, fasta, fai, dict, circular_target, mitochondrion ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion ]
                            }

    ch_indexmapper_for_reference = ch_fasta_for_mapperindex.skip.mix(ch_indexed_formix)

    emit:
    reference = ch_indexmapper_for_reference // [ meta, fasta, fai, dict, mapindex, circular_target, mitochondrion ]
    versions = ch_versions
}
