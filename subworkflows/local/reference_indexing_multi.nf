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
                                                    def meta                   = [:]
                                                    meta.id                    = row["reference_name"]
                                                    def fasta                  = file(row["fasta"], checkIfExists: true) // mandatory parameter!
                                                    def fai                    = row["fai"] != "" ? file(row["fai"], checkIfExists: true) : ""
                                                    def dict                   = row["dict"] != "" ? file(row["dict"], checkIfExists: true) : ""
                                                    def mapper_index           = row["mapper_index"] != "" ? file(row["mapper_index"], checkIfExists: true) : ""
                                                    def circular_target        = row["circular_target"]
                                                    def mitochondrion          = row["mitochondrion_header"]
                                                    def capture_bed            = row["snpcapture_bed"] != "" ? file(row["snpcapture_bed"], checkIfExists: true) : ""
                                                    def pileupcaller_bed       = row["pileupcaller_bedfile"] != "" ? file(row["pileupcaller_bedfile"], checkIfExists: true) : ""
                                                    def pileupcaller_snp       = row["pileupcaller_snpfile"] != "" ? file(row["pileupcaller_snpfile"], checkIfExists: true) : ""
                                                    def hapmap                 = row["hapmap_file"] != "" ? file(row["hapmap_file"], checkIfExists: true) : ""
                                                    def pmd_mask               = row["pmdtools_masked_fasta"] != "" ? file(row["pmdtools_masked_fasta"], checkIfExists: true) : ""
                                                    def sexdet_bed             = row["sexdeterrmine_snp_bed"] != "" ? file(row["sexdeterrmine_snp_bed"], checkIfExists: true) : ""
                                                    def bedtools_feature       = row["bedtools_feature_file"] != "" ? file(row["bedtools_feature_file"], checkIfExists: true) : ""
                                                    def genotyping_gatk_ploidy = row["genotyping_gatk_ploidy"] != "" ? row["genotyping_gatk_ploidy"] : ""
                                                    def genotyping_gatk_dbsnp  = row["genotyping_gatk_dbsnp"] != "" ? file(row["genotyping_gatk_dbsnp"], checkIfExists: true) : ""
                                                    [ meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion, capture_bed, pileupcaller_bed, pileupcaller_snp, hapmap, pmd_mask, sexdet_bed, bedtools_feature, genotyping_gatk_ploidy, genotyping_gatk_dbsnp ]
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

ch_input_from_referencesheet = ch_splitreferencesheet_for_branch
                            .multiMap {
                                meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion, capture_bed, pileupcaller_bed, pileupcaller_snp, hapmap, pmd_mask, sexdet_bed, bedtools_feature, genotyping_gatk_ploidy, genotyping_gatk_dbsnp ->
                                generated:            [ meta, fasta, fai, dict, mapper_index, circular_target ]
                                mitochondrion_header: [ meta, mitochondrion ]
                                angsd_hapmap:         [ meta, hapmap ]
                                pmd_mask:             [ meta, pmd_mask, capture_bed ]
                                snp_bed:              [ meta, capture_bed ]
                                pileupcaller_snp:     [ meta, pileupcaller_bed, pileupcaller_snp ]
                                sexdeterrmine_bed:    [ meta, sexdet_bed ]
                                bedtools_feature:     [ meta, bedtools_feature ]
                                dbsnp:                [ meta + [ ploidy: genotyping_gatk_ploidy ], genotyping_gatk_dbsnp ] // Include ploidy of the reference in dbsnp meta.
                            }

    // Detect if fasta is gzipped or not
    ch_fasta_for_gunzip = ch_input_from_referencesheet.generated
                            .branch {
                                meta, fasta, fai, dict, mapper_index, circular_target ->
                                    forgunzip: fasta.extension == "gz"
                                    skip: true
                            }

    // Pull out name/file to match cardinality for GUNZIP module
    ch_gunzip_input = ch_fasta_for_gunzip.forgunzip
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                gunzip:    [ meta, fasta ]
                remainder: [ meta, fai, dict, mapper_index, circular_target ]
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
            meta, fasta, fai, dict, mapper_index, circular_target ->
                forfaidx: fai == ""
                skip: true
        }

    // Split channel to ensure cardinality matching
    ch_faidx_input = ch_fasta_for_faidx
        .forfaidx
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                faidx:      [ meta, fasta ]
                remainder:  [ meta, fasta, dict, mapper_index, circular_target ] // we drop fai here as we are going to make it
        }

    SAMTOOLS_FAIDX ( ch_faidx_input.faidx )
    ch_version = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )

    // Rejoin output channel with main reference indicies channel elements
    ch_faidxed_formix =  SAMTOOLS_FAIDX.out.fai
                            .join( ch_faidx_input.remainder, failOnMismatch: true )
                            .map {
                                meta, fai, fasta, dict, mapper_index, circular_target ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target ]
                            }

    // Mix back newly faidx'd references with the pre-indexed ones
    ch_fasta_for_dictindexing = ch_fasta_for_faidx.skip.mix(ch_faidxed_formix)

    //
    // INDEXING: dict
    //

    ch_fasta_for_dict = ch_fasta_for_dictindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                fordict: dict == ""
                skip: true
        }

    ch_dict_input = ch_fasta_for_dict
        .fordict
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                dict:      [ meta, fasta ]
                remainder:  [ meta, fasta, fai, mapper_index, circular_target ]
        }

    PICARD_CREATESEQUENCEDICTIONARY ( ch_dict_input.dict )
    ch_version = ch_versions.mix( PICARD_CREATESEQUENCEDICTIONARY.out.versions )

    ch_dicted_formix =  PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
                            .join( ch_dict_input.remainder, failOnMismatch: true )
                            .map {
                                meta, dict, fasta, fai, mapper_index, circular_target ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target ]
                            }

    ch_dict_formapperindexing = ch_fasta_for_dict.skip.mix(ch_dicted_formix)

    //
    // INDEXING: Mapping indicies
    //

    // Generate mapper indicies if not supplied, and if supplied generate meta

    ch_fasta_for_mapperindex = ch_dict_formapperindexing
        .branch {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                forindex: mapper_index == ""
                skip: true
        }

    ch_mapindex_input = ch_fasta_for_mapperindex
        .forindex
        .multiMap {
            meta, fasta, fai, dict, mapper_index, circular_target ->
                index:      [ meta, fasta ]
                remainder:  [ meta, fasta, fai, dict, circular_target ]
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
                                meta, mapper_index, fasta, fai, dict, circular_target ->

                                [ meta, fasta, fai, dict, mapper_index, circular_target ]
                            }

    ch_indexmapper_for_reference = ch_fasta_for_mapperindex.skip.mix(ch_indexed_formix)

    emit:
    reference            = ch_indexmapper_for_reference                      // [ meta, fasta, fai, dict, mapindex, circular_target ]
    mitochondrion_header = ch_input_from_referencesheet.mitochondrion_header // [ meta, mitochondrion ]
    hapmap               = ch_input_from_referencesheet.angsd_hapmap         // [ meta, hapmap ]
    pmd_mask             = ch_input_from_referencesheet.pmd_mask             // [ meta, pmd_mask, capture_bed ]
    snp_capture_bed      = ch_input_from_referencesheet.snp_bed              // [ meta, capture_bed ]
    pileupcaller_snp     = ch_input_from_referencesheet.pileupcaller_snp     // [ meta, pileupcaller_snp, pileupcaller_bed ]
    sexdeterrmine_bed    = ch_input_from_referencesheet.sexdeterrmine_bed    // [ meta, sexdet_bed ]
    bedtools_feature     = ch_input_from_referencesheet.bedtools_feature     // [ meta, bedtools_feature ]
    dbsnp                = ch_input_from_referencesheet.dbsnp                // [ meta + [ ploidy: genotyping_gatk_ploidy ], genotyping_gatk_dbsnp ]
    versions             = ch_versions
}
