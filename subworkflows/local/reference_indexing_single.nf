//
// Index input reference as required
//
include { grabUngzippedExtension          } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

include { GUNZIP                          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { MAPAD_INDEX                     } from '../../modules/nf-core/mapad/index/main'

workflow REFERENCE_INDEXING_SINGLE {

    take:
    fasta                // file: /path/to/name.{fasta,fa,fna,fas,fasta.gz,fa.gz,fna.gz,fas.gz}
    fasta_fai            // file: /path/to/name.{fasta,fa,fna,fas}.fai
    fasta_dict           // file: /path/to/fasta.dict
    fasta_mapperindexdir // file: /path/to/

    main:

    ch_versions = Channel.empty()

    def fasta_ext = grabUngzippedExtension(fasta)
    def clean_name = fasta.name.toString() - fasta_ext

    // Detect if fasta is gzipped or not, unzip if necessary, and generate meta ID by sanitizing file
    if ( fasta.extension == 'gz' ) {
        ch_gz_ref = Channel.fromPath(fasta).map{[[], it]}
        GUNZIP ( ch_gz_ref )
        ch_ungz_ref = GUNZIP.out.gunzip.map{[[id: clean_name], it[1] ]}
        ch_versions = ch_versions.mix( GUNZIP.out.versions.first())
    } else {
        ch_ungz_ref = Channel.fromPath(fasta).map{[[id: clean_name], it ]}
    }

    // Generate FAI if not supplied, and if supplied generate meta ID
    if ( !fasta_fai ) {
        ch_fasta_fai = SAMTOOLS_FAIDX ( ch_ungz_ref, [ [], [] ] ).fai.map{[ [id: clean_name - '.fai'], it[1] ] }
        ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions.first())
    } else {
        ch_fasta_fai = Channel.fromPath(fasta_fai).map{[[id: clean_name], it ]}
    }

    // Generate DICT if not supplied, and if supplied generate meta
    if ( !fasta_dict ) {
        ch_fasta_dict = PICARD_CREATESEQUENCEDICTIONARY ( ch_ungz_ref ).reference_dict.map{[ [id: clean_name - '.fai'], it[1] ] }
        ch_versions = ch_versions.mix( PICARD_CREATESEQUENCEDICTIONARY.out.versions.first())
    } else {
        ch_fasta_dict = Channel.fromPath(fasta_dict).map{[[id: clean_name], it ]}
    }

    // Generate mapper indicies if not supplied, and if supplied generate meta
    if ( params.mapping_tool == 'bwaaln' || params.mapping_tool == 'bwamem' || params.mapping_tool == 'circularmapper' ){

        if ( !fasta_mapperindexdir ) {
            ch_fasta_mapperindexdir = BWA_INDEX ( ch_ungz_ref ).index
            ch_versions = ch_versions.mix( BWA_INDEX.out.versions.first())
        } else {
            ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
        }

    } else if ( params.mapping_tool == "bowtie2" ) {

        if ( !fasta_mapperindexdir ) {
            ch_fasta_mapperindexdir = BOWTIE2_BUILD ( ch_ungz_ref ).index
            ch_versions = ch_versions.mix( BOWTIE2_BUILD.out.versions.first())
        } else {
            ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
        }

    } else if ( params.mapping_tool == "mapad" ) {

        if ( !fasta_mapperindexdir ) {
            ch_fasta_mapperindexdir = MAPAD_INDEX ( ch_ungz_ref ).index
            ch_versions = ch_versions.mix( MAPAD_INDEX.out.versions.first())
        } else {
            ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
        }
    }

    // Join all together into a single map. failOnMismatch allows check if
    // a user supplies indicies with different 'base' names.
    ch_reference_for_mapping = ch_ungz_ref
                                .join(ch_fasta_fai, failOnMismatch: true)
                                .join(ch_fasta_dict, failOnMismatch: true)
                                .join(ch_fasta_mapperindexdir, failOnMismatch: true)
                                .map{
                                    meta, fasta, fai, dict, mapper_index ->
                                    def contamination_estimation_angsd_hapmap = params.contamination_estimation_angsd_hapmap != null ? file( params.contamination_estimation_angsd_hapmap, checkIfExists: true ) : ""
                                    def pmd_masked_fasta                      = params.damage_manipulation_pmdtools_masked_reference != null ? file(params.damage_manipulation_pmdtools_masked_reference, checkIfExists: true ) : ""
                                    def pmd_bed_for_masking                   = params.damage_manipulation_pmdtools_reference_mask != null ? file(params.damage_manipulation_pmdtools_reference_mask, checkIfExists: true ) : ""
                                    def capture_bed                           = params.snpcapture_bed != null ? file(params.snpcapture_bed, checkIfExists: true ) : ""
                                    def pileupcaller_bed                      = params.genotyping_pileupcaller_bedfile != null ? file(params.genotyping_pileupcaller_bedfile, checkIfExists: true ) : ""
                                    def pileupcaller_snp                      = params.genotyping_pileupcaller_snpfile != null ? file(params.genotyping_pileupcaller_snpfile, checkIfExists: true ) : ""
                                    def sexdet_bed                            = params.sexdeterrmine_bedfile != null ? file(params.sexdeterrmine_bedfile, checkIfExists: true ) : ""
                                    def bedtools_feature                      = params.mapstats_bedtools_featurefile != null ? file(params.mapstats_bedtools_featurefile, checkIfExists: true ) : ""
                                    def genotyping_reference_ploidy           = params.genotyping_reference_ploidy
                                    def genotyping_gatk_dbsnp                 = params.genotyping_gatk_dbsnp != null ? file(params.genotyping_gatk_dbsnp, checkIfExists: true ) : ""
                                    def circularmapper_elongated_fasta        = params.fasta_circularmapper_elongatedfasta != null ? file( params.fasta_circularmapper_elongatedfasta, checkIfExists: true ) : ""
                                    def circularmapper_elongated_index        = params.fasta_circularmapper_elongatedindex != null ? file( params.fasta_circularmapper_elongatedindex, checkIfExists: true ) : ""
                                    [ meta + [ ploidy: genotyping_reference_ploidy ], fasta, fai, dict, mapper_index, params.fasta_circular_target, params.mitochondrion_header, contamination_estimation_angsd_hapmap, pmd_masked_fasta, pmd_bed_for_masking, capture_bed, pileupcaller_bed, pileupcaller_snp, sexdet_bed, bedtools_feature, genotyping_gatk_dbsnp, circularmapper_elongated_fasta, circularmapper_elongated_index ]
                                }

    ch_ref_index_single = ch_reference_for_mapping
                                .multiMap{
                                    meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion_header, contamination_estimation_angsd_hapmap, pmd_masked_fasta, pmd_bed_for_masking, capture_bed, pileupcaller_bed, pileupcaller_snp, sexdet_bed, bedtools_feature, genotyping_gatk_dbsnp, circularmapper_elongated_fasta, circularmapper_elongated_index ->
                                    reference:              [ meta, fasta, fai, dict, mapper_index ]
                                    circularmapper:         [ meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ]
                                    mito_header:            [ meta, mitochondrion_header ]
                                    hapmap:                 [ meta, contamination_estimation_angsd_hapmap ]
                                    pmd_masked_fasta:       [ meta, pmd_masked_fasta ]
                                    pmd_bed_for_masking:    [ meta, pmd_bed_for_masking ]
                                    snp_bed:                [ meta, capture_bed ]
                                    pileupcaller_bed_snp:   [ meta, pileupcaller_bed, pileupcaller_snp ]
                                    sexdeterrmine_bed:      [ meta, sexdet_bed ]
                                    bedtools_feature:       [ meta, bedtools_feature ]
                                    dbsnp:                  [ meta, genotyping_gatk_dbsnp ]
                                }

    emit:
    reference            = ch_ref_index_single.reference             // [ meta, fasta, fai, dict, mapindex ]
    elongated_reference  = ch_ref_index_single.circularmapper        // [ meta, circular_target, circularmapper_elongated_fasta, circularmapper_elongated_index ]
    mitochondrion_header = ch_ref_index_single.mito_header           // [ meta, mito_header ]
    hapmap               = ch_ref_index_single.hapmap                // [ meta, hapmap ]
    pmd_masked_fasta     = ch_ref_index_single.pmd_masked_fasta      // [ meta, pmd_masked_fasta ]
    pmd_bed_for_masking  = ch_ref_index_single.pmd_bed_for_masking   // [ meta, pmd_bed_for_masking ]
    snp_capture_bed      = ch_ref_index_single.snp_bed               // [ meta, capture_bed ]
    pileupcaller_bed_snp = ch_ref_index_single.pileupcaller_bed_snp  // [ meta, pileupcaller_bed, pileupcaller_snp ]
    sexdeterrmine_bed    = ch_ref_index_single.sexdeterrmine_bed     // [ meta, sexdet_bed ]
    bedtools_feature     = ch_ref_index_single.bedtools_feature      // [ meta, bedtools_feature ]
    dbsnp                = ch_ref_index_single.dbsnp                 // [ meta, genotyping_gatk_dbsnp ]
    versions             = ch_versions

}
