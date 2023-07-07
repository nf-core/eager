//
// Index input reference as required
//

include { GUNZIP                          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
// TODO missing: circulargeneraotr?

workflow REFERENCE_INDEXING_SINGLE {

    take:
    fasta                // file: /path/to/name.{fasta,fa,fna,fas,fasta.gz,fa.gz,fna.gz,fas.gz}
    fasta_fai            // file: /path/to/name.{fasta,fa,fna,fas}.fai
    fasta_dict           // file: /path/to/fasta.dict
    fasta_mapperindexdir // file: /path/to/

    main:

    ch_versions = Channel.empty()

    def fasta_ext = WorkflowEager.grabUngzippedExtension(fasta)
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
        ch_fasta_fai = SAMTOOLS_FAIDX ( ch_ungz_ref ).fai.map{[ [id: clean_name - '.fai'], it[1] ] }
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
    if ( params.mapping_tool == 'bwaaln' || params.mapping_tool == 'bwamem' ){

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

    }

    // Join all together into a single map. failOnMismatch allows check if
    // a user supplies indicies with different 'base' names.
    ch_reference_for_mapping = ch_ungz_ref
                                .join(ch_fasta_fai, failOnMismatch: true)
                                .join(ch_fasta_dict, failOnMismatch: true)
                                .join(ch_fasta_mapperindexdir, failOnMismatch: true)
                                .map{
                                    meta, fasta, fai, dict, mapper_index ->
                                    //TODO: add params for snpcapturebed, snpeigenstrat and sexdetbed once implemented
                                    def contamination_estimation_angsd_hapmap = params.contamination_estimation_angsd_hapmap != null ? file( params.contamination_estimation_angsd_hapmap, checkIfExists: true ) : ""
                                    def pmd_mask                              = params.damage_manipulation_pmdtools_reference_mask != null ? file(params.damage_manipulation_pmdtools_reference_mask, checkIfExists: true ) : ""
                                    def capture_bed                           = ""
                                    def snp_eigenstrat                        = ""
                                    def sexdet_bed                            = ""
                                    [ meta, fasta, fai, dict, mapper_index, params.fasta_circular_target, params.mitochondrion_header, contamination_estimation_angsd_hapmap, pmd_mask, capture_bed, snp_eigenstrat, sexdet_bed ]
                                }

    ch_ref_index_single = ch_reference_for_mapping
                                .multiMap{
                                    meta, fasta, fai, dict, mapper_index, circular_target, mitochondrion_header, contamination_estimation_angsd_hapmap, pmd_mask, capture_bed, snp_eigenstrat, sexdet_bed ->
                                    reference:         [ meta, fasta, fai, dict, mapper_index, circular_target ]
                                    mito_header:       [ meta, mitochondrion_header ]
                                    hapmap:            [ meta, contamination_estimation_angsd_hapmap ]
                                    pmd_mask:          [ meta, pmd_mask, capture_bed ]
                                    snp_bed:           [ meta, capture_bed ]
                                    snp_eigenstrat:    [ meta, snp_eigenstrat ]
                                    sexdeterrmine_bed: [ meta, sexdet_bed ]
                                }

    emit:
    reference            = ch_ref_index_single.reference          // [ meta, fasta, fai, dict, mapindex, circular_target ]
    mitochondrion_header = ch_ref_index_single.mito_header        // [ meta, mito_header ]
    hapmap               = ch_ref_index_single.hapmap             // [ meta, hapmap ]
    pmd_mask             = ch_ref_index_single.pmd_mask           // [ meta, pmd_mask, capture_bed ]
    snp_capture_bed      = ch_ref_index_single.snp_bed            // [ meta, capture_bed ]
    snp_eigenstrat       = ch_ref_index_single.snp_eigenstrat     // [ meta, snp_eigenstrat ]
    sexdeterrmine_bed    = ch_ref_index_single.sexdeterrmine_bed  // [ meta, sexdet_bed ]
    versions             = ch_versions

}
