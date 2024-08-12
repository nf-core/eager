//
// Prepare reference indexing for downstream
//

include { REFERENCE_INDEXING_SINGLE        } from '../../subworkflows/local/reference_indexing_single.nf'
include { REFERENCE_INDEXING_MULTI         } from '../../subworkflows/local/reference_indexing_multi.nf'
include { GUNZIP as GUNZIP_PMDBED          } from '../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_PMDFASTA        } from '../../modules/nf-core/gunzip/main.nf'
include { GUNZIP as GUNZIP_SNPBED          } from '../../modules/nf-core/gunzip/main.nf'
include { ELONGATE_REFERENCE               } from '../../subworkflows/local/elongate_reference.nf'

workflow REFERENCE_INDEXING {
    take:
    fasta // file: /path/to/samplesheet.csv
    fasta_fai
    fasta_dict
    fasta_mapperindexdir

    main:
    ch_versions = Channel.empty()

    // Warn user if they've given a reference sheet that already includes fai/dict/mapper index etc.
    if ( ( fasta.extension == 'csv' || fasta.extension == 'tsv' ) && ( fasta_fai || fasta_dict || fasta_mapperindexdir )) log.warn("A TSV or CSV has been supplied to `--fasta_sheet` as well as e.g. `--fasta_fai`. --fasta_sheet CSV/TSV takes priority and --fasta_* parameters will be ignored.")
    if ( ( fasta.extension == 'csv' || fasta.extension == 'tsv' ) && ( params.mitochondrion_header || params.contamination_estimation_angsd_hapmap || params.damage_manipulation_pmdtools_reference_mask || params.damage_manipulation_pmdtools_reference_mask || params.snpcapture_bed || params.genotyping_pileupcaller_bedfile || params.genotyping_pileupcaller_snpfile || params.sexdeterrmine_bedfile || params.mapstats_bedtools_featurefile || params.genotyping_reference_ploidy || params.genotyping_gatk_dbsnp || params.fasta_circular_target || params.circularmapper_elongated_fasta || params.circularmapper_elongated_fai )) log.warn("A TSV or CSV has been supplied to `--fasta_sheet` as well as individual reference-specific input files, e.g. `--contamination_estimation_angsd_hapmap`. Input files specified in the --fasta_sheet CSV/TSV take priority and other input parameters will be ignored.")

    if ( fasta.extension == 'csv' || fasta.extension == 'tsv' ) {
        // If input (multi-)reference sheet supplied
        REFERENCE_INDEXING_MULTI ( fasta )
        ch_reference_for_mapping = REFERENCE_INDEXING_MULTI.out.reference
        ch_reference_to_elongate   = REFERENCE_INDEXING_MULTI.out.elongated_reference
        ch_mitochondrion_header  = REFERENCE_INDEXING_MULTI.out.mitochondrion_header
        ch_hapmap                = REFERENCE_INDEXING_MULTI.out.hapmap
        ch_pmd_masked_fasta      = REFERENCE_INDEXING_MULTI.out.pmd_masked_fasta
        ch_pmd_bed_for_masking   = REFERENCE_INDEXING_MULTI.out.pmd_bed_for_masking
        ch_snp_capture_bed       = REFERENCE_INDEXING_MULTI.out.snp_capture_bed
        ch_pileupcaller_bed_snp  = REFERENCE_INDEXING_MULTI.out.pileupcaller_bed_snp
        ch_sexdeterrmine_bed     = REFERENCE_INDEXING_MULTI.out.sexdeterrmine_bed
        ch_bedtools_feature      = REFERENCE_INDEXING_MULTI.out.bedtools_feature
        ch_dbsnp                 = REFERENCE_INDEXING_MULTI.out.dbsnp
        ch_versions = ch_versions.mix( REFERENCE_INDEXING_MULTI.out.versions )
    } else {
        // If input FASTA and/or indicies supplied
        REFERENCE_INDEXING_SINGLE ( fasta, fasta_fai, fasta_dict, fasta_mapperindexdir )
        ch_reference_to_elongate   = REFERENCE_INDEXING_SINGLE.out.elongated_reference
        ch_mitochondrion_header  = REFERENCE_INDEXING_SINGLE.out.mitochondrion_header
        ch_hapmap                = REFERENCE_INDEXING_SINGLE.out.hapmap
        ch_pmd_masked_fasta      = REFERENCE_INDEXING_SINGLE.out.pmd_masked_fasta
        ch_pmd_bed_for_masking   = REFERENCE_INDEXING_SINGLE.out.pmd_bed_for_masking
        ch_snp_capture_bed       = REFERENCE_INDEXING_SINGLE.out.snp_capture_bed
        ch_pileupcaller_bed_snp  = REFERENCE_INDEXING_SINGLE.out.pileupcaller_bed_snp
        ch_sexdeterrmine_bed     = REFERENCE_INDEXING_SINGLE.out.sexdeterrmine_bed
        ch_bedtools_feature      = REFERENCE_INDEXING_SINGLE.out.bedtools_feature
        ch_reference_for_mapping = REFERENCE_INDEXING_SINGLE.out.reference
        ch_dbsnp                 = REFERENCE_INDEXING_SINGLE.out.dbsnp
        ch_versions = ch_versions.mix( REFERENCE_INDEXING_SINGLE.out.versions )
    }

    // Filter out input options that are not provided and unzip if necessary
    ch_mitochondrion_header = ch_mitochondrion_header
                    .filter{ it[1] != "" }

    ch_hapmap = ch_hapmap
                    .filter{ it[1] != "" }

    ch_pmd_masked_fasta = ch_pmd_masked_fasta
                    .branch {
                        meta, pmd_masked_fasta ->
                        input: pmd_masked_fasta != ""
                        skip: true
                    }
    ch_pmd_masked_fasta_gunzip = ch_pmd_masked_fasta.input
                    .branch {
                        meta, pmd_masked_fasta ->
                        forgunzip: pmd_masked_fasta.extension == "gz"
                        skip: true
                    }
    GUNZIP_PMDFASTA( ch_pmd_masked_fasta_gunzip.forgunzip )
    ch_pmd_masked_fasta = ch_pmd_masked_fasta_gunzip.skip.mix( GUNZIP_PMDFASTA.out.gunzip ).mix( ch_pmd_masked_fasta.skip )
    ch_version = ch_versions.mix( GUNZIP_PMDFASTA.out.versions.first() )

    ch_pmd_bed_for_masking = ch_pmd_bed_for_masking
                    .branch {
                        meta, pmd_bed_for_masking ->
                        input: pmd_bed_for_masking != ""
                        skip: true
                    }
    ch_pmd_bed_for_masking_gunzip = ch_pmd_bed_for_masking.input
                    .branch {
                        meta, pmd_bed_for_masking ->
                        forgunzip: pmd_bed_for_masking.extension == "gz"
                        skip: true
                    }
    GUNZIP_PMDBED( ch_pmd_bed_for_masking_gunzip.forgunzip )
    ch_pmd_bed_for_masking = ch_pmd_bed_for_masking_gunzip.skip.mix( GUNZIP_PMDBED.out.gunzip ).mix( ch_pmd_bed_for_masking.skip )
    ch_version = ch_versions.mix( GUNZIP_PMDBED.out.versions.first() )

    ch_pmd_masking = ch_pmd_masked_fasta
                    .combine( by: 0, ch_pmd_bed_for_masking )

    ch_capture_bed = ch_snp_capture_bed //bed file input is optional, so no filtering
                    .branch {
                        meta, capture_bed ->
                        input: capture_bed != ""
                        skip: true
                    }
    ch_capture_bed_gunzip = ch_capture_bed.input //unzip
                    .branch {
                        meta, capture_bed ->
                        forgunzip: capture_bed.extension == "gz"
                        skip: true
                    }
    GUNZIP_SNPBED( ch_capture_bed_gunzip.forgunzip )
    ch_capture_bed = GUNZIP_SNPBED.out.gunzip.mix( ch_capture_bed_gunzip.skip ).mix( ch_capture_bed.skip )
    ch_version = ch_versions.mix( GUNZIP_SNPBED.out.versions.first() )

    ch_pileupcaller_bed_snp = ch_pileupcaller_bed_snp
                    .filter { it[1] != "" || it[2] != "" } // They go together or not at all.
                    // Check if the channel is empty, and throw an error. Will only trigger for tsv fasta input. Single reference gets validated immediately.
                    .ifEmpty { if(params.run_genotyping && params.genotyping_tool == 'pileupcaller') { error "[nf-core/eager] ERROR: Genotyping with pileupcaller requires that both '--genotyping_pileupcaller_bedfile' AND '--genotyping_pileupcaller_snpfile' are provided for at least one reference genome." } }
                    .filter{ it != null } // Remove null channel which arises if empty cause error returns null.

    ch_sexdeterrmine_bed = ch_sexdeterrmine_bed
                    .filter { it[1] != "" }

    ch_bedtools_feature = ch_bedtools_feature
                    .filter { it[1] != "" }

    ch_dbsnp = ch_dbsnp
        .filter { it[1] != "" }

    // Elongate reference for circularmapper if requested
    if ( params.mapping_tool == "circularmapper" ) {
        // Throw errors if required parameters are missing
        //  A circular target is required even when an elongated reference has been provided.
        ch_elongated_for_gunzip = ch_reference_to_elongate
                        .filter{
                            meta, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex ->
                            circular_target != ""
                        }
                        .ifEmpty{ error "[nf-core/eager] ERROR: Mapping with circularmapper requires either a circular target for at least one reference." }

        // This ELONGATE_REFERENCE subworkflow also checks if the provided reference is gzipped, and unzips it if necessary.
        ELONGATE_REFERENCE( ch_reference_for_mapping, ch_reference_to_elongate )
        ch_version = ch_versions.mix( ELONGATE_REFERENCE.out.versions )
        ch_elongated_indexed_reference = ELONGATE_REFERENCE.out.circular_reference
        ch_elongated_chr_list = ELONGATE_REFERENCE.out.elongated_chr_list

    } else {
        ch_elongated_indexed_reference = ch_reference_to_elongate
        ch_elongated_chr_list = Channel.empty()
    }

    emit:
    reference            = ch_reference_for_mapping       // [ meta, fasta, fai, dict, mapindex ]
    elongated_reference  = ch_elongated_indexed_reference // [ meta, circularmapper_elongated_fasta, circularmapper_elongated_index ]
    elongated_chr_list   = ch_elongated_chr_list          // [ meta, elongated_chr_list ]
    mitochondrion_header = ch_mitochondrion_header        // [ meta, mitochondrion_header ]
    hapmap               = ch_hapmap                      // [ meta, hapmap ]
    pmd_masking          = ch_pmd_masking                 // [ meta, pmd_masked_fasta, pmd_bed_for_masking ]
    pmd_bed_for_masking  = ch_pmd_bed_for_masking         // [ meta, pmd_bed_for_masking ]
    snp_capture_bed      = ch_capture_bed                 // [ meta, capture_bed ]
    pileupcaller_bed_snp = ch_pileupcaller_bed_snp        // [ meta, pileupcaller_bed, pileupcaller_snp ]
    sexdeterrmine_bed    = ch_sexdeterrmine_bed           // [ meta, sexdet_bed ]
    bedtools_feature     = ch_bedtools_feature            // [ meta, bedtools_feature ]
    dbsnp                = ch_dbsnp                       // [ meta, dbsnp ]
    versions             = ch_versions

}
