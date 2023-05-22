//
// Calculate PMD scores, trim, or rescale DNA damage from mapped reads.
//

include { MAPDAMAGE2                                              } from '../../modules/nf-core/mapdamage2/main'
include { PMDTOOLS_FILTER                                         } from '../../modules/nf-core/pmdtools/filter/main'
include { BAMUTIL_TRIMBAM                                         } from '../../modules/nf-core/bamutil/trimbam/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_DAMAGE_RESCALED     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_DAMAGE_FILTERED     } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX    as SAMTOOLS_INDEX_DAMAGE_TRIMMED      } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_DAMAGE_FILTERED  } from '../../modules/nf-core/samtools/flagstat/main'

// TODO: Add required channels and channel manipulations for reference-dependent bed masking before pmdtools. Requires multi-ref support before implementation.
workflow MANIPULATE_DAMAGE {
    take:
    ch_bam_bai  // [ [ meta ], bam , bai ]
    ch_fasta    // [ [ meta ], fasta ]

    main:
    ch_versions              = Channel.empty()
    ch_multiqc_files         = Channel.empty()
    ch_rescaled_bams         = Channel.empty()
    ch_pmd_filtered_bams     = Channel.empty()
    ch_trimmed_bams          = Channel.empty()
    ch_pmd_filtered_flagstat = Channel.empty() // Only run flagstat on pmd filtered bam, since rescaling and trimming does not change the number of reads

    // Ensure correct reference is associated with each bam_bai pair
    ch_refs = ch_fasta
        .map {
            // Create additional map containing only meta.reference for combining samples and reference fastas
            meta, fasta ->
                meta2 = [:]
                meta2.reference = meta.id
            [ meta2, meta, fasta ]
        }

    ch_input_for_damage_manipulation = ch_bam_bai
        .map {
            meta, bam, bai ->
                meta2 = [:]
                meta2.reference = meta.reference
            [ meta2, meta, bam, bai ]
        }
        .combine(ch_refs, by: 0 ) // [ [combine_meta], [meta], bam, bai, [ref_meta], fasta ]

    if ( params.run_mapdamage_rescaling ) {
        ch_mapdamage_prep = ch_input_for_damage_manipulation
            .branch {
                ignore_me, meta, bam, bai, ref_meta, fasta ->
                skip:    meta.damage_treatment == 'full'
                no_skip: true
            }

        ch_skip_rescale = ch_mapdamage_prep.skip
            .map {
                ignore_me, meta, bam, bai, ref_meta, fasta ->
                [ meta, bam, bai ]
            }

        ch_mapdamage_input = ch_mapdamage_prep.no_skip
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta ->
                    bam: [ meta, bam ]
                    fasta: fasta
            }

        MAPDAMAGE2( ch_mapdamage_input.bam, ch_mapdamage_input.fasta )
        ch_versions       = ch_versions.mix( MAPDAMAGE2.out.versions.first() )

        SAMTOOLS_INDEX_DAMAGE_RESCALED( MAPDAMAGE2.out.rescaled )
        ch_versions       = ch_versions.mix( SAMTOOLS_INDEX_DAMAGE_RESCALED.out.versions.first() )
        ch_rescaled_index = params.fasta_largeref ? SAMTOOLS_INDEX_DAMAGE_RESCALED.out.csi : SAMTOOLS_INDEX_DAMAGE_RESCALED.out.bai

        // TODO When adding library-level data merging pre-genotyping, make sure that rescaled bams are not merged in any way as the underlying damage model could differ between libraries
        ch_rescaled_bams  = MAPDAMAGE2.out.rescaled.join(ch_rescaled_index)
            .mix(ch_skip_rescale) // Should these be mixed actually, or excluded? Might not make sense to take rescaled and non-rescaled bams togetehr for anything downstream.
    }

    if ( params.run_pmd_filtering ) {
        // TODO Add module to produce Masked reference from given references and bed file (with meta specifying the reference it matches)?
        // if ( params.pmdtools_reference_mask) {
        //     MASK_REFERENCE_BY_BED()
        // }

        ch_pmdtools_input = ch_input_for_damage_manipulation
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta ->
                    bam: [ meta, bam, bai ]
                    fasta: fasta
            }

        PMDTOOLS_FILTER( ch_pmdtools_input.bam, params.damage_manipulation_pmdtools_threshold, ch_pmdtools_input.fasta )
        ch_versions       = ch_versions.mix( PMDTOOLS_FILTER.out.versions.first() )

        SAMTOOLS_INDEX_DAMAGE_FILTERED( PMDTOOLS_FILTER.out.bam )
        ch_versions       = ch_versions.mix( SAMTOOLS_INDEX_DAMAGE_FILTERED.out.versions.first() )
        ch_filtered_index = params.fasta_largeref ? SAMTOOLS_INDEX_DAMAGE_FILTERED.out.csi : SAMTOOLS_INDEX_DAMAGE_FILTERED.out.bai

        ch_pmd_filtered_bams = PMDTOOLS_FILTER.out.bam.join( ch_filtered_index )

        SAMTOOLS_FLAGSTAT_DAMAGE_FILTERED( ch_pmd_filtered_bams )
        ch_pmd_filtered_flagstat = SAMTOOLS_FLAGSTAT_DAMAGE_FILTERED.out.flagstat
    }

    if ( params.run_trim_bam ) {
        if ( params.run_pmd_filtering ) {
            ch_to_trim = ch_pmd_filtered_bams
                .map{
                    meta, bam, bai ->
                    [ meta, bam ]
                }
        } else {
            ch_to_trim = ch_bam_bai
                .map {
                    meta, bam, bai ->
                    [ meta, bam ]
                }
        }

        ch_trimbam_input = ch_to_trim
            .map {
                meta, bam ->
                    trim_left  = meta.strandedness == 'single' ? ( meta.damage_treatment == 'none' ? params.damage_manipulation_bamutils_clip_single_stranded_none_udg_left  : meta.damage_treatment == 'half' ? params.damage_manipulation_bamutils_clip_single_stranded_half_udg_left  : 0 ) : ( meta.damage_treatment == 'none' ? params.damage_manipulation_bamutils_clip_double_stranded_none_udg_left  : meta.damage_treatment == 'half' ? params.damage_manipulation_bamutils_clip_double_stranded_half_udg_left  : 0 )
                    trim_right = meta.strandedness == 'single' ? ( meta.damage_treatment == 'none' ? params.damage_manipulation_bamutils_clip_single_stranded_none_udg_right : meta.damage_treatment == 'half' ? params.damage_manipulation_bamutils_clip_single_stranded_half_udg_right : 0 ) : ( meta.damage_treatment == 'none' ? params.damage_manipulation_bamutils_clip_double_stranded_none_udg_right : meta.damage_treatment == 'half' ? params.damage_manipulation_bamutils_clip_double_stranded_half_udg_right : 0 )
                [ meta, bam, trim_left, trim_right ]
            }

        BAMUTIL_TRIMBAM( ch_trimbam_input )
        ch_versions      = ch_versions.mix( BAMUTIL_TRIMBAM.out.versions.first() )

        SAMTOOLS_INDEX_DAMAGE_TRIMMED( BAMUTIL_TRIMBAM.out.bam )
        ch_versions      = ch_versions.mix( SAMTOOLS_INDEX_DAMAGE_TRIMMED.out.versions.first() )
        ch_trimmed_index = params.fasta_largeref ? SAMTOOLS_INDEX_DAMAGE_TRIMMED.out.csi : SAMTOOLS_INDEX_DAMAGE_TRIMMED.out.bai

        ch_trimmed_bams  = BAMUTIL_TRIMBAM.out.bam.join( ch_trimmed_index )
    }

    emit:
    rescaled = ch_rescaled_bams         // [ meta, bam, bai ]
    filtered = ch_pmd_filtered_bams     // [ meta, bam, bai ]
    trimmed  = ch_trimmed_bams          // [ meta, bam, bai ]
    flagstat = ch_pmd_filtered_flagstat // [ meta, flagstat ]
    versions = ch_versions
}
