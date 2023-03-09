//
// Calculate PMD scores, trim, or rescale DNA damage from mapped reads.
//

include { MAPDAMAGE2                                          } from '../../modules/nf-core/mapdamage2/main'
include { PMDTOOLS_FILTER                                     } from '../../modules/nf-core/pmdtools/filter/main'
include { BAMUTIL_TRIMBAM                                     } from '../../modules/nf-core/bamutil/trimbam/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_DAMAGE_MANIPULATED } from '../../modules/nf-core/samtools/index/main'


workflow MANIPULATE_DAMAGE {
    take:
    ch_bam_bai  // [ [ meta ], bam , bai ]
    ch_fasta    // [ [ meta ], fasta ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    // Ensure correct reference is associated with each bam_bai pair
    ch_refs = ch_fasta.map {
        // Create additional map containing only meta.reference for combining samples and reference fastas
        meta, fasta ->
            meta2 = [:]
            meta2.reference = meta.id
        [ meta2, meta, fasta ]
    }

    ch_input_for_damage_manipulation = ch_bam_bai.map {
        meta, bam, bai ->
            meta2 = [:]
            meta2.reference = meta.reference
        [ meta2, meta, bam, bai ]
    }
    .combine(ch_refs, by: 0 )

    if ( params.run_mapdamage_rescaling ) {
        ch_mapdamage_input = ch_input_for_damage_manipulation.multiMap {
            ignore_me, meta, bam, bai, ref_meta, fasta ->
                bam: [ meta, bam ]
                fasta: fasta
        }

        MAPDAMAGE2( ch_mapdamage_input.bam, ch_mapdamage_input.fasta )
        ch_versions = ch_versions.mix( MAPDAMAGE2.out.versions.first() )

        ch_rescaled_bams = MAPDAMAGE2.out.rescaled
    }

    if ( params.run_pmd_filtering ) {
        // TODO Add local module to produce Masked reference from given references and bed file (with meta specifying the reference it matches)?
        // if ( params.pmdtools_reference_mask) {
        //     MASK_REFERENCE_BY_BED()
        // }

        // TODO no PMD fitlering for full_UDG libraries
        ch_pmdtools_input = ch_input_for_damage_manipulation.multiMap {
            ignore_me, meta, bam, bai, ref_meta, fasta ->
                bam: [ meta, bam, bai ]
                fasta: fasta
        }


        PMDTOOLS_FILTER(ch_pmdtools_input.bam, params.pmdtools_threshold, ch_pmdtools_input.fasta)
        ch_versions = ch_versions.mix( PMDTOOLS_FILTER.out.versions.first() )

        ch_damage_filtered_bams = PMDTOOLS_FILTER.out.bam
    }

    if ( params.run_trim_bam ) {
        if ( params.run_pmd_filtering ) {
            ch_to_trim = ch_damage_filtered_bams
        } else {
            ch_to_trim = ch_bam_bai.map {
                meta, bam, bai ->
                [ meta, bam ]
            }
        }

        ch_trimbam_input = ch_to_trim.map {
            meta, bam ->
            // TODO this is a mess and throws syntax error
                trim_left  = meta.strandedness == 'single' ? ( meta.damage_treatment == 'none' ? ( params.bamutils_clip_single_stranded_none_udg_left  : meta.damage_treatment == 'half' ? bamutils_clip_single_stranded_half_udg_left  : 0  ) ) : ( meta.damage_treatment == 'none' ? ( params.bamutils_clip_double_stranded_none_udg_left  : meta.damage_treatment == 'half' ? bamutils_clip_double_stranded_half_udg_left  : 0 ) )
                trim_right = meta.strandedness == 'single' ? ( meta.damage_treatment == 'none' ? ( params.bamutils_clip_single_stranded_none_udg_right : meta.damage_treatment == 'half' ? bamutils_clip_single_stranded_half_udg_right : 0  ) ) : ( meta.damage_treatment == 'none' ? ( params.bamutils_clip_double_stranded_none_udg_right : meta.damage_treatment == 'half' ? bamutils_clip_double_stranded_half_udg_right : 0 ) )
            [ meta, bam, trim_left, trim_right ]
        }

        BAMUTIL_TRIMBAM( ch_trimbam_input )
        ch_versions = ch_versions.mix( BAMUTIL_TRIMBAM.out.versions.first() )

        ch_trimmed_bams = BAMUTIL_TRIMBAM.out.bam
    }

    ch_damage_manipulated_bam = params.run_mapdamage_rescaling ? ch_rescaled_bams : ( params.run_trim_bam ? ch_trimmed_bams : ch_damage_filtered_bams )

    // INDEX DAMAGE MANIPULATED BAM
    SAMTOOLS_INDEX_DAMAGE_MANIPULATED( ch_damage_manipulated_bam )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX_DAMAGE_MANIPULATED.out.versions.first() )
    ch_damage_manipulated_index =  params.fasta_largeref ? SAMTOOLS_INDEX_DAMAGE_MANIPULATED.out.csi : SAMTOOLS_INDEX_DAMAGE_MANIPULATED.out.bai

    emit:
    bam      = ch_damage_manipulated_bam    // [ meta, bam ]
    bai      = ch_damage_manipulated_index  // [ meta, bai ]
    versions = ch_versions
}
