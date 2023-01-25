//
//  Carry out per-chromosome deduplication
//

include { BUILD_INTERVALS        } from '../../modules/local/build_intervals'
include { BAM_SPLIT_BY_REGION    } from '../../subworkflows/nf-core/bam_split_by_region/main'
include { PICARD_MARKDUPLICATES  } from '../../modules/nf-core/picard/markduplicates/main'
include { DEDUP                  } from '../../modules/nf-core/dedup/main'
include { SAMTOOLS_MERGE         } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT          } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX         } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT      } from '../../modules/nf-core/samtools/flagstat/main'

workflow DEDUPLICATE {
    take:
    ch_bam_bai  // [ [ meta ], bam , bai ]
    fasta       // [ [ meta ], fasta ]
    fasta_fai   // [ [ meta ], fasta_fai ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files  = Channel.empty()
    ch_refs = fasta.join(fasta_fai)
        .map{
            // Create additional map containing only meta.reference for combining samples and intervals
            meta, fasta, fai ->
            meta2 = [:]
            meta2.reference = meta.id
            [ meta2, meta, fasta, fai ]
        }

    // Create genomic regions file for splitting the bam before deduplication
    BUILD_INTERVALS(fasta_fai)
    ch_versions   = ch_versions.mix( BUILD_INTERVALS.out.versions )

    // Ensure input bam matches the regions file
    ch_bam_for_split = ch_bam_bai
        .map{
            // Create additional map containing only meta.reference for combining with intervals for reference
            meta, bam, bai ->
            meta2 = [:]
            meta2.reference = meta.reference
            [ meta2, meta, bam, bai ]
        }
        // .dump(tag:"bams_extra_meta", pretty: true)
        .combine(
            by: 0,
            BUILD_INTERVALS.out.bed.map{
                // Rename reference meta.id to meta.reference to allow combining with bams of specific reference
                meta, bed ->
                meta2 = [:]
                meta2.reference = meta.id
                [ meta2, bed ]
            }
        )
        .map{
            ignore_me, meta, bam, bai, regions ->
            [ meta, bam, bai, regions ]
        }

    //Split input bam by region
    BAM_SPLIT_BY_REGION(ch_bam_for_split)
    ch_versions   = ch_versions.mix( BAM_SPLIT_BY_REGION.out.versions )

    if ( params.deduplication_tool == 'markduplicates' ) {

        ch_markduplicates_input = BAM_SPLIT_BY_REGION.out.bam_bai
            .map{
                // Create additional map containing only meta.reference for combining with intervals for reference
                meta, bam, bai ->
                meta2 = [:]
                meta2.reference = meta.reference
                [ meta2, meta, bam, bai ]
            }
            .combine(
                by:0,
                ch_refs
            )
            .multiMap{
                ignore_me, meta, bam, bai, meta2, fasta, fasta_fai ->
                bam: [ meta, bam ]
                fasta: fasta
                fasta_fai: fasta_fai
            }
            // Dedup each bam
            PICARD_MARKDUPLICATES(
                ch_markduplicates_input.bam,
                ch_markduplicates_input.fasta,
                ch_markduplicates_input.fasta_fai
            )
            ch_versions   = ch_versions.mix( PICARD_MARKDUPLICATES.out.versions.first() )

            // ch_dedup_multiqc_files        = ch_multiqc_files.mix( PICARD_MARKDUPLICATES.out.metrics )
            ch_dedupped_region_bam        = PICARD_MARKDUPLICATES.out.bam

    } else if ( params.deduplication_tool == "dedup" ) {
        // Other deduppers go here.
        ch_dedup_input = BAM_SPLIT_BY_REGION.out.bam_bai
            .map{
                meta, bam, bai ->
                [ meta, bam ]
            }

        DEDUP( ch_dedup_input )
        ch_versions = ch_versions.mix( DEDUP.out.versions.first() )

        ch_dedupped_region_bam = DEDUP.out.bam

    }

ch_input_for_samtools_merge = ch_dedupped_region_bam
    .map{
        meta, bam ->
        meta2 = meta.clone().findAll{ it.key != 'genomic_region' }
        [ meta2, bam ]
    }
    .groupTuple()
    .map{
        meta, bam ->
        ref_meta = [:]
        ref_meta.reference = meta.reference
        [ ref_meta, meta, bam ]
    }
    .combine(
        by:0,
        ch_refs
    )
    .multiMap{
        // bam here is a list of bams
        ignore_me, meta, bam, meta2, fasta, fasta_fai ->
        bam:        [ meta, bam ]
        fasta:      fasta
        fasta_fai:  fasta_fai
    }

    // Merge the bams for each region into one bam
    SAMTOOLS_MERGE(
        ch_input_for_samtools_merge.bam,
        ch_input_for_samtools_merge.fasta,
        ch_input_for_samtools_merge.fasta_fai
    )
    ch_versions   = ch_versions.mix( SAMTOOLS_MERGE.out.versions )


    // Sort the merged bam and index
    SAMTOOLS_SORT ( SAMTOOLS_MERGE.out.bam )
    ch_dedup_bam  = SAMTOOLS_SORT.out.bam
    ch_versions   = ch_versions.mix( SAMTOOLS_SORT.out.versions )

    SAMTOOLS_INDEX ( ch_dedup_bam )
    ch_dedup_bai  = SAMTOOLS_INDEX.out.bai
    ch_versions   = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    // Finally run flagstat on the dedupped bam
    SAMTOOLS_FLAGSTAT(
        ch_dedup_bam.join( ch_dedup_bai )
    )
    ch_dedup_flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
    ch_versions       = ch_versions.mix( SAMTOOLS_FLAGSTAT.out.versions )

    emit:
    bam        = ch_dedup_bam
    bai        = ch_dedup_bai
    versions   = ch_versions
    // mqc        = ch_dedup_multiqc_files // The deduplication rate for the overall bam will come from endorspy, so nothing to pass along here.
}
