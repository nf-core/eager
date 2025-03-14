//
// Filter BAMs for mapping quality, length, unmapped etc.
//

include { FILTER_BAM_FRAGMENT_LENGTH                             } from '../../modules/local/filter_bam_fragment_length'
include { SAMTOOLS_VIEW  as SAMTOOLS_VIEW_BAM_FILTERING          } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_INDEX as SAMTOOLS_LENGTH_FILTER_INDEX         } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_FILTER_INDEX                } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_FILTERED        } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_METAGENOMICS          } from '../../modules/nf-core/samtools/fastq/main'
include { CAT_FASTQ as CAT_FASTQ_METAGENOMICS                    } from '../../modules/nf-core/cat/fastq'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_SAVEBAMFILTERINGREADS } from '../../modules/nf-core/samtools/fastq/main'


workflow FILTER_BAM {

    take:
    bam // [ [meta], [bam], [bai/csi] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()
    ch_flagstats_file = Channel.empty()

    //
    // GENOMICS BAM GENERATION
    //

    // Solution to the Wang-Schiffels Problem: Generate length filtered BAM of mapped reads
    // - While it would be slightly more efficient to do this prior length filter
    //   as samtools will be faster, so would reduce data size into slower script,
    //   length_filter doesn't include indexing

    if ( params.bamfiltering_minreadlength != 0  ) {
        FILTER_BAM_FRAGMENT_LENGTH ( bam )
        ch_versions = ch_versions.mix( FILTER_BAM_FRAGMENT_LENGTH.out.versions.first() )

        SAMTOOLS_LENGTH_FILTER_INDEX ( FILTER_BAM_FRAGMENT_LENGTH.out.bam )
        ch_length_filtered_index = params.fasta_largeref ? SAMTOOLS_LENGTH_FILTER_INDEX.out.csi : SAMTOOLS_LENGTH_FILTER_INDEX.out.bai
        ch_versions = ch_versions.mix( SAMTOOLS_LENGTH_FILTER_INDEX.out.versions.first() )

        ch_bam_for_qualityfilter = FILTER_BAM_FRAGMENT_LENGTH.out.bam.join( ch_length_filtered_index )

    } else {
        ch_bam_for_qualityfilter = bam
    }

    // Generate BAM file of quality filtered and mapped-only reads,
    // optionally retaining unmapped reads, defined in modules.config

    SAMTOOLS_VIEW_BAM_FILTERING ( ch_bam_for_qualityfilter, [[], []], [] ) // fasta isn't needed until we support CRAM
    ch_versions = ch_versions.mix( SAMTOOLS_VIEW_BAM_FILTERING.out.versions.first() )

    SAMTOOLS_FILTER_INDEX ( SAMTOOLS_VIEW_BAM_FILTERING.out.bam )
    ch_filtered_bam_index = params.fasta_largeref ? SAMTOOLS_FILTER_INDEX.out.csi : SAMTOOLS_FILTER_INDEX.out.bai
    ch_versions = ch_versions.mix( SAMTOOLS_FILTER_INDEX.out.versions.first() )

    ch_bam_for_genomics = SAMTOOLS_VIEW_BAM_FILTERING.out.bam.join( ch_filtered_bam_index )

    // Only run if we actually remove mapped reads
    if ( params.bamfiltering_mappingquality != 0 || params.bamfiltering_minreadlength != 0  ) {
        SAMTOOLS_FLAGSTAT_FILTERED ( ch_bam_for_genomics )
        ch_versions      = ch_versions.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.flagstat )
        ch_flagstats_file = ch_flagstats_file.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.flagstat )
    }

    // FASTQ generation for saving reads (INDEPENDENT of metagenomics)
    // Output of reads is determined parameters are set by parameters for bamfiltering above (see SAMTOOLS_VIEW_BAM_FILTERING module config)
    // Extension based on #945 possible here duirng hackathon
    // Only possible with run_bamfiltering parameter set
    if ( params.bamfiltering_generatefastq ) {
        SAMTOOLS_FASTQ_SAVEBAMFILTERINGREADS ( SAMTOOLS_VIEW_BAM_FILTERING.out.bam, false )
        ch_versions = ch_versions.mix ( SAMTOOLS_FASTQ_SAVEBAMFILTERINGREADS.out.versions.first() )
    }

    // FASTQ generation for metagenomics (OR just generation for saving mapped/unmapped reads)
    // - FASTQ generation is now separate from BAM filtering -
    //    No length/quality filtering applies to metagenomic bam files (could be extension)
    //    All bam -> fastq filtering options (-F 4, -f 4 or none) will be dealt with within modules.config

    if ( params.run_metagenomics )  {
        // Execute fastq generation on original bam mapping (independent of above bamfiltering)
        SAMTOOLS_FASTQ_METAGENOMICS ( bam.map{[ it[0], it[1] ]}, false )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_METAGENOMICS.out.versions.first() )

        if ( params.preprocessing_skippairmerging ) {
            // Splitting of paired vs single end data
            ch_paired_fastq_for_cat_metagenomics = SAMTOOLS_FASTQ_METAGENOMICS.out.fastq.filter { !it[0].single_end }
            ch_single_fastq_for_cat_metagenomics = SAMTOOLS_FASTQ_METAGENOMICS.out.fastq
                                        .mix(SAMTOOLS_FASTQ_METAGENOMICS.out.singleton)
                                        .mix(SAMTOOLS_FASTQ_METAGENOMICS.out.other)
                                        .groupTuple()
                                        .filter{ it[0].single_end }
            CAT_FASTQ_METAGENOMICS ( ch_single_fastq_for_cat_metagenomics )
            ch_fastq_for_metagenomics = CAT_FASTQ_METAGENOMICS.out.reads.mix( ch_paired_fastq_for_cat_metagenomics )
            ch_versions = ch_versions.mix( CAT_FASTQ_METAGENOMICS.out.versions.first() )
        }
        else {
            ch_fastq_for_metagenomics = SAMTOOLS_FASTQ_METAGENOMICS.out.other
        }
    } else {
        ch_fastq_for_metagenomics = Channel.empty()
    }

    emit:
    genomics         = ch_bam_for_genomics
    metagenomics     = ch_fastq_for_metagenomics
    versions         = ch_versions
    flagstat         = ch_flagstats_file    // [ [ meta ], stats ]
    mqc              = ch_multiqc_files

}
