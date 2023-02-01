//
// Filter BAMs for mapping quality, length, unmapped etc.
//

include { FILTER_BAM_FRAGMENT_LENGTH                      } from '../../modules/local/filter_bam_fragment_length'
include { SAMTOOLS_VIEW                                   } from '../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_FILTERED } from '../../modules/nf-core/samtools/flagstat/main'
// include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED       } from '../../modules/nf-core/samtools/fastq/main'
// include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../../modules/nf-core/samtools/fastq/main'

workflow FILTER_BAM {

    take:
    bam // [ [meta], [bam], [bai] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    //
    // GENOMICS BAM GENERATION
    //

   // Solution to the Wang-Schiffels Problem: Generate length filtered BAM of mapped reads
    // - While it would be slightly more efficient to do this prior length filter
    //   as samtools will be faster, so would reduce data size into slower script,
    //   length_filter doesn't include indexing

    if ( params.bamfiltering_minreadlength != 0  ) {
        FILTER_BAM_FRAGMENT_LENGTH ( bam )
        ch_bam_for_qualityfilter = FILTER_BAM_FRAGMENT_LENGTH.out.bam
        ch_versions = ch_versions.mix( FILTER_BAM_FRAGMENT_LENGTH.out.versions.first() )
    } else {
        ch_bam_for_qualityfilter = bam
    }

    // Generate BAM file of quality filtered and mapped-only reads,
    // optionally retaining unmapped reads, defined in modules.config
    // TODO: add reference fasta in second channel?
    // TODO: document by unmapped reads always removed
    SAMTOOLS_VIEW ( ch_bam_for_qualityfilter, [], [] )
    ch_bam_for_genomics = SAMTOOLS_VIEW.out.bam
    ch_versions = ch_versions.mix( SAMTOOLS_VIEW.out.versions.first() )

    // Only run if we actually remove mapped reads
    if ( params.bamfiltering_mappingquality != 0 || params.bamfiltering_minreadlength != 0  ) {
        SAMTOOLS_FLAGSTAT_FILTERED ( ch_bam_for_genomics )
        ch_versions      = ch_versions.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.flagstat )
    }

//     //
//     // Metagenomics FASTQ generation for metagenomics (or just generation)
//     // - FASTQ generation is now separate from BAM filtering -
//     //    no length/quality filtering applies to metagenomic bam
//     //

//     // Generate unmapped bam (no additional filtering) if the unmapped bam OR unmapped for metagneomics selected
//     if ( params.bamfiltering_generateunmappedfastq || ( run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) ) {
//         SAMTOOLS_FASTQ_UNMAPPED ( bam )
//         ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_UNMAPPED.out.versions.first() )
//     }

//     // Solution to the Andrades Valtueña-Light Problem: mapped bam for metagenomics (with options for quality- and length filtered)
//     // TODO: document that BAMtoFASTQs never go through read length filtering! Only genomic BAM
//     if ( params.bamfiltering_generatemappedfastq ||  ( params.run_metagenomicscreening && ( params.metagenomicscreening_input == 'mapped' || params.metagenomicscreening_input == 'all' ) ) ) {
//         SAMTOOLS_FASTQ_MAPPED ( bam )
//         ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_MAPPED.out.versions.first() )
//     }

//     // Routing for metagenomic screeninng
//     if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) {
//         ch_bam_for_metagenomics = SAMTOOLS_FASTQ_UNMAPPED.out.fastq
//     } else if ( params.run_metagenomicscreening && ( params.metagenomicscreening_input == 'mapped' : params.metagenomicscreening_input == 'all' )) {
//         ch_bam_for_metagenomics = SAMTOOLS_FASTQ_MAPPED.out.fastq
//     } else if ( !params.run_metagenomicscreening ) {
//         ch_bamfiltered_for_metagenomics = Channel.empty()
//     }

    emit:
        genomics         = Channel.empty() // ch_bam_for_genomics
        metagenomics     = Channel.empty() // ch_bam_for_metagenomics
        versions         = ch_versions
        mqc              = ch_multiqc_files

}
