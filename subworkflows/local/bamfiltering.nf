//
// Filter BAMs for mapping quality, length, unmapped etc.
//

include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_UNMAPPED       } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_VIEW                                   } from '../../modules/nf-core/samtools/view/main'
include { LENGTH_FILTER                                   } from '../../modules/local/length_filter'
include { SAMTOOLS_FASTQ as SAMTOOLS_FASTQ_MAPPED         } from '../../modules/nf-core/samtools/fastq/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_FILTERED } from '../../modules/nf-core/samtools/flagstat/main'

workflow BAM_FILTER {

    take:
    bam // [ [meta], [bai] ]

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

    if ( params.bamfiltering_readlength != 0  ) {
        ch_bam_for_lengthfilter = LENGTH_FILTER ( ch_bam_for_lengthfilter )
        ch_versions = ch_versions.mix( LENGTH_FILTER.out.versions )
    } else {
        ch_bam_for_lengthfilter = ch_bam_for_lengthfilter
    }

    // Generate BAM file of quality filtered reads
    if ( params.bamfiltering_mappingquality != 0 ) {
        SAMTOOLS_VIEW ( bam )
        ch_bam_for_genomics = SAMTOOLS_VIEW.out.bam
        ch_versions = ch_versions.mix( SAMTOOLS_VIEW.out.versions )
    } else {
        ch_bam_for_genomics = bam
    }

    if ( params.bamfiltering_mappingquality != 0 || bamfiltering_minlength != 0  ) {
        SAMTOOLS_FLAGSTAT_FILTERED ( ch_bam_for_genomics )
        ch_versions      = ch_versions.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.versions )
        ch_multiqc_files = ch_multiqc.mix( SAMTOOLS_FLAGSTAT_FILTERED.out.flagstat )
    }

    //
    // Metagenomics FASTQ generation for metagenomics (or just generation)
    // - FASTQ generation is now separate from BAM filtering -
    //    no length/quality filtering applies to metagenomic bam
    //

    // Generate unmapped bam (no additional filtering) if the unmapped bam OR unmapped for metagneomics selected
    if ( params.bamfiltering_generateunmappedfastq || ( run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) ) {
        SAMTOOLS_FASTQ_UNMAPPED ( bam )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_UNMAPPED.out.versions )
    }

    // Solution to the Andrades Valtueña-Light Problem: mapped bam for metagenomics (with options for quality- and length filtered)
    // TODO: document that BAMtoFASTQs never go through read length filtering! Only genomic BAM
    if ( params.bamfiltering_generatemappedfastq ||  ( params.run_metagenomicscreening && ( params.metagenomicscreening_input == 'mapped' || params.metagenomicscreening_input == 'all' ) ) ) {
        SAMTOOLS_FASTQ_MAPPED ( bam )
        ch_versions = ch_versions.mix( SAMTOOLS_FASTQ_MAPPED.out.versions )
    }

    // Routing for metagenomic screeninng
    if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) {
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ_UNMAPPED.out.fastq
    } else if ( params.run_metagenomicscreening && ( params.metagenomicscreening_input == 'mapped' : params.metagenomicscreening_input == 'all' )) {
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ_MAPPED.out.fastq
    } else if ( !params.run_metagenomicscreening ) {
        ch_bamfiltered_for_metagenomics = Channel.empty()
    }

    emit:
    genomics_bam     = ch_bam_for_genomics
    metagenomics     = ch_bam_for_metagenomics

}
