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

    /* possible conditions


    ✅ filtering mapped bam
    ✅ filtering mapped fastq (a.k.a. metagenomic mapped)
    ✅ filtering mapped both (Confiugred in publishDir !)

    filtering filtering && unmapped bam
    filtering filtering && unmapped fastq (a.k.a. metagenomic unmapped)
    filtering filtering && unmapped both

    filtering mapped/unmapped fastq (a.k.a. metagenomic unmapped)
    filtering mapped/unmapped bam
    filtering mapped/unmapped both

    withfiltering length filter
    nofiltering length filter

    */

    // TODO: add a check in eager.nf: if ( params.run_metagenomicscreening && !params.run_bamfiltering ) { exit 1, "If you wish to run metagenomic screening you must turn on --run_bamfiltering" }
    // TODO  add check - if want unmapped bam/fastq or metagenomic used,  mustrestain umapped
    // TODO: add a check in eager.nf: if ( params.run_metagenomicscreening && !params.bamfiltering_generatefastq ) { exit 1, "If you wish to run metagenomic screening yo umust include --bamfiltering_generatefastq" }
    // TODO: add check in eager.nf: if ( params.run_metagenomicscreening && --metagenomicscreening_input 'unmapped' && !params.bamfiltering_generateunmapped ) { exit 1, "If you want unmapped input, you need to also generate this with `--bamfiltering_generatedunmmaped`  }
    // TODO: write length filter module!
    // TODO: add all params to schema



    /* ATTEMPT 3

        Qs:
            - Do you want to run quality filtering or not?
            - Do you want to do length filtering or not?

            - Do you want mapped reads only?
            - Do you want unmapped reads only?
            - Do you want unmapped + mapped reads separately?
            - Do you want unmapped + mapped reads together

            - Do you want FASTQ files of your mapped BAM file
            - Do you want FASTQ files of your unmapped BAM file

            - Do you want to run metagenomic screening on your mapped reads (fastq only)
            - Do you want to run metagenomic screening on your unmapped reads (fastq only)

    Order: mapped (+ quality filtering?) / unmapped separateion first ; length filtering ;

    SAMTOOLS_FASTQ -> if/else mapped and/or unmapped

    MAPPED ONLY: length filtering
    MAPPED ONLY: quality filtering

    -> issue what if people want quality filter the mapped reads to metagenomics -> flip view and fastq?


    */

    // GENERAL TODO: VERSIONS! MULTIQC

    // Generate unmapped bam (no additional filtering) if the unmapped bam OR unmapped for metagneomics selected
    if ( params.bamfiltering_generateunmappedfastq || ( run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) ) {
        SAMTOOLS_FASTQ_UNMAPPED ( bam )
    }

    // Solution to the Wang-Schiffels Problem: Generate length filtered BAM of mapped reads
    if ( params.bamfiltering_readlength != 0  ) {
        ch_bam_for_qualityfilter = LENGTH_FILTER ( ch_bam_for_lengthfilter )
        ch_versions = ch_versions.mix( LENGTH_FILTER.out.versions )
    } else {
        ch_bam_for_qualityfilter = SAMTOOLS_MAPPED.out.bam
    }

    // Generate BAM file of quality filtered reads
    if ( params.bamfiltering_mappingquality != 0 ) {
        ch_bam_for_genomics = SAMTOOLS_VIEW ( ch_bam_for_qualityfilter )
    } else {
        ch_bam_for_genomics = ch_bam_for_qualityfilter
    }

    // Solution to the Andrades-Valtueña-Light Problem: mapped bam for metagenomics (with options for quality- and length filtered)
    // TODO: document cannot send 'raw' mapped bam to metagenomics, must be same as genomics
    // TODO: make optional in modules config if we remove unmapped reads from the i.e, bamfiltering_retainunmappedinbam... TODO: ADD FOLLOWING CONDITOIN TO MODULES.CONF
    if ( params.bamfiltering_generatemappedfastq ||  ( params.run_metagenomicscreening && ( params.metagenomicscreening_input == 'mapped' || params.metagenomicscreening_input == 'all' ) ) {
        SAMTOOLS_FASTQ_MAPPED ( bam )
    }

    // Routing for metagenomic screeninng
    if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'unmapped' ) {
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ_UNMAPPED.out.fastq
    } else if ( params.run_metagenomicscreening && params.metagenomicscreening_input == 'mapped' ) {
        ch_bam_for_metagenomics = SAMTOOLS_FASTQ_MAPPED.out.fastq
    } else if ( !params.run_metagenomicscreening ) {
        ch_bamfiltered_for_metagenomics = Channel.empty()
    }

    if ( params.bamfiltering_mappingquality != 0 || bamfiltering_minlength != 0  ) {
        SAMTOOLS_STATS_FILTERED ( ch_bam_for_genomics )
    }


    emit:
    genomics_bam     = ch_bam_for_genomics
    genomics_bai     = ch_bai_for_genomics // TODO Join here or later?
    metagenomics     = ch_bam_for_metagenomics

}
