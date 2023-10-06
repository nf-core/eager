//
// Run sex determine
//

include { SAMTOOLS_DEPTH } from '../modules/nf-core/samtools/depth/main'
include { SEXDETERRMINE } from '../modules/nf-core/sexdeterrmine/main'

workflow RUN_SEXDETERMINE {

    take:
    samtoolsdepth_input // channel: [ val(meta1), [ bam ], [ meta2 ] [ intervals_bed ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.run_sexdeterrmine ) {

        // Generate sex_determine input

        SAMTOOLS_DEPTH(samtoolsdepth_input).out.tsv
        ch_sex_determine_input = SAMTOOLS_DEPTH.out.tsv
        ch_multiqc_files       = ch_multiqc_files.mix( SAMTOOLS_DEPTH.out.mqc )
        ch_versions            = ch_versions.mix( SAMTOOLS_DEPTH.out.versions )

        // Run sex determination with samtools depth input
        SEXDETERRMINE(ch_sex_determine_input)
        ch_coverages = SEXDETERRMINE.out.tsv
        ch_multiqc_files      = ch_multiqc_files.mix( SEXDETERRMINE.out.json )
        ch_versions           = ch_versions.mix( SEXDETERRMINE.out.versions )
    }

    emit:
    covergages         = ch_coverages
    mqc                = ch_multiqc_files
    versions           = ch_versions
}
