//
// Run sex determine
//

include { SAMTOOLS_DEPTH } from '../../modules/nf-core/samtools/depth/main'
include { SEXDETERRMINE } from '../../modules/nf-core/sexdeterrmine/main'

workflow RUN_SEXDETERMINE {

    take:
    // samtoolsdepth_input // channel: [ val(meta1), [ bam ], [ meta2 ] [ intervals_bed ] ]
    sexdeterrmine_bam // channel: [ val(meta1), [ bam ] ]
    sexdeterrmine_bed // channel: [ [ meta2 ] [ intervals_bed ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()


    // Generate sex_determine input
    // Break here the input not in the workflow

    ch_samtoolsdepth_input = sexdeterrmine_bam
                            .combine(sexdeterrmine_bed)
                            .multiMap{
                            meta, bam, meta2, bed ->
                            bam       : [ meta, bam ]
                            intervals : [ meta2, bed ]
                                    }

    SAMTOOLS_DEPTH(ch_samtoolsdepth_input)
    ch_sex_determine_input = SAMTOOLS_DEPTH.out.tsv
    ch_versions            = ch_versions.mix( SAMTOOLS_DEPTH.out.versions )

    // Run sex determination with samtools depth input
    // Does it need a sample list file to run?
    SEXDETERRMINE(ch_sex_determine_input, sexdeterrmine_list)
    ch_coverages          = SEXDETERRMINE.out.tsv
    ch_multiqc_files      = ch_multiqc_files.mix( SEXDETERRMINE.out.json )
    ch_versions           = ch_versions.mix( SEXDETERRMINE.out.versions )

    emit:
    coverages          = ch_coverages // channel: [ val(meta), path("*.tsv") ]
    mqc                = ch_multiqc_files // channel: [ val(meta), path("*.json") ]
    versions           = ch_versions // channel: path(versions.yml)
}
