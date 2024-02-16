//
// Run sex determine
//

include { SAMTOOLS_DEPTH as SAMTOOLS_DEPTH_SEXDETERMINE } from '../../modules/nf-core/samtools/depth/main'
include { SEXDETERRMINE                                 } from '../../modules/nf-core/sexdeterrmine/main'

workflow RUN_SEXDETERMINE {

    take:
    // samtoolsdepth_input // channel: [ val(meta1), [ bam ], [ meta2 ] [ intervals_bed ] ]
    sexdeterrmine_bam // channel: [ val(meta1), [ bam ] ]
    sexdeterrmine_bed // channel: [ [ meta2 ] [ intervals_bed ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.run_sexdeterrmine ) {
        // Generate sex_determine input
        // Break here the input not in the workflow

        ch_bed = sexdeterrmine_bed
            .map {
                // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }

        ch_input_sexdetermine = sexdeterrmine_bam
            .map {
                // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_bed, by: 0 ) // [ [combine_meta], [meta], bam, bai, [ref_meta], bed ]

        ch_samtoolsdepth_input = ch_input_sexdetermine
                .multiMap {
                    ignore_me, meta, bam, bai, ref_meta, bed ->
                        bam: [ meta, bam ]
                        bed: [ ref_meta, bed ]
                }

        SAMTOOLS_DEPTH_SEXDETERMINE(ch_samtoolsdepth_input.bam, ch_samtoolsdepth_input.bed)
        ch_sex_determine_input = SAMTOOLS_DEPTH_SEXDETERMINE.out.tsv
        ch_versions            = ch_versions.mix( SAMTOOLS_DEPTH_SEXDETERMINE.out.versions )

        // Run sex determination with samtools depth input
        SEXDETERRMINE(ch_sex_determine_input, [])
        ch_coverages          = SEXDETERRMINE.out.tsv
        ch_multiqc_files      = ch_multiqc_files.mix( SEXDETERRMINE.out.json )
        ch_versions           = ch_versions.mix( SEXDETERRMINE.out.versions )

        }

    emit:
    coverages          = ch_coverages // channel: [ val(meta), path("*.tsv") ]
    mqc                = ch_multiqc_files // channel: [ val(meta), path("*.json") ]
    versions           = ch_versions // channel: path(versions.yml)
}
