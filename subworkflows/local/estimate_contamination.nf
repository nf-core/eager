//
// Estimate contamination
//

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../../subworkflows/nf-core/bam_docounts_contamination_angsd/main'
include { PRINT_CONTAMINATION_ANGSD        } from '../../modules/local/print_contamination_angsd'

workflow ESTIMATE_CONTAMINATION {

    take:
    contamination_input // channel: [ val(meta), [ bam ], [ bai ] ]
    hapmap_input        // channel: [ val(meta), [ hapmap_file ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

// fix unexpected input/syntax error
    if ( params.run_contamination_estimation_angsd ) {
        angsd_input_hapmap = hapmap_input
            .map {
                // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }
        angsd_input = contamination_input
            .map {
                // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine(
                by: 0,
                angsd_input_hapmap
            )
            .multiMap {
                ignore_me, meta, bam, bai, meta2, hapmap ->
                bam:    [ meta, bam ]
                bai:    [ meta, bai ]
                hapmap: [ meta, hapmap ]
            }

        BAM_DOCOUNTS_CONTAMINATION_ANGSD( angsd_input.bam, angsd_input.bai, angsd_input.hapmap )

        ch_versions     = ch_versions.mix( BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.versions )
        ch_angsd_contam = BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.txt

        ch_angsd_output = ch_angsd_contam.collect()

        PRINT_CONTAMINATION_ANGSD( ch_angsd_output )
        ch_contam_angsd_print = PRINT_CONTAMINATION_ANGSD.out.txt
        ch_multiqc_files      = ch_multiqc_files.mix( PRINT_CONTAMINATION_ANGSD.out.json )
        ch_versions           = ch_versions.mix( PRINT_CONTAMINATION_ANGSD.out.versions )
    }

    emit:
    contam_angsd_print = ch_contam_angsd_print
    mqc                = ch_multiqc_files
    versions           = ch_versions
}
