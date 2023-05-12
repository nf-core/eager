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
    if ( params.run_contamination_angsd ) {
        angsd_input_hapmap = hapmap_input
        .map {
            // Create additional map containing only meta.id for combining samples and hapmap
            meta, hapmap ->
                meta2 = [:]
                meta2.reference = meta.id
            [ meta2, meta, hapmap ]
        }
        angsd_input = contamination_input
        .map {
            // Create additional map containing only meta.reference for combining samples and hapmap
            meta, bam, bai ->
                meta2 = [:]
                meta2.reference = meta.reference
            [ meta2, meta, bam, bai ]
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

        PRINT_CONTAMINATION_ANGSD( ch_angsd_output)
        ch_contam_angsd_print = PRINT_CONTAMINATION_ANGSD.out.contam_angsd_txt
        ch_contam_mqc         = PRINT_CONTAMINATION_ANGSD.out.contam_angsd_multiqc
    }

    emit:
    contam_angsd_print = ch_contam_angsd_print
    contam_mqc         = ch_contam_mqc
    versions           = ch_versions

}
