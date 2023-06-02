//
// Estimate contamination
//

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../../subworkflows/nf-core/bam_docounts_contamination_angsd/main'
include { PRINT_CONTAMINATION_ANGSD        } from '../../modules/local/print_contamination_angsd'
include { AUTHENTICT_DEAM2CONT             } from '../../modules/nf-core/authentict/deam2cont/main'

workflow ESTIMATE_CONTAMINATION {

    take:
    contamination_input // channel: [ val(meta), [ bam ], [ bai ] ]
    hapmap_input        // channel: [ val(meta), [ hapmap_file ] ]
    position_input      // channel: [ val(meta), [ position_file ] ]

    main:
    ch_versions           = Channel.empty()
    ch_multiqc_files      = Channel.empty()
    ch_contam_angsd_print = Channel.empty()
    ch_contam_authentict  = Channel.empty()


// fix unexpected input/syntax error
    if ( params.run_contamination_estimation_angsd ) {
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
                hapmap: [ meta2, hapmap ]
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

    if ( params.run_contamination_estimation_authentict ) {
        // TODO: only run on single-stranded and non-UDG data
        if ( position_input ) {
            authnetict_input_position = position_input
                .map {
                    // Create additional map containing only meta.id for combining samples and hapmap
                    meta, position ->
                        meta2 = [:]
                        meta2.reference = meta.id
                    [ meta2, meta, position ]
                }
            authentict_input = contamination_input
                .map {
                // Create additional map containing only meta.reference for combining samples and hapmap
                    meta, bam, bai ->
                        meta2 = [:]
                        meta2.reference = meta.reference
                    [ meta2, meta, bam, bai ]
                }
                .combine(
                    by: 0,
                    authentict_input_position
                )
                .multiMap {
                    ignore_me, meta, bam, bai, meta2, position ->
                    bam:    [ meta, bam ]
                    position: [ meta2, position ]
                }
            authentict_input_bam = authentict_input.bam
            authentict_input_pos = authentict_input.position
        } else {
            authentict_input_bam = contamination_input
                .map {
                    meta, bam, bai ->
                    [ meta, bam ]
                }
            authentict_input_pos = Channel.empty()
        }

        AUTHENTICT_DEAM2CONT( authentict_input_bam, [[], []], authentict_input_pos.ifEmpty( [[], []] ) )

        ch_contam_authentict  = AUTHENTICT_DEAM2CONT.out.txt
        //ch_multiqc_files      = ch_multiqc_files.mix( AUTHENTICT_DEAM2CONT.out.json )
        ch_versions           = ch_versions.mix( AUTHENTICT_DEAM2CONT.out.versions )
    }

    emit:
    contam_angsd_print = ch_contam_angsd_print
    contam_authentict  = ch_contam_authentict.dump()
    mqc                = ch_multiqc_files
    versions           = ch_versions
}
