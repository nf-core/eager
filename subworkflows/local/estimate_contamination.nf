//
// Estimate contamination
//

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../../subworkflows/nf-core/bam_docounts_contamination_angsd/main'

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
            // Create additional map containing only meta.reference for combining samples and hapmap
            meta, hapmap ->
                meta2 = [:]
                meta2.reference = meta.id
            [ meta2, meta, hapmap ]
        }// works fine [DUMP] [['reference':'hs37d5_chr21'], ['id':'hs37d5_chr21'], /nf-core/test-datasets/blob/9a51c46c2ee9e46c906a3b136275a06155f3cdbc/data/delete_me/angsd/HapMapChrX.gz]

        angsd_input = contamination_input
        .map {
            // Create additional map containing only meta.reference for combining samples and hapmap
            meta, bam, bai ->
                meta2 = [:]
                meta2.reference = meta.reference
            [ meta2, meta, bam, bai ]
        }.dump() // no reference in meta [DUMP] [['reference':null], ['id':'JK2067_JK2067', 'sample_id':'JK2067', 'library_id':'JK2067', 'strandedness':'double', 'damage_treatment':'full'], /Users/carlhoff/Documents/git/eager/test/work/c3/43a2c63895c0f069d6a6ffd489ef68/JK2067_JK2067_JK2067_filtered.bam, /Users/carlhoff/Documents/git/eager/test/work/71/f02947a9b4a7a4f1eb62e7927decdd/JK2067_JK2067_JK2067_filtered.bam.bai]
        .combine(
            by: 0,
            angsd_input_hapmap
        )
        .multiMap {
            meta, bam, bai, hapmap ->
            bam:    [ meta, bam ]
            bai:    [ meta, bai ]
            hapmap: [ meta, hapmap ]
        }

        BAM_DOCOUNTS_CONTAMINATION_ANGSD( angsd_input.bam, angsd_input.bai, angsd_input.hapmap )

        ch_versions     = ch_versions.mix( BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.versions.first() )
        ch_angsd_contam = BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.txt
    }

    emit:
    angsd_contam      = ch_angsd_contam
    versions          = ch_versions
    mqc               = ch_multiqc_files

}
