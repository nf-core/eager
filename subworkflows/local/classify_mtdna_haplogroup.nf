//
// Run classify mtdna haplogroup
//

include { addNewMetaFromAttributes                          } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

include { HAPLOGREP3_CLASSIFY as HAPLOGREP3_CLASSIFY_MTDNA  } from '../../modules/nf-core/haplogrep3/classify/main'

workflow CLASSIFY_MTDNA_HAPLOGROUP {

    take:
    mtdna_haplogroup_bam // channel: [ val(meta1), [ bam ], [ bai ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()
    ch_haplogroups    = Channel.empty()

    if ( params.run_mtdna_haplogroup ) {

        ch_input_haplogrep3 = mtdna_haplogroup_bam
            .map {
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .map { meta, bam, bai ->
                [meta, bam]
            }

        HAPLOGREP3_CLASSIFY_MTDNA(ch_input_haplogrep3)
        ch_haplogroups      = HAPLOGREP3_CLASSIFY_MTDNA.out.txt
        ch_versions         = ch_versions.mix(HAPLOGREP3_CLASSIFY_MTDNA.out.versions)
    }

    emit:
    haplogroups        = ch_haplogroups    // channel: [ val(meta), path("*.txt") ]
    versions           = ch_versions       // channel: path(versions.yml)
}
