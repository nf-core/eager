//
// Calculate damage profile
//

include { DAMAGEPROFILER } from '../../modules/nf-core/damageprofiler/main'

workflow CALCULATE_DAMAGE {
    take:
    ch_bam_bai //  [ [ meta ], [ bam ], [ bai ] ]
    fasta       // [ [ meta ], fasta ]
    fasta_fai   // [ [ meta ], fasta_fai ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    //From deduplicate.nf
    ch_refs = fasta.join(fasta_fai)
        .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
            WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
        }
    // Ensure input bam matches the regions file
    ch_damageprofiler_input = ch_bam_bai
        .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
            WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
        }
        .combine(
                by:0,
                ch_refs
            )
            .multiMap{
                ignore_me, meta, bam, bai, meta2, fasta, fasta_fai ->
                bam: [ meta, bam ]
                fasta: fasta
                fasta_fai: fasta_fai
            }
    // Calculate damage
    DAMAGEPROFILER(
                ch_damageprofiler_input.bam,
                ch_damageprofiler_input.fasta,
                ch_damageprofiler_input.fasta_fai,
                []
            )
            ch_versions       = ch_versions.mix( DAMAGEPROFILER.out.versions.first() )
            ch_multiqc_files  = ch_multiqc_files.mix( DAMAGEPROFILER.out.results )

    emit:
    versions = ch_versions          // channel: [ versions.yml ]
    mqc      = ch_multiqc_files
}
