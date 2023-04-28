//
// Contamination estimation
//

include { BAM_DOCOUNTS_CONTAMINATION_ANGSD } from '../subworkflows/nf-core/bam_docounts_contamination_angsd/main'

workflow CONTAMINATION_ESTIMATION {

    take:
    dedup_bam
    dedup_bai
    hapmap_file

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    BAM_DOCOUNTS_CONTAMINATION_ANGSD ( ch_dedup_bam, ch_dedup_bai, ch_hapmap_file )
    ch_versions     = ch_versions.mix ( BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.versions.first() )
    ch_angsd_contam = BAM_DOCOUNTS_CONTAMINATION_ANGSD.out.txt

    emit:
    angsd_contam      = ch_angsd_contam
    versions          = ch_versions
    mqc               = ch_multiqc_files

}
