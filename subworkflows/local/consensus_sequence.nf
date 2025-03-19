//
// Produce consensus sequences from reads aligned to a reference some targetted to aDNA
//

include { MULTIVCFANALYZER         } from '../../modules/nf-core/multivcfanalyzer'
include { addNewMetaFromAttributes } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

workflow CONSENSUS_SEQUENCE {
    //
    // Produce a snp alignment and statistics files from genotyping results allowing extra checks for the genotyping geared to aDNA
    //

    //Need to ensure that it only runs when gatk3 unifiedgenotyper is run

    take:
        ch_genotypes_vcf
        ch_mva_files // [meta, additional_vcf, reference_gff, reference_gff_exclude ]
        fasta // [ meta, fasta, fai  ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.consensus_tool == 'multivcfanalyzer' ) {

        write_allele_frequencies = params.consensus_multivcfanalyzer_write_allele_frequencies ? "T" : "F"
        ch_mva_input = ch_mva_files
                            .multiMap{
                                meta, additional_vcf, reference_gff, reference_gff_exclude, reference_snpeff_results ->
                                vcfs:                      [ meta, additional_vcf ]
                                reference_gff:             [ meta, reference_gff ]
                                reference_gff_exclude:     [ meta, reference_gff_exclude ]
                                reference_snpeff_results : [ meta, reference_snpeff_results ]
                            }

        //Mix in the vcf from the additional vcf channel
        ch_mva_vcf  = ch_genotypes_vcf
                        .mix(
                            ch_mva_input.vcfs //Mix additional vcfs
                        )

        MULTIVCFANALYZER ( ch_mva_vcf,
                    fasta,
                    ch_mva_input.reference_snpeff_results,
                    ch_mva_input.reference_gff,
                    write_allele_frequencies,
                    params.consensus_multivcfanalyzer_min_genotype_quality,
                    params.consensus_multivcfanalyzer_min_base_coverage,
                    params.consensus_multivcfanalyzer_allele_freq_hom,
                    params.consensus_multivcfanalyzer_allele_freq_het,
                    ch_mva_input.reference_gff_exclude
        )

        ch_full_alignment_mva                = MULTIVCFANALYZER.out.full_alignment
        ch_info_mva                          = MULTIVCFANALYZER.out.info_txt
        ch_snp_alignment_mva                 = MULTIVCFANALYZER.out.snp_alignment
        ch_snp_genome_alignment_mva           = MULTIVCFANALYZER.out.snp_genome_alignment
        ch_snp_statistics_mva                = MULTIVCFANALYZER.out.snpstatistics
        ch_snp_table_mva                     = MULTIVCFANALYZER.out.snptable
        ch_snp_table_snpeff_mva              = MULTIVCFANALYZER.out.snptable_snpeff
        ch_snp_table_uncertainty_mva         = MULTIVCFANALYZER.out.snptable_uncertainty
        ch_structure_genotypes_mva           = MULTIVCFANALYZER.out.structure_genotypes
        ch_structure_genotypes_nomissing_mva = MULTIVCFANALYZER.out.structure_genotypes_nomissing
        ch_versions                          = ch_versions.mix( MULTIVCFANALYZER.out.versions )
        ch_multiqc_files                     = ch_multiqc_files.mix( MULTIVCFANALYZER.out.json )

    }

    emit:
        full_alignment_mva                = ch_full_alignment_mva
        info_mva                          = ch_info_mva
        snp_alignment_mva                 = ch_snp_alignment_mva
        snp_genome_alignment_mva          = ch_snp_genome_alignment_mva
        snpstatistics_mva                 = ch_snp_statistics_mva
        snptable_mva                      = ch_snp_table_mva
        snptable_snpeff_mva               = ch_snp_table_snpeff_mva
        snptable_uncertainty_mva          = ch_snp_table_uncertainty_mva
        structure_genotypes_mva           = ch_structure_genotypes_mva
        structure_genotypes_nomissing_mva = ch_structure_genotypes_nomissing_mva
        versions                          = ch_versions // channel: path(versions.yml)
        mqc                               = ch_multiqc_files // channel: [ val(meta), path("*.json") ]
}
