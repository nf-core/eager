//
// Produce consensus sequences from reads aligned to a reference some targetted to aDNA
//

include { MULTIVCFANALYZER         } from '../../modules/nf-core/multivcfanalyzer'
include { addNewMetaFromAttributes } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'
include { TABIX_BGZIP as UG_BGZIP  } from '../../modules/nf-core/tabix/bgzip'
include { GUNZIP as REF_MVA_GUNZIP } from '../../modules/nf-core/gunzip'

workflow CONSENSUS_SEQUENCE {
    take:
    ch_genotypes_vcf    // [ meta, genotyped_vcf, genotyped_vcf_index ]
    ch_samplesheet_vcfs // [ meta, additional_vcf]
    ch_mva_files        // [ meta, reference_gff, reference_gff_exclude, reference_snpeff_results ]
    ch_fasta            // [ meta, fasta ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    if (params.consensus_tool == 'multivcfanalyzer') {
        // Prepare VCF sets for MVCFA input.
        ch_mvcfa_vcf_input = ch_genotypes_vcf
            .map {
                meta, vcfs, vcf_index ->
                [meta, vcfs]
            }
            .map {
                addNewMetaFromAttributes(it, "reference", "reference", true)
            }
            .mix(
                // Mix in any additional VCFs
                ch_samplesheet_vcfs
                    .map {
                        addNewMetaFromAttributes(it, "vcf_reference_id", "reference", true)
                    }
            )
            // Filter out lines with missing VCF. Serves to remove the ifEmpty input of ch_samplesheet_vcfs.
            .filter{
                merge_meta, vcfs ->
                vcfs != []
            }
            .groupTuple() // [ meta, [genotyping_vcfs + additional_vcfs] ]

        ch_fasta_final = ch_fasta
            .map {
                addNewMetaFromAttributes(it, "id", "reference", true)
            }
            .map { meta, fasta, fai, dict, mapindex ->
                [meta, fasta]
            }

        ch_mva_input = ch_mva_files
            .map{
                addNewMetaFromAttributes(it, "id", "reference", true)
            }
            .join(ch_mvcfa_vcf_input)
            .join(ch_fasta_final)
            .multiMap{ meta, reference_gff, reference_gff_exclude, reference_snpeff_results, vcf_inputs, fasta ->
                vcfs: [meta, vcf_inputs]
                reference_gff: [meta, reference_gff ?: []]
                reference_gff_exclude: [meta, reference_gff_exclude ?: []]
                reference_snpeff_results: [meta, reference_snpeff_results ?: []]
                reference_fasta: [meta, fasta]
            }

        MULTIVCFANALYZER(
            ch_mva_input.vcfs,
            ch_mva_input.reference_fasta,
            ch_mva_input.reference_snpeff_results,
            ch_mva_input.reference_gff,
            params.consensus_multivcfanalyzer_write_allele_frequencies,
            params.consensus_multivcfanalyzer_min_genotype_quality,
            params.consensus_multivcfanalyzer_min_base_coverage,
            params.consensus_multivcfanalyzer_allele_freq_hom,
            params.consensus_multivcfanalyzer_allele_freq_het,
            ch_mva_input.reference_gff_exclude,
        )

        ch_full_alignment_mva                = MULTIVCFANALYZER.out.full_alignment
        ch_info_mva                          = MULTIVCFANALYZER.out.info_txt
        ch_snp_alignment_mva                 = MULTIVCFANALYZER.out.snp_alignment
        ch_snp_genome_alignment_mva          = MULTIVCFANALYZER.out.snp_genome_alignment
        ch_snp_statistics_mva                = MULTIVCFANALYZER.out.snpstatistics
        ch_snp_table_mva                     = MULTIVCFANALYZER.out.snptable
        ch_snp_table_snpeff_mva              = MULTIVCFANALYZER.out.snptable_snpeff
        ch_snp_table_uncertainty_mva         = MULTIVCFANALYZER.out.snptable_uncertainty
        ch_structure_genotypes_mva           = MULTIVCFANALYZER.out.structure_genotypes
        ch_structure_genotypes_nomissing_mva = MULTIVCFANALYZER.out.structure_genotypes_nomissing
        ch_versions                          = ch_versions.mix(MULTIVCFANALYZER.out.versions_multivcfanalyzer, MULTIVCFANALYZER.out.versions_tabix)
        ch_multiqc_files                     = ch_multiqc_files.mix(MULTIVCFANALYZER.out.json)
    }

    emit:
    full_alignment_mva                = ch_full_alignment_mva                // channel: [ val(meta), path("fullAlignment.fasta.gz") ]
    info_mva                          = ch_info_mva                          // channel: [ val(meta), path("info.txt") ]
    snp_alignment_mva                 = ch_snp_alignment_mva                 // channel: [ val(meta), path("snpAlignment.fasta.gz") ]
    snp_genome_alignment_mva          = ch_snp_genome_alignment_mva          // channel: [ val(meta), path("snpAlignmentIncludingRefGenome.fasta.gz") ]
    snp_statistics_mva                = ch_snp_statistics_mva                // channel: [ val(meta), path("snpStatistics.tsv") ]
    snp_table_mva                     = ch_snp_table_mva                     // channel: [ val(meta), path("snpTable.tsv") ]
    snp_table_snpeff_mva              = ch_snp_table_snpeff_mva              // channel: [ val(meta), path("snpTableForSnpEff.tsv") ]
    snp_table_uncertainty_mva         = ch_snp_table_uncertainty_mva         // channel: [ val(meta), path("snpTableWithUncertaintyCalls.tsv") ]
    structure_genotypes_mva           = ch_structure_genotypes_mva           // channel: [ val(meta), path("structureGenotypes.tsv") ]
    structure_genotypes_nomissing_mva = ch_structure_genotypes_nomissing_mva // channel: [ val(meta), path("structureGenotypes_noMissingData-Columns.tsv") ]
    versions                          = ch_versions                          // channel: [path(versions.yml)]
    mqc                               = ch_multiqc_files                     // channel: [ val(meta), path("MultiVCFAnalyzer.json") ]
}
