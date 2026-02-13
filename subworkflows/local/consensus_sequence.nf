//
// Produce consensus sequences from reads aligned to a reference some targetted to aDNA
//

include { MULTIVCFANALYZER         } from '../../modules/nf-core/multivcfanalyzer'
include { addNewMetaFromAttributes } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'
include { TABIX_BGZIP as UG_BGZIP  } from '../../modules/nf-core/tabix/bgzip'
include { GUNZIP as REF_MVA_GUNZIP } from '../../modules/nf-core/gunzip'

workflow CONSENSUS_SEQUENCE {
    take:
    ch_genotypes_vcf
    ch_samplesheet_vcfs // [meta, additional_vcf]
    ch_mva_files        // [meta, reference_gff, reference_gff_exclude,  ]
    ch_fasta            // [ meta, fasta ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    if (params.consensus_tool == 'multivcfanalyzer') {

        write_allele_frequencies = params.consensus_multivcfanalyzer_write_allele_frequencies ? "T" : "F"

        ch_genotypes_unzip = ch_genotypes_vcf.map { meta, vcfs, vcf_index ->
            [meta, vcfs]
        }

        UG_BGZIP(ch_genotypes_unzip)

        ch_genotypes_vcf_final = UG_BGZIP.out.output
            .map {
                addNewMetaFromAttributes(it, "reference", "reference", false)
            }
            .groupTuple()
            .map { metaref, meta, vcfs ->
                [metaref, vcfs]
            }
            .dump(tag: "consensus_genotyped_vcfs")

        REF_MVA_GUNZIP(ch_samplesheet_vcfs)
        
        ch_additional_vcfs = REF_MVA_GUNZIP.out.gunzip
            .dump(tag: "additional_vcfs")

        ch_fasta_final = ch_fasta
            .map { meta, fasta, fai, dict, mapindex ->
                def new_meta = meta.subMap(['id'])
                [[reference: new_meta.id], fasta]
            }
            .dump(tag: "consensus_fasta")

        ch_mva_input = ch_mva_files
            .dump(tag: "consensus_ref_related_files")
            .map { meta, reference_gff, reference_gff_exclude, reference_snpeff_results ->
                def new_meta = meta.subMap(['id'])
                [[reference: new_meta.id], reference_gff, reference_gff_exclude, reference_snpeff_results]
            }
            .join(ch_genotypes_vcf_final)
            .dump(tag: "consensus_postjoin")
            .join(ch_additional_vcfs)
            .join(ch_fasta_final)
            .dump(tag: "consensus_postjoin2")
//            .multiMap { meta, reference_gff, reference_gff_exclude, reference_snpeff_results, ug_vcfs, additional_vcf, fasta ->
//                vcfs: [meta, additional_vcf + ug_vcfs]
//                reference_gff: [meta, reference_gff ?: []]
//                reference_gff_exclude: [meta, reference_gff_exclude ?: []]
//                reference_snpeff_results: [meta, reference_snpeff_results ?: []]
//                reference_fasta: [meta, fasta]
//            }
//            .dump(tag: "consensus_sequence_final")

//        MULTIVCFANALYZER(
//            ch_mva_input.vcfs,
//            ch_mva_input.reference_fasta,
//            ch_mva_input.reference_snpeff_results,
//            ch_mva_input.reference_gff,
//            write_allele_frequencies,
//            params.consensus_multivcfanalyzer_min_genotype_quality,
//            params.consensus_multivcfanalyzer_min_base_coverage,
//            params.consensus_multivcfanalyzer_allele_freq_hom,
//            params.consensus_multivcfanalyzer_allele_freq_het,
//            ch_mva_input.reference_gff_exclude,
//        )

//      ch_full_alignment_mva = MULTIVCFANALYZER.out.full_alignment
//      ch_info_mva = MULTIVCFANALYZER.out.info_txt
//      ch_snp_alignment_mva = MULTIVCFANALYZER.out.snp_alignment
//      ch_snp_genome_alignment_mva = MULTIVCFANALYZER.out.snp_genome_alignment
//      ch_snp_statistics_mva = MULTIVCFANALYZER.out.snpstatistics
//      ch_snp_table_mva = MULTIVCFANALYZER.out.snptable
//      ch_snp_table_snpeff_mva = MULTIVCFANALYZER.out.snptable_snpeff
//      ch_snp_table_uncertainty_mva = MULTIVCFANALYZER.out.snptable_uncertainty
//      ch_structure_genotypes_mva = MULTIVCFANALYZER.out.structure_genotypes
//      ch_structure_genotypes_nomissing_mva = MULTIVCFANALYZER.out.structure_genotypes_nomissing
//      ch_versions = ch_versions.mix(MULTIVCFANALYZER.out.versions)
//      ch_multiqc_files = ch_multiqc_files.mix(MULTIVCFANALYZER.out.json)
    }

    emit:
    versions = ch_versions // channel: path(versions.yml)
    mqc      = ch_multiqc_files // channel: [ val(meta), path("*.json") ]
}
