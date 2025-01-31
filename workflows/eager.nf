/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap                                    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_eager_pipeline'
include { addNewMetaFromAttributes                            } from '../subworkflows/local/utils_nfcore_eager_pipeline/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// TODO rename to active: index_reference, filter_bam etc.
include { REFERENCE_INDEXING                                  } from '../subworkflows/local/reference_indexing'
include { PREPROCESSING                                       } from '../subworkflows/local/preprocessing'
include { MAP                                                 } from '../subworkflows/local/map'
include { MERGE_LANES_INPUTBAM                                } from '../subworkflows/local/merge_lanes_inputbam'
include { FILTER_BAM                                          } from '../subworkflows/local/bamfiltering.nf'
include { DEDUPLICATE                                         } from '../subworkflows/local/deduplicate'
include { MANIPULATE_DAMAGE                                   } from '../subworkflows/local/manipulate_damage'
include { METAGENOMICS_COMPLEXITYFILTER                       } from '../subworkflows/local/metagenomics_complexityfilter'
include { ESTIMATE_CONTAMINATION                              } from '../subworkflows/local/estimate_contamination'
include { CALCULATE_DAMAGE                                    } from '../subworkflows/local/calculate_damage'
include { RUN_SEXDETERRMINE                                   } from '../subworkflows/local/run_sex_determination'
include { MERGE_LIBRARIES                                     } from '../subworkflows/local/merge_libraries'
include { MERGE_LIBRARIES as MERGE_LIBRARIES_GENOTYPING       } from '../subworkflows/local/merge_libraries'
include { GENOTYPE                                            } from '../subworkflows/local/genotype'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                                              } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { SAMTOOLS_COLLATEFASTQ as SAMTOOLS_CONVERT_BAM_INPUT } from '../modules/nf-core/samtools/collatefastq/main'
include { CAT_FASTQ as CAT_FASTQ_CONVERTED_BAM                } from '../modules/nf-core/cat/fastq/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BAM_INPUT          } from '../modules/nf-core/samtools/index/main'
include { PRESEQ_CCURVE                                       } from '../modules/nf-core/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP                                     } from '../modules/nf-core/preseq/lcextrap/main'
include { FALCO                                               } from '../modules/nf-core/falco/main'
include { MTNUCRATIO                                          } from '../modules/nf-core/mtnucratio/main'
include { HOST_REMOVAL                                        } from '../modules/local/host_removal'
include { ENDORSPY                                            } from '../modules/nf-core/endorspy/main'
include { BEDTOOLS_COVERAGE as BEDTOOLS_COVERAGE_DEPTH        } from '../modules/nf-core/bedtools/coverage/main'
include { BEDTOOLS_COVERAGE as BEDTOOLS_COVERAGE_BREADTH      } from '../modules/nf-core/bedtools/coverage/main'
include { SAMTOOLS_VIEW_GENOME                                } from '../modules/local/samtools_view_genome.nf'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_NOBED              } from '../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_WITHBED            } from '../modules/nf-core/qualimap/bamqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow EAGER {
    take:
    ch_samplesheet_fastqs // channel: samplesheet FASTQ entries read in from --input
    ch_samplesheet_bams   // channel: samplesheet BAM entries read in from --input

    main:

    log.info("Schaffa, Schaffa, Genome Baua!")

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Reference
    fasta_fn = params.fasta ? file(params.fasta, checkIfExists: true) : params.fasta_sheet ? file(params.fasta_sheet, checkIfExists: true) : []
    fasta_fai = params.fasta_fai ? file(params.fasta_fai, checkIfExists: true) : []
    fasta_dict = params.fasta_dict ? file(params.fasta_dict, checkIfExists: true) : []
    fasta_mapperindexdir = params.fasta_mapperindexdir ? file(params.fasta_mapperindexdir, checkIfExists: true) : []

    // Preprocessing
    adapterlist = params.preprocessing_skipadaptertrim ? [] : params.preprocessing_adapterlist ? file(params.preprocessing_adapterlist, checkIfExists: true) : []

    if (params.preprocessing_adapterlist && !params.preprocessing_skipadaptertrim) {
        if (params.preprocessing_tool == 'adapterremoval' && !(adapterlist.extension == 'txt')) {
            error("[nf-core/eager] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --preprocessing_adapterlist ${params.preprocessing_adapterlist}")
        }
        if (params.preprocessing_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)")) {
            error("[nf-core/eager] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --preprocessing_adapterlist ${params.preprocessing_adapterlist}")
        }
    }

    //
    // MODULE: Convert input BAMs back to FastQ
    //

    if (params.convert_inputbam) {
        // Convert input BAMs back to FastQ with non-interleaved output.
        SAMTOOLS_CONVERT_BAM_INPUT(ch_samplesheet_bams, [[], []], false)

        // if BAM is single-end, pull R1 output as well as 'other' output and merge (in case collapsed reads have their R1 and R2 flags both set to 0 or 1)
        ch_single_end_reads = SAMTOOLS_CONVERT_BAM_INPUT.out.fastq
            .filter { meta, reads ->
                meta.single_end
            }
            .join(SAMTOOLS_CONVERT_BAM_INPUT.out.fastq_other)
            .map { meta, read1, fastq_other ->
                [meta, [read1, fastq_other]]
            }

        // Put all the converted FASTQs with single-end reads back together again
        CAT_FASTQ_CONVERTED_BAM(ch_single_end_reads)

        //if BAM is paired-end, pull R1 and R2 outputs, discarding 'other' output and singletons
        ch_paired_end_reads = SAMTOOLS_CONVERT_BAM_INPUT.out.fastq.filter { meta, reads ->
            !meta.single_end
        }

        ch_fastqs_from_converted_bams = CAT_FASTQ_CONVERTED_BAM.out.reads
            .mix(ch_paired_end_reads)
            .map { meta, reads ->
                [meta - meta.subMap('reference', 'id_index'), reads]
            }

        // Mix the converted fastqs with the original fastqs
        ch_fastqs_for_preprocessing = ch_fastqs_from_converted_bams.mix(ch_samplesheet_fastqs)
    }
    else {
        // If BAM conversion is not activated , just use the original fastqs
        ch_fastqs_for_preprocessing = ch_samplesheet_fastqs
    }

    //
    // SUBWORKFLOW: Indexing of reference files
    //

    REFERENCE_INDEXING(fasta_fn, fasta_fai, fasta_dict, fasta_mapperindexdir)
    ch_versions = ch_versions.mix(REFERENCE_INDEXING.out.versions)

    //
    // MODULE: Run FastQC or Falco
    //

    if (params.sequencing_qc_tool == "falco") {
        FALCO(ch_fastqs_for_preprocessing)
        ch_versions = ch_versions.mix(FALCO.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FALCO.out.txt.collect { it[1] }.ifEmpty([]))
    }
    else {
        FASTQC(ch_fastqs_for_preprocessing)
        ch_versions = ch_versions.mix(FASTQC.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: Read preprocessing (clipping, merging, fastq trimming etc. )
    //

    if (!params.skip_preprocessing) {
        PREPROCESSING(ch_fastqs_for_preprocessing, adapterlist)
        ch_reads_for_mapping = PREPROCESSING.out.reads
        ch_versions = ch_versions.mix(PREPROCESSING.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(PREPROCESSING.out.mqc.collect { it[1] }.ifEmpty([]))
    }
    else {
        ch_reads_for_mapping = ch_fastqs_for_preprocessing
    }

    //
    // SUBWORKFLOW: Reference mapping
    //
    ch_reference_for_mapping = REFERENCE_INDEXING.out.reference.map { meta, fasta, fai, dict, index ->
        [meta, index, fasta]
    }

    MAP(ch_reads_for_mapping, ch_reference_for_mapping, REFERENCE_INDEXING.out.elongated_reference, REFERENCE_INDEXING.out.elongated_chr_list)

    ch_versions = ch_versions.mix(MAP.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(MAP.out.mqc.collect { it[1] }.ifEmpty([]))

    //
    //  MODULE: indexing of user supplied unconverted input BAMs
    //

    if (!params.convert_inputbam) {
        SAMTOOLS_INDEX_BAM_INPUT(ch_samplesheet_bams)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX_BAM_INPUT.out.versions)

        if (params.fasta_largeref) {
            ch_bams_from_input = ch_samplesheet_bams.join(SAMTOOLS_INDEX_BAM_INPUT.out.csi)
        }
        else {
            ch_bams_from_input = ch_samplesheet_bams.join(SAMTOOLS_INDEX_BAM_INPUT.out.bai)
        }

        // SUBWORKFLOW: Merging lanes for ch_bams_from_input

        MERGE_LANES_INPUTBAM(ch_bams_from_input)
        ch_bams_from_input_lanemerged = MERGE_LANES_INPUTBAM.out.bam
                                            .join(MERGE_LANES_INPUTBAM.out.bai)

    } else {
        ch_bams_from_input    = Channel.empty()
        ch_flagstat_input_bam = Channel.empty()
    }


    //
    // SUBWORKFLOW: bam filtering (length, mapped/unmapped, quality etc.)
    //

    if (params.run_bamfiltering || params.run_metagenomics) {

        ch_mapped_for_bamfilter = MAP.out.bam
                                    .join(MAP.out.bai)
                                    .mix(ch_bams_from_input_lanemerged)
        FILTER_BAM(ch_mapped_for_bamfilter)
        ch_bamfiltered_for_deduplication = FILTER_BAM.out.genomics
        ch_bamfiltered_for_metagenomics = FILTER_BAM.out.metagenomics
        ch_versions = ch_versions.mix(FILTER_BAM.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(FILTER_BAM.out.mqc.collect { it[1] }.ifEmpty([]))
    }
    else {
        ch_bamfiltered_for_deduplication = MAP.out.bam
                                                .join(MAP.out.bai)
                                                .mix(ch_bams_from_input_lanemerged)
    }

    ch_reads_for_deduplication = ch_bamfiltered_for_deduplication

    //
    // SUBWORKFLOW: genomic BAM deduplication
    //

    ch_fasta_for_deduplication = REFERENCE_INDEXING.out.reference.multiMap { meta, fasta, fai, dict, index ->
        fasta: [meta, fasta]
        fasta_fai: [meta, fai]
    }

    if (!params.skip_deduplication) {
        DEDUPLICATE(ch_reads_for_deduplication, ch_fasta_for_deduplication.fasta, ch_fasta_for_deduplication.fasta_fai)
        ch_dedupped_bams = DEDUPLICATE.out.bam.join(DEDUPLICATE.out.bai)
        ch_dedupped_flagstat = DEDUPLICATE.out.flagstat
        ch_versions = ch_versions.mix(DEDUPLICATE.out.versions)
    }
    else {
        ch_dedupped_bams = ch_reads_for_deduplication
        ch_dedupped_flagstat = Channel.empty()
    }

    //
    // SUBWORKFLOW: Merge libraries per sample
    //

    MERGE_LIBRARIES(ch_dedupped_bams)
    ch_versions = ch_versions.mix(MERGE_LIBRARIES.out.versions)
    ch_merged_dedup_bams = MERGE_LIBRARIES.out.bam_bai
    ch_multiqc_files = ch_multiqc_files.mix(MERGE_LIBRARIES.out.mqc.collect { it[1] }.ifEmpty([]))

    //
    // MODULE QUALIMAP
    //

    if (!params.skip_qualimap) {
        ch_snp_capture_bed = REFERENCE_INDEXING.out.snp_capture_bed.map {
            addNewMetaFromAttributes(it, "id", "reference", false)
        }
        ch_qualimap_input = ch_merged_dedup_bams
            .map { meta, bam, bai ->
                [meta, bam]
            }
            .map {
                addNewMetaFromAttributes(it, "reference", "reference", false)
            }
            .combine(
                ch_snp_capture_bed,
                by: 0
            )
            .branch { ignore_meta, meta, bam, meta2, snp_capture_bed ->
                withbed: snp_capture_bed != ""
                nobed: true
            }
        ch_qualimap_input_with = ch_qualimap_input.withbed.multiMap { ignore_meta, meta, bam, meta2, snp_capture_bed ->
            bam: [meta, bam]
            snp_capture_bed: [snp_capture_bed]
        }

        QUALIMAP_BAMQC_WITHBED(ch_qualimap_input_with.bam, ch_qualimap_input_with.snp_capture_bed)
        ch_qualimap_input_without = ch_qualimap_input.nobed.map { ignore_meta, meta, bam, meta2, snp_capture_bed ->
            [meta, bam]
        }

        QUALIMAP_BAMQC_NOBED(ch_qualimap_input_without, [])
        ch_qualimap_output = QUALIMAP_BAMQC_WITHBED.out.results.mix(QUALIMAP_BAMQC_NOBED.out.results)
        ch_versions = ch_versions.mix(QUALIMAP_BAMQC_NOBED.out.versions).mix(QUALIMAP_BAMQC_WITHBED.out.versions)
    }

    //
    // MODULE: remove reads mapping to the host from the raw fastq
    //

    if (params.run_host_removal) {
        // Preparing bam channel for host removal to be combined with the input fastq channel
        // The bam channel consist of [meta, bam, bai] and in the meta we have in addition 'single_end' always set as TRUE and 'reference' set
        // To be able to join it with fastq channel, we need to remove them from the meta (done in map) and stored in new_meta
        ch_bam_for_host_removal = MAP.out.bam
            .join(MAP.out.bai)
            .map { meta, bam, bai ->
                new_meta = meta.clone().findAll { it.key !in ['single_end', 'reference'] }
                [new_meta, meta, bam, bai]
            }
        // Preparing fastq channel for host removal to be combined with the bam channel
        // The meta of the fastq channel contains additional fields when compared to the meta from the bam channel: lane, colour_chemistry,
        // and not necessarily matching single_end. Those fields are dropped of the meta in the map and stored in new_meta
        ch_fastqs_for_host_removal = ch_fastqs_for_preprocessing.map { meta, fastqs ->
            new_meta = meta.clone().findAll { it.key !in ['lane', 'colour_chemistry', 'single_end'] }
            [new_meta, meta, fastqs]
        }
        // We join the bam and fastq channel with now matching metas (new_meta) referred as meta_join
        // and remove the meta_join from the final channel, keeping the original metas for the bam and the fastqs
        ch_input_for_host_removal = ch_bam_for_host_removal
            .join(ch_fastqs_for_host_removal)
            .map { meta_join, meta_bam, bam, bai, meta_fastq, fastqs ->
                [meta_bam, bam, bai, meta_fastq, fastqs]
            }

        HOST_REMOVAL(ch_input_for_host_removal)

        ch_versions = ch_versions.mix(HOST_REMOVAL.out.versions)
    }

    //
    // Section: Metagenomics
    //

    if (params.run_metagenomics) {

        ch_database = Channel.fromPath(params.metagenomics_profiling_database)

        // this is for MALT
        ch_tax_list = Channel.empty()
        ch_ncbi_dir = Channel.empty()

        if (params.metagenomics_run_postprocessing && params.metagenomics_profiling_tool == 'malt') {
            ch_tax_list = Channel.fromPath(params.metagenomics_maltextract_taxonlist, checkIfExists: true)
            ch_ncbi_dir = Channel.fromPath(params.metagenomics_maltextract_ncbidir, checkIfExists: true)
        }

        METAGENOMICS(ch_bamfiltered_for_metagenomics, ch_database, ch_tax_list, ch_ncbi_dir)
        ch_versions = ch_versions.mix(METAGENOMICS.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(METAGENOMICS.out.ch_multiqc_files)
    }

    //
    // MODULE: MTNUCRATIO
    //

    if (params.run_mtnucratio) {
        ch_mito_header = REFERENCE_INDEXING.out.mitochondrion_header.map {
            addNewMetaFromAttributes(it, "id", "reference", false)
        }
        mtnucratio_input = ch_dedupped_bams
            .map {
                addNewMetaFromAttributes(it, "reference", "reference", false)
            }
            .combine(
                ch_mito_header,
                by: 0
            )
            .multiMap { ignore_meta, meta, bam, bai, meta2, mito_header ->
                bam: [meta, bam]
                mito_header: [meta2, mito_header]
            }

        MTNUCRATIO(mtnucratio_input.bam, mtnucratio_input.mito_header.map { it[1] })
        ch_multiqc_files = ch_multiqc_files.mix(MTNUCRATIO.out.mtnucratio.collect { it[1] }.ifEmpty([]))
        ch_versions = ch_versions.mix(MTNUCRATIO.out.versions)
    }

    //
    // MODULE: ENDORSPY (raw, filtered, deduplicated)
    //

    ch_flagstat_for_endorspy_raw    = MAP.out.flagstat
                                            .mix( MERGE_LANES_INPUTBAM.out.flagstat )

    if (params.run_bamfiltering & !params.skip_deduplication) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
            .join(FILTER_BAM.out.flagstat)
            .join(DEDUPLICATE.out.flagstat)
    }
    else if (params.run_bamfiltering & params.skip_deduplication) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
            .join(FILTER_BAM.out.flagstat)
            .map { meta, flags_raw, flags_filtered ->
                [meta, flags_raw, flags_filtered, []]
            }
    }
    else if (!params.run_bamfiltering & !params.skip_deduplication) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
            .join(DEDUPLICATE.out.flagstat)
            .map { meta, flags_raw, flags_dedup ->
                [meta, flags_raw, [], flags_dedup]
            }
    }
    else {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw.map { meta, flags_raw ->
            [meta, flags_raw, [], []]
        }
    }

    ENDORSPY(ch_for_endorspy)

    ch_versions = ch_versions.mix(ENDORSPY.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(ENDORSPY.out.json.collect { it[1] }.ifEmpty([]))

    //
    // MODULE: PreSeq
    //

    if (!params.mapstats_skip_preseq && params.mapstats_preseq_mode == 'c_curve') {
        PRESEQ_CCURVE(ch_reads_for_deduplication.map { [it[0], it[1]] })
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_CCURVE.out.c_curve.collect { it[1] }.ifEmpty([]))
        ch_versions = ch_versions.mix(PRESEQ_CCURVE.out.versions)
    }
    else {
        (!params.mapstats_skip_preseq && params.mapstats_preseq_mode == 'lc_extrap').call {
            PRESEQ_LCEXTRAP(ch_reads_for_deduplication.map { [it[0], it[1]] })
            ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect { it[1] }.ifEmpty([]))
            ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions)
        }
    }


    //
    // MODULE: Bedtools coverage
    //

    if (params.run_bedtools_coverage) {

        ch_bedtools_feature = REFERENCE_INDEXING.out.bedtools_feature.map {
            addNewMetaFromAttributes(it, "id", "reference", false)
        }

        ch_bedtools_prep = ch_merged_dedup_bams
            .map {
                addNewMetaFromAttributes(it, "reference", "reference", false)
            }
            .combine(
                ch_bedtools_feature,
                by: 0
            )
            .map { ignore_meta, meta, bam, bai, meta2, bedtools_feature ->
                [meta, bedtools_feature, bam, bai]
            }
            .branch { meta, bedtools_feature, bam, bai ->
                withfeature: bedtools_feature != ""
                nobed: true
            }

        // Running samtools view to get header
        ch_bedtools_input = ch_bedtools_prep.withfeature.multiMap { meta, bedtools_feature, bam, bai ->
            bam: [meta, bam, bai]
            withfeature: [meta, bedtools_feature, bam]
        }

        SAMTOOLS_VIEW_GENOME(ch_bedtools_input.bam)

        ch_genome_for_bedtools = SAMTOOLS_VIEW_GENOME.out.genome

        BEDTOOLS_COVERAGE_DEPTH(ch_bedtools_input.withfeature, ch_genome_for_bedtools)

        ch_versions = ch_versions.mix(SAMTOOLS_VIEW_GENOME.out.versions)
        //ch_versions = ch_versions.mix( BEDTOOLS_COVERAGE_BREADTH.out.versions )
        ch_versions = ch_versions.mix(BEDTOOLS_COVERAGE_DEPTH.out.versions)
    }

    //
    // SUBWORKFLOW: Calculate Damage
    //

    ch_fasta_for_damagecalculation = REFERENCE_INDEXING.out.reference.multiMap { meta, fasta, fai, dict, index ->
        fasta: [meta, fasta]
        fasta_fai: [meta, fai]
    }

    if (!params.skip_damagecalculation) {
        CALCULATE_DAMAGE(ch_dedupped_bams, ch_fasta_for_damagecalculation.fasta, ch_fasta_for_damagecalculation.fasta_fai)
        ch_versions = ch_versions.mix(CALCULATE_DAMAGE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_DAMAGE.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: Run Sex Determination
    //

    if (params.run_sexdeterrmine) {
        ch_sexdeterrmine_input = ch_merged_dedup_bams

        RUN_SEXDETERRMINE(ch_sexdeterrmine_input, REFERENCE_INDEXING.out.sexdeterrmine_bed)
        ch_versions = ch_versions.mix(RUN_SEXDETERRMINE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(RUN_SEXDETERRMINE.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: Contamination estimation
    //

    if (params.run_contamination_estimation_angsd) {
        contamination_input = ch_dedupped_bams

        ESTIMATE_CONTAMINATION(contamination_input, REFERENCE_INDEXING.out.hapmap)
        ch_versions = ch_versions.mix(ESTIMATE_CONTAMINATION.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(ESTIMATE_CONTAMINATION.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    //
    // SUBWORKFLOW: aDNA Damage Manipulation
    //

    if (params.run_mapdamage_rescaling || params.run_pmd_filtering || params.run_trim_bam) {
        MANIPULATE_DAMAGE(ch_dedupped_bams, ch_fasta_for_deduplication.fasta, REFERENCE_INDEXING.out.pmd_masking)
        ch_multiqc_files = ch_multiqc_files.mix(MANIPULATE_DAMAGE.out.flagstat.collect { it[1] }.ifEmpty([]))
        ch_versions = ch_versions.mix(MANIPULATE_DAMAGE.out.versions)
        ch_bams_for_library_merge = params.genotyping_source == 'rescaled' ? MANIPULATE_DAMAGE.out.rescaled : params.genotyping_source == 'pmd' ? MANIPULATE_DAMAGE.out.filtered : params.genotyping_source == 'trimmed' ? MANIPULATE_DAMAGE.out.trimmed : ch_merged_dedup_bams

        // SUBWORKFLOW: merge libraries for genotyping
        MERGE_LIBRARIES_GENOTYPING(ch_bams_for_library_merge)
        ch_versions = ch_versions.mix(MERGE_LIBRARIES_GENOTYPING.out.versions)
        ch_bams_for_genotyping = MERGE_LIBRARIES_GENOTYPING.out.bam_bai
        ch_multiqc_files = ch_multiqc_files.mix(MERGE_LIBRARIES_GENOTYPING.out.mqc.collect { it[1] }.ifEmpty([]))
    }
    else {
        ch_bams_for_genotyping = ch_merged_dedup_bams
    }

    //
    // SUBWORKFLOW: Genotyping
    //

    if (params.run_genotyping) {
        ch_reference_for_genotyping = REFERENCE_INDEXING.out.reference.map { meta, fasta, fai, dict, mapindex ->
            [meta, fasta, fai, dict]
        }
        GENOTYPE(
            ch_bams_for_genotyping,
            ch_reference_for_genotyping,
            REFERENCE_INDEXING.out.pileupcaller_bed_snp.ifEmpty([[], [], []]),
            REFERENCE_INDEXING.out.dbsnp.ifEmpty([[], []])
        )

        ch_versions = ch_versions.mix(GENOTYPE.out.versions)
        ch_multiqc_files = ch_multiqc_files.mix(GENOTYPE.out.mqc.collect { it[1] }.ifEmpty([]))
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'eager_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    if (!params.skip_qualimap) {
        ch_multiqc_files = ch_multiqc_files.mix(ch_qualimap_output.collect { it[1] }.ifEmpty([]))
    }

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
