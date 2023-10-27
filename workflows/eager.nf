/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowEager.initialise(params, log)

// Check failing parameter combinations
if ( params.bamfiltering_retainunmappedgenomicbam && params.bamfiltering_mappingquality > 0  ) { exit 1, ("[nf-core/eager] ERROR: You cannot both retain unmapped reads and perform quality filtering, as unmapped reads have a mapping quality of 0. Pick one or the other functionality.") }
if ( params.genotyping_source == 'trimmed'        && ! params.run_trim_bam                   ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'trimmed' unless BAM trimming is turned on with `--run_trim_bam`.") }
if ( params.genotyping_source == 'pmd'            && ! params.run_pmd_filtering              ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'pmd' unless PMD-filtering is ran.") }
if ( params.genotyping_source == 'rescaled'       && ! params.run_mapdamage_rescaling        ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'rescaled' unless aDNA damage rescaling is ran.") }
if ( params.metagenomics_complexity_tool == 'prinseq' && params.metagenomics_prinseq_mode == 'dust' && params.metagenomics_complexity_entropy != 0.3 ) {
    // entropy score was set but dust method picked. If no dust-score provided, assume it was an error and fail
    if (params.metagenomics_prinseq_dustscore == 0.5) {
            exit 1, ("[nf-core/eager] ERROR: Metagenomics: You picked PRINSEQ++ with 'dust' mode but provided an entropy score. Please specify a dust filter threshold using the --metagenomics_prinseq_dustscore flag")
    }
}
if ( params.metagenomics_complexity_tool == 'prinseq' && params.metagenomics_prinseq_mode == 'entropy' && params.metagenomics_prinseq_dustscore != 0.5 ) {
    // dust score was set but entropy method picked. If no entropy-score provided, assume it was an error and fail
    if (params.metagenomics_complexity_entropy == 0.3) {
            exit 1, ("[nf-core/eager] ERROR: Metagenomics: You picked PRINSEQ++ with 'entropy' mode but provided a dust score. Please specify an entropy filter threshold using the --metagenomics_complexity_entropy flag")
    }
}

// TODO What to do when params.preprocessing_excludeunmerged is provided but the data is SE?
if ( params.deduplication_tool == 'dedup' && ! params.preprocessing_excludeunmerged ) { exit 1, "[nf-core/eager] ERROR: Dedup can only be used on collapsed (i.e. merged) PE reads. For all other cases, please set --deduplication_tool to 'markduplicates'."}

// Report possible warnings
if ( params.preprocessing_skipadaptertrim && params.preprocessing_adapterlist ) log.warn("[nf-core/eager] --preprocessing_skipadaptertrim will override --preprocessing_adapterlist. Adapter trimming will be skipped!")

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// TODO rename to active: index_reference, filter_bam etc.
include { INPUT_CHECK                   } from '../subworkflows/local/input_check'
include { REFERENCE_INDEXING            } from '../subworkflows/local/reference_indexing'
include { PREPROCESSING                 } from '../subworkflows/local/preprocessing'
include { MAP                           } from '../subworkflows/local/map'
include { FILTER_BAM                    } from '../subworkflows/local/bamfiltering.nf'
include { DEDUPLICATE                   } from '../subworkflows/local/deduplicate'
include { MANIPULATE_DAMAGE             } from '../subworkflows/local/manipulate_damage'
include { METAGENOMICS_COMPLEXITYFILTER } from '../subworkflows/local/metagenomics_complexityfilter'
include { ESTIMATE_CONTAMINATION        } from '../subworkflows/local/estimate_contamination'
include { CALCULATE_DAMAGE              } from '../subworkflows/local/calculate_damage'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { FASTQC                                            } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                           } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_INDEX                                    } from '../modules/nf-core/samtools/index/main'
include { PRESEQ_CCURVE                                     } from '../modules/nf-core/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP                                   } from '../modules/nf-core/preseq/lcextrap/main'
include { FALCO                                             } from '../modules/nf-core/falco/main'
include { MTNUCRATIO                                        } from '../modules/nf-core/mtnucratio/main'
include { HOST_REMOVAL                                      } from '../modules/local/host_removal'
include { ENDORSPY                                          } from '../modules/nf-core/endorspy/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTATS_BAM_INPUT } from '../modules/nf-core/samtools/flagstat/main'
include { BEDTOOLS_COVERAGE as BEDTOOLS_COVERAGE_DEPTH ; BEDTOOLS_COVERAGE as BEDTOOLS_COVERAGE_BREADTH } from '../modules/nf-core/bedtools/coverage/main'
include { SAMTOOLS_VIEW_GENOME                              } from '../modules/local/samtools_view_genome.nf'
include { QUALIMAP_BAMQC as QUALIMAP_BAMQC_NOBED ; QUALIMAP_BAMQC as QUALIMAP_BAMQC_WITHBED             } from '../modules/nf-core/qualimap/bamqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EAGER {

    log.info "Schaffa, Schaffa, Genome Baua!"

    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    //
    // Input file checks
    //

    // Reference
    fasta                = file(params.fasta, checkIfExists: true)
    fasta_fai            = params.fasta_fai ? file(params.fasta_fai, checkIfExists: true) : []
    fasta_dict           = params.fasta_dict ? file(params.fasta_dict, checkIfExists: true) : []
    fasta_mapperindexdir = params.fasta_mapperindexdir ? file(params.fasta_mapperindexdir, checkIfExists: true) : []

    // Preprocessing
    adapterlist          = params.preprocessing_skipadaptertrim ? [] : params.preprocessing_adapterlist ? file(params.preprocessing_adapterlist, checkIfExists: true) : []


    if ( params.preprocessing_adapterlist && !params.preprocessing_skipadaptertrim ) {
        if ( params.preprocessing_tool == 'adapterremoval' && !(adapterlist.extension == 'txt') ) error "[nf-core/eager] ERROR: AdapterRemoval2 adapter list requires a `.txt` format and extension. Check input: --preprocessing_adapterlist ${params.preprocessing_adapterlist}"
        if ( params.preprocessing_tool == 'fastp' && !adapterlist.extension.matches(".*(fa|fasta|fna|fas)") ) error "[nf-core/eager] ERROR: fastp adapter list requires a `.fasta` format and extension (or fa, fas, fna). Check input: --preprocessing_adapterlist ${params.preprocessing_adapterlist}"
    }


    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //

    INPUT_CHECK (
        file(params.input)
    )
    ch_versions = ch_versions.mix( INPUT_CHECK.out.versions )
    // TODO: OPTIONAL, you can use nf-validation plugin to create an input channel from the samplesheet with Channel.fromSamplesheet("input")
    // See the documentation https://nextflow-io.github.io/nf-validation/samplesheets/fromSamplesheet/
    // ! There is currently no tooling to help you write a sample sheet schema

    //
    // SUBWORKFLOW: Indexing of reference files
    //

    REFERENCE_INDEXING ( fasta, fasta_fai, fasta_dict, fasta_mapperindexdir )
    ch_versions = ch_versions.mix( REFERENCE_INDEXING.out.versions )

    //
    // MODULE: Run FastQC or Falco
    //

    if ( params.sequencing_qc_tool == "falco" ) {
        FALCO ( INPUT_CHECK.out.fastqs )
        ch_versions = ch_versions.mix( FALCO.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( FALCO.out.txt.collect{it[1]}.ifEmpty([]) )
    } else {
        FASTQC ( INPUT_CHECK.out.fastqs )
        ch_versions = ch_versions.mix( FASTQC.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( FASTQC.out.zip.collect{it[1]}.ifEmpty([]) )
    }

    //
    // SUBWORKFLOW: Read preprocessing (clipping, merging, fastq trimming etc. )
    //

    if ( !params.skip_preprocessing ) {
        PREPROCESSING ( INPUT_CHECK.out.fastqs, adapterlist )
        ch_reads_for_mapping = PREPROCESSING.out.reads
        ch_versions          = ch_versions.mix( PREPROCESSING.out.versions )
        ch_multiqc_files     = ch_multiqc_files.mix( PREPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    } else {
        ch_reads_for_mapping = INPUT_CHECK.out.fastqs
    }

    //
    // SUBWORKFLOW: Reference mapping
    //
    ch_reference_for_mapping = REFERENCE_INDEXING.out.reference
            .map{
                meta, fasta, fai, dict, index, circular_target ->
                [ meta, index ]
            }

    MAP ( ch_reads_for_mapping, ch_reference_for_mapping )

    ch_versions       = ch_versions.mix( MAP.out.versions )
    ch_multiqc_files  = ch_multiqc_files.mix( MAP.out.mqc.collect{it[1]}.ifEmpty([]) )

    //
    //  MODULE: indexing of user supplied input BAMs
    //

    SAMTOOLS_INDEX ( INPUT_CHECK.out.bams )
    ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions )

    if ( params.fasta_largeref )
        ch_bams_from_input = INPUT_CHECK.out.bams.join( SAMTOOLS_INDEX.out.csi )
    else {
        ch_bams_from_input = INPUT_CHECK.out.bams.join( SAMTOOLS_INDEX.out.bai )
    }


    //
    // MODULE: flagstats of user supplied input BAMs
    //
    ch_bam_bai_input = INPUT_CHECK.out.bams
                            .join(SAMTOOLS_INDEX.out.bai)

    SAMTOOLS_FLAGSTATS_BAM_INPUT ( ch_bam_bai_input )
    ch_versions = ch_versions.mix( SAMTOOLS_FLAGSTATS_BAM_INPUT.out.versions )

    //
    // SUBWORKFLOW: bam filtering (length, mapped/unmapped, quality etc.)
    //

    if ( params.run_bamfiltering || params.run_metagenomicscreening ) {

        ch_mapped_for_bamfilter = MAP.out.bam
                                    .join(MAP.out.bai)
                                    .mix(ch_bams_from_input)

        FILTER_BAM ( ch_mapped_for_bamfilter )
        ch_bamfiltered_for_deduplication = FILTER_BAM.out.genomics
        ch_bamfiltered_for_metagenomics  = FILTER_BAM.out.metagenomics
        ch_versions                      = ch_versions.mix( FILTER_BAM.out.versions )
        ch_multiqc_files                 = ch_multiqc_files.mix( FILTER_BAM.out.mqc.collect{it[1]}.ifEmpty([]) )

    } else {
        ch_bamfiltered_for_deduplication = MAP.out.bam
                                                .join(MAP.out.bai)
                                                .mix(ch_bams_from_input)
    }

    ch_reads_for_deduplication = ch_bamfiltered_for_deduplication

    //
    // SUBWORKFLOW: genomic BAM deduplication
    //

    ch_fasta_for_deduplication = REFERENCE_INDEXING.out.reference
        .multiMap{
            meta, fasta, fai, dict, index, circular_target ->
            fasta:      [ meta, fasta ]
            fasta_fai:  [ meta, fai ]
        }

    if ( !params.skip_deduplication ) {
        DEDUPLICATE( ch_reads_for_deduplication, ch_fasta_for_deduplication.fasta, ch_fasta_for_deduplication.fasta_fai )
        ch_dedupped_bams = DEDUPLICATE.out.bam
            .join( DEDUPLICATE.out.bai )
        ch_dedupped_flagstat = DEDUPLICATE.out.flagstat
        ch_versions                   = ch_versions.mix( DEDUPLICATE.out.versions )

    } else {
        ch_dedupped_bams     = ch_reads_for_deduplication
        ch_dedupped_flagstat = Channel.empty()
    }

    //
    // MODULE QUALIMAP
    //

    if ( !params.skip_qualimap ) {
        ch_snp_capture_bed = REFERENCE_INDEXING.out.snp_capture_bed
            .map{
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }
        ch_qualimap_input = ch_dedupped_bams
            .map {
            meta, bam, bai ->
                [ meta, bam ]
            }
            .map {
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine(
                by: 0,
                ch_snp_capture_bed
            )
            .branch {
                ignore_meta, meta, bam, meta2, snp_capture_bed ->
                withbed: snp_capture_bed != ""
                nobed: true
            }
        ch_qualimap_input_with = ch_qualimap_input.withbed
            .multiMap{
                ignore_meta, meta, bam, meta2, snp_capture_bed ->
                bam:             [ meta, bam ]
                snp_capture_bed: [ snp_capture_bed ]
            }

        QUALIMAP_BAMQC_WITHBED( ch_qualimap_input_with.bam, ch_qualimap_input_with.snp_capture_bed )
        ch_qualimap_input_without = ch_qualimap_input.nobed
            .map{
                ignore_meta, meta, bam, meta2, snp_capture_bed ->
                [ meta, bam ]
            }

        QUALIMAP_BAMQC_NOBED( ch_qualimap_input_without, [] )
        ch_qualimap_output = QUALIMAP_BAMQC_WITHBED.out.results.mix( QUALIMAP_BAMQC_NOBED.out.results )
        ch_versions = ch_versions.mix( QUALIMAP_BAMQC_NOBED.out.versions ).mix( QUALIMAP_BAMQC_WITHBED.out.versions )
    }

    //
    // MODULE: remove reads mapping to the host from the raw fastq
    //
    if ( params.run_host_removal ) {
        // Preparing bam channel for host removal to be combined with the input fastq channel
        // The bam channel consist of [meta, bam, bai] and in the meta we have in addition 'single_end' always set as TRUE and 'reference' set
        // To be able to join it with fastq channel, we need to remove them from the meta (done in map) and stored in new_meta
        ch_bam_for_host_removal= MAP.out.bam.join(MAP.out.bai)
                                            .map{
                                                meta, bam, bai ->
                                                new_meta = meta.clone().findAll{ it.key !in [ 'single_end', 'reference' ] }
                                                [ new_meta, meta, bam, bai ]
                                                }
        // Preparing fastq channel for host removal to be combined with the bam channel
        // The meta of the fastq channel contains additional fields when compared to the meta from the bam channel: lane, colour_chemistry,
        // and not necessarily matching single_end. Those fields are dropped of the meta in the map and stored in new_meta
        ch_fastqs_for_host_removal= INPUT_CHECK.out.fastqs.map{
                                                        meta, fastqs ->
                                                        new_meta = meta.clone().findAll{ it.key !in [ 'lane', 'colour_chemistry', 'single_end' ] }
                                                        [ new_meta, meta, fastqs ]
                                                    }
        // We join the bam and fastq channel with now matching metas (new_meta) referred as meta_join
        // and remove the meta_join from the final channel, keeping the original metas for the bam and the fastqs
        ch_input_for_host_removal = ch_bam_for_host_removal.join(ch_fastqs_for_host_removal)
                                                    .map{
                                                        meta_join, meta_bam, bam, bai, meta_fastq, fastqs ->
                                                        [ meta_bam, bam, bai, meta_fastq, fastqs]
                                                    }

        HOST_REMOVAL ( ch_input_for_host_removal )

        ch_versions = ch_versions.mix( HOST_REMOVAL.out.versions )
    }

    //
    // Section: Metagenomics screening
    //

    if( params.run_metagenomicscreening ) {
        ch_bamfiltered_for_metagenomics = ch_bamfiltered_for_metagenomics
            .map{ meta, fastq ->
                [meta+['single_end':true], fastq]
            }

        // Check if a complexity filter is wanted?
        if ( params.run_metagenomics_complexityfiltering ) {
            METAGENOMICS_COMPLEXITYFILTER( ch_bamfiltered_for_metagenomics )
            ch_reads_for_metagenomics = METAGENOMICS_COMPLEXITYFILTER.out.fastq
            ch_versions = ch_versions.mix(METAGENOMICS_COMPLEXITYFILTER.out.versions.first())
            ch_multiqc_files = ch_multiqc_files.mix(METAGENOMICS_COMPLEXITYFILTER.out.fastq.collect{it[1]}.ifEmpty([]))
        } else {
            ch_reads_for_metagenomics = ch_bamfiltered_for_metagenomics
        }
    }

    //
    // MODULE: MTNUCRATIO
    //

    if ( params.run_mtnucratio ) {
        ch_mito_header = REFERENCE_INDEXING.out.mitochondrion_header
            .map{
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }
        mtnucratio_input = ch_dedupped_bams
            .map {
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine(
                by: 0,
                ch_mito_header
            )
            .multiMap{
                ignore_meta, meta, bam, bai, meta2, mito_header ->
                bam:         [ meta, bam ]
                mito_header: [ meta2, mito_header ]
            }

        MTNUCRATIO( mtnucratio_input.bam, mtnucratio_input.mito_header.map{ it[1] } )
        ch_multiqc_files = ch_multiqc_files.mix(MTNUCRATIO.out.mtnucratio.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix( MTNUCRATIO.out.versions )
    }

    //
    // MODULE: ENDORSPY (raw, filtered, deduplicated)
    //

    ch_flagstat_for_endorspy_raw    = MAP.out.flagstat
                                            .mix( SAMTOOLS_FLAGSTATS_BAM_INPUT.out.flagstat )

    if ( params.run_bamfiltering & !params.skip_deduplication ) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
                                                .join (FILTER_BAM.out.flagstat)
                                                .join (DEDUPLICATE.out.flagstat)
    } else if ( params.run_bamfiltering & params.skip_deduplication ) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
                                                        .join (FILTER_BAM.out.flagstat)
                                                        .map{
                                                            meta, flags_raw, flags_filtered ->
                                                            [ meta, flags_raw, flags_filtered, [] ]
                                                            }
    } else if ( !params.run_bamfiltering & !params.skip_deduplication) {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
                                                        .join (DEDUPLICATE.out.flagstat)
                                                        . map{
                                                            meta, flags_raw, flags_dedup ->
                                                            [ meta, flags_raw, [], flags_dedup ]
                                                            }
    } else {
        ch_for_endorspy = ch_flagstat_for_endorspy_raw
                                                        .map {
                                                            meta, flags_raw ->
                                                            [ meta, flags_raw, [], [] ]
                                                        }
    }

    ENDORSPY ( ch_for_endorspy )

    ch_versions       = ch_versions.mix( ENDORSPY.out.versions )
    ch_multiqc_files  = ch_multiqc_files.mix( ENDORSPY.out.json.collect{it[1]}.ifEmpty([]) )

    //
    // MODULE: PreSeq
    //

    if ( !params.mapstats_skip_preseq && params.mapstats_preseq_mode == 'c_curve') {
        PRESEQ_CCURVE(ch_reads_for_deduplication.map{[it[0],it[1]]})
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_CCURVE.out.c_curve.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix( PRESEQ_CCURVE.out.versions )
    } else ( !params.mapstats_skip_preseq && params.mapstats_preseq_mode == 'lc_extrap') {
        PRESEQ_LCEXTRAP(ch_reads_for_deduplication.map{[it[0],it[1]]})
        ch_multiqc_files = ch_multiqc_files.mix(PRESEQ_LCEXTRAP.out.lc_extrap.collect{it[1]}.ifEmpty([]))
        ch_versions = ch_versions.mix( PRESEQ_LCEXTRAP.out.versions )
    }


    //
    // MODULE: Bedtools coverage
    //

    if ( params.run_bedtools_coverage ) {

        ch_bedtools_feature = REFERENCE_INDEXING.out.bedtools_feature
                                .map{
                                    WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
                                }

        ch_bedtools_prep = ch_dedupped_bams
                    .map {
                        WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
                    }
                    .combine(
                        by: 0,
                        ch_bedtools_feature
                    )
                    .map{
                        ignore_meta, meta, bam, bai, meta2, bedtools_feature ->
                        [ meta, bedtools_feature, bam, bai ]
                    }
                    .branch{
                        meta, bedtools_feature, bam, bai ->
                        withfeature: bedtools_feature != ""
                        nobed: true
                    }

        // Running samtools view to get header
        ch_bedtools_input = ch_bedtools_prep.withfeature
                        .multiMap{
                            meta, bedtools_feature, bam, bai ->
                            bam: [ meta, bam, bai ]
                            withfeature: [ meta, bedtools_feature, bam ]
                        }

        SAMTOOLS_VIEW_GENOME( ch_bedtools_input.bam )

        ch_genome_for_bedtools = SAMTOOLS_VIEW_GENOME.out.genome

        BEDTOOLS_COVERAGE_BREADTH(ch_bedtools_input.withfeature, ch_genome_for_bedtools)
        BEDTOOLS_COVERAGE_DEPTH(ch_bedtools_input.withfeature, ch_genome_for_bedtools)

        ch_versions = ch_versions.mix( SAMTOOLS_VIEW_GENOME.out.versions )
        ch_versions = ch_versions.mix( BEDTOOLS_COVERAGE_BREADTH.out.versions )
        ch_versions = ch_versions.mix( BEDTOOLS_COVERAGE_DEPTH.out.versions )
    }

    //
    // SUBWORKFLOW: Calculate Damage
    //

    ch_fasta_for_damagecalculation = REFERENCE_INDEXING.out.reference
        .multiMap{
            meta, fasta, fai, dict, index, circular_target ->
            fasta:      [ meta, fasta ]
            fasta_fai:  [ meta, fai ]
        }

    if ( !params.skip_damage_calculation ) {
        CALCULATE_DAMAGE( ch_dedupped_bams, ch_fasta_for_damagecalculation.fasta, ch_fasta_for_damagecalculation.fasta_fai )
        ch_versions      = ch_versions.mix( CALCULATE_DAMAGE.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix(CALCULATE_DAMAGE.out.mqc.collect{it[1]}.ifEmpty([]))

    }

    //
    // SUBWORKFLOW: Contamination estimation
    //

    if ( params.run_contamination_estimation_angsd ) {
        contamination_input = ch_dedupped_bams

        ESTIMATE_CONTAMINATION( contamination_input, REFERENCE_INDEXING.out.hapmap )
        ch_versions      = ch_versions.mix( ESTIMATE_CONTAMINATION.out.versions )
        ch_multiqc_files = ch_multiqc_files.mix( ESTIMATE_CONTAMINATION.out.mqc.collect{it[1]}.ifEmpty([]) )
    }

    //
    // SUBWORKFLOW: aDNA Damage Manipulation

    if ( params.run_mapdamage_rescaling || params.run_pmd_filtering || params.run_trim_bam ) {
        MANIPULATE_DAMAGE( ch_dedupped_bams, ch_fasta_for_deduplication.fasta, REFERENCE_INDEXING.out.pmd_masked_fasta, REFERENCE_INDEXING.pmd_bed_for_masking )
        ch_multiqc_files       = ch_multiqc_files.mix( MANIPULATE_DAMAGE.out.flagstat.collect{it[1]}.ifEmpty([]) )
        ch_versions            = ch_versions.mix( MANIPULATE_DAMAGE.out.versions )
        ch_bams_for_genotyping = params.genotyping_source == 'rescaled' ? MANIPULATE_DAMAGE.out.rescaled : params.genotyping_source == 'pmd' ? MANIPULATE_DAMAGE.out.filtered : params.genotyping_source == 'trimmed' ? MANIPULATE_DAMAGE.out.trimmed : ch_dedupped_bams
    } else {
        ch_bams_for_genotyping = ch_dedupped_bams
    }

    //
    // MODULE: MultiQC
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    workflow_summary    = WorkflowEager.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowEager.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([])) // Replaced with custom mixing

    if ( !params.skip_qualimap ) {
        ch_multiqc_files = ch_multiqc_files.mix( ch_qualimap_output.collect{it[1]}.ifEmpty([]) )
    }

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
