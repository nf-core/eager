/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEager.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check failing parameter combinations
if ( params.bamfiltering_retainunmappedgenomicbam && params.bamfiltering_mappingquality > 0  ) { exit 1, ("[nf-core/eager] ERROR: You cannot both retain unmapped reads and perform quality filtering, as unmapped reads have a mapping quality of 0. Pick one or the other functionality.") }
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
include { INPUT_CHECK        } from '../subworkflows/local/input_check'
include { REFERENCE_INDEXING } from '../subworkflows/local/reference_indexing'
include { PREPROCESSING      } from '../subworkflows/local/preprocessing'
include { MAP                } from '../subworkflows/local/map'
include { FILTER_BAM         } from '../subworkflows/local/bamfiltering.nf'
include { DEDUPLICATE        } from '../subworkflows/local/deduplicate'
include { METAGENOMICS_COMPLEXITYFILTER } from '../subworkflows/local/metagenomics_complexityfilter'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SAMTOOLS_INDEX              } from '../modules/nf-core/samtools/index/main'
include { PRESEQ_CCURVE               } from '../modules/nf-core/preseq/ccurve/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { FALCO                       } from '../modules/nf-core/falco/main'
include { MTNUCRATIO                  } from '../modules/nf-core/mtnucratio/main'
include { ENDORSPY                    } from '../modules/nf-core/endorspy/main'

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
        ch_input
    )
    ch_versions = ch_versions.mix( INPUT_CHECK.out.versions )

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

    MAP ( ch_reads_for_mapping, REFERENCE_INDEXING.out.reference.map{meta, fasta, fai, dict, index -> [meta, index]} )
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
            meta, fasta, fai, dict, index ->
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
        mtnucratio_input = ch_dedupped_bams
        .map {
            meta, bam, bai ->
            [ meta, bam ]
        }

        MTNUCRATIO( mtnucratio_input, params.mitochondrion_header )
        ch_multiqc_files = ch_multiqc_files.mix(MTNUCRATIO.out.mtnucratio.collect{it[1]}.ifEmpty([]))
        ch_versions      = ch_versions.mix( MTNUCRATIO.out.versions )
    }

    //
    // MODULE: ENDORSPY (raw, filtered, deduplicated)
    //

    //  TODO: what if not mapping executed due to bam input, figure it out later

    ch_for_endorspy_rawbam    = MAP.out.flagstat //CHECK IF INPUT BAM, need to add a flagstat for the input bam (DANGEROUS if the provided bam file does not contain all reads --> wrong endogenous calculations!), else MAP.out.flagstats
    ch_for_endorspy_filterbam = params.run_bamfiltering ? FILTER_BAM.out.flagstat : Channel.empty()
    ch_for_endorspy_dedupped  = !params.skip_deduplication ? DEDUPLICATE.out.flagstat : Channel.empty()

    ch_for_endorspy = ch_for_endorspy_rawbam
                                    .join(ch_for_endorspy_filterbam.ifEmpty([ [], [] ]))
                                    .join(ch_for_endorspy_dedupped.ifEmpty([ [], [] ]))

    ENDORSPY ( ch_for_endorspy )
    ENDORSPY.out.json.dump()

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
    // MODULE: MultiQC
    //

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    workflow_summary    = WorkflowEager.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowEager.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
