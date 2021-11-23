/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEager.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [
        params.input, params.multiqc_config, params.fasta,
        params.bwa_index_dir, params.bt2_index_dir, params.fasta_fai,
        params.fasta_dict
    ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/nf-core_eager_dummy.txt", checkIfExists: true)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//

// Set additional custom parameters from modules config
def gunzip_options                           = modules['fasta_gunzip']
def bwa_index_options                        = modules['bwa_index_options']
def bowtie2_build_options                    = modules['bowtie2_build_options']
def samtools_faidx_options                   = modules['samtools_faidx_options']
def gatk4_createsequencedictionary_options   = modules['gatk4_createsequencedictionary_options']

println "B4" + gunzip_options

if ( params.save_reference ) {
    gunzip_options['publish_files']                         = true
    bwa_index_options['publish_files']                      = true
    bowtie2_build_options['publish_files']                  = true
    samtools_faidx_options['publish_files']                 = true
    gatk4_createsequencedictionary_options['publish_files'] = true
}

println "A4" + gunzip_options

include { INPUT_CHECK  } from '../subworkflows/local/input_check' addParams( options: [:] )
include { INDEX_FASTA  } from '../subworkflows/local/index_fasta' addParams( gunzip_options: gunzip_options, bwa_index_options: bwa_index_options, bowtie2_build_options: bowtie2_build_options, samtools_faidx_options: samtools_faidx_options, gatk4_createsequencedictionary_options: gatk4_createsequencedictionary_options )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'  addParams( options: [publish_files : ['_versions.yml':'']] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow EAGER {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    INDEX_FASTA (
        ch_dummy_file
    )

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowEager.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
