//
// Subworkflow with functionality specific to the nf-core/eager pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input samplesheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //
    ch_samplesheet = Channel.fromSamplesheet("input")
                                    .map {
                                        meta, r1, r2, bam ->
                                            meta.single_end = meta.pairment == "single" ? true : false
                                            meta.id = meta.sample_id
                                        [ meta, r1, r2, bam ]
                                    }

    ch_samplesheet_for_branch = ch_samplesheet
                                    .branch {
                                        meta, r1, r2, bam ->
                                            bam: bam.toString().endsWith(".bam")
                                            fastq: true
                                    }

    ch_samplesheet_fastqs = ch_samplesheet_for_branch.fastq
                                .map {
                                    meta, r1, r2, bam ->
                                        reads = meta.single_end ? [ r1 ] : [ r1, r2 ]
                                    [ meta - meta.subMap('pairment', 'bam_reference_id'), reads ]
                                }

    ch_samplesheet_bams = ch_samplesheet_for_branch.bam
                            .map {
                                meta, r1, r2, bam ->
                                    meta.reference = meta.bam_reference_id
                                    meta.id_index = meta.bam_reference_id
                                [ meta - meta.subMap('pairment', 'bam_reference_id'), bam ]
                            }

    // Extra validation
    // - Only paired end specified when R2 provided
    // - No single-ended data allowed when using dedup
    ch_samplesheet_for_branch.fastq
        .map {
            meta, r1, r2, bam ->
                if ( meta.pairment == "single" && r2 != [] ) {
                    exit 1, "[nf-core] ERROR: Validation of 'input' file failed. Reads 2 cannot be provided when sequencing pairment is set to 'single'."
                }
                if ( meta.pairment == "paired" && r2 == [] ) {
                    exit 1, "[nf-core] ERROR: Validation of 'input' file failed. Reads 2 have to be provided when sequencing pairment is set to 'paired'."
                }
                if ( meta.pairment == "single" && params.deduplication_tool == "dedup" ) {
                    exit 1, "[nf-core] ERROR: Invalid input/parameter combination. '--deduplication_tool' cannot be 'dedup' on runs that include SE data. Use 'markduplicates' for runs with both SE and PE data or separate SE and PE data into separate runs."
                }
            [ meta, r1, r2, bam ]
        }

    // - Only single-ended specified for BAM files
    ch_samplesheet_for_branch.bam
        .map {
            meta, r1, r2, bam ->
                if ( meta.pairment == "paired" && bam != [] ) {
                    exit 1, "[nf-core] ERROR: Validation of 'input' file failed. Sequencing pairment has to be 'single' when BAM files are provided."
                }
            [ meta, r1, r2, bam ]
        }

    // - No single- and double-stranded libraries with same sample ID
    ch_samplesheet_test = ch_samplesheet
                            .map {
                                meta, r1, r2, bam ->
                                [ meta.subMap('sample_id'), meta.subMap('strandedness') ]
                            }
                            .groupTuple()
                            .map { meta, singlestrand ->
                                    if ( singlestrand.toList().unique().size() > 1 ) {
                                        exit 1, "[nf-core] ERROR: Validation of 'input' file failed. Sample IDs can only be identical if library strandedness is identical."
                                    }
                                [ meta, singlestrand ]
                            }

    emit:
    samplesheet_fastqs = ch_samplesheet_fastqs
    samplesheet_bams   = ch_samplesheet_bams
    versions           = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
    if ( !params.fasta                                     && !params.fasta_sheet ) { exit 1, "[nf-core/eager] ERROR: Neither FASTA file --fasta nor reference sheet --fasta_sheet have been provided."}
    if ( params.fasta                                      && params.fasta_sheet ) { exit 1, "[nf-core/eager] ERROR: A FASTA file --fasta and a reference sheet --fasta_sheet have been provided."}
    if ( params.preprocessing_adapterlist                  && params.preprocessing_skipadaptertrim ) { log.warn("[nf-core/eager] WARNING: --preprocessing_skipadaptertrim will override --preprocessing_adapterlist. Adapter trimming will be skipped!") }
    if ( params.deduplication_tool == 'dedup'              && ! params.preprocessing_excludeunmerged ) { exit 1, "[nf-core/eager] ERROR: Dedup can only be used on collapsed (i.e. merged) PE reads. For all other cases, please set --deduplication_tool to 'markduplicates'."}
    if ( params.bamfiltering_retainunmappedgenomicbam      && params.bamfiltering_mappingquality > 0  ) { exit 1, ("[nf-core/eager] ERROR: You cannot both retain unmapped reads and perform quality filtering, as unmapped reads have a mapping quality of 0. Pick one or the other functionality.") }
    if ( params.genotyping_source == 'trimmed'             && ! params.run_trim_bam                   ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'trimmed' unless BAM trimming is turned on with `--run_trim_bam`.") }
    if ( params.genotyping_source == 'pmd'                 && ! params.run_pmd_filtering              ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'pmd' unless PMD-filtering is ran.") }
    if ( params.genotyping_source == 'rescaled'            && ! params.run_mapdamage_rescaling        ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'rescaled' unless aDNA damage rescaling is ran.") }
    if ( params.metagenomics_complexity_tool == 'prinseq'  && params.metagenomics_prinseq_mode == 'dust' && params.metagenomics_complexity_entropy != 0.3 ) {
        if (params.metagenomics_prinseq_dustscore == 0.5) { exit 1, ("[nf-core/eager] ERROR: Metagenomics: You picked PRINSEQ++ with 'dust' mode but provided an entropy score. Please specify a dust filter threshold using the --metagenomics_prinseq_dustscore flag") }
    }
    if ( params.metagenomics_complexity_tool == 'prinseq'  && params.metagenomics_prinseq_mode == 'entropy' && params.metagenomics_prinseq_dustscore != 0.5 ) {
        if (params.metagenomics_complexity_entropy == 0.3) { exit 1, ("[nf-core/eager] ERROR: Metagenomics: You picked PRINSEQ++ with 'entropy' mode but provided a dust score. Please specify an entropy filter threshold using the --metagenomics_complexity_entropy flag") }
    }
    if ( params.run_genotyping                        && ! params.genotyping_tool                ) { exit 1, ("[nf-core/eager] ERROR: --run_genotyping was specified, but no --genotyping_tool was specified.") }
    if ( params.run_genotyping                        && ! params.genotyping_source              ) { exit 1, ("[nf-core/eager] ERROR: --run_genotyping was specified, but no --genotyping_source was specified.") }
    if ( params.genotyping_source == 'trimmed'        && ! params.run_trim_bam                   ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'trimmed' unless BAM trimming is turned on with `--run_trim_bam`.") }
    if ( params.genotyping_source == 'pmd'            && ! params.run_pmd_filtering              ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'pmd' unless PMD-filtering is ran.") }
    if ( params.genotyping_source == 'rescaled'       && ! params.run_mapdamage_rescaling        ) { exit 1, ("[nf-core/eager] ERROR: --genotyping_source cannot be 'rescaled' unless aDNA damage rescaling is ran.") }
    if ( params.fasta && params.run_genotyping && params.genotyping_tool == 'pileupcaller' && ! (params.genotyping_pileupcaller_bedfile || params.genotyping_pileupcaller_snpfile ) ) { exit 1, ("[nf-core/eager] ERROR: Genotyping with pileupcaller requires both '--genotyping_pileupcaller_bedfile' AND '--genotyping_pileupcaller_snpfile' to be provided.") }

}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
//
// Get attribute from genome config file e.g. fasta
//
def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

def grabUngzippedExtension(infile) {

    def split_name = infile.toString().tokenize('.')
    def output = split_name.reverse().first() == 'gz' ? split_name.reverse()[1,0].join('.') : split_name.reverse()[0]

    return '.' + output

}

/*
This function can be applied to the contents of a channel row-by-row with a .map operator.
It assumes that the first element of the channel row is a map of metadata. It will then create a new metadata map
that consists of the source_attributes of the original metadata map, but named after the corresponding target_attributes.
The new metadata map is then prepended to the channel row, becoming the new first element.
If the remove flag is set to true, then the original metadata map is removed from the channel row.

Example:
ch_my_channel=Channel.of( [ [id:'id', sample_id:'sample_id', single_end:true, reference:"hs37d5" ], "bam", "bai" ] )

ch_my_channel.map{ row -> addNewMetaFromAttributes(row, "id", "new_attribute", false) }
// This will create a new channel with the following rows:
[ [ new_attribute: 'id' ], [id:'id', sample_id:'sample_id', single_end:true, reference:"hs37d5" ], "bam", "bai" ]

ch_my_channel.map{ row -> addNewMetaFromAttributes(row, ["id", "single_end"], ["new_attribute", "endedness"] , true) }
// This will create a new channel with the following rows:
[ [ new_attribute: 'id' , endedness: true ], "bam", "bai" ]

*/
def addNewMetaFromAttributes( ArrayList row, Object source_attributes, Object target_attributes, boolean remove = false) {
    def meta = row[0]
    def meta2 = [:]

    // Read in target and source attributes and create a mapping between them
    // Option A: both attributes are Strings
    if ((source_attributes instanceof String) && (target_attributes instanceof String)) {
            meta2[target_attributes] = meta[source_attributes]

    } else if ((source_attributes instanceof List) && (target_attributes instanceof List)) {
        if (source_attributes.size() == target_attributes.size()) {
            for (int i = 0; i < source_attributes.size(); i++) {
                // Option B: Both are lists of same size
                meta2[target_attributes[i]] = meta[source_attributes[i]]
            }
        } else {
            // Option C: Both lists, but uneven. Error.
            throw new IllegalArgumentException("Error: The target_attributes and source_attributes lists do not have the same size.")
        }
    } else {
        // Option D: Not both the same type or acceptable types. Error.
        throw new IllegalArgumentException("Error: target_attributes and source_attributes must be of same type (both Strings or both Lists).")

    }

    def new_row = [ meta2 ] + row

    // If replace is true, then remove the old meta
    if (remove && new_row.size() > 1) {
        new_row.remove(1)
    }

    return new_row
}
