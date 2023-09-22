//
// This file holds several functions specific to the workflow/eager.nf in the nf-core/eager pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine

class WorkflowEager {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {

        genomeExistsError(params, log)

        if (!params.fasta) {
            Nextflow.error "Genome fasta file not specified with e.g. '--fasta genome.fa' or via a detectable config file."
        }
    }

    //
    // Get workflow summary for MultiQC
    //
    public static String paramsSummaryMultiqc(workflow, summary) {
        String summary_section = ''
        for (group in summary.keySet()) {
            def group_params = summary.get(group)  // This gets the parameters of that particular group
            if (group_params) {
                summary_section += "    <p style=\"font-size:110%\"><b>$group</b></p>\n"
                summary_section += "    <dl class=\"dl-horizontal\">\n"
                for (param in group_params.keySet()) {
                    summary_section += "        <dt>$param</dt><dd><samp>${group_params.get(param) ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>\n"
                }
                summary_section += '    </dl>\n'
            }
        }

        String yaml_file_text  = "id: '${workflow.manifest.name.replace('/', '-')}-summary'\n"
        yaml_file_text        += "description: ' - this information is collected when the pipeline is started.'\n"
        yaml_file_text        += "section_name: '${workflow.manifest.name} Workflow Summary'\n"
        yaml_file_text        += "section_href: 'https://github.com/${workflow.manifest.name}'\n"
        yaml_file_text        += "plot_type: 'html'\n"
        yaml_file_text        += 'data: |\n'
        yaml_file_text        += "${summary_section}"
        return yaml_file_text
    }

    //
    // Generate methods description for MultiQC
    //

    public static String toolCitationText(params) {

        // TODO Optionally add in-text citation tools to this list.
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

    public static String toolBibliographyText(params) {

        // TODO Optionally add bibliographic entries to this list.
        // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
        // Uncomment function in methodsDescriptionText to render in MultiQC report
        def reference_text = [
                "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
                "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
            ].join(' ').trim()

        return reference_text
    }

    public static String methodsDescriptionText(run_workflow, mqc_methods_yaml, params) {
        // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
        def meta = [:]
        meta.workflow = run_workflow.toMap()
        meta['manifest_map'] = run_workflow.manifest.toMap()

        // Pipeline DOI
        meta["doi_text"] = meta.manifest_map.doi ? "(doi: <a href=\'https://doi.org/${meta.manifest_map.doi}\'>${meta.manifest_map.doi}</a>)" : ""
        meta["nodoi_text"] = meta.manifest_map.doi ? "": "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"
        // Tool references
        meta["tool_citations"] = ""
        meta["tool_bibliography"] = ""

        // TODO Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
        //meta["tool_citations"] = toolCitationText(params).replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
        //meta["tool_bibliography"] = toolBibliographyText(params)


        def methods_text = mqc_methods_yaml.text

        def engine =  new SimpleTemplateEngine()
        def description_html = engine.createTemplate(methods_text).make(meta)

        return description_html
    }

    //
    // Exit pipeline if incorrect --genome key provided
    //
    private static void genomeExistsError(params, log) {
        if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
            def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
                "  Currently, the available genome keys are:\n" +
                "  ${params.genomes.keySet().join(", ")}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
            Nextflow.error(error_string)
        }
    }

    def public static String grabUngzippedExtension(infile) {

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
    def public static ArrayList addNewMetaFromAttributes( ArrayList row, Object source_attributes, Object target_attributes, boolean remove = false) {
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
}
