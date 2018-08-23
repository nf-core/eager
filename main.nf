#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/EAGER2
========================================================================================
 EAGER2 Analysis Pipeline. Started 2018-06-05
 #### Homepage / Documentation
 https://github.com/nf-core/eager
 #### Authors
 Alexander Peltzer apeltzer <alex.peltzer@gmail.com> - https://github.com/apeltzer>
 James A. Fellows Yates <jfy133@gmail.com> - https://github.com/jfy133
 Stephen Clayton <clayton@shh.mpg.de> - https://github.com/sc13-bioinf
----------------------------------------------------------------------------------------
*/


/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.bwa_index = false
params.seq_dict = false
params.fasta_index = false
params.saveReference = false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
Channel.fromPath("${params.fasta}")
    .ifEmpty { exit 1, "No genome specified! Please specify one with --fasta or --bwa_index"}
    .into {ch_fasta_for_bwa_indexing;ch_fasta_for_faidx_indexing;ch_fasta_for_dict_indexing}

//AWSBatch sanity checking
if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { read_files_fastqc; read_files_trimming }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { ch_read_files_clip; ch_read_files_fastqc }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { ch_read_files_clip; ch_read_files_fastqc }
}



// Header log info
log.info "========================================="
log.info " nf-core/eager v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container']    = workflow.container
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


/*
 * Parse software version numbers: TODO testing this
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into ch_software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    AdapterRemoval -version > v_adapterremoval.txt
    picard MarkDuplicates -version &> v_markduplicates.txt  || true
    dedup -h > v_dedup.txt
    #angsd > v_angsd.txt
    #realignsamfile > v_circularmapper.txt
    #schmutzi > v_schmutzi.txt
    gatk BaseRecalibrator --version &> v_gatk.txt
    qualimap --version > v_qualimap.txt
    vcf2genome > v_vcf2genome.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}

/* 
* Create BWA indices if they are not present
*/ 

if(!params.bwa_index && params.fasta && params.aligner == 'bwa'){
    process makeBWAIndex {
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_bwa_indexing

        output:
        file "*.{amb,ann,bwt,pac,sa,fasta,fa}" into ch_bwa_index

        script:
        """
        bwa index $fasta
        """
    }
}

/*
 * PREPROCESSING - Index Fasta file
 */
if(!params.fasta_index && params.fasta && params.aligner == 'bwa'){
    process makeFastaIndex {
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_faidx_indexing

        output:
        file "${fasta}.fai" into ch_fasta_faidx_index

        script:
        """
        samtools faidx $fasta
        """
    }
}

/*
 * PREPROCESSING - Create Sequence Dictionary for FastA
 */
if(!params.seq_dict && params.fasta){
    process makeSeqDict {
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_dict_indexing

        output:
        file "*.dict" into ch_seq_dict

        script:
        """
        picard CreateSequenceDictionary R=$fasta O="${fasta.baseName}.dict"
        """
    }
}


/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/01-FastQC", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from ch_read_files_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q $reads
    """
}

/*
 * STEP 2 - Adapter Clipping / Read Merging
 */


// process adapter_removal {
//     tag "$name"
//     publishDir "${params.outdir}/02-Merging", mode: 'copy'

//     input:
//     set val(name), file(reads) from ch_read_files_clip

//     output:
//     file "*.combined.fq.gz" into ch_clipped_reads

//     script:
//     prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
//     """
//     AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --baseName ${prefix} --gzip --threads ${process.cpus} --trimns --trimqualities --adapter1 ${params.clip.forward_adaptor} --adapter2 ${params.clip.reverse_adaptor} --minlength ${params.clip.readlength} --minquality ${params.clip.min_read_quality} --minadapteroverlap ${params.min_adap_overlap} --collapse
//     #Fix Prefixes
//     AdapterRemovalFixPrefix  TODO
//     #Combine files
//     zcat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > ${prefix}.combined.fq.gz
//     """
// }

// process adapter_removal_fixprefix {
//       tag "$name"
//       publishDir "${params.outdir}/02-Merging", mode: 'copy'

//       input:
//       set val(name), file(reads) from ch_clipped_reads

//       output:
//       file "*.fastq.prefixed.gz" into ch_mappable_reads

//       script:
//       '''
//       AdapterRemovalFixPrefix ${reads} ${reads}.fastq.prefixed.gz
//       '''
// }



/*
Step 3: Mapping with BWA, CircularMapper
Step 4: Conversion to BAM; sorting
Step 5: Keep unmapped/remove unmapped reads
Step 5.1: Preseq
Step 5.2: DMG Assessment
Step 5.3: Qualimap (before or after Dedup?)
Step 6: DeDup / MarkDuplicates
Step 7: angsd
Step 7: GATK
Step 8: vcf2genome

*/


/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from ch_fastqc_results.collect()
    file ('software_versions/*') from ch_software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    """
}



/*
 * STEP 3 - Output Description HTML
 */
process output_documentation {
    tag "$prefix"
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}



/**
Useful functionality, e.g. help messages etc
* 
*/ 


def helpMessage() {
    log.info"""
    =========================================
    EAGER2 v${params.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/EAGER2 --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference
      --bwa_index                   Path to BWA index

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}
