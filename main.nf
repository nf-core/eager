#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/eager
========================================================================================
 EAGER Analysis Pipeline. Started 2018-06-05
 #### Homepage / Documentation
 https://github.com/nf-core/eager
 #### Authors
 Alexander Peltzer apeltzer <alex.peltzer@gmail.com> - https://github.com/apeltzer>
 James A. Fellows Yates <jfy133@gmail.com> - https://github.com/jfy133
 Stephen Clayton <clayton@shh.mpg.de> - https://github.com/sc13-bioinf
========================================================================================
*/

def helpMessage() {
    log.info"""
    =========================================
    eager v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/eager --reads '*_R{1,2}.fastq.gz' -profile docker

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
/*
 * SET UP CONFIGURATION VARIABLES
 */

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Configurable variables
params.name = false
params.singleEnd = false
params.genome = "Custom"
params.snpcapture = false
params.bedfile = ''
params.fasta = false
params.bwa_index = false
params.seq_dict = false
params.fasta_index = false
params.saveReference = false
params.udg = false 
params.udg_type = 'Half'

params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

// Skipping parts of the pipeline for impatient users
params.skip_preseq = false
params.skip_damage_calculation = false
params.skip_qualimap = false
params.skip_deduplication = false

//Complexity filtering reads
params.complexity_filter = false
params.complexity_filter_poly_g_min = 10

//Read clipping and merging parameters
params.clip_forward_adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
params.clip_reverse_adaptor = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
params.clip_readlength = 30
params.clip_min_read_quality = 20
params.min_adap_overlap = 1

//Read mapping parameters (default = BWA aln default)
params.bwaalnn = 0.04
params.bwaalnk = 2
params.bwaalnl = 32

//BAM Filtering steps (default = keep mapped and unmapped in BAM file)
params.bam_keep_mapped_only = false
params.bam_keep_all = true
params.bam_filter_reads = false
params.bam_mapping_quality_threshold = 0

//DamageProfiler settings
params.damageprofiler_length = 100
params.damageprofiler_threshold = 15

//DeDuplication settings
params.dedupper = 'dedup' //default value dedup

//Preseq settings
params.preseq_step_size = 1000

//PMDTools settings
params.run_pmdtools = false
params.pmdtools_range = 10
params.pmdtools_threshold = 3
params.pmdtools_reference_mask = ''
params.pmdtools_max_reads = 10000



multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

// Validate inputs
Channel.fromPath("${params.fasta}")
    .ifEmpty { exit 1, "No genome specified! Please specify one with --fasta or --bwa_index"}
    .into {ch_fasta_for_bwa_indexing;ch_fasta_for_faidx_indexing;ch_fasta_for_dict_indexing; ch_fasta_for_bwa_mapping; ch_fasta_for_damageprofiler; ch_fasta_for_qualimap; ch_fasta_for_pmdtools}

//AWSBatch sanity checking

if(workflow.profile == 'awsbatch'){
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    if (!workflow.workDir.startsWith('s3') || !params.outdir.startsWith('s3')) exit 1, "Specify S3 URLs for workDir and outdir parameters on AWSBatch!"
}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */

if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .dump(tag:'input')
            .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filtering }
            
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .dump(tag:'input')
            .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filtering }
            
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .dump(tag:'input')
        .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filtering }
        
}



// Header log info
log.info "========================================="
log.info " nf-core/eager v${params.pipelineVersion}"
log.info "========================================="
def summary = [:]
summary['Pipeline Name']  = 'nf-core/eager'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output dir']   = params.outdir
summary['Working dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Working dir']    = workflow.workDir
summary['Output dir']     = params.outdir
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region'] = params.awsregion
   summary['AWS Queue'] = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="


def create_workflow_summary(summary) {

    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-eager-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/eager Workflow Summary'
    section_href: 'https://github.com/nf-core/eager'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version &> v_pipeline.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    AdapterRemoval -version  &> v_adapterremoval.txt 2>&1 || true
    picard MarkDuplicates --version &> v_markduplicates.txt  2>&1 || true
    dedup -v &> v_dedup.txt 2>&1 || true
    preseq &> v_preseq.txt 2>&1 || true
    gatk BaseRecalibrator --version 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    vcf2genome &> v_vcf2genome.txt 2>&1 || true
    fastp --version &> v_fastp.txt 2>&1 || true
    bam --version &> v_bamutil.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/* 
* Create BWA indices if they are not present
*/ 

if(!params.bwa_index && params.fasta && params.aligner == 'bwa'){
    process makeBWAIndex {
        tag {fasta}
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
        tag {fasta}
        publishDir path: "${params.outdir}/reference_genome", saveAs: { params.saveReference ? it : null }, mode: 'copy'

        input:
        file fasta from ch_fasta_for_faidx_indexing

        output:
        file "${fasta}.fai" into ch_fasta_faidx_index
        file "${fasta}"

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
        tag {fasta}
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


/* STEP 2.0 - FastP
* Optional poly-G complexity filtering step before read merging/adapter clipping etc
* Note: Clipping, Merging, Quality Trimning are turned off here - we leave this to adapter removal itself!
*/

process fastp {
    tag "$name"
    publishDir "${params.outdir}/01-FastP", mode: 'copy'

    when: params.complexity_filter

    input:
    set val(name), file(reads) from ch_read_files_complexity_filtering

    output:
    set val(name), file("*pG.fq.gz") into ch_clipped_reads_complexity_filtered
    file("*.json") into ch_fastp_for_multiqc

    script:
    if(params.singleEnd){
    """
    fastp --in1 ${reads[0]} --out1 "${reads[0].baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json "${reads[0].baseName}"_fastp.json 
    """
    } else {
    """
    fastp --in1 ${reads[0]} --in2 ${reads[1]} --out1 "${reads[0].baseName}.pG.fq.gz" --out2 "${reads[1].baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json "${reads[0].baseName}"_fastp.json 
    """
    }
}


/*
 * STEP 2 - Adapter Clipping / Read Merging
 */


process adapter_removal {
    tag "$name"
    publishDir "${params.outdir}/02-Merging", mode: 'copy'

    input:
    set val(name), file(reads) from ( params.complexity_filter ? ch_clipped_reads_complexity_filtered : ch_read_files_clip )

    output:
    file "*.combined*.gz" into (ch_clipped_reads, ch_clipped_reads_for_fastqc)
    file "*.settings" into ch_adapterremoval_logs

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    //Readprefixing only required for PE data with merging
    fixprefix = (params.singleEnd) ? "" : "AdapterRemovalFixPrefix ${prefix}.combined.fq.gz ${prefix}.combined.prefixed.fq.gz"
    
    if( !params.singleEnd ){
    """
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${prefix} --gzip --threads ${task.cpus} --trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap} --collapse
    #Combine files
    zcat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > ${prefix}.combined.fq.gz
    ${fixprefix}
    rm ${prefix}.combined.fq.gz
    """
    } else {
    """
    AdapterRemoval --file1 ${reads[0]} --basename ${prefix} --gzip --threads ${task.cpus} --trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} 
    # Pseudo-Combine
    mv *.truncated.gz ${prefix}.combined.fq.gz
    """
    }
}

/*
 * STEP 2.1 - FastQC after clipping/merging (if applied!)
 */
process fastqc_after_clipping {
    tag "${reads[0].baseName}"
    publishDir "${params.outdir}/01-FastQC/after_clipping", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    file(reads) from ch_clipped_reads_for_fastqc

    output:
    file "*_fastqc.{zip,html}" optional true into ch_fastqc_after_clipping

    script:
    """
    fastqc -q $reads
    """
}

//Close that channel
ch_fastqc_after_clipping.close()

/*
Step 3: Mapping with BWA, SAM to BAM, Sort BAM
*/

process bwa {
    tag "$prefix"
    publishDir "${params.outdir}/03-Mapping", mode: 'copy'

    input:
    file(reads) from ch_clipped_reads
    file "*" from ch_bwa_index
    file fasta from ch_fasta_for_bwa_mapping

    output:
    file "*.sorted.bam" into ch_mapped_reads_idxstats,ch_mapped_reads_filter,ch_mapped_reads_preseq, ch_mapped_reads_damageprofiler
    

    script:
    prefix = reads[0].toString() - ~/(_R1)?(\.combined\.)?(prefixed)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """ 
    bwa aln -t ${task.cpus} $fasta $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f "${reads.baseName}.sai"
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta "${reads.baseName}".sai $reads | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    """
}

//TODO Multi Lane BAM merging here (!)


/*
* Step 4 - IDXStats
*/

process samtools_idxstats {
    tag "$prefix"
    publishDir "${params.outdir}/04-Samtools", mode: 'copy'

    input:
    file(bam) from ch_mapped_reads_idxstats

    output:
    file "*.stats" into ch_idxstats_for_multiqc

    script:
    prefix = "$bam" - ~/(\.bam)?$/
    """
    samtools flagstat $bam > ${prefix}.stats
    """
}


/*
* Step 5: Keep unmapped/remove unmapped reads
*/

process samtools_filter {
    tag "$prefix"
    publishDir "${params.outdir}/04-Samtools", mode: 'copy', pattern: "*.bam"

    input: 
    file bam from ch_mapped_reads_filter

    output:
    file "*filtered.bam" into ch_bam_filtered_qualimap, ch_bam_filtered_dedup, ch_bam_filtered_markdup, ch_bam_filtered_pmdtools, ch_bam_filtered_angsd, ch_bam_filtered_gatk

    when: "${params.bam_filter_reads}"

    script:
    prefix="$bam" - ~/(\.bam)?/

    if("${params.bam_keep_mapped_only}"){
    """
    samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
    samtools fastq -tn "${prefix}.unmapped.bam" | gzip > "${prefix}.unmapped.fq.gz"
    """
    } else {
    """
    samtools view -h $bam | samtools view - -@ ${task.cpus} -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
    """
    }  
}


/*
Step 6: DeDup / MarkDuplicates
*/ 

process dedup{
    tag "${bam.baseName}"
    publishDir "${params.outdir}/5-DeDup"

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    file bam from ch_bam_filtered_dedup

    output:
    file "*.hist" into ch_hist_for_preseq
    file "*.log" into ch_dedup_results_for_multiqc
    file "*_rmdup.bam" into ch_dedup_bam

    script:
    if(params.singleEnd) {
    """
    dedup -i $bam -m -o . -u 
    """  
    } else {
    """
    dedup -i $bam -o . -u 
    """  
    }
}

/*
Step 5.1: Preseq
*/

process preseq {
    tag "${input.baseName}"
    publishDir "${params.outdir}/08-Preseq", mode: 'copy'

    when:
    !params.skip_preseq

    input:
    file input from (params.skip_deduplication ? ch_mapped_reads_preseq : ch_hist_for_preseq )

    output:
    file "${input.baseName}.ccurve" into ch_preseq_results

    script:
    if(!params.skip_deduplication){
    """
    preseq c_curve -s ${params.preseq_step_size} -o ${input.baseName}.ccurve -H $input
    """

    } else {
    """
    preseq c_curve -s ${params.preseq_step_size} -o ${input.baseName}.ccurve -B $input
    """
    }
}

/*
Step 5.2: DMG Assessment
*/ 

process damageprofiler {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/09-DamageProfiler", mode: 'copy'

    when:
    !params.skip_damage_calculation

    input:
    file bam from ch_mapped_reads_damageprofiler
    file fasta from ch_fasta_for_damageprofiler

    output:
    file "*" into ch_damageprofiler_results

    script:
    """
    damageprofiler -i $bam -r $fasta -l ${params.damageprofiler_length} -t ${params.damageprofiler_threshold} -o . 
    """
}

/* 
Step 5.3: Qualimap
*/

process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/06-QualiMap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    file bam from ch_bam_filtered_qualimap
    file fasta from ch_fasta_for_qualimap

    output:
    file "*" into ch_qualimap_results

    script:
    snpcap = ''
    if(params.snpcapture) snpcap = "-gff ${params.bedfile}"
    """
    qualimap bamqc -bam $bam -nt ${task.cpus} -outdir . -outformat "HTML" ${snpcap}
    """
}



/*
 Step 6: MarkDuplicates
 */

process markDup{
    tag "${bam.baseName}"
    publishDir "${params.outdir}/5-DeDup"

    when:
    !params.skip_deduplication && params.dedupper != 'dedup'

    input:
    file bam from ch_bam_filtered_markdup

    output:
    file "*.metrics" into ch_markdup_results_for_multiqc
    file "*.markDup.bam" into ch_markdup_bam

    script:
    prefix = "${bam.baseName}"
    """
    picard MarkDuplicates INPUT=$bam OUTPUT=${prefix}.markDup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE=${prefix}.markdup.metrics" VALIDATION_STRINGENCY=SILENT
    """
}

//If no deduplication runs, the input is mixed directly from samtools filter, if it runs either markdup or dedup is used thus mixed from these two channels
ch_dedup_for_pmdtools = Channel.create()

if(!params.skip_deduplication){
    ch_dedup_for_pmdtools.mix(ch_markdup_bam,ch_dedup_bam).set {ch_for_pmdtools}
} else {
    ch_dedup_for_pmdtools.mix(ch_markdup_bam,ch_dedup_bam,ch_bam_filtered_pmdtools).set {ch_for_pmdtools}
}

if(!params.run_pmdtools){
    ch_dedup_for_pmdtools.close()
}

process pmdtools {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/04-Samtools", mode: 'copy'

    when: params.run_pmdtools

    input: 
    file bam from ch_for_pmdtools
    file fasta from ch_fasta_for_pmdtools

    output:
    file "*.bam" into ch_bam_after_pmdfiltering
    file "*.cpg.range*.txt"

    script:
    //Check which treatment for the libraries was used
    def treatment = params.udg ? (params.udg_type =='half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(params.snpcapture){
        snpcap = (params.pmdtools_reference_mask != '') ? "--refseq ${params.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    """
    #Run Filtering step 
    samtools fillmd -b $bam $fasta | pmdtools --threshold ${params.pmdtools_threshold} $treatment $snpcap --header | samtools view -@ ${task.cpus} -Sb - > "${bam.baseName}".pmd.bam
    #Run Calc Range step
    #samtools fillmd -b $bam $fasta | pmdtools --deamination --range ${params.pmdtools_range} $treatment $snpcap -n ${params.pmdtools_max_reads} > "${bam.baseName}".cpg.range."${params.pmdtools_range}".txt 
    """
}

//Close Channel for pmdtools as it 




/*
Processing missing:
- pmdtools

Genotyping tools:
- angsd
- gatk (if even suitable anymore?)
- snpAD
- sequenceTools

Downstream VCF tools:
- vcf2genome
- gencons
- READ/mcMLKin
- popGen output? PLINK? 
*/




/*
 * STEP 2 - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from ch_fastqc_results.collect().ifEmpty([])
    file ('software_versions/*') from software_versions_yaml.collect().ifEmpty([])
    file ('idxstats/*') from ch_idxstats_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from ch_preseq_results.collect().ifEmpty([])
    file ('damageprofiler/*') from ch_damageprofiler_results.collect().ifEmpty([])
    file ('qualimap/*') from ch_qualimap_results.collect().ifEmpty([])
    file ('markdup/*') from ch_markdup_results_for_multiqc.collect().ifEmpty([])
    file ('dedup/*') from ch_dedup_results_for_multiqc.collect().ifEmpty([])
    file ('fastp/*') from ch_fastp_for_multiqc.collect().ifEmpty([])

    file workflow_summary from create_workflow_summary(summary)

    output:
    file "*multiqc_report.html" into ch_multiqc_report
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




/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/eager] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/eager] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = workflow.manifest.version
    email_fields['runName'] = custom_runName ?: workflow.runName
    email_fields['success'] = workflow.success
    email_fields['dateComplete'] = workflow.complete
    email_fields['duration'] = workflow.duration
    email_fields['exitStatus'] = workflow.exitStatus
    email_fields['errorMessage'] = (workflow.errorMessage ?: 'None')
    email_fields['errorReport'] = (workflow.errorReport ?: 'None')
    email_fields['commandLine'] = workflow.commandLine
    email_fields['projectDir'] = workflow.projectDir
    email_fields['summary'] = summary
    email_fields['summary']['Date Started'] = workflow.start
    email_fields['summary']['Date Completed'] = workflow.complete
    email_fields['summary']['Pipeline script file path'] = workflow.scriptFile
    email_fields['summary']['Pipeline script hash ID'] = workflow.scriptId
    if(workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if(workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if(workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$baseDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$baseDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir" ]
    def sf = new File("$baseDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (params.email) {
        try {
          if( params.plaintext_email ){ throw GroovyException('Send plaintext e-mail, not HTML') }
          // Try to send HTML e-mail using sendmail
          [ 'sendmail', '-t' ].execute() << sendmail_html
          log.info "[nf-core/eager] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/eager] Sent summary e-mail to $params.email (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File( "${params.outdir}/Documentation/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    log.info "[nf-core/eager] Pipeline Complete"
}
