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
      -profile                      Hardware config to use (standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentatoin.
      --singleEnd                   Specifies that the input is single end reads (required if not pairedEnd)
      --pairedEnd                   Specifies that the input is paired end reads (required if not singleend)
      --fasta                       Path to Fasta reference (required if not iGenome reference)
      --genome                      Name of iGenomes reference (required if not fasta reference)

    Input Data Additional Options:
      --snpcapture                  Runs in SNPCapture mode (specify a BED file if you do this!)
      --udg                         Specify that your libraries are treated with UDG
      --udg_type                    Specify here if you have UDG half treated libraries, Set to 'Half' in that case

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --bwa_index                   Path to BWA index
      --bedfile                     Path to BED file for SNPCapture methods
      --seq_dict                    Path to sequence dictionary file
      --fasta_index                 Path to FastA index 
      --saveReference               Saves reference genome indices for later reusage

    Skipping                        Skip any of the mentioned steps
      --skip_preseq
      --skip_damage_calculation
      --skip_qualimap
      --skip_deduplication
    
    Complexity Filtering 
      --complexity_filtering            Run complexity filtering on FastQ files
      --complexity_filter_poly_g_min    Specify poly-g min filter (default: 10) for filtering
    
    Clipping / Merging
      --clip_forward_adaptor        Specify adapter sequence to be clipped off (forward)
      --clip_reverse_adaptor        Specify adapter sequence to be clipped off (reverse)
      --clip_readlength             Specify read minimum length to be kept for downstream analysis
      --clip_min_read_quality       Specify minimum base quality for not trimming off bases
      --min_adap_overlap            Specify minimum adapter overlap
    
    BWA Mapping
      --bwaalnn                     Specify the -n parameter for BWA aln
      --bwaalnk                     Specify the -k parameter for BWA aln
      --bwaalnl                     Specify the -l parameter for BWA aln
    
    CircularMapper
      --circularmapper              Turn on CircularMapper (CM)
      --circularextension           Specify the number of bases to extend
      --circulartarget              Specify the target chromosome for CM
      --circularfilter              Specify to filter off-target reads
    
    BWA Mem Mapping
      --bwamem                      Turn on BWA Mem instead of CM/BWA aln for mapping
    
    BAM Filtering
      --bam_keep_mapped_only            Only consider mapped reads for downstream analysis. Unmapped reads are extracted to separate output.
      --bam_filter_reads                Keep all reads in BAM file for downstream analysis
      --bam_mapping_quality_threshold   Minimum mapping quality for reads filter
    
    DeDuplication
      --dedupper                    Deduplication method to use
      --dedup_all_merged            Treat all reads as merged reads
    
    Library Complexity Estimation
      --preseq_step_size            Specify the step size of Preseq
    
    (aDNA) Damage Analysis
      --damageprofiler_length       Specify length filter for DamageProfiler
      --damageprofiler_threshold    Specify number of bases to consider for damageProfiler
      --run_pmdtools                Turn on PMDtools
      --pmdtools_range              Specify range of bases for PMDTools
      --pmdtools_threshold          Specify PMDScore threshold for PMDTools
      --pmdtools_reference_mask     Specify a reference mask for PMDTools
      --pmdtools_max_reads          Specify the max. number of reads to consider for metrics generation
    
    BAM Trimming
      --trim_bam                    Turn on BAM trimming for UDG(+ or 1/2) protocols
      --bamutils_clip_left / --bamutils_clip_right  Specify the number of bases to clip off reads
      --bamutils_softclip           Use softclip instead of hard masking


    For a full list and more information of available parameters, consider the documentation.


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
params.pairedEnd = false
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

//Mapper to use, by default BWA aln will be used
params.circularmapper = false
params.circularextension = 500
params.circulartarget = 'MT'
params.circularfilter = false

//BWAMem Specific Settings 
params.bwamem = false

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
params.dedup_all_merged = false

//Preseq settings
params.preseq_step_size = 1000

//PMDTools settings
params.run_pmdtools = false
params.pmdtools_range = 10
params.pmdtools_threshold = 3
params.pmdtools_reference_mask = ''
params.pmdtools_max_reads = 10000

//bamUtils trimbam settings
params.trim_bam = false 
params.bamutils_clip_left = 1 
params.bamutils_clip_right = 1 
params.bamutils_softclip = false 




multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

// Validate inputs
Channel.fromPath("${params.fasta}")
    .ifEmpty { exit 1, "No genome specified! Please specify one with --fasta or --bwa_index"}
    .into {ch_fasta_for_bwa_indexing;ch_fasta_for_faidx_indexing;ch_fasta_for_dict_indexing; ch_fasta_for_bwa_mapping; ch_fasta_for_damageprofiler; ch_fasta_for_qualimap; ch_fasta_for_pmdtools; ch_fasta_for_circularmapper; ch_fasta_for_circularmapper_index;ch_fasta_for_bwamem_mapping}

//Validate that either pairedEnd or singleEnd has been specified by the user!
if( params.singleEnd || params.pairedEnd ){
} else {
    exit 1, "Please specify either --singleEnd or --pairedEnd to execute the pipeline!"
}


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
log.info " nf-core/eager v${workflow.manifest.version}"
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
process makeBWAIndex {
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/bwa_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: !params.bwa_index && params.fasta && params.aligner == 'bwa'

    input:
    file fasta from ch_fasta_for_bwa_indexing
    file wherearemyfiles

    output:
    file "*.{amb,ann,bwt,pac,sa,fasta,fa}" into (ch_bwa_index,ch_bwa_index_bwamem)
    file "where_are_my_files.txt"

    script:
    """
    bwa index $fasta
    """
}


/*
 * PREPROCESSING - Index Fasta file if not specified on CLI 
 */
process makeFastaIndex {
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/fasta_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }
    when: !params.fasta_index && params.fasta && params.aligner == 'bwa'

    input:
    file fasta from ch_fasta_for_faidx_indexing
    file wherearemyfiles

    output:
    file "${fasta}.fai" into ch_fasta_faidx_index
    file "${fasta}"
    file "where_are_my_files.txt"

    script:
    """
    samtools faidx $fasta
    """
}


/*
 * PREPROCESSING - Create Sequence Dictionary for FastA if not specified on CLI 
 */

process makeSeqDict {
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/seq_dict", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: !params.seq_dict && params.fasta

    input:
    file fasta from ch_fasta_for_dict_indexing
    file wherearemyfiles

    output:
    file "*.dict" into ch_seq_dict
    file "where_are_my_files.txt"

    script:
    """
    picard CreateSequenceDictionary R=$fasta O="${fasta.baseName}.dict"
    """
}



/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/FastQC", mode: 'copy',
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
    publishDir "${params.outdir}/FastP", mode: 'copy'

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
    publishDir "${params.outdir}/read_merging", mode: 'copy'

    input:
    set val(name), file(reads) from ( params.complexity_filter ? ch_clipped_reads_complexity_filtered : ch_read_files_clip )

    output:
    file "*.combined*.gz" into (ch_clipped_reads, ch_clipped_reads_for_fastqc,ch_clipped_reads_circularmapper,ch_clipped_reads_bwamem)
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
    publishDir "${params.outdir}/FastQC/after_clipping", mode: 'copy',
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

/*
Step 3: Mapping with BWA, SAM to BAM, Sort BAM
*/

process bwa {
    tag "$prefix"
    publishDir "${params.outdir}/mapping/bwa", mode: 'copy'

    when: !params.circularmapper && !params.bwamem

    input:
    file(reads) from ch_clipped_reads
    file "*" from ch_bwa_index
    file fasta from ch_fasta_for_bwa_mapping

    output:
    file "*.sorted.bam" into ch_mapped_reads_idxstats,ch_mapped_reads_filter,ch_mapped_reads_preseq, ch_mapped_reads_damageprofiler
    file "*.bai" 
    

    script:
    prefix = reads[0].toString() - ~/(_R1)?(\.combined\.)?(prefixed)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """ 
    bwa aln -t ${task.cpus} $fasta $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f "${reads.baseName}.sai"
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta "${reads.baseName}".sai $reads | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    samtools index -@ ${task.cpus} "${prefix}".sorted.bam
    """
}

process circulargenerator{
    tag "$prefix"
    publishDir "${params.outdir}/reference_genome/circularmapper_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: params.circularmapper

    input:
    file fasta from ch_fasta_for_circularmapper_index

    output:
    file "*.fasta*" into ch_circularmapper_indices

    script:
    """
    circulargenerator -e ${params.circularextension} -i $fasta -s ${params.circulartarget}
    bwa index "${fasta.baseName}_${params.circularextension}.fasta"
    """

}


process circularmapper{
    tag "$prefix"
    publishDir "${params.outdir}/mapping/circularmapper", mode: 'copy'

    when: params.circularmapper

    input:
    file reads from ch_clipped_reads_circularmapper
    file fasta from ch_fasta_for_circularmapper
    file "*" from ch_circularmapper_indices

    output:
    file "*.sorted.bam" into ch_mapped_reads_idxstats_cm,ch_mapped_reads_filter_cm,ch_mapped_reads_preseq_cm, ch_mapped_reads_damageprofiler_cm
    file "*.bai" 
    
    script:
    filter = "${params.circularfilter}" ? '' : '-f true -x false'
    prefix = reads[0].toString() - ~/(_R1)?(\.combined\.)?(prefixed)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """ 
    bwa aln -t ${task.cpus} "${fasta.baseName}_${params.circularextension}.fasta" $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f "${reads.baseName}.sai"
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" "${fasta.baseName}_${params.circularextension}.fasta" "${reads.baseName}".sai $reads > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > "${prefix}".sorted.bam
    samtools index -@ ${task.cpus} "${prefix}".sorted.bam
    """
}

process bwamem {
    tag "$prefix"
    publishDir "${params.outdir}/mapping/bwamem", mode: 'copy'

    when: params.bwamem && !params.circularmapper

    input:
    file(reads) from ch_clipped_reads_bwamem
    file "*" from ch_bwa_index_bwamem
    file fasta from ch_fasta_for_bwamem_mapping

    output:
    file "*.sorted.bam" into ch_bwamem_mapped_reads_idxstats,ch_bwamem_mapped_reads_filter,ch_bwamem_mapped_reads_preseq, ch_bwamem_mapped_reads_damageprofiler
    file "*.bai" 
    

    script:
    prefix = reads[0].toString() - ~/(_R1)?(\.combined\.)?(prefixed)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    """
    bwa mem -t ${task.cpus} $fasta $reads -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    samtools index -@ ${task.cpus} "${prefix}".sorted.bam
    """
}

/*
* Step 4 - IDXStats
*/

process samtools_idxstats {
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    input:
    file(bam) from ch_mapped_reads_idxstats.mix(ch_mapped_reads_idxstats_cm,ch_bwamem_mapped_reads_idxstats)

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
    publishDir "${params.outdir}/samtools/filter", mode: 'copy',
    saveAs: {filename ->
            if (filename.indexOf(".fq.gz") > 0) "unmapped/$filename"
            else if (filename.indexOf(".unmapped.bam") > 0) "unmapped/$filename"
            else if (filename.indexOf(".filtered.bam")) filename
            else null
    }

    input: 
    file bam from ch_mapped_reads_filter.mix(ch_mapped_reads_filter_cm,ch_bwamem_mapped_reads_filter)

    output:
    file "*filtered.bam" into ch_bam_filtered_qualimap, ch_bam_filtered_dedup, ch_bam_filtered_markdup, ch_bam_filtered_pmdtools, ch_bam_filtered_angsd, ch_bam_filtered_gatk
    file "*.fq.gz" optional true
    file "*.unmapped.bam" optional true
    file "*.bai"

    when: "${params.bam_filter_reads}"

    script:
    prefix="$bam" - ~/(\.bam)?/

    if("${params.bam_keep_mapped_only}"){
    """
    samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
    samtools fastq -tn "${prefix}.unmapped.bam" | gzip > "${prefix}.unmapped.fq.gz"
    samtools index -@ ${task.cpus} ${prefix}.filtered.bam
    """
    } else {
    """
    samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
    samtools index -@ ${task.cpus} ${prefix}.filtered.bam
    """
    }  
}


/*
Step 6: DeDup / MarkDuplicates
*/ 

process dedup{
    tag "${bam.baseName}"
    publishDir "${params.outdir}/deduplication/dedup"

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    file bam from ch_bam_filtered_dedup

    output:
    file "*.hist" into ch_hist_for_preseq
    file "*.log" into ch_dedup_results_for_multiqc
    file "${prefix}.sorted.bam" into ch_dedup_bam
    file "*.bai"

    script:
    prefix="${bam.baseName}"
    treat_merged="${params.dedup_all_merged}" ? '-m' : ''

    if(params.singleEnd) {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix".sorted.bam
    samtools index -@ ${task.cpus} "$prefix".sorted.bam
    """  
    } else {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix".sorted.bam
    samtools index -@ ${task.cpus} "$prefix".sorted.bam
    """  
    }
}

/*
Step 5.1: Preseq
*/

process preseq {
    tag "${input.baseName}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_preseq

    input:
    file input from (params.skip_deduplication ? ch_mapped_reads_preseq.mix(ch_mapped_reads_preseq_cm,ch_bwamem_mapped_reads_preseq) : ch_hist_for_preseq )

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
    publishDir "${params.outdir}/damageprofiler", mode: 'copy'

    when:
    !params.skip_damage_calculation

    input:
    file bam from ch_mapped_reads_damageprofiler.mix(ch_mapped_reads_damageprofiler_cm,ch_bwamem_mapped_reads_damageprofiler)
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
    publishDir "${params.outdir}/qualimap", mode: 'copy'

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
    publishDir "${params.outdir}/deduplication/markdup"

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

//Bamutils TrimBam Channel
ch_for_bamutils = Channel.create()

if(!params.skip_deduplication){
    ch_dedup_for_pmdtools.mix(ch_markdup_bam,ch_dedup_bam).into {ch_for_pmdtools;ch_for_bamutils}
} else {
    ch_dedup_for_pmdtools.mix(ch_markdup_bam,ch_dedup_bam,ch_bam_filtered_pmdtools).into {ch_for_pmdtools;ch_for_bamutils}
}

if(!params.run_pmdtools){
    ch_dedup_for_pmdtools.close()
}

process pmdtools {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/pmdtools", mode: 'copy'

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

/*
* Optional BAM Trimming step using bamUtils 
* Can be used for UDGhalf protocols to clip off -n bases of each read
*/

process bam_trim {
    tag "${prefix}" 
    publishDir "${params.outdir}/trimmed_bam", mode: 'copy'
 
    when: params.trim_bam

    input:
    file bam from ch_for_bamutils  

    output: 
    file "*.trimmed.bam" into ch_trimmed_bam_for_genotyping
    file "*.bai"

    script:
    prefix="${bam.baseName}"
    softclip = "${params.bamutils_softclip}" ? '-c' : '' 
    """
    bam trimBam $bam tmp.bam -L ${params.bamutils_clip_left} -R ${params.bamutils_clip_right} ${softclip}
    samtools sort -@ ${task.cpus} tmp.bam -o ${prefix}.trimmed.bam 
    samtools index ${prefix}.trimmed.bam
    """
}




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
    file ('adapter_removal/*') from ch_adapterremoval_logs.collect().ifEmpty([])
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
