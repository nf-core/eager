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
 Maxime Borry <borry@shh.mpg.de.de> - https://github.com/maxibor
========================================================================================
*/

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    =========================================
    eager v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/eager --reads '*_R{1,2}.fastq.gz' -profile docker

    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      -profile                      Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation.
      --singleEnd                   Specifies that the input is single end reads (required if not pairedEnd)
      --pairedEnd                   Specifies that the input is paired end reads (required if not singleEnd)
      --bam                         Specifies that the input is in BAM format
      --fasta                       Path and name of FASTA reference file (required if not iGenome reference). File suffixes can be: '.fa', '.fn', '.fna', '.fasta'
      --genome                      Name of iGenomes reference (required if not fasta reference)

    Input Data Additional Options:
      --snpcapture                  Runs in SNPCapture mode (specify a BED file if you do this!)

    References                      If not specified in the configuration file, or you wish to overwrite any of the references.
      --bwa_index                   Path and name of a a bwa indexed FASTA reference file with index suffixes (i.e. everything before the endings '.amb' '.ann' '.bwt'. Most likely the same value supplied with the --fasta option)
      --bedfile                     Path to BED file for SNPCapture methods
      --seq_dict                    Path to picard sequence dictionary file (typically ending in '.dict')
      --fasta_index                 Path to samtools FASTA index (typically ending in '.fai')
      --saveReference               Saves reference genome indices for later reusage

    Skipping                        Skip any of the mentioned steps
      --skip_fastqc                 Skips both pre- and post-Adapter Removal FastQC steps.
      --skip_adapterremoval         
      --skip_preseq
      --skip_damage_calculation
      --skip_qualimap
      --skip_deduplication
    
    Complexity Filtering 
      --complexity_filter_poly_g        Run poly-G removal on FASTQ files
      --complexity_filter_poly_g_min    Specify length of poly-g min for clipping to be performed (default: 10)
    
    Clipping / Merging
      --clip_forward_adaptor        Specify adapter sequence to be clipped off (forward)
      --clip_reverse_adaptor        Specify adapter sequence to be clipped off (reverse)
      --clip_readlength             Specify read minimum length to be kept for downstream analysis
      --clip_min_read_quality       Specify minimum base quality for not trimming off bases
      --min_adap_overlap            Specify minimum adapter overlap
      --skip_collapse               Skip merging forward and reverse reads together. (Only for PE samples)
      --skip_trim                   Skip adaptor and quality trimming
    
    BWA Mapping
      --bwaalnn                     Specify the -n parameter for BWA aln.
      --bwaalnk                     Specify the -k parameter for BWA aln
      --bwaalnl                     Specify the -l parameter for BWA aln
 
    Stripping
      --strip_input_fastq           Create pre-Adapter Removal FASTQ files without reads that mapped to reference (e.g. for public upload of privacy sensitive non-host data)
      --strip_mode                  Stripping mode. Remove mapped reads completely from FASTQ (strip) or just mask mapped reads sequence by N (replace)
    
    CircularMapper
      --circularmapper              Turn on CircularMapper (CM)
      --circularextension           Specify the number of bases to extend reference by
      --circulartarget              Specify the target chromosome for CM
      --circularfilter              Specify to filter off-target reads
    
    BWA Mem Mapping
      --bwamem                      Turn on BWA Mem instead of BWA aln for mapping
    
    BAM Filtering
      --run_bam_filtering           Turn on mapping quality or unmapped reads filtering using samtools filter.
      --bam_mapping_quality_threshold   Minimum mapping quality for reads filter, default 0.
      --bam_discard_unmapped        Discards unmapped reads in either FASTQ or BAM format, depending on choice (see --bam_unmapped_type).
      --bam_unmapped_type           Defines whether to discard all unmapped reads, keep only bam and/or keep only fastq format (options: discard, bam, fastq, both).
    
    DeDuplication
      --dedupper                    Deduplication method to use (options: dedup, markduplicates). Default: dedup
      --dedup_all_merged            Treat all reads as merged reads
    
    Library Complexity Estimation
      --preseq_step_size            Specify the step size of Preseq
    
    (aDNA) Damage Analysis
      --damageprofiler_length       Specify length filter for DamageProfiler
      --damageprofiler_threshold    Specify number of bases to consider for damageProfiler
      --run_pmdtools                Turn on PMDtools
      --udg_type                    Specify here if you have UDG half treated libraries, Set to 'half' in that case, or 'full' for UDG+. If not set, libraries are set to UDG-.
      --pmdtools_range              Specify range of bases for PMDTools
      --pmdtools_threshold          Specify PMDScore threshold for PMDTools
      --pmdtools_reference_mask     Specify a reference mask for PMDTools
      --pmdtools_max_reads          Specify the max. number of reads to consider for metrics generation
    
    BAM Trimming
      --trim_bam                    Turn on BAM trimming for UDG(+ or 1/2) protocols
      --bamutils_clip_left / --bamutils_clip_right  Specify the number of bases to clip off reads
      --bamutils_softclip           Use softclip instead of hard masking

    Other options:     
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --plaintext_email             Receive plain text emails rather than HTML
      --maxMultiqcEmailFileSize     Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --max_memory                  Memory limit for each step of pipeline. Should be in form e.g. --max_memory '8.GB'
      --max_time                    Time limit for each step of the pipeline. Should be in form e.g. --max_memory '2.h'
      --max_cpus                    Maximum number of CPUs to use for each step of the pipeline. Should be in form e.g. --max_cpus 1
      
    For a full list and more information of available parameters, consider the documentation (https://github.com/nf-core/eager/).

      
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
params.fasta = ''
params.seq_dict = false
params.fasta_index = false
params.saveReference = false
params.pmd_udg_type = 'half'

params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

// Skipping parts of the pipeline for impatient users
params.skip_fastqc = false 
params.skip_adapterremoval = false
params.skip_preseq = false
params.skip_damage_calculation = false
params.skip_qualimap = false
params.skip_deduplication = false

//Complexity filtering reads
params.complexity_filter_poly_g = false
params.complexity_filter_poly_g_min = 10

//Read clipping and merging parameters
params.clip_forward_adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
params.clip_reverse_adaptor = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"
params.clip_readlength = 30
params.clip_min_read_quality = 20
params.min_adap_overlap = 1
params.skip_collapse = false
params.skip_trim = false

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
params.run_bam_filtering = false
params.bam_discard_unmapped = false
params.bam_unmapped_type = ''

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

//unmap

params.strip_input_fastq = false
params.strip_mode = 'strip'


multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
where_are_my_files = file("$baseDir/assets/where_are_my_files.txt")
// Validate inputs
if ( params.fasta.isEmpty () ){
    exit 1, "Please specify --fasta with the path to your reference"
} else if("${params.fasta}".endsWith(".gz")){
    //Put the zip into a channel, then unzip it and forward to downstream processes. DONT unzip in all steps, this is inefficient as NXF links the files anyways from work to work dir
    zipped_fasta = file("${params.fasta}")

    rm_gz = params.fasta - '.gz'
    lastPath = rm_gz.lastIndexOf(File.separator)
    bwa_base = rm_gz.substring(lastPath+1)

    process unzip_reference{
        tag "${zipped_fasta}"

        input:
        file zipped_fasta

        output:
        file "*.{fa,fn,fna,fasta}" into fasta_for_indexing

        script:
        rm_zip = zipped_fasta - '.gz'
        """
        pigz -f -d -p ${task.cpus} $zipped_fasta
        """
        }
       
    } else {
    fasta_for_indexing = file("${params.fasta}")
    lastPath = params.fasta.lastIndexOf(File.separator)
    bwa_base = params.fasta.substring(lastPath+1)
}

//Index files provided? Then check whether they are correct and complete
if (params.aligner != 'bwa' && !params.circularmapper && !params.bwamem){
    exit 1, "Invalid aligner option. Default is bwa, but specify --circularmapper or --bwamem to use these."
}
if( params.bwa_index && (params.aligner == 'bwa' | params.bwamem)){
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir =  params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)

    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .ifEmpty { exit 1, "BWA index directory not found: ${bwa_dir}" }
        .into {bwa_index; bwa_index_bwamem}
}



//Validate that either pairedEnd or singleEnd has been specified by the user!
if( params.singleEnd || params.pairedEnd || params.bam){
} else {
    exit 1, "Please specify either --singleEnd, --pairedEnd to execute the pipeline on FastQ files and --bam for previously processed BAM files!"
}

//Validate that skip_collapse is only set to True for pairedEnd reads!
if (params.skip_collapse  && params.singleEnd){
    exit 1, "--skip_collapse can only be set for pairedEnd samples!"
}

//Strip mode sanity checking
if (params.strip_input_fastq){
    if (!(['strip','replace'].contains(params.strip_mode))) {
        exit 1, "--strip_mode can only be set to strip or replace"
    }
}



// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

if( workflow.profile == 'awsbatch') {
  // AWSBatch sanity checking
  if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
  // Check outdir paths to be S3 buckets if running on AWSBatch
  // related: https://github.com/nextflow-io/nextflow/issues/813
  if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
  // Prevent trace files to be stored on S3 since S3 does not support rolling files.
  if (workflow.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}
/*damageprofiler, ch_bam_filtered_dedup, ch_bam_filtered_markdup, ch_bam_filtered_pmdtools, ch_bam_filtered_angsd, ch_bam_filtered_gatk
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */

if( params.readPaths ){
    if( params.singleEnd && !params.bam) {
        Channel
            .from( params.readPaths )
            .map { row -> [ row[0], [ file( row[1][0] ) ] ] }
            .ifEmpty { exit 1, "params.readPaths or params.bams was empty - no input files supplied!" }
            .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filter_poly_g ; ch_read_unmap}
            ch_bam_to_fastq_convert = Channel.empty()
    } else if (!params.bam){
        Channel
            .from( params.readPaths )
            .map { row -> [ row[0], [ file( row[1][0] ), file( row[1][1] ) ] ] }
            .ifEmpty { exit 1, "params.readPaths or params.bams was empty - no input files supplied!" }
            .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filter_poly_g; ch_read_unmap }
            ch_bam_to_fastq_convert = Channel.empty()
    } else {
        Channel
            .from( params.readPaths )
            .map { row -> [ file( row )  ] }
            .ifEmpty { exit 1, "params.readPaths or params.bams was empty - no input files supplied!" }
            .dump()
            .set { ch_bam_to_fastq_convert }

            //Set up clean channels
            ch_read_files_fastqc = Channel.empty()
            ch_read_files_complexity_filter_poly_g = Channel.empty()
            ch_read_files_clip = Channel.empty()
            ch_read_unmap = Channel.empty()
    }
} else if (!params.bam){
     Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs" +
            "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { ch_read_files_clip; ch_read_files_fastqc; ch_read_files_complexity_filter_poly_g; ch_read_unmap }
        ch_bam_to_fastq_convert = Channel.empty()
} else {
     Channel
        .fromPath( params.reads )
        .map { row -> [  file( row )  ] }
        .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nNB: Path needs" +
            "to be enclosed in quotes!\n" }
        .dump() //For debugging purposes
        .set { ch_bam_to_fastq_convert }

        //Set up clean channels
        ch_read_files_fastqc = Channel.empty()
        ch_read_files_complexity_filter_poly_g = Channel.empty()
        ch_read_files_clip = Channel.empty()
        ch_read_unmap = Channel.empty()

}

// Header log info
log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/eager'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['BAM Index Type'] = (params.large_ref == "") ? 'BAI' : 'CSI'
if(params.bwa_index) summary['BWA Index'] = params.bwa_index
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Skip Collapsing'] = params.skip_collapse ? 'Yes' : 'No'
summary['Skip Trimming']  = params.skip_trim  ? 'Yes' : 'No' 
summary['Output stripped fastq'] = params.strip_input_fastq ? 'Yes' : 'No'
if (params.strip_input_fastq){
    summary['Strip mode'] = params.strip_mode
}
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
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
summary['Config Profile'] = workflow.profile
if(params.config_profile_description) summary['Config Description'] = params.config_profile_description
if(params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if(params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if(params.email) {
  summary['E-mail Address']  = params.email
  summary['MultiQC maxsize'] = params.maxMultiqcEmailFileSize
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "\033[2m----------------------------------------------------\033[0m"

// Check the hostnames against configured profiles
checkHostname()

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
* PREPROCESSING - Create BWA indices if they are not present
*/ 

if(!params.bwa_index && !params.fasta.isEmpty() && (params.aligner == 'bwa' || params.bwamem)){
process makeBWAIndex {
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/bwa_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    file fasta from fasta_for_indexing
    file where_are_my_files

    output:
    file "BWAIndex" into (bwa_index, bwa_index_bwamem)
    file "where_are_my_files.txt"

    script:
    """
    bwa index $fasta
    mkdir BWAIndex && mv ${fasta}* BWAIndex
    """
    }
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
    when: !params.fasta_index && !params.fasta.isEmpty() && params.aligner == 'bwa'

    input:
    file fasta from fasta_for_indexing
    file where_are_my_files

    output:
    file "*.fai" into ch_fasta_faidx_index
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
    
    when: !params.seq_dict && !params.fasta.isEmpty()

    input:
    file fasta from fasta_for_indexing
    file where_are_my_files

    output:
    file "*.dict" into ch_seq_dict
    file "where_are_my_files.txt"

    script:
    """
    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M CreateSequenceDictionary R=$fasta O="${fasta.baseName}.dict"
    """
}

/*
* PREPROCESSING - Convert BAM to FastQ if BAM input is specified instead of FastQ file(s)
*
*/ 

process convertBam {
    tag "$bam"
    
    when: params.bam

    input: 
    file bam from ch_bam_to_fastq_convert

    output:
    set val("${base}"), file("*.fastq.gz") into (ch_read_files_converted_fastqc, ch_read_files_converted_fastp, ch_read_files_converted_mapping_bwa, ch_read_files_converted_mapping_cm, ch_read_files_converted_mapping_bwamem,
    ch_read_unmap_convertBam)

    script:
    base = "${bam.baseName}"
    """
    samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.fastq.gz
    """ 
}



/*
 * STEP 1a - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/FastQC", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: !params.skip_fastqc

    input:
    set val(name), file(reads) from ch_read_files_fastqc.mix(ch_read_files_converted_fastqc)

    output:
    file "*_fastqc.{zip,html}" into ch_fastqc_results

    script:
    """
    fastqc -q $reads
    rename 's/_fastqc\\.zip\$/_raw_fastqc.zip/' *_fastqc.zip
    rename 's/_fastqc\\.html\$/_raw_fastqc.html/' *_fastqc.html
    """
}


/* STEP 1b - FastP
* Optional poly-G complexity filtering step before read merging/adapter clipping etc
* Note: Clipping, Merging, Quality Trimning are turned off here - we leave this to adapter removal itself!
*/

process fastp {
    tag "$name"
    publishDir "${params.outdir}/FastP", mode: 'copy'

    when: params.complexity_filter_poly_g

    input:
    set val(name), file(reads) from ch_read_files_complexity_filter_poly_g.mix(ch_read_files_converted_fastp)

    output:
    set val(name), file("*pG.fq.gz") into ch_clipped_reads_complexity_filtered_poly_g
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
//Initialize empty channel if we skip adapterremoval entirely
if(params.skip_adapterremoval) {
    //No logs if no AR is run
    ch_adapterremoval_logs = Channel.empty()
    //Either coming from complexity filtering, or directly use reads normally directed to clipping first and push them through to the other channels downstream! 
    ch_clipped_reads_complexity_filtered_poly_g.mix(ch_read_files_clip).into { ch_clipped_reads;ch_clipped_reads_for_fastqc;ch_clipped_reads_circularmapper;ch_clipped_reads_bwamem }
} else {
process adapter_removal {
    tag "$name"
    publishDir "${params.outdir}/read_merging", mode: 'copy'

    when: !params.bam && !params.skip_adapterremoval

    input:
    set val(name), file(reads) from ( params.complexity_filter_poly_g ? ch_clipped_reads_complexity_filtered_poly_g : ch_read_files_clip )

    output:
    set val(base), file("output/*.gz") into (ch_clipped_reads,ch_clipped_reads_for_fastqc,ch_clipped_reads_circularmapper,ch_clipped_reads_bwamem)
    file("*.settings") into ch_adapterremoval_logs

    script:
    base = reads[0].baseName
    //This checks whether we skip trimming and defines a variable respectively
    trim_me = params.skip_trim ? '' : "--trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}"
    collapse_me = params.skip_collapse ? '' : '--collapse'
    
    //PE mode, dependent on trim_me and collapse_me the respective procedure is run or not :-) 
    if (!params.singleEnd && !params.skip_collapse && !params.skip_trim){
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} ${trim_me} --gzip --threads ${task.cpus} ${collapse_me}
    #Combine files
    zcat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > output/${base}.combined.fq.gz
    """
    //PE, don't collapse, but trim reads
    } else if (!params.singleEnd && params.skip_collapse && !params.skip_trim) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} --gzip --threads ${task.cpus} ${trim_me} ${collapse_me}
    mv ${base}.pair*.truncated.gz output/
    """
    //PE, collapse, but don't trim reads
    } else if (!params.singleEnd && !params.skip_collapse && params.skip_trim) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} --gzip --threads ${task.cpus} --basename ${base} ${collapse_me} ${trim_me}
    
    mv ${base}.pair*.truncated.gz output/
    """
    } else {
    //SE, collapse not possible, trim reads
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --basename ${base} --gzip --threads ${task.cpus} ${trim_me}
    
    mv *.truncated.gz output/
    """
    }
}
}



/*
* STEP 2b - FastQC after clipping/merging (if applied!)
*/
process fastqc_after_clipping {
    tag "${name}"
    publishDir "${params.outdir}/FastQC/after_clipping", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: !params.bam && !params.skip_adapterremoval && !params.skip_fastqc

    input:
    set val(name), file(reads) from ch_clipped_reads_for_fastqc

    output:
    file "*_fastqc.{zip,html}" optional true into ch_fastqc_after_clipping

    script:
    """
    fastqc -q $reads
    """
}

/*
Step 3a  - Mapping with BWA, SAM to BAM, Sort BAM
*/

process bwa {
    tag "${name}"
    publishDir "${params.outdir}/mapping/bwa", mode: 'copy'

    when: !params.circularmapper && !params.bwamem

    input:
    set val(name), file(reads) from ch_clipped_reads.mix(ch_read_files_converted_mapping_bwa)
    file index from bwa_index.collect()


    output:
    file "*.sorted.bam" into ch_mapped_reads_flagstat,ch_mapped_reads_filter,ch_mapped_reads_preseq,ch_bwa_mapped_reads_strip
    file "*.{bai,csi}"
    

    script:
    size = "${params.large_ref}" ? '-c' : ''
    fasta = "${index}/${bwa_base}"

    //PE data without merging, PE data without any AR applied
    if (!params.singleEnd && (params.skip_collapse || params.skip_adapterremoval || params.skip_trim)){
    prefix = "${reads[0].baseName}"
    """
    bwa aln -t ${task.cpus} $fasta ${reads[0]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r1.sai
    bwa aln -t ${task.cpus} $fasta ${reads[1]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta ${prefix}.r1.sai ${prefix}.r2.sai ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -O bam - > ${prefix}.sorted.bam
    samtools index "${size}" "${prefix}".sorted.bam
    """
    } else {
    //PE collapsed, or SE data 
    prefix = "${reads.baseName}"
    """
    bwa aln -t ${task.cpus} $fasta $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta ${prefix}.sai $reads | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    samtools index "${size}" "${prefix}".sorted.bam
    """
    }
    
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
    file fasta from fasta_for_indexing

    output:
    file "${prefix}.{amb,ann,bwt,sa,pac}" into ch_circularmapper_indices

    script:
    prefix = "${fasta.baseName}_${params.circularextension}.fasta"
    """
    circulargenerator -e ${params.circularextension} -i $fasta -s ${params.circulartarget}
    bwa index $prefix
    """

}


process circularmapper{
    tag "$prefix"
    publishDir "${params.outdir}/mapping/circularmapper", mode: 'copy'

    when: params.circularmapper

    input:
    set val(name), file(reads) from ch_clipped_reads_circularmapper.mix(ch_read_files_converted_mapping_cm)
    file index from ch_circularmapper_indices.collect()
    file fasta from fasta_for_indexing

    output:
    file "*.sorted.bam" into ch_mapped_reads_flagstat_cm,ch_mapped_reads_filter_cm,ch_mapped_reads_preseq_cm,ch_circular_mapped_reads_strip
    file "*.{bai,csi}"
    
    script:
    filter = "${params.circularfilter}" ? '' : '-f true -x false'
    elongated_root = "${fasta.baseName}_${params.circularextension}.fasta"

    
    size = "${params.large_ref}" ? '-c' : ''

    if (!params.singleEnd && params.skip_collapse ){
    prefix = reads[0].toString().tokenize('.')[0]
    """ 
    bwa aln -t ${task.cpus} $elongated_root ${reads[0]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r1.sai
    bwa aln -t ${task.cpus} $elongated_root ${reads[1]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $elongated_root ${prefix}.r1.sai ${prefix}.r2.sai ${reads[0]} ${reads[1]} > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > ${prefix}.sorted.bam
    samtools index "${size}" ${prefix}.sorted.bam
    """
    } else {
    prefix = reads[0].toString().tokenize('.')[0]
    """ 
    bwa aln -t ${task.cpus} $elongated_root $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $elongated_root ${prefix}.sai $reads > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > "${prefix}".sorted.bam
    samtools index "${size}" "${prefix}".sorted.bam
    """
    }
    
}

process bwamem {
    tag "$prefix"
    publishDir "${params.outdir}/mapping/bwamem", mode: 'copy'

    when: params.bwamem && !params.circularmapper

    input:
    set val(name), file(reads) from ch_clipped_reads_bwamem.mix(ch_read_files_converted_mapping_bwamem)
    file index from bwa_index_bwamem.collect()

    output:
    file "*.sorted.bam" into ch_bwamem_mapped_reads_flagstat,ch_bwamem_mapped_reads_filter,ch_bwamem_mapped_reads_preseq,ch_bwamem_mapped_reads_strip
    file "*.{bai,csi}"
    

    script:
    fasta = "${index}/${bwa_base}"
    prefix = reads[0].toString() - ~/(_R1)?(\.combined\.)?(prefixed)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    size = "${params.large_ref}" ? '-c' : ''

    if (!params.singleEnd && params.skip_collapse){
    """
    bwa mem -t ${task.cpus} $fasta ${reads[0]} ${reads[1]} -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    samtools index "${size}" -@ ${task.cpus} "${prefix}".sorted.bam
    """
    } else {
    """
    bwa mem -t ${task.cpus} $fasta $reads -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".sorted.bam
    samtools index "${size}" -@ ${task.cpus} "${prefix}".sorted.bam
    """
    }
    
}

/*
* Step 3b - flagstat
*/
process samtools_flagstat {
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'
    input:
    file(bam) from ch_mapped_reads_flagstat.mix(ch_mapped_reads_flagstat_cm,ch_bwamem_mapped_reads_flagstat)
    output:
    file "*.stats" into ch_flagstat_for_multiqc
    script:
    prefix = "$bam" - ~/(\.bam)?$/
    """
    samtools flagstat $bam > ${prefix}.stats
    """
}
// If BAM filtering not run, directly send mapped reads to channels for DeDuping
if (!params.run_bam_filtering) {
ch_bam_filtered_dedup = Channel.empty()
ch_bam_filtered_markdup = Channel.empty()
ch_filtered_bam_for_downstream = Channel.empty()

  ch_mapped_reads_filter.mix(ch_mapped_reads_filter_cm,ch_bwamem_mapped_reads_filter)
    .into{ch_bam_filtered_dedup;ch_bam_filtered_markdup;ch_filtered_bam_for_downstream}

 ch_mapped_reads_filter.close()
 ch_mapped_reads_filter_cm.close()
 ch_bwamem_mapped_reads_filter.close()

}

/*
* Step 4a - Keep unmapped/remove unmapped reads
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
    
    when:
    params.run_bam_filtering
    
    input: 
    file bam from ch_mapped_reads_filter.mix(ch_mapped_reads_filter_cm,ch_bwamem_mapped_reads_filter)
    
    output:
    file "*filtered.bam" into ch_bam_filtered_flagstat,ch_bam_filtered_dedup,ch_bam_filtered_markdup,ch_filtered_bam_for_downstream 
    file "*.fastq.gz" optional true
    file "*.unmapped.bam" optional true
    file "*.{bai,csi}" into ch_bam_index_filtered_qualimap
    
    script:
    prefix="$bam" - ~/(\.bam)?/
    size = "${params.large_ref}" ? '-c' : ''
    
    if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "discard"){
        """
        samtools view -h -b $bam -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
        samtools index "${size}" ${prefix}.filtered.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "bam"){
        """
        samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
        samtools index "${size}" ${prefix}.filtered.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "fastq"){
        """
        samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
        samtools index "${size}" ${prefix}.filtered.bam
        samtools fastq -tn ${prefix}.unmapped.bam | pigz -p ${task.cpus} > ${prefix}.unmapped.fastq.gz
        rm ${prefix}.unmapped.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "both"){
        """
        samtools view -h $bam | tee >(samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam) >(samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
        samtools index "${size}" ${prefix}.filtered.bam
        samtools fastq -tn ${prefix}.unmapped.bam | pigz -p ${task.cpus} > ${prefix}.unmapped.fastq.gz
        """
    } else { //Only apply quality filtering, default
        """
        samtools view -h -b $bam -@ ${task.cpus} -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
        samtools index "${size}" ${prefix}.filtered.bam
        """
    }  
}
process strip_input_fastq {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/samtools/stripped_fastq", mode: 'copy'
    when: params.strip_input_fastq
    input: 
    set val(name), file(fq) from ch_read_unmap.mix(ch_read_unmap_convertBam)
    file bam from ch_bwa_mapped_reads_strip.mix(ch_circular_mapped_reads_strip, ch_bwamem_mapped_reads_strip)
    output:
    file "*.fq.gz" into unmapped_fq_ch
    script:
    if (params.singleEnd) {
        out_fwd = bam.baseName+'.stripped.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${fq[0]} -m ${params.strip_mode} -of $out_fwd -p ${task.cpus}
        """
    } else {
        out_fwd = bam.baseName+'.stripped.fwd.fq.gz'
        out_rev = bam.baseName+'.stripped.rev.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${fq[0]} -rev ${fq[1]} -m  ${params.strip_mode} -of $out_fwd -or $out_rev -p ${task.cpus}
        """ 
    }
    
}
/*
* Step 4b: Keep unmapped/remove unmapped reads flagstat
*/
process samtools_flagstat_after_filter {
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'
    
    when params.run_bam_filtering
    input:
    file(bam) from ch_bam_filtered_flagstat
    output:
    file "*.stats" into ch_bam_filtered_flagstat_for_multiqc
    script:
    prefix = "$bam" - ~/(\.bam)?$/
    """
    samtools flagstat $bam > ${prefix}.stats
    """
}


/*
Step 5a: DeDup / MarkDuplicates
*/ 

process dedup{
    tag "${bam.baseName}"
    publishDir "${params.outdir}/deduplication/dedup", mode: 'copy',
        saveAs: {filename -> "${prefix}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    file bam from ch_bam_filtered_dedup

    output:
    file "*.hist" into ch_hist_for_preseq
    file "*.log" into ch_dedup_results_for_multiqc
    file "${prefix}_rmdup.sorted.bam" into ch_dedup_bam_for_downstream
    file "*.{bai,csi}" into ch_dedup_bam_index_for_genotyping_ug

    script:
    prefix="${bam.baseName}"
    treat_merged="${params.dedup_all_merged}" ? '-m' : ''
    size = "${params.large_ref}" ? '-c' : ''
    
    if(params.singleEnd) {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix"_rmdup.sorted.bam
    samtools index "${size}" "$prefix"_rmdup.sorted.bam
    """  
    } else {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix"_rmdup.sorted.bam
    samtools index "${size}" "$prefix"_rmdup.sorted.bam
    """  
    }
}

/*
 Step 5b: MarkDuplicates
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
    file "*_rmdup.sorted.bam" into ch_markdup_bam,ch_markdup_bam_for_downstream
    file "*.{bai,csi}"


    script:
    prefix = "${bam.baseName}"
    size = "${params.large_ref}" ? '-c' : ''
    """
    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M MarkDuplicates INPUT=$bam OUTPUT=${prefix}_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="${prefix}_rmdup.metrics" VALIDATION_STRINGENCY=SILENT
    samtools sort -@ ${task.cpus} "$prefix"_rmdup -o "$prefix"_rmdup.bam.sorted.bam
    samtools index "${size}" "$prefix"_rmdup.sorted.bam
    """
}

/*
Step 6: Preseq
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
Switch for sending filtered or rmdup'd BAMs downstream
*/

//If no deduplication runs, the input is mixed directly from samtools filter, else selects the deduplicated BAMs 
ch_bams_for_damageprofiler = Channel.empty()
ch_bams_for_qualimap = Channel.empty()
ch_bams_for_bamutils = Channel.empty()
ch_bams_for_pmdtools = Channel.empty()


if (!params.skip_deduplication) {
  ch_filtered_bam_for_downstream
    .mix(ch_dedup_bam_for_downstream ,ch_markdup_bam_for_downstream)
    .filter { it =~/.*_rmdup.sorted.bam/ }
    .into{ ch_bams_for_damageprofiler; ch_bams_for_qualimap; ch_bams_for_bamutils; ch_bams_for_pmdtools; ch_bams_for_angsd; ch_bams_for_gatk } 
} else {
  ch_filtered_bam_for_downstream
    .mix(ch_dedup_bam_for_downstream ,ch_markdup_bam_for_downstream)
    .filter { it !=~/.*_rmdup.sorted.bam/ }
    .into{ ch_bams_for_damageprofiler; ch_bams_for_qualimap; ch_bams_for_bamutils; ch_bams_for_pmdtools; ch_bams_for_angsd; ch_bams_for_gatk } 
}

if(params.skip_damage_calculation){
  ch_bams_for_damageprofiler.close()
}

if(params.skip_qualimap){
  ch_bams_for_qualimap.close()
}

if(!params.run_pmdtools){
    ch_bams_for_pmdtools.close()
}

if(!params.trim_bam){
    ch_bams_for_bamutils.close()
}

/*
Step 7a: DMG Assessment
*/ 

process damageprofiler {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/damageprofiler", mode: 'copy'

    when:
    !params.skip_damage_calculation

    input:
    file bam from ch_bams_for_damageprofiler
    file fasta from fasta_for_indexing    

    output:
    file "*"
    file "${base}/dmgprof.json" into ch_damageprofiler_results, ch_damageprofiler_for_software_versions

    script:
    base = "${bam.baseName}"
    """
    damageprofiler -i $bam -r $fasta -l ${params.damageprofiler_length} -t ${params.damageprofiler_threshold} -o . 
    """
}

/* 
Step 8: Qualimap
*/

process qualimap {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    file bam from ch_bams_for_qualimap
    file fasta from fasta_for_indexing

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
 Step 9: PMDtools
 */

process pmdtools {
    tag "${bam.baseName}"
    publishDir "${params.outdir}/pmdtools", mode: 'copy'

    when: params.run_pmdtools

    input: 
    file bam from ch_bams_for_pmdtools
    file fasta from fasta_for_indexing

    output:
    file "*.bam" into ch_bam_after_pmdfiltering
    file "*.cpg.range*.txt"

    script:
    //Check which treatment for the libraries was used
    def treatment = params.pmd_udg_type ? (params.pmd_udg_type =='half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(params.snpcapture){
        snpcap = (params.pmdtools_reference_mask != '') ? "--refseq ${params.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    """
    #Run Filtering step 
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --threshold ${params.pmdtools_threshold} $treatment $snpcap --header | samtools view -@ ${task.cpus} -Sb - > "${bam.baseName}".pmd.bam
    #Run Calc Range step
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --deamination --range ${params.pmdtools_range} $treatment $snpcap -n ${params.pmdtools_max_reads} > "${bam.baseName}".cpg.range."${params.pmdtools_range}".txt 
    """
}

/*
* Step 10 - BAM Trimming step using bamUtils 
* Can be used for UDGhalf protocols to clip off -n bases of each read
*/

process bam_trim {
    tag "${prefix}" 
    publishDir "${params.outdir}/trimmed_bam", mode: 'copy'
 
    when: params.trim_bam

    input:
    file bam from ch_bams_for_bamutils  

    output: 
    file "*.trimmed.bam" into ch_trimmed_bam_for_genotyping
    file "*.{bai,csi}"

    script:
    prefix="${bam.baseName}"
    softclip = "${params.bamutils_softclip}" ? '-c' : '' 
    size = "${params.large_ref}" ? '-c' : ''
    """
    bam trimBam $bam tmp.bam -L ${params.bamutils_clip_left} -R ${params.bamutils_clip_right} ${softclip}
    samtools sort -@ ${task.cpus} tmp.bam -o ${prefix}.trimmed.bam 
    samtools index "${size}" ${prefix}.trimmed.bam
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
 * Step 11a - Output Description HTML
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
 * Step 11b - Parse software version numbers
 */
process get_software_versions {

    input:
    file json from ch_damageprofiler_for_software_versions

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
    cat $json &> v_damageprofiler.txt 2>&1 || true 
    
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * Step 11c - MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc_raw/*') from ch_fastqc_results.collect().ifEmpty([])
    file('fastqc/*') from ch_fastqc_after_clipping.collect().ifEmpty([])
    file ('software_versions/software_versions_mqc*') from software_versions_yaml.collect().ifEmpty([])
    file ('adapter_removal/*') from ch_adapterremoval_logs.collect().ifEmpty([])
    file ('flagstat/*') from ch_flagstat_for_multiqc.collect().ifEmpty([])
    file ('flagstat_filtered/*') from ch_bam_filtered_flagstat_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from ch_preseq_results.collect().ifEmpty([])
    file ('damageprofiler/dmgprof*/*') from ch_damageprofiler_results.collect().ifEmpty([])
    file ('qualimap/qualimap*/*') from ch_qualimap_results.collect().ifEmpty([])
    file ('markdup/*') from ch_markdup_results_for_multiqc.collect().ifEmpty([])
    file ('dedup*/*') from ch_dedup_results_for_multiqc.collect().ifEmpty([])
    file ('fastp/*') from ch_fastp_for_multiqc.collect().ifEmpty([])

    file workflow_summary from create_workflow_summary(summary)

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
 * Step 11d - Completion e-mail notification
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
    if(workflow.container) email_fields['summary']['Docker image'] = workflow.container
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList){
                log.warn "[nf-core/eager] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
    }

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
    def smail_fields = [ email: params.email, subject: subject, email_txt: email_txt, email_html: email_html, baseDir: "$baseDir", mqcFile: mqc_report, mqcMaxSize: params.maxMultiqcEmailFileSize.toBytes() ]
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
    def output_d = new File( "${params.outdir}/pipeline_info/" )
    if( !output_d.exists() ) {
      output_d.mkdirs()
    }
    def output_hf = new File( output_d, "pipeline_report.html" )
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File( output_d, "pipeline_report.txt" )
    output_tf.withWriter { w -> w << email_txt }

    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";

    if (workflow.stats.ignoredCountFmt > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCountFmt} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCountFmt} ${c_reset}"
    }

    if(workflow.success){
        log.info "${c_purple}[nf-core/eager]${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}[nf-core/eager]${c_red} Pipeline completed with errors${c_reset}"
    }

}


def nfcoreHeader(){
    // Log colors ANSI codes
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";

    return """    ${c_dim}----------------------------------------------------${c_reset}
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/eager v${workflow.manifest.version}${c_reset}
    ${c_dim}----------------------------------------------------${c_reset}
    """.stripIndent()
}

def checkHostname(){
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if(params.hostnames){
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if(hostname.contains(hname) && !workflow.profile.contains(prof)){
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}
