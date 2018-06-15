#!/usr/bin/env nextflow
/*
========================================================================================
                         nf-core/EAGER2
========================================================================================
 EAGER2 Analysis Pipeline. Started 2018-06-05
 #### Homepage / Documentation
 https://github.com/nf-core/EAGER2
 #### Authors
 Alexander Peltzer apeltzer <alex.peltzer@gmail.com> - https://github.com/apeltzer>
 James A. Fellows Yates <jfy133@gmail.com> - https://github.com/jfy133
----------------------------------------------------------------------------------------
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
      --lanes                       If your fastq files are per lane then you must specify a pattern to extract the lane from the fastq filename.
      --reads                       Path to input data (must be surrounded with quotes)
      TODO: Maybe convert the above into a pattern to find the _R1_, _R2_ in the fastq file name ?
      --indir                       Path to input data dir containing a directory for each libraries fastqs
      --genome                      Name of iGenomes reference
      -profile                      Hardware config to use. docker / aws

    Options:
      --singleEnd                   Specifies that the input is single end reads

    References                      If not specified in the configuration file or you wish to overwrite any of the references.
      --fasta                       Path to Fasta reference

    Other options:
      --outdir                      The output directory where the results will be saved
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */
params.help = false
params.name = false

// Configurable analysis variables
params.lanes = "^.*_(L[0-9]+)_.*\\.fastq\\.gz\$"
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.multiqc_config = "$baseDir/conf/multiqc_config.yaml"
params.email = false
params.plaintext_email = false

lane_pattern = java.util.regex.Pattern.compile (params.lanes)
multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")

params.clip_forward_adaptor = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
params.clip_reverse_adaptor = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA"

params.clip_3pclip = 0
params.clip_5p_clip = 0

params.clip_min_length = 30
params.clip_min_quality = 20
params.clip_min_overlap = 1

params.merge_method ="AdapterRemoval"

params.mapq = 30
params.bwaalnn = 0.04
params.bwaalnl = 32
params.bwaalnk = 2

// Show help emssage
if (params.help){
    helpMessage()
    exit 0
}

// Validate inputs
if ( params.fasta ){
    fasta = file(params.fasta)
    if( !fasta.exists() ) exit 1, "Fasta file not found: ${params.fasta}"
}
//
// NOTE - THIS IS NOT USED IN THIS PIPELINE, EXAMPLE ONLY
// If you want to use the above in a process, define the following:
//   input:
//   file fasta from fasta
//


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}

/*
 * Create a channel for input read files
 */
//Channel
//    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
//    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
//    .into { ch_read_files_clip, ch_read_files_fastqc }

/*
 * Setup channels
 */

ch_fastq_files = Channel.fromFilePairs ("${params.indir}/*/*.fastq.gz",size: -1) { file -> def m = lane_pattern.matcher (file.name); if ( params.lanes.isEmpty () ) { "" } else if ( m.matches () ) { file.parent.name + "_" + m.group (1) } else { throw new Exception ("You specified a lane pattern '" + params.lanes + "' but no match could be found in file name '" + file.name + "'. Set the lane pattern using --lanes (use \"\" if your input files have no lane)"); } }

ch_fastq_files_single = Channel.create ()
ch_fastq_files_paired = Channel.create ()

ch_fastqc_merge_group_single = Channel.create ()
ch_fastqc_merge_group_paired = Channel.create ()

ch_fastqc_results = Channel.create ()

ch_adapter_clip_log_slurp_single = Channel.create ()
ch_adapter_clip_log_slurp_paired = Channel.create ()

ch_adapter_clip_log_single = Channel.create ()
ch_adapter_clip_log_paired = Channel.create ()

ch_adapter_clip_log_parsed_single = Channel.create ()
ch_adapter_clip_log_parsed_paired = Channel.create ()

ch_adapters_clipped_single = Channel.create ()
ch_adapters_clipped_paired = Channel.create ()

ch_bwa_align_single = Channel.create ()
ch_bwa_align_paired = Channel.create ()

ch_bwa_align = ch_bwa_align_single.mix (ch_bwa_align_paired)
ch_bwa_samse = Channel.create ()

ch_samtools_merge_group = Channel.create ()
ch_samtools_view = Channel.create ()
ch_samtools_flagstat_default = Channel.create ()
ch_samtools_flagstat_extract_mapped_sorted = Channel.create ()
ch_samtools_flagstat_dedup_sorted = Channel.create ()
ch_samtools_extract_mapped_sorted = Channel.create ()
ch_samtools_index = Channel.create ()

ch_picard_clean = Channel.create ()
ch_dedup = Channel.create ()

ch_fastq_files.view().choice (ch_fastq_files_single, ch_fastq_files_paired) { t -> t[1].size () == 1 ? 0 : 1 }

//ch_fastq_files_single.tap (ch_fastqc_merge_group_single)
//ch_fastq_files_paired.tap (ch_fastqc_merge_group_paired)

ch_fastqc_merge_single = ch_fastqc_merge_group_single.map { t -> tuple( String.join ("_", Arrays.asList (t[0].split ("_")).subList (0,t[0].split ("_").size () -1)), t[1]) }.groupTuple ()

ch_adapter_clip_log_merge_single = ch_adapter_clip_log_parsed_single.map { t -> tuple( String.join ("_", Arrays.asList (t[0].split ("_")).subList (0,t[0].split ("_").size () -1)), t[1]) }.groupTuple ()
ch_adapter_clip_log_merge_paired = ch_adapter_clip_log_parsed_paired.map { t -> tuple( String.join ("_", Arrays.asList (t[0].split ("_")).subList (0,t[0].split ("_").size () -1)), t[1]) }.groupTuple ()
ch_adapter_clip_log_merge = ch_adapter_clip_log_merge_single.mix (ch_adapter_clip_log_merge_paired)

ch_samtools_merge = ch_samtools_merge_group.map { t -> tuple( String.join ("_", Arrays.asList (t[0].split ("_")).subList (0,t[0].split ("_").size () -1)), t[1]) }.groupTuple ()



// Header log info
log.info "========================================="
log.info " nf-core/EAGER2 v${params.version}"
log.info "========================================="
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Fasta Ref']    = params.fasta
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Input dir']    = params.indir
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


// Check that Nextflow version is up to date enough
// try / throw / catch works for NF versions < 0.25 when this was implemented
try {
    if( ! nextflow.version.matches(">= $params.nf_required_version") ){
        throw GroovyException('Nextflow version too old')
    }
} catch (all) {
    log.error "====================================================\n" +
              "  Nextflow version $params.nf_required_version required! You are running v$workflow.nextflow.version.\n" +
              "  Pipeline execution will continue, but things may break.\n" +
              "  Please run `nextflow self-update` to update Nextflow.\n" +
              "============================================================"
}


/*
 * Parse software version numbers: TODO testing this
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $params.version > v_pipeline.txt
    echo $workflow.nextflow.version > v_nextflow.txt
    fastqc --version > v_fastqc.txt
    echo \$(bwa 2>&1) > v_bwa.txt
    samtools --version > v_samtools.txt
    AdapterRemoval -version > v_adapterremoval.txt
    echo \$(ClipAndMerge 2>&1) v_clipandmerge.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    dedup -h > v_dedup.txt
    #angsd > v_angsd.txt
    #realignsamfile > v_circularmapper.txt
    #schmutzi > v_schmutzi.txt
    gatk --version > v_gatk.txt
    qualimap --version > v_qualimap.txt
    vcf2genome > v_vcf2genome.txt
    multiqc --version > v_multiqc.txt
    scrape_software_versions.py > software_versions_mqc.yaml
    """
}



/*
 * STEP 1 - FastQC
 */
process fastqc_merge_single {
    tag "$name"
    publishDir "${params.outdir}/01-FastQC", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
      set pair_id, file(fastq_file) from ch_fastqc_merge_single.view ().map { t -> tuple ( t[0], t[1].flatten () ) }

    output:
      file "fastqc/*_fastqc.{zip,html}" into ch_fastqc_results

    shell:
    """\
#!/usr/bin/env bash
mkdir fastqc
fastqc -q --extract -f fastq -o fastqc ${fastq_file}
"""
}

/*
 * STEP 2 - Adapter Clipping
 */

if ( params.merge_method.equals ("AdapterRemoval") ) {

 process adapter_removal_single {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set pair_id, fastq_file_pair from ch_fastq_files_single.tap (ch_fastqc_merge_group_single)

 	output:
 		set pair_id, file("adapters-removed.settings") into ch_adapter_clip_log_slurp_single
 		set pair_id, file("adapters-removed.truncated.gz") into ch_adapters_clipped_single

 	shell:
 	"""\
 #!/usr/bin/env bash

 AdapterRemoval --file1 ${fastq_file_pair[0]} --basename adapters-removed --gzip --threads 4 --trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_min_length} --minquality ${params.clip_min_quality} --minadapteroverlap ${params.clip_min_overlap}
 """
 }

 } else {

 process clip_merge_single {

   executor "slurm"
   cpus = 5
   clusterOptions = "--partition \"short\" --mem 6000"

 input:
   set pair_id, fastq_file_pair from ch_fastq_files_single

 output:
   set pair_id, file("ClipAndMergeStats.log") into ch_adapter_clip_log_slurp_single
   set pair_id, file("adapters-removed.truncated.gz") into ch_adapters_clipped_single

   shell:
   """\
 #!/usr/bin/env bash

 ClipAndMerge -in1 ${reads[0]} -f ${params.clip_forward_adaptor} -r ${params.clip_reverse_adaptor} -trim3p ${params.clip_3p_clip} -trim5p ${params.clip_five_p_clip} -l ${params.clip_min_length} -m ${params.clip_min_overlap} -qt -q ${params.clip_min_quality} -log "ClipAndMergeStats.log"

   """
 }

 }

 process adapter_removal_paired {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set pair_id, fastq_file_pair from ch_fastq_files_paired.tap (ch_fastqc_merge_group_paired)

 	output:
 		set pair_id, file("adapters-removed.settings") into ch_adapter_clip_log_slurp_paired
 		set pair_id, file("adapters-removed.collapsed.gz"), file("adapters-removed.collapsed.truncated.gz") into ch_adapters_clipped_paired

 	shell:
 	"""\
 #!/usr/bin/env bash

 AdapterRemoval --file1 ${fastq_file_pair[0]} --file2 ${fastq_file_pair[1]} --basename adapters-removed --gzip --threads 4 --trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_min_length} --minquality ${params.clip_min_quality} --minadapteroverlap ${params.clip_min_overlap} --collapse
 """
 }

 process adapter_clip_log_slurp_single {

 	input:
 		set pair_id, file(settings) from ch_adapter_clip_log_slurp_single

 	output:
 		set pair_id, stdout into ch_adapter_clip_log_single

 	shell:
 	"""\
 #!/usr/bin/env bash

 cat ${settings}
 """
 }

 process adapter_clip_log_slurp_paired {

 	input:
 		set pair_id, file(settings) from ch_adapter_clip_log_slurp_paired

 	output:
 		set pair_id, stdout into ch_adapter_clip_log_paired

 	shell:
 	"""\
 #!/usr/bin/env bash

 cat ${settings}
 """
 }

 process adapter_clip_log_single {

 	def m = new HashMap ();

 	input:
 		set pair_id, val(settings_content) from ch_adapter_clip_log_single

 	output:
 		set pair_id, val { m.get (pair_id) } into ch_adapter_clip_log_parsed_single

 	exec:
 		def stack = new ArrayDeque ()
 		def lm = new HashMap ()
 		settings_content.split ("\\R", -1).each { if ( it.startsWith ("[") && it.endsWith ("]") ) { stack.push (it) } else if ( it.isEmpty () && stack.size () > 0 ) { stack.pop () } else if ( stack.size () > 0 && stack.peekLast ().equals ("[Trimming statistics]") ) { lm.put (it.split (":")[0],it.split (":")[1].trim ()); } }
 		m.put (pair_id, lm)
 }

 process adapter_clip_log_paired {

 	def m = new HashMap ();

 	input:
 		set pair_id, val(settings_content) from ch_adapter_clip_log_paired

 	output:
 		set pair_id, val { m.get (pair_id) } into ch_adapter_clip_log_parsed_paired

 	exec:
 		def stack = new ArrayDeque ()
 		def lm = new HashMap ()
 		settings_content.split ("\\R", -1).each { if ( it.startsWith ("[") && it.endsWith ("]") ) { stack.push (it) } else if ( it.isEmpty () && stack.size () > 0 ) { stack.pop () } else if ( stack.size () > 0 && stack.peekLast ().equals ("[Trimming statistics]") ) { lm.put (it.split (":")[0],it.split (":")[1].trim ()); } }
 		m.put (pair_id, lm)
 }

 process adapter_clip_log_merge {

 	input:
 		set val(library_id), settings_list from ch_adapter_clip_log_merge

 	exec:
 		println "settings list";
 		println settings_list;
 		def m = new HashMap ();
 		settings_list.each { it.each { i -> if ( i.key.equals ("Average length of retained reads") ) { m.put (i.key, m.getOrDefault (i.key, 0) + Float.valueOf (i.value) * Long.valueOf (it.get ("Number of retained reads") ) ) } else {  m.put (i.key, m.getOrDefault (i.key, 0) + Long.valueOf (i.value) ) } } };
 		m.put("Average length of retained reads", m.get ("Average length of retained reads") /  m.get("Number of retained reads") );
 		println "m";
 		println m;

 		def fnb = new StringBuilder ();
 		fnb.append (params.outdir);
 		fnb.append (File.separator);
 		fnb.append (library_id);
 		fnb.append (File.separator);
 		fnb.append ("adapter-removal.settings");
 		def ofile = new File ( fnb.toString () );
 		ofile.getParentFile ().mkdirs ();
 		ofile.withWriter { w -> w.println "[Trimming statistics]"; m.each { w.println it.key + ": " + it.value; } };
 }

 process prefix_fastq_single {

 	executor "slurm"
 	cpus = 4
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set pair_id, file(truncated_file) from ch_adapters_clipped_single

 	output:
 		set pair_id, file("adapters-removed.prefixed.fastq.gz") into ch_bwa_align_single

 	shell:
 	"""\

 #!/usr/bin/env bash
 AdapterRemovalFixPrefix ${truncated_file} adapters-removed.prefixed.fastq.gz
 """
 }

 process prefix_fastq_paired {

 	executor "slurm"
 	cpus = 4
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set pair_id, file(collapsed_file), file(collapsed_truncated_file) from ch_adapters_clipped_paired

 	output:
 		set pair_id, file("adapters-removed.combined.prefixed.fastq.gz") into ch_bwa_align_paired

 	shell:
 	"""\
 #!/usr/bin/env bash

 echo "${pair_id}"
 zcat ${collapsed_file} ${collapsed_truncated_file} | gzip > adapters-removed.combined.fastq.gz
 AdapterRemovalFixPrefix adapters-removed.combined.fastq.gz adapters-removed.combined.prefixed.fastq.gz
 """
 }

 /*
  * STEP 2 - Adapter Clipping with Clip and Merge
  */



 /*
  * STEP 2 - Mapping
  */

 process bwa_align {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"medium\" --time \"2-0\" --mem 16000"

 	input:
 		set library_lane_id, file(combined_prefixed_file) from ch_bwa_align

 	output:
 		set library_lane_id, file(combined_prefixed_file), file("adapters-removed.combined.prefixed.sai") into ch_bwa_samse

 	shell:
 	"""\
 #!/usr/bin/env bash

 bwa aln -t 4 ${params.fasta} ${combined_prefixed_file} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk}  -f adapters-removed.combined.prefixed.sai
 """
 }

 process bwa_samse {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"medium\" --time \"2-0\" --mem 16000"

 	input:
 		set library_lane_id, file(combined_prefixed_file), file(combined_prefixed_sai_file) from ch_bwa_samse

 	output:
 		set library_lane_id, file("adapters-removed.combined.prefixed.sam") into ch_samtools_view

 	shell:
 	"""\
 #!/usr/bin/env bash

 bwa samse -r "@RG\tID:ILLUMINA-${library_lane_id}\tSM:${library_lane_id}\tPL:illumina" ${reference_path} ${combined_prefixed_sai_file} ${combined_prefixed_file} -f adapters-removed.combined.prefixed.sam
 """
 }

 process samtools_view {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 2000"

 	input:
 		set library_lane_id, file(combined_prefixed_sam_file) from ch_samtools_view

 	output:
 		set library_lane_id, file("adapters-removed.combined.prefixed.bam") into ch_samtools_merge_group

 	shell:
 	"""\
 #!/usr/bin/env bash

 samtools view -@ 4 -bS ${combined_prefixed_sam_file} -o adapters-removed.combined.prefixed.bam
 """
 }

 process samtools_merge {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set library_id, file("bams") from ch_samtools_merge

 	output:
 		set library_id, file("${library_id}.merged.bam") into ch_samtools_flagstat_default
 		set library_id, file("${library_id}.merged.bam") into ch_samtools_extract_mapped_sorted

 	shell:
 	"""\
 #!/usr/bin/env bash

 samtools merge ${library_id}.merged.bam ${bams}
 """
 }

 process samtools_flagstat_default {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 2000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_flagstat_default

 	shell:
 	"""\
 #!/usr/bin/env bash

 samtools flagstat ${library_bam} > ${library_bam}.stats
 mkdir -p ${params.outdir}/${library_id}
 cp ${library_bam}.stats ${params.outdir}/${library_id}/
 """
 }

 process samtools_extract_mapped_sorted {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 40000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_extract_mapped_sorted

 	output:
 		set library_id, file("${library_id}.merged.sorted.qf.bam") into ch_samtools_flagstat_extract_mapped_sorted
 		set library_id, file("${library_id}.merged.sorted.qf.bam") into ch_samtools_index

 	shell:
 	"""\
 #!/usr/bin/env bash
 samtools view -F4 -q ${params.mapq} -u ${library_bam} | samtools sort -@ 4 -m 8G - -o ${library_id}.merged.sorted.qf.bam
 mkdir -p ${params.outdir}/${library_id}
 cp ${library_id}.merged.sorted.qf.bam ${params.outdir}/${library_id}/
 """

 }

 process samtools_flagstat_extract_mapped_sorted {

 	executor "slurm"
 	cpus = 5
 	clusterOptions = "--partition \"short\" --mem 2000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_flagstat_extract_mapped_sorted

 	shell:
 	"""\
 #!/usr/bin/env bash

 samtools flagstat ${library_bam} > ${library_bam}.stats
 cp ${library_bam}.stats ${params.outdir}/${library_id}/
 """
 }


 process samtools_index {

 	executor "slurm"
 	cpus = 2
 	clusterOptions = "--partition \"short\" --mem 2000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_index

 	output:
 		set library_id, file(library_bam), file("${library_bam}.bai") into ch_picard_clean

 	shell:
 	"""\
 #!/usr/bin/env bash
 samtools index ${library_bam}
 mkdir -p ${params.outdir}/${library_id}
 cp ${library_bam}.bai ${params.outdir}/${library_id}/
 """
 }

 process picard_clean {

 	executor "slurm"
 	cpus = 2
 	clusterOptions = "--partition \"short\" --mem 6000"

 	input:
 		set library_id, file(library_bam), file(library_bai) from ch_picard_clean

 	output:
 		set library_id, file("${library_id}.merged.sorted.qf.clean.bam") into ch_dedup

 	shell:
 	"""\
 #!/usr/bin/env bash
 picard CleanSam INPUT=${library_bam} OUTPUT=${library_id}.merged.sorted.qf.clean.bam
 """
 }

 process dedup {

 	executor "slurm"
 	cpus = 2
 	clusterOptions = "--partition \"medium\" --mem 32000"

 	input:
 		set library_id, file(library_bam) from ch_dedup

 	output:
 		set library_id, file("${library_id}.merged.sorted.qf.clean_rmdup.bam") into ch_samtools_sort_dedup

 	shell:
 	"""\
 #!/usr/bin/env bash
 dedup -i ${library_bam} -o .
 """
 }

 process samtools_sort_dedup {

 	executor "slurm"
 	cpus = 5
         clusterOptions = "--partition \"short\" --mem 40000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_sort_dedup

 	output:
 		set library_id, file("${library_id}.merged.sorted.qf.clean_rmdup.sorted.bam") into ch_samtools_index_dedup_sorted

 	shell:
 	"""\
 samtools sort -@ 4 -m 8G ${library_bam} -o ${library_id}.merged.sorted.qf.clean_rmdup.sorted.bam
 mkdir -p ${params.outdir}/${library_id}
 cp ${library_id}.merged.sorted.qf.clean_rmdup.sorted.bam ${params.outdir}/${library_id}/
 """
 }

 process samtools_index_dedup_sorted {

 	executor "slurm"
 	cpus = 2
 	clusterOptions = "--partition \"short\" --mem 2000"

 	input:
 		set library_id, file(library_bam) from ch_samtools_index_dedup_sorted

 	output:
 		set library_id, file(library_bam) into ch_samtools_flagstat_dedup_sorted

 	shell:
 	"""\
 samtools index ${library_bam}
 mkdir -p ${params.outdir}/${library_id}
 cp ${library_bam}.bai ${params.outdir}/${library_id}/
 """
 }

 process samtools_flagstat_dedup_sorted {

 	input:
 		set library_id, file(bam_file) from ch_samtools_flagstat_dedup_sorted

 	shell:
 	"""\
 #!/usr/bin/env bash

 samtools flagstat ${bam_file} > ${bam_file}.stats
 mkdir -p ${params.outdir}/${library_id}
 cp ${bam_file}.stats ${params.outdir}/${library_id}/
 """
 }

/*
Step 2.1: AdapterRemovalfixprefix
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
    file ('software_versions/*') from software_versions_yaml

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc -f $rtitle $rfilename --config $multiqc_config .
    touch foo_multiqc_report.html
    touch foo_data
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

    shell:
    """\
#!/usr/bin/env bash

#markdown_to_html.r ${output_docs} results_description.html
touch results_description.html
"""
}



/*
 * Completion e-mail notification
 */
workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/EAGER2] Successful: $workflow.runName"
    if(!workflow.success){
      subject = "[nf-core/EAGER2] FAILED: $workflow.runName"
    }
    def email_fields = [:]
    email_fields['version'] = params.version
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
    email_fields['software_versions'] = software_versions
    email_fields['software_versions']['Nextflow Build'] = workflow.nextflow.build
    email_fields['software_versions']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

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
          log.info "[nf-core/EAGER2] Sent summary e-mail to $params.email (sendmail)"
        } catch (all) {
          // Catch failures and try with plaintext
          [ 'mail', '-s', subject, params.email ].execute() << email_txt
          log.info "[nf-core/EAGER2] Sent summary e-mail to $params.email (mail)"
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

    log.info "[nf-core/EAGER2] Pipeline Complete"

}
