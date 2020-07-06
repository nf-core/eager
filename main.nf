#!/usr/bin/env nextflow
/*
============================================================================================================
                         nf-core/eager
============================================================================================================
 EAGER Analysis Pipeline. Started 2018-06-05
 #### Homepage / Documentation
 https://github.com/nf-core/eager
 #### Authors
 For a list of authors and contributors, see: https://github.com/nf-core/eager/tree/dev#authors-alphabetical
============================================================================================================
*/
///////////////////////////////////////////////////////////////////////////////
/* --                SET UP CONFIGURATION VARIABLES                       -- */
///////////////////////////////////////////////////////////////////////////////

//Parse JSON schema and provide params with usage
def paramsWithUsage = CommonFunctions.readParamsFromJsonSettings("$baseDir/nextflow_schema.json")

// Show help message
params.help = false
if (params.help){
    log.info CommonFunctions.nfcoreHeader(params, workflow)
    log.info CommonFunctions.helpMessage(paramsWithUsage, workflow)
    exit 0
}

// Small console separator to make it easier to read errors after launch
println ""

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// Validate reference inputs
if ( params.fasta.isEmpty () ){
    exit 1, "[nf-core/eager] error: please specify --fasta with the path to your reference"
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
        file "*.{fa,fn,fna,fasta}" into ch_fasta_for_bwaindex,ch_fasta_for_bt2index,ch_fasta_for_faidx,ch_fasta_for_seqdict,ch_fasta_for_circulargenerator,ch_fasta_for_circularmapper,ch_fasta_for_damageprofiler,ch_fasta_for_qualimap,ch_fasta_for_pmdtools,ch_fasta_for_genotyping_ug,ch_fasta_for_genotyping_hc,ch_fasta_for_genotyping_freebayes,ch_fasta_for_genotyping_pileupcaller,ch_fasta_for_vcf2genome,ch_fasta_for_multivcfanalyzer,ch_fasta_for_genotyping_angsd

        script:
        rm_zip = zipped_fasta - '.gz'
        """
        pigz -f -d -p ${task.cpus} $zipped_fasta
        """
        }
       
    } else {
    fasta_for_indexing = Channel
    .fromPath("${params.fasta}", checkIfExists: true)
    .into{ ch_fasta_for_bwaindex; ch_fasta_for_bt2index; ch_fasta_for_faidx; ch_fasta_for_seqdict; ch_fasta_for_circulargenerator; ch_fasta_for_circularmapper; ch_fasta_for_damageprofiler; ch_fasta_for_qualimap; ch_fasta_for_pmdtools; ch_fasta_for_genotyping_ug; ch_fasta__for_genotyping_hc; ch_fasta_for_genotyping_hc; ch_fasta_for_genotyping_freebayes; ch_fasta_for_genotyping_pileupcaller; ch_fasta_for_vcf2genome; ch_fasta_for_multivcfanalyzer;ch_fasta_for_genotyping_angsd }
    
    lastPath = params.fasta.lastIndexOf(File.separator)
    bwa_base = params.fasta.substring(lastPath+1)
    bt2_base = params.fasta.substring(lastPath+1)
}

// Check that fasta index file path ends in '.fai'
if (params.fasta_index && !params.fasta_index.endsWith(".fai")) {
    exit 1, "The specified fasta index file (${params.fasta_index}) is not valid. Fasta index files should end in '.fai'."
}

// Check if genome exists in the config file. params.genomes is from igenomes.conf, params.genome specified by user
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "[nf-core/eager] error: the provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}."
}

// Mapper sanity checking
if (params.mapper != 'bwaaln' && !params.mapper == 'circularmapper' && !params.mapper == 'bwamem' && !params.mapper == "bowtie2"){
    exit 1, "[nf-core/eager] error: invalid mapper option. Options are: 'bwaaln', 'bwamem', 'circularmapper', 'bowtie2'. Default: 'bwaaln'. You gave: --mapper '${params.mapper}'."
}

if (params.mapper == 'bowtie2' && params.bt2_alignmode != 'local' && params.bt2_alignmode != 'end-to-end' ) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 alignment mode. Options: 'local', 'end-to-end'. You gave: -bt2_alignmode '${params.bt2_alignmode}'"
}

if (params.mapper == 'bowtie2' && params.bt2_sensitivity != 'no-preset' && params.bt2_sensitivity != 'very-fast' && params.bt2_sensitivity != 'fast' && params.bt2_sensitivity != 'sensitive' && params.bt2_sensitivity != 'very-sensitive' ) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 sensitivity mode. Options: 'no-preset', 'very-fast', 'fast', 'sensitive', 'very-sensitive'. Options are for both alignmodes You gave: --bt2_sensitivity '${params.bt2_sensitivity}'."
}

if (params.bt2n != 0 && params.bt2n != 1) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 --bt2n (-N) parameter. Options: 0, 1. You gave: --bt2n ${params.bt2n}."

}

// Index files provided? Then check whether they are correct and complete
if( params.bwa_index != '' && (params.mapper == 'bwaaln' | params.mapper == 'bwamem')){
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir =  params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)

    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .ifEmpty { exit 1, "[nf-core/eager] error: bwa indicies not found in: ${bwa_dir}." }
        .into {bwa_index; bwa_index_bwamem}

    bt2_index = ''
}

if( params.bt2_index != '' && params.mapper == 'bowtie2' ){
    lastPath = params.bt2_index.lastIndexOf(File.separator)
    bt2_dir =  params.bt2_index.substring(0,lastPath+1)
    bt2_base = params.bt2_index.substring(lastPath+1)

    Channel
        .fromPath(bt2_dir, checkIfExists: true)
        .ifEmpty { exit 1, "[nf-core/eager] error: bowtie2 indicies not found in: ${bt2_dir}." }
        .into {bt2_index; bt2_index_bwamem}

    bwa_index = ''
}

// Validate BAM input isn't set to paired_end
if ( params.bam && !params.single_end ) {
  exit 1, "[nf-core/eager] error: bams can only be specified with --single_end. Please check input command."
}

// Validate that skip_collapse is only set to True for paired_end reads!
if (!has_extension(params.input, "tsv") && params.skip_collapse  && params.single_end){
    exit 1, "[nf-core/eager] error: --skip_collapse can only be set for paired_end samples."
}

// Strip mode sanity checking
if (params.strip_input_fastq){
    if (!(['strip','replace'].contains(params.strip_mode))) {
        exit 1, "[nf-core/eager] error: --strip_mode can only be set to strip or replace."
    }
}

if (params.bam_discard_unmapped && params.bam_unmapped_type == '') {
    exit 1, "[nf-core/eager] error: please specify valid unmapped read output format. Options: 'discard', 'bam', 'fastq', 'both'. You gave --bam_unmapped_type '${params.bam_unmapped_type}'."
}

// Bedtools sanity checking
if(params.run_bedtools_coverage && params.anno_file == ''){
  exit 1, "[nf-core/eager] error: you have turned on bedtools coverage, but not specified a BED or GFF file with --anno_file. Please validate your parameters."
}

// BAM filtering sanity checking
if (!params.run_bam_filtering && params.bam_mapping_quality_threshold != 0) {
  exit 1, "[nf-core/eager] error: please turn on BAM filtering if you want to perform mapping quality filtering! Give --run_bam_filtering."
}

if (!params.run_bam_filtering && params.bam_discard_unmapped) {
  exit 1, "[nf-core/eager] error: please turn on BAM filtering before trying to indicate how to deal with unmapped reads! Give --run_bam_filtering."
}

if (params.run_bam_filtering && params.bam_discard_unmapped && params.bam_unmapped_type == '') {
  exit 1, "[nf-core/eager] error: please specify how to deal with unmapped reads. Options: 'discard', 'bam', 'fastq', 'both'."
}

if (params.run_bam_filtering && !params.bam_discard_unmapped && params.bam_unmapped_type != 'discard') {
  exit 1, "[nf-core/eager] error: Please turned on unmapped read discarding, if you have specifed a different unmapped type. Give: --bam_discard_unmapped."
}

// Deduplication sanity checking
if (params.dedupper != 'dedup' && params.dedupper != 'markduplicates') {
  exit 1, "[nf-core/eager] error: Selected deduplication tool is not recognised. Options: 'dedup' or 'markduplicates'. You gave: --dedupper '${params.dedupper}'."
}

// Genotyping sanity checking
if (params.run_genotyping){
  if (params.genotyping_tool != 'ug' && params.genotyping_tool != 'hc' && params.genotyping_tool != 'freebayes' && params.genotyping_tool != 'pileupcaller' && params.genotyping_tool != 'angsd' ) {
  exit 1, "[nf-core/eager] error: please specify a genotyper. Options: 'ug', 'hc', 'freebayes', 'pileupcaller'. You gave: --genotyping_tool '${params.genotyping_tool}'."
  }

  if (params.genotyping_tool == 'ug' && params.gatk_ug_jar == '') {
  exit 1, "[nf-core/eager] error: please specify path to a GATK 3.5 .jar file with --gatk_ug_jar."
  }

  if (params.genotyping_tool == 'ug' && !params.gatk_ug_jar.endsWith('.jar') ) {
    exit 1, "[nf-core/eager] error: please specify path with --gatk_ug_jar to a valid GATK 3.5 binary that ends with .jar!. You gave: --gatk_ug_jar '${params.gatk_ug_jar}'."
  }
  
  if (params.gatk_ug_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_ug_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_ug_out_mode != 'EMIT_ALL_SITES') {
  exit 1, "[nf-core/eager] error: please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. You gave: --gatk_ug_out_mode '${params.gatk_out_mode}'."
  }

  if (params.gatk_hc_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_hc_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_hc_out_mode != 'EMIT_ALL_ACTIVE_SITES') {
  exit 1, "[nf-core/eager] error: please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. You gave: --gatk_out_mode '${params.gatk_out_mode}'."
  }
  
  if (params.genotyping_tool == 'ug' && (params.gatk_ug_genotype_model != 'SNP' && params.gatk_ug_genotype_model != 'INDEL' && params.gatk_ug_genotype_model != 'BOTH' && params.gatk_ug_genotype_model != 'GENERALPLOIDYSNP' && params.gatk_ug_genotype_model != 'GENERALPLOIDYINDEL')) {
    exit 1, "[nf-core/eager] error: please check your UnifiedGenotyper genotype model. Options: 'SNP', 'INDEL', 'BOTH', 'GENERALPLOIDYSNP', 'GENERALPLOIDYINDEL'. You gave: --gatk_ug_genotype_model '${params.gatk_ug_genotype_model}'."
  }

  if (params.genotyping_tool == 'hc' && (params.gatk_hc_emitrefconf != 'NONE' && params.gatk_hc_emitrefconf != 'GVCF' && params.gatk_hc_emitrefconf != 'BP_RESOLUTION')) {
    exit 1, "[nf-core/eager] error: please check your HaplotyperCaller reference confidence parameter. Options: 'NONE', 'GVCF', 'BP_RESOLUTION'. You gave: --gatk_hc_emitrefconf '${params.gatk_hc_emitrefconf}'."
  }

  if (params.genotyping_tool == 'pileupcaller' && ! ( params.pileupcaller_method == 'randomHaploid' || params.pileupcaller_method == 'randomDiploid' || params.pileupcaller_method == 'majorityCall' ) ) {
	exit 1, "[nf-core/eager] error: please check your pileupCaller method parameter. Options: 'randomHaploid', 'randomDiploid', 'majorityCall'. You gave: --pileupcaller_method '${params.pileupcaller_method}'."
  }

  if (params.genotyping_tool == 'pileupcaller' && ( params.pileupcaller_bedfile == '' || params.pileupcaller_snpfile == '' ) ) {
    exit 1, "[nf-core/eager] error: please check your pileupCaller bed file and snp file parameters. You must supply a bed file and a snp file."
  }

  if (params.genotyping_tool == 'angsd' && ! ( params.angsd_glmodel == 'samtools' || params.angsd_glmodel == 'gatk' || params.angsd_glmodel == 'soapsnp' || params.angsd_glmodel == 'syk' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD genotyping model! Options: 'samtools', 'gatk', 'soapsnp', 'syk'. You gave: --angsd_glmodel' ${params.angsd_glmodel}'."
  }

  if (params.genotyping_tool == 'angsd' && ! ( params.angsd_glformat == 'text' || params.angsd_glformat == 'binary' || params.angsd_glformat == 'binary_three' || params.angsd_glformat == 'beagle' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD genotyping model! Options: 'text', 'binary', 'binary_three', 'beagle'. You gave: --angsd_glmodel '${params.angsd_glmodel}'."
  }

  if ( !params.angsd_createfasta && params.angsd_fastamethod != 'random' ) {
    exit 1, "[nf-core/eager] error: to output a ANGSD FASTA file, please turn on FASTA creation with --angsd_createfasta."
  }

  if ( params.angsd_createfasta && !( params.angsd_fastamethod == 'random' || params.angsd_fastamethod == 'common' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD FASTA file creation method. Options: 'random', 'common'. You gave: --angsd_fastamethod '${params.angsd_fastamethod}'."
  }

}

// Consensus sequence generation sanity checking
if (params.run_vcf2genome) {
    if (!params.run_genotyping) {
      exit 1, "[nf-core/eager] error: consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Please check your genotyping parameters."
    }

    if (params.genotyping_tool != 'ug') {
      exit 1, "[nf-core/eager] error: consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Please check your genotyping parameters."
    }
}

// MultiVCFAnalyzer sanity checking
if (params.run_multivcfanalyzer) {
  if (!params.run_genotyping) {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer requires genotyping to be turned on with the parameter --run_genotyping. Please check your genotyping parameters."
  }

  if (params.genotyping_tool != "ug") {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer only accepts VCF files from GATK UnifiedGenotyper. Please check your genotyping parameters."
  }

  if (params.gatk_ploidy != '2') {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer only accepts VCF files generated with a GATK ploidy set to 2. Please check your genotyping parameters."
  }
}

// Metagenomic sanity checking
if (params.run_metagenomic_screening) {
  if ( !params.bam_discard_unmapped ) {
  exit 1, "[nf-core/eager] error: metagenomic classification can only run on unmapped reads. Please supply --bam_discard_unmapped and --bam_unmapped_type 'fastq'."
  }

  if (params.bam_discard_unmapped && params.bam_unmapped_type != 'fastq' ) {
  exit 1, "[nf-core/eager] error: metagenomic classification can only run on unmapped reads in FASTSQ format. Please supply --bam_unmapped_type 'fastq'. You gave: --bam_unmapped_type '${params.bam_unmapped_type}'."
  }

  if (params.metagenomic_tool != 'malt' &&  params.metagenomic_tool != 'kraken') {
    exit 1, "[nf-core/eager] error: metagenomic classification can currently only be run with 'malt' or 'kraken' (kraken2). Please check your classifer. You gave: --metagenomic_tool '${params.metagenomic_tool}'."
  }

  if (params.database == '' ) {
    exit 1, "[nf-core/eager] error: metagenomic classification requires a path to a database directory. Please specify one with --database '/path/to/database/'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_mode != 'BlastN' && params.malt_mode != 'BlastP' && params.malt_mode != 'BlastX') {
    exit 1, "[nf-core/eager] error: unknown MALT mode specified. Options: 'BlastN', 'BlastP', 'BlastX'. You gave: --malt_mode '${params.malt_mode}'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_alignment_mode != 'Local' && params.malt_alignment_mode != 'SemiGlobal') {
    exit 1, "[nf-core/eager] error: unknown MALT alignment mode specified. Options: 'Local', 'SemiGlobal'. You gave: --malt_alignment_mode '${params.malt_alignment_mode}'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_min_support_mode == 'percent' && params.metagenomic_min_support_reads != 1) {
    exit 1, "[nf-core/eager] error: incompatible MALT min support configuration. Percent can only be used with --malt_min_support_percent. You modified --metagenomic_min_support_reads."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_min_support_mode == 'reads' && params.malt_min_support_percent != 0.01) {
    exit 1, "[nf-core/eager] error: incompatible MALT min support configuration. Reads can only be used with --malt_min_supportreads. You modified --malt_min_support_percent."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_memory_mode != 'load' && params.malt_memory_mode != 'page' && params.malt_memory_mode != 'map') {
    exit 1, "[nf-core/eager] error: unknown MALT memory mode specified. Options: 'load', 'page', 'map'. You gave: --malt_memory_mode '${params.malt_memory_mode}'."
  }

  if (!params.metagenomic_min_support_reads.toString().isInteger()){
    exit 1, "[nf-core/eager] error: incompatible min_support_reads configuration. min_support_reads can only be used with integers. --metagenomic_min_support_reads You gave: ${params.metagenomic_min_support_reads}."
  }
}

// MaltExtract Sanity checking
if (params.run_maltextract) {

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "[nf-core/eager] error: MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'."
  }

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "[nf-core/eager] error: MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'."
  }

  if (params.maltextract_taxon_list == '') {
    exit 1, "[nf-core/eager] error: MaltExtract requires a taxon list specify target taxa of interest. Specify the file with --params.maltextract_taxon_list."
  }

  if (params.maltextract_filter != 'def_anc' && params.maltextract_filter != 'default' && params.maltextract_filter != 'ancient' && params.maltextract_filter != 'scan' && params.maltextract_filter != 'crawl' && params.maltextract_filter != 'srna') {
    exit 1, "[nf-core/eager] error: unknown MaltExtract filter specified. Options are: 'def_anc', 'default', 'ancient', 'scan', 'crawl', 'srna'. You gave: --maltextract_filter '${params.maltextract_filter}'."
  }

}

// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

if (workflow.profile.contains('awsbatch')) {
    CommonFunctions.checkAWSbatch(params, workflow)
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config = file("$baseDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_eager_logo = file("$baseDir/docs/images/nf-core_eager_logo.png")
ch_output_docs = file("$baseDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$baseDir/docs/images/", checkIfExists: true)
where_are_my_files = file("$baseDir/assets/where_are_my_files.txt")

///////////////////////////////////////////////////
/* --    INPUT FILE LOADING AND VALIDATING    -- */
///////////////////////////////////////////////////

// check we have valid --reads or --input
if (params.input == null) {
  exit 1, "[nf-core/eager] error: --input was not supplied! Please see --help and documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
tsv_path = null
if (params.input && (has_extension(params.input, "tsv"))) tsv_path = params.input

ch_input_sample = Channel.empty()
if (tsv_path) {

    tsv_file = file(tsv_path)
    ch_input_sample = extract_data(tsv_file)

} else if (params.input && !has_extension(params.input, "tsv")) {

    log.info ""
    log.info "No TSV file provided - creating TSV from supplied directory."
    log.info "Reading path(s): ${params.input}\n"
    inputSample = retrieve_input_paths(params.input, params.colour_chemistry, params.single_end, params.single_stranded, params.udg_type, params.bam)
    ch_input_sample = inputSample

} else exit 1, "[nf-core/eager] error: --input file(s) not correctly not supplied or improperly defined, see --help and documentation under 'running the pipeline' for details."

ch_input_sample
  .into { ch_input_sample_downstream; ch_input_sample_check }

///////////////////////////////////////////////////
/* --         INPUT CHANNEL CREATION          -- */
///////////////////////////////////////////////////

// Check we don't have any duplicate file names
ch_input_sample_check
    .map {
      it ->
        def r1 = file(it[8]).getName()
        def r2 = file(it[9]).getName()
        def bam = file(it[10]).getName()

      [r1, r2, bam]

    }
    .collect()
    .map{
       file -> 
       filenames = file
       filenames -= 'NA'
       
       if( filenames.size() != filenames.unique().size() )
           exit 1, "[nf-core/eager] error: You have duplicate input FASTQ and/or BAM file names! All files must have unique names, different directories are not sufficent. Please check your input."
    }

// Drop samples with R1/R2 to fastQ channel, BAM samples to other channel
ch_branched_input = ch_input_sample_downstream.branch{
    fastq: it[8] != 'NA' //These are all fastqs
    bam: it[10] != 'NA' //These are all BAMs
}

//Removing BAM/BAI in case of a FASTQ input
ch_fastq_channel = ch_branched_input.fastq.map {
  samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ->
    [samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2]
}

//Removing R1/R2 in case of BAM input
ch_bam_channel = ch_branched_input.bam.map {
  samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ->
    [samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, bam]
}

// Prepare starting channels, here we go
ch_input_for_convertbam = Channel.empty()

ch_bam_channel
  .into { ch_input_for_convertbam; ch_input_for_indexbam; }

// Also need to send raw files for lane merging, if we want to strip fastq
ch_fastq_channel
  .into { ch_input_for_skipconvertbam; ch_input_for_lanemerge_stripfastq }

///////////////////////////////////////////////////
/* --             HEADER LOG INFO             -- */
///////////////////////////////////////////////////

log.info CommonFunctions.nfcoreHeader(params, workflow)
def summary = [:]
summary['Pipeline Name']  = 'nf-core/eager'
summary['Pipeline Version'] = workflow.manifest.version
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Input']        = params.input
summary['Convert input BAM?'] = params.run_convertinputbam ? 'Yes' : 'No'
summary['Fasta Ref']    = params.fasta
summary['BAM Index Type'] = (params.large_ref == "") ? 'BAI' : 'CSI'
if(params.bwa_index || params.bt2_index ) summary['BWA Index'] = "Yes"
summary['Skipping FASTQC?'] = params.skip_fastqc ? 'Yes' : 'No'
summary['Skipping AdapterRemoval?'] = params.skip_adapterremoval ? 'Yes' : 'No'
if (!params.skip_adapterremoval) {
  summary['Skip Read Merging'] = params.skip_collapse ? 'Yes' : 'No'
  summary['Skip Adapter Trimming']  = params.skip_trim  ? 'Yes' : 'No' 
}
summary['Running BAM filtering'] = params.run_bam_filtering ? 'Yes' : 'No'
if (params.run_bam_filtering) {
  summary['Skip Read Merging'] = params.bam_discard_unmapped ? 'Yes' : 'No'
  summary['Skip Read Merging'] = params.bam_unmapped_type
}
summary['Run Fastq Stripping'] = params.strip_input_fastq ? 'Yes' : 'No'
if (params.strip_input_fastq){
    summary['Strip mode'] = params.strip_mode
}
summary['Skipping Preseq?'] = params.skip_preseq ? 'Yes' : 'No'
summary['Skipping Deduplication?'] = params.skip_deduplication ? 'Yes' : 'No'
summary['Skipping DamageProfiler?'] = params.skip_damage_calculation ? 'Yes' : 'No'
summary['Skipping Qualimap?'] = params.skip_qualimap ? 'Yes' : 'No'
summary['Run BAM Trimming?'] = params.run_trim_bam ? 'Yes' : 'No'
summary['Run PMDtools?'] = params.run_pmdtools ? 'Yes' : 'No'
summary['Run Genotyping?'] = params.run_genotyping ? 'Yes' : 'No'
if (params.run_genotyping){
  summary['Genotyping Tool?'] = params.genotyping_tool
  summary['Genotyping BAM Input?'] = params.genotyping_source
}
summary['Run MultiVCFAnalyzer'] = params.run_multivcfanalyzer ? 'Yes' : 'No'
summary['Run VCF2Genome'] = params.run_vcf2genome ? 'Yes' : 'No'
summary['Run SexDetErrMine'] = params.run_sexdeterrmine ? 'Yes' : 'No'
summary['Run Nuclear Contamination Estimation'] = params.run_nuclear_contamination ? 'Yes' : 'No'
summary['Run Bedtools Coverage'] = params.run_bedtools_coverage ? 'Yes' : 'No'
summary['Run Metagenomic Binning'] = params.run_metagenomic_screening ? 'Yes' : 'No'
if (params.run_metagenomic_screening) {
  summary['Metagenomic Tool'] = params.metagenomic_tool
  summary['Run MaltExtract'] = params.run_maltextract ? 'Yes' : 'No'
}
summary['Max Memory']   = params.max_memory
summary['Max CPUs']     = params.max_cpus
summary['Max Time']     = params.max_time
summary['Output Dir']   = params.outdir
summary['Working Dir']  = workflow.workDir
summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['Current Home']   = workflow.homeDir
summary['Current User']   = workflow.userName
summary['Working Dir']    = workflow.workDir
summary['Output Dir']     = params.outdir
summary['Script Dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
if(workflow.profile == 'awsbatch'){
   summary['AWS Region']    = params.awsregion
   summary['AWS Queue']     = params.awsqueue
}
if(params.email) summary['E-mail Address'] = params.email
summary['Config Profile'] = workflow.profile
if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url
if (params.email || params.email_on_fail) {
    summary['E-mail Address']    = params.email
    summary['E-mail on failure'] = params.email_on_fail
    summary['MultiQC maxsize']   = params.max_multiqc_email_size
}
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
log.info "-\033[2m--------------------------------------------------\033[0m-"

// Check the hostnames against configured profiles
CommonFunctions.checkHostname(params)

Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'nf-core-eager-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/eager Workflow Summary'
    section_href: 'https://github.com/nf-core/eager'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

///////////////////////////////////////////////////
/* --          REFERENCE FASTA INDEXING       -- */
///////////////////////////////////////////////////

// BWA Index
if( params.bwa_index == '' && !params.fasta.isEmpty() && (params.mapper == 'bwaaln' || params.mapper == 'bwamem' || params.mapper == 'circularmapper')){
  process makeBWAIndex {
    label 'sc_medium'
    tag "${fasta}"
    publishDir path: "${params.outdir}/reference_genome/bwa_index", mode: 'copy', saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    file fasta from ch_fasta_for_bwaindex
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
    bt2_index = 'none'
}

// bowtie2 Index
if(params.bt2_index == '' && !params.fasta.isEmpty() && params.mapper == "bowtie2"){
  process makeBT2Index {
    label 'sc_medium'
    tag "${fasta}"
    publishDir path: "${params.outdir}/reference_genome/bt2_index", mode: 'copy', saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    file fasta from ch_fasta_for_bt2index
    file where_are_my_files

    output:
    file "BT2Index" into (bt2_index)
    file "where_are_my_files.txt"

    script:
    """
    bowtie2-build $fasta $fasta
    mkdir BT2Index && mv ${fasta}* BT2Index
    """
    }

  bwa_index = 'none'
  bwa_index_bwamem = 'none'

}

// FASTA Index (FAI)
if (params.fasta_index != '') {
  Channel
    .fromPath( params.fasta_index )
    .set { ch_fai_for_skipfastaindexing }
} else {
  Channel
    .empty()
    .set { ch_fai_for_skipfastaindexing }
}

process makeFastaIndex {
    label 'sc_small'
    tag "${fasta}"
    publishDir path: "${params.outdir}/reference_genome/fasta_index", mode: 'copy', saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.fasta_index == '' && !params.fasta.isEmpty() && ( params.mapper == 'bwaaln' || params.mapper == 'bwamem' || params.mapper == 'circularmapper')

    input:
    file fasta from ch_fasta_for_faidx
    file where_are_my_files

    output:
    file "*.fai" into ch_fasta_faidx_index
    file "where_are_my_files.txt"

    script:
    """
    samtools faidx $fasta
    """
}

ch_fai_for_skipfastaindexing.mix(ch_fasta_faidx_index) 
  .into { ch_fai_for_ug; ch_fai_for_hc; ch_fai_for_freebayes; ch_fai_for_pileupcaller; ch_fai_for_angsd }

// Stage dict index file if supplied, else load it into the channel

if (params.seq_dict != '') {
  Channel
    .fromPath( params.seq_dict )
    .set { ch_dict_for_skipdict }
} else {
  Channel
    .empty()
    .set { ch_dict_for_skipdict }
}

process makeSeqDict {
    label 'sc_medium'
    tag "${fasta}"
    publishDir path: "${params.outdir}/reference_genome/seq_dict", mode: 'copy', saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.seq_dict == '' && !params.fasta.isEmpty()

    input:
    file fasta from ch_fasta_for_seqdict
    file where_are_my_files

    output:
    file "*.dict" into ch_seq_dict
    file "where_are_my_files.txt"

    script:
    """
    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M CreateSequenceDictionary R=$fasta O="${fasta.baseName}.dict"
    """
}

ch_dict_for_skipdict.mix(ch_seq_dict)
  .into { ch_dict_for_ug; ch_dict_for_hc; ch_dict_for_freebayes; ch_dict_for_pileupcaller; ch_dict_for_angsd }

//////////////////////////////////////////////////
/* --         BAM INPUT PREPROCESSING        -- */
//////////////////////////////////////////////////

// Convert to FASTQ if re-mapping is requested
process convertBam {
    label 'mc_small'
    tag "$libraryid"
    
    when: 
    params.run_convertinputbam

    input: 
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(bam) from ch_input_for_convertbam 

    output:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file("*fastq.gz"), val('NA') into ch_output_from_convertbam

    script:
    base = "${bam.baseName}"
    """
    samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.converted.fastq.gz
    """ 
}

// If not converted to FASTQ generate pipeline compatible BAM index file (i.e. with correct samtools version) 
process indexinputbam {
  label 'sc_small'
  tag "$libraryid"

  when: 
  bam != 'NA' && !params.run_convertinputbam

  input:
  tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(bam) from ch_input_for_indexbam 

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file("*.{bai,csi}")  into ch_indexbam_for_filtering

  script:
  size = "${params.large_ref}" ? '-c' : ''
  """
  samtools index "${size}" ${bam}
  """
}

// convertbam bypass
    ch_input_for_skipconvertbam.mix(ch_output_from_convertbam)
        .into { ch_convertbam_for_fastp; ch_convertbam_for_fastqc } 

//////////////////////////////////////////////////
/* -- SEQUENCING QC AND FASTQ PREPROCESSING  -- */
//////////////////////////////////////////////////

// Raw sequencing QC - allow user evaluate if sequencing any good?

process fastqc {
    label 'sc_small'
    tag "${libraryid}_L${lane}"
    publishDir "${params.outdir}/FastQC/input_fastq", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: 
    !params.skip_fastqc

    input:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_convertbam_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_prefastqc_for_multiqc

    script:
    if ( seqtype == 'PE' ) {
    """
    fastqc -t ${task.cpus} -q $r1 $r2
    rename 's/_fastqc\\.zip\$/_raw_fastqc.zip/' *_fastqc.zip
    rename 's/_fastqc\\.html\$/_raw_fastqc.html/' *_fastqc.html
    """
    } else {
    """
    fastqc -q $r1
    rename 's/_fastqc\\.zip\$/_raw_fastqc.zip/' *_fastqc.zip
    rename 's/_fastqc\\.html\$/_raw_fastqc.html/' *_fastqc.html
    """
    }
}

// Poly-G clipping for 2-colour chemistry sequencers, to reduce erroenous mapping of sequencing artefacts

if (params.complexity_filter_poly_g) {
  ch_input_for_fastp = ch_convertbam_for_fastp.branch{
    twocol: it[3] == '2' // Nextseq/Novaseq data with possible sequencing artefact
    fourcol: it[3] == '4'  // HiSeq/MiSeq data where polyGs would be true
  }

} else {
  ch_input_for_fastp = ch_convertbam_for_fastp.branch{
    twocol: it[3] == "dummy" // seq/Novaseq data with possible sequencing artefact
    fourcol: it[3] == '4' || it[3] == '2'  // HiSeq/MiSeq data where polyGs would be true
  }

}

process fastp {
    label 'mc_small'
    tag "${libraryid}_L${lane}"
    publishDir "${params.outdir}/FastP", mode: 'copy'

    when: 
    params.complexity_filter_poly_g

    input:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_input_for_fastp.twocol

    output:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file("*.pG.fq.gz") into ch_output_from_fastp
    file("*.json") into ch_fastp_for_multiqc

    script:
    if( seqtype == 'SE' ){
    """
    fastp --in1 ${r1} --out1 "${r1.baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json "${r1.baseName}"_L${lane}_fastp.json 
    """
    } else {
    """
    fastp --in1 ${r1} --in2 ${r2} --out1 "${r1.baseName}.pG.fq.gz" --out2 "${r2.baseName}.pG.fq.gz" -A -g --poly_g_min_len "${params.complexity_filter_poly_g_min}" -Q -L -w ${task.cpus} --json "${libraryid}"_L${lane}_fastp.json 
    """
    }
}

// Colour column only useful for fastp, so dropping now to reduce complexity downstream
ch_input_for_fastp.fourcol
  .map {
      def samplename = it[0]
      def libraryid  = it[1]
      def lane = it[2]
      def seqtype = it[4]
      def organism = it[5]
      def strandedness = it[6]
      def udg = it[7]
      def r1 = it[8]
      def r2 = seqtype == "PE" ? it[9] : 'NA'
      
      [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

    }
 .set { ch_skipfastp_for_merge }

ch_output_from_fastp
  .map{
    def samplename = it[0]
    def libraryid  = it[1]
    def lane = it[2]
    def seqtype = it[4]
    def organism = it[5]
    def strandedness = it[6]
    def udg = it[7]
    def r1 = it[8].getClass() == ArrayList ? it[8].sort()[0] : it[8]
    def r2 = seqtype == "PE" ? it[8].sort()[1] : 'NA'

    [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

  }
  .set{ ch_fastp_for_merge }

ch_skipfastp_for_merge.mix(ch_fastp_for_merge)
  .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval } 

// Sequencing adapter clipping and optional paired-end merging in preparation for mapping

process adapter_removal {
    label 'mc_small'
    tag "${libraryid}_L${lane}"
    publishDir "${params.outdir}/AdapterRemoval", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_fastp_for_adapterremoval

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("output/*{combined.fq,.se.truncated,pair1.truncated}.gz") into ch_output_from_adapterremoval_r1
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("output/*pair2.truncated.gz") optional true into ch_output_from_adapterremoval_r2
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("output/*.settings") into ch_adapterremoval_logs

    when: 
    !params.skip_adapterremoval

    script:
    base = "${r1.baseName}_L${lane}"
    //This checks whether we skip trimming and defines a variable respectively
    trim_me = params.skip_trim ? '' : "--trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}"
    collapse_me = params.skip_collapse ? '' : '--collapse'
    preserve5p = params.preserve5p ? '--preserve5p' : ''
    mergedonly = params.mergedonly ? "Y" : "N"
    
    //PE mode, dependent on trim_me and collapse_me the respective procedure is run or not :-) 
    if ( seqtype == 'PE'  && !params.skip_collapse && !params.skip_trim ){
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe ${trim_me} --gzip --threads ${task.cpus} ${collapse_me} ${preserve5p}
    
    #Combine files
    if [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "N" ]; then 
      cat *.collapsed.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.fq.gz
    elif [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "Y" ] ; then
      cat *.collapsed.gz > output/${base}.pe.combined.fq.gz
    elif [ ${mergedonly} = "Y" ] ; then
      cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.fq.gz
    else
      cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.fq.gz
    fi
   
    mv *.settings output/
    """
    //PE, don't collapse, but trim reads
    } else if ( seqtype == 'PE' && params.skip_collapse && !params.skip_trim ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} ${trim_me} ${collapse_me} ${preserve5p}
    mv *.settings ${base}.pe.pair*.truncated.gz output/
    """
    //PE, collapse, but don't trim reads
    } else if ( seqtype == 'PE' && !params.skip_collapse && params.skip_trim ) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe --gzip --threads ${task.cpus} ${collapse_me} ${trim_me}
    
    if [ ${mergedonly} = "Y" ]; then
      cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.fq.gz
    else
      cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz  > output/${base}.pe.combined.fq.gz
    fi

    mv *.settings output/
    """
    } else if ( seqtype != 'PE' ) {
    //SE, collapse not possible, trim reads
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --basename ${base}.se --gzip --threads ${task.cpus} ${trim_me} ${preserve5p}
    
    mv *.settings *.se.truncated.gz output/
    """
    }
}

// When not collapsing paired-end data, re-merge the R1 and R2 files into single map. Otherwise if SE or collapsed PE, R2 now becomes NA
// Sort to make sure we get consistent R1 and R2 ordered when using `-resume`, even if not needed for FastQC
if ( params.skip_collapse ){
  ch_output_from_adapterremoval_r1
    .mix(ch_output_from_adapterremoval_r2)
    .groupTuple(by: [0,1,2,3,4,5,6])
    .map{
      it -> 
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = file(it[7].sort()[0])
        def r2 = seqtype == "PE" ? file(it[7].sort()[1]) : 'NA'

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

    }
    .into { ch_output_from_adapterremoval; ch_adapterremoval_for_postfastqc }
} else {
  ch_output_from_adapterremoval_r1
    .map{
      it -> 
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = file(it[7])
        def r2 = 'NA'

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]
    }
    .into { ch_output_from_adapterremoval; ch_adapterremoval_for_postfastqc }
}

// AdapterRemoval bypass when not running it
if (!params.skip_adapterremoval) {
    ch_output_from_adapterremoval.mix(ch_fastp_for_skipadapterremoval)
        .filter { it =~/.*combined.fq.gz|.*truncated.gz/ }
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_lanemerge; } 
} else {
    ch_fastp_for_skipadapterremoval
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_lanemerge; } 
}

// Lane merging for libraries sequenced over multiple lanes (e.g. NextSeq)

ch_branched_for_lanemerge = ch_adapterremoval_for_lanemerge
  .groupTuple(by: [0,1,3,4,5,6])
  .branch {
    skip_merge: it[7].size() == 1 // Can skip merging if only single lanes
    merge_me: it[7].size() > 1
  }

process lanemerge {
  label 'mc_tiny'
  tag "${libraryid}"
  publishDir "${params.outdir}/lanemerging", mode: 'copy'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_branched_for_lanemerge.merge_me

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.fq.gz") into ch_lanemerge_for_mapping

  script:
  if ( seqtype == 'PE' && ( params.skip_collapse || params.skip_adapterremoval ) ){
  lane = 0
  """
  cat ${r1} > "${libraryid}"_R1_lanemerged.fq.gz
  cat ${r2} > "${libraryid}"_R2_lanemerged.fq.gz
  """
  } else {
  """
  cat ${r1} > "${libraryid}"_lanemerged.fq.gz
  """
  }

}

// TODO check lane merged skipped non-collapse'd R2
ch_lanemerge_for_mapping
  .map {
      def samplename = it[0]
      def libraryid  = it[1]
      def lane = it[2]
      def seqtype = it[3]
      def organism = it[4]
      def strandedness = it[5]
      def udg = it[6]
      def reads = arrayify(it[7])
      def r1 = it[7].getClass() == ArrayList ? reads[0] : it[7]
      def r2 = reads[1] ? reads[1] : "NA"      

      [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

  }
  .mix(ch_branched_for_lanemerge.skip_merge)
  .into { ch_lanemerge_for_skipmap; ch_lanemerge_for_bwa; ch_lanemerge_for_cm; ch_lanemerge_for_bwamem; ch_lanemerge_for_bt2 } 

// ENA upload doesn't do separate lanes, so merge raw FASTQs for mapped-reads stripping 

// Per-library lane grouping done within process
process lanemerge_stripfastq {
  label 'mc_tiny'
  tag "${libraryid}"

  when: 
  params.strip_input_fastq

  input:
  tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_input_for_lanemerge_stripfastq.groupTuple(by: [0,1,3,4,5,6,7])

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.fq.gz") into ch_fastqlanemerge_for_stripfastq

  script:
  if ( seqtype == 'PE' ){
  lane = 0
  """
  cat ${r1} > "${libraryid}"_R1_lanemerged.fq.gz
  cat ${r2} > "${libraryid}"_R2_lanemerged.fq.gz
  """
  } else {
  """
  cat ${r1} > "${libraryid}"_R1_lanemerged.fq.gz
  """
  }

}

// Post-preprocessing QC to help user check pre-processing removed all sequencing artefacts

process fastqc_after_clipping {
    label 'sc_small'
    tag "${libraryid}_L${lane}"
    publishDir "${params.outdir}/FastQC/after_clipping", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: !params.skip_adapterremoval && !params.skip_fastqc

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_adapterremoval_for_fastqc_after_clipping

    output:
    file("*_fastqc.{zip,html}") into ch_fastqc_after_clipping

    script:
    if ( params.skip_collapse && seqtype == "PE" ) {
    """
    fastqc -t ${task.cpus} -q ${r1} ${r2}
    """
    } else {
    """
    fastqc -t ${task.cpus} -q ${r1}
    """
    }

}

//////////////////////////////////////////////////
/* --    READ MAPPING AND POSTPROCESSING     -- */
//////////////////////////////////////////////////

// bwa aln as standard aDNA mapper

process bwa {
    label 'mc_medium'
    tag "${libraryid}"
    publishDir "${params.outdir}/mapping/bwa", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_bwa
    file index from bwa_index.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mapped.bam"), file("*.{bai,csi}") into ch_output_from_bwa   

    when: 
    params.mapper == 'bwaaln'

    script:
    size = "${params.large_ref}" ? '-c' : ''
    fasta = "${index}/${bwa_base}"

    //PE data without merging, PE data without any AR applied
    if ( seqtype == 'PE' && ( params.skip_collapse || params.skip_adapterremoval ) ){
    """
    bwa aln -t ${task.cpus} $fasta ${r1} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r1.sai
    bwa aln -t ${task.cpus} $fasta ${r2} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $fasta ${libraryid}.r1.sai ${libraryid}.r2.sai ${r1} ${r2} | samtools sort -@ ${task.cpus} -O bam - > ${libraryid}.mapped.bam
    samtools index "${size}" "${libraryid}".mapped.bam
    """
    } else {
    //PE collapsed, or SE data 
    """
    bwa aln -t ${task.cpus} ${fasta} ${r1} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $fasta ${libraryid}.sai $r1 | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}".mapped.bam
    samtools index "${size}" "${libraryid}".mapped.bam
    """
    }
    
}

// bwa mem for more complex or for modern data mapping

process bwamem {
    label 'mc_medium'
    tag "$libraryid"
    publishDir "${params.outdir}/mapping/bwamem", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_bwamem
    file index from bwa_index_bwamem.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mapped.bam"), file("*.{bai,csi}") into ch_output_from_bwamem

    when: 
    params.mapper == 'bwamem'

    script:
    fasta = "${index}/${bwa_base}"
    size = "${params.large_ref}" ? '-c' : ''

    if (!params.single_end && params.skip_collapse){
    """
    bwa mem -t ${task.cpus} $fasta $r1 $r2 -R "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}".mapped.bam
    samtools index "${size}" -@ ${task.cpus} "${libraryid}".mapped.bam
    """
    } else {
    """
    bwa mem -t ${task.cpus} $fasta $r1 -R "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}".mapped.bam
    samtools index "${size}" -@ ${task.cpus} "${libraryid}".mapped.bam
    """
    }
    
}

// CircularMapper reference preparation and mapping for circular genomes e.g. mtDNA

process circulargenerator{
    label 'sc_tiny'
    tag "$prefix"
    publishDir "${params.outdir}/reference_genome/circularmapper_index", mode: 'copy', saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    file fasta from ch_fasta_for_circulargenerator

    output:
    file "${prefix}.{amb,ann,bwt,sa,pac}" into ch_circularmapper_indices

    when: 
    params.mapper == 'circularmapper'

    script:
    prefix = "${fasta.baseName}_${params.circularextension}.fasta"
    """
    circulargenerator -e ${params.circularextension} -i $fasta -s ${params.circulartarget}
    bwa index $prefix
    """

}

process circularmapper{
    label 'mc_medium'
    tag "$libraryid"
    publishDir "${params.outdir}/mapping/circularmapper", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_cm
    file index from ch_circularmapper_indices.collect()
    file fasta from ch_fasta_for_circularmapper.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mapped.bam"), file("*.{bai,csi}") into ch_output_from_cm, ch_outputindex_from_cm

    when: 
    params.mapper == 'circularmapper'

    script:
    filter = "${params.circularfilter}" ? '' : '-f true -x false'
    elongated_root = "${fasta.baseName}_${params.circularextension}.fasta"

    size = "${params.large_ref}" ? '-c' : ''

    if (!params.single_end && params.skip_collapse ){
    """
    bwa aln -t ${task.cpus} $elongated_root $r1 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r1.sai
    bwa aln -t ${task.cpus} $elongated_root $r2 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $elongated_root ${libraryid}.r1.sai ${libraryid}.r2.sai $r1 $r2 > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > ${libraryid}.mapped.bam
    samtools index "${size}" ${libraryid}.mapped.bam
    """
    } else {
    """ 
    bwa aln -t ${task.cpus} $elongated_root $r1 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $elongated_root ${libraryid}.sai $r1 > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > "${libraryid}".mapped.bam
    samtools index "${size}" "${libraryid}".mapped.bam
    """
    }
    
}

process bowtie2 {
    label 'mc_medium'
    tag "${libraryid}"
    publishDir "${params.outdir}/mapping/bt2", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_bt2
    file index from bt2_index.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mapped.bam"), file("*.{bai,csi}") into ch_output_from_bt2
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*_bt2.log") into ch_bt2_for_multiqc

    when: 
    params.mapper == 'bowtie2'

    script:
    size = "${params.large_ref}" ? '-c' : ''
    fasta = "${index}/${bt2_base}"
    trim5 = "${params.bt2_trim5}" != 0 ? "--trim5 ${params.bt2_trim5}" : ""
    trim3 = "${params.bt2_trim3}" != 0 ? "--trim3 ${params.bt2_trim3}" : ""
    bt2n = "${params.bt2n}" != 0 ? "-N ${params.bt2n}" : ""
    bt2l = "${params.bt2l}" != 0 ? "-L ${params.bt2l}" : ""

    if ( "${params.bt2_alignmode}" == "end-to-end"  ) {
      switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast"; break
        case "fast":
        sensitivity = "--fast"; break
        case "sensitive":
        sensitivity = "--sensitive"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive"; break
        default:
        sensitivity = ""; break
        }
      } else if ("${params.bt2_alignmode}" == "local") {
      switch ( "${params.bt2_sensitivity}" ) {
        case "no-preset":
        sensitivity = ""; break
        case "very-fast":
        sensitivity = "--very-fast-local"; break
        case "fast":
        sensitivity = "--fast-local"; break
        case "sensitive":
        sensitivity = "--sensitive-local"; break
        case "very-sensitive":
        sensitivity = "--very-sensitive-local"; break
        default:
        sensitivity = ""; break

        }
      }

    //PE data without merging, PE data without any AR applied
    if ( seqtype == 'PE' && ( params.skip_collapse || params.skip_adapterremoval ) ){
    """
    bowtie2 -x ${fasta} -1 ${r1} -2 ${r2} -p ${task.cpus} ${sensitivity} ${bt2n} ${bt2l} ${trim5} ${trim3} 2> "${libraryid}"_bt2.log | samtools sort -@ ${task.cpus} -O bam > "${libraryid}".mapped.bam
    samtools index "${size}" "${libraryid}".mapped.bam
    """
    } else {
    //PE collapsed, or SE data 
    """
    bowtie2 -x ${fasta} -U ${r1} -p ${task.cpus} ${sensitivity} ${bt2n} ${bt2l} ${trim5} ${trim3} 2> "${libraryid}"_bt2.log | samtools sort -@ ${task.cpus} -O bam > "${libraryid}".mapped.bam
    samtools index "${size}" "${libraryid}".mapped.bam
    """
    }
    
}

// Gather all mapped BAMs from all possible mappers into common channels to send downstream
ch_output_from_bwa.mix(ch_output_from_bwamem, ch_output_from_cm, ch_indexbam_for_filtering, ch_output_from_bt2)
  .into { ch_mapping_for_skipfiltering; ch_mapping_for_filtering;  ch_mapping_for_samtools_flagstat }

// Post-mapping QC

process samtools_flagstat {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_mapping_for_samtools_flagstat

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*stats") into ch_flagstat_for_multiqc,ch_flagstat_for_endorspy

    script:
    """
    samtools flagstat $bam > ${libraryid}_flagstat.stats
    """
}

// BAM filtering e.g. to extract unmapped reads for downstream or stricter mapping quality

process samtools_filter {
    label 'mc_medium'
    tag "$libraryid"
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
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_mapping_for_filtering

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*filtered.bam"), file("*.{bai,csi}") into ch_output_from_filtering,ch_outputindex_from_filtering
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.unmapped.fastq.gz") optional true into ch_bam_filtering_for_metagenomic
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.unmapped.bam") optional true

    script:
    size = "${params.large_ref}" ? '-c' : ''
    
    if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "discard"){
        """
        samtools view -h -b $bam -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.filtered.bam
        samtools index "${size}" ${libraryid}.filtered.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "bam"){
        """
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.unmapped.bam
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.filtered.bam
        samtools index "${size}" ${libraryid}.filtered.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "fastq"){
        """
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.unmapped.bam
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.filtered.bam
        samtools index "${size}" ${libraryid}.filtered.bam
        samtools fastq -tn ${libraryid}.unmapped.bam | pigz -p ${task.cpus} > ${libraryid}.unmapped.fastq.gz
        rm ${libraryid}.unmapped.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "both"){
        """
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.unmapped.bam
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.filtered.bam
        samtools index "${size}" ${libraryid}.filtered.bam
        samtools fastq -tn ${libraryid}.unmapped.bam | pigz -p ${task.cpus} > ${libraryid}.unmapped.fastq.gz
        """
    } else { //Only apply quality filtering, default
        """
        samtools view -h -b $bam -@ ${task.cpus} -q ${params.bam_mapping_quality_threshold} -o ${libraryid}.filtered.bam
        samtools index "${size}" ${libraryid}.filtered.bam
        """
    }  
}

// samtools_filter bypass in case not run
if (params.run_bam_filtering) {
    ch_mapping_for_skipfiltering.mix(ch_output_from_filtering)
        .filter { it =~/.*filtered.bam/ }
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_stripfastq; ch_filtering_for_flagstat } 

} else {
    ch_mapping_for_skipfiltering
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_stripfastq; ch_filtering_for_flagstat } 

}

// Synchronise the mapped input FASTQ and input non-remapped BAM channels
ch_fastqlanemerge_for_stripfastq
    .map {
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def reads = arrayify(it[7])
        def r1 = it[7].getClass() == ArrayList ? reads[0] : it[7]
        def r2 = it[7].getClass() == ArrayList ? reads[1] : "NA"      

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

    }
    .mix(ch_filtering_for_stripfastq)
    .groupTuple(by: [0,1,3,4,5,6])
    .map {
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = it[7][0]
        def r2 = it[8][0]
        def bam = it[7][1]
        def bai = it[8][1]

       [ samplename, libraryid, seqtype, organism, strandedness, udg, r1, r2, bam, bai ]

    }
    .filter{ it[8] != null }
    .set { ch_synced_for_stripfastq }

// Remove mapped reads from original (lane merged) input FASTQ e.g. for sensitive host data when running metagenomic data

process strip_input_fastq {
    label 'mc_medium'
    tag "${libraryid}"
    publishDir "${params.outdir}/stripped_fastq", mode: 'copy'

    when: 
    params.strip_input_fastq

    input: 
    tuple samplename, libraryid, seqtype, organism, strandedness, udg, file(r1), file(r2), file(bam), file(bai) from ch_synced_for_stripfastq

    output:
    tuple samplename, libraryid, seqtype, organism, strandedness, udg, file("*.fq.gz") into ch_output_from_stripfastq

    script:
    if ( seqtype == 'SE' ) {
        out_fwd = bam.baseName+'.stripped.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${r1} -m ${params.strip_mode} -of $out_fwd -p ${task.cpus}
        """
    } else {
        out_fwd = bam.baseName+'.stripped.fwd.fq.gz'
        out_rev = bam.baseName+'.stripped.rev.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${r1} -rev ${r2} -m  ${params.strip_mode} -of $out_fwd -or $out_rev -p ${task.cpus}
        """ 
    }
    
}

// Post filtering mapping QC - particularly to help see how much was removed from mapping quality filtering

process samtools_flagstat_after_filter {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/samtools/filtered_stats", mode: 'copy'

    when:
    params.run_bam_filtering

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_filtering_for_flagstat

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.stats") into ch_bam_filtered_flagstat_for_multiqc, ch_bam_filtered_flagstat_for_endorspy

    script:
    """
    samtools flagstat $bam > ${libraryid}_postfilterflagstat.stats
    """
}

if (params.run_bam_filtering) {
  ch_flagstat_for_endorspy
    .join(ch_bam_filtered_flagstat_for_endorspy, by: [0,1,2,3,4,5,6])
    .set{ ch_allflagstats_for_endorspy }

} else {
  // Add a file entry to match expected no. tuple elements for endorS.py even if not giving second file
  ch_flagstat_for_endorspy
    .map { it -> 
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]     
        def stats = file(it[7])
        def poststats = file("$baseDir/assets/dummy_postfilterflagstat.stats")

      [samplename, libraryid, lane, seqtype, organism, strandedness, udg, stats, poststats ] }
    .set{ ch_allflagstats_for_endorspy }
}

// Endogenous DNA calculator to say how much of a library contained 'on-target' DNA

process endorSpy {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/endorSpy", mode: 'copy'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(stats), file(poststats) from ch_allflagstats_for_endorspy

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.json") into ch_endorspy_for_multiqc

    script:

    if (params.run_bam_filtering) {
      """
      endorS.py -o json -n ${libraryid} ${stats} ${poststats}
      """
    } else {
      """
      endorS.py -o json -n ${libraryid} ${stats}
      """
    }
}

// Post-mapping PCR amplicon removal because these lab artefacts inflate coverage statistics

process dedup{
    label 'mc_small'
    tag "${libraryid}"
    publishDir "${params.outdir}/deduplication/", mode: 'copy',
        saveAs: {filename -> "${libraryid}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_filtering_for_dedup

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.hist") into ch_hist_for_preseq
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.json") into ch_dedup_results_for_multiqc
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${libraryid}_rmdup.bam"), file("*.{bai,csi}") into ch_output_from_dedup

    script:
    outname = "${bam.baseName}"
    treat_merged="${params.dedup_all_merged}" ? '-m' : ''
    size = "${params.large_ref}" ? '-c' : ''
    
    if(seqtype == 'SE') {
    """
    ## To make sure direct BAMs have a clean name
    if [[ "${bam}" != "${libraryid}.bam" ]]; then
      mv ${bam} ${libraryid}.bam
    fi
    
    dedup -i ${libraryid}.bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "${libraryid}"_rmdup.bam -o "${libraryid}"_rmdup.bam
    samtools index "${size}" "${libraryid}"_rmdup.bam
    """  
    } else {
    """
    ## To make sure direct BAMs have a clean name
    if [[ "${bam}" != "${libraryid}.bam" ]]; then 
      mv ${bam} ${libraryid}.bam
    fi
    
    dedup -i ${libraryid}.bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "${libraryid}"_rmdup.bam -o "${libraryid}"_rmdup.bam
    samtools index "${size}" "${libraryid}"_rmdup.bam
    """  
    }
}

process markDup{
    label 'mc_small'
    tag "${outname}"
    publishDir "${params.outdir}/deduplication/", mode: 'copy',
        saveAs: {filename -> "${libraryid}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'markduplicates'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_filtering_for_markdup

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.metrics") into ch_markdup_results_for_multiqc
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${libraryid}_rmdup.bam"), file("*.{bai,csi}") into ch_output_from_markdup

    script:
    outname = "${bam.baseName}"
    size = "${params.large_ref}" ? '-c' : ''
    """
    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M MarkDuplicates INPUT=$bam OUTPUT=${libraryid}_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="${libraryid}_rmdup.metrics" VALIDATION_STRINGENCY=SILENT
    samtools index "${size}" ${libraryid}_rmdup.bam
    """
}

// Merge independent libraries sequenced but with same treatment (often done to improve complexity). Different strand/UDG libs not merged because bamtrim/pmdtools needs UDG info

// Step one: work out which are single libraries (from skipping rmdup and both dedups) that do not need merging and pass to a skipping
if ( params.skip_deduplication ) {
  ch_input_for_librarymerging = ch_filtering_for_skiprmdup
    .groupTuple(by:[0,4,5,6])
    .branch{
      clean_libraryid: it[7].size() == 1
      merge_me: it[7].size() > 1
    }
} else {
    ch_input_for_librarymerging = ch_output_from_dedup.mix(ch_output_from_markdup)
    .groupTuple(by:[0,4,5,6])
    .branch{
      clean_libraryid: it[7].size() == 1
      merge_me: it[7].size() > 1
    }
}

// For non-merging libraries, fix group libraryIDs into single values. 
// This is a bit hacky as theoretically could have different, but this should
// rarely be the case.

ch_input_for_librarymerging.clean_libraryid
  .map{
    it ->
      def libraryid = it[1][0]
      [it[0], libraryid, it[2], it[3], it[4], it[5], it[6], it[7], it[8] ]
    }
  .set { ch_input_for_skiplibrarymerging }

process library_merge {
  label 'mc_tiny'
  tag "${samplename}"
  publishDir "${params.outdir}/merged_bams/initial", mode: 'copy'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_input_for_librarymerging.merge_me

  output:
  tuple samplename, val("${samplename}_libmerged"), lane, seqtype, organism, strandedness, udg, file("*_libmerged_rg_rmdup.bam"), file("*_libmerged_rg_rmdup.bam.{bai,csi}") into ch_output_from_librarymerging

  script:
  size = "${params.large_ref}" ? '-c' : ''
  """
  samtools merge ${samplename}_libmerged_rmdup.bam ${bam}
  picard AddOrReplaceReadGroups I=${samplename}_libmerged_rmdup.bam O=${samplename}_libmerged_rg_rmdup.bam RGID=1 RGLB="${samplename}_merged" RGPL=illumina RGPU=4410 RGSM="${samplename}_merged"
  samtools index "${size}" ${samplename}_libmerged_rg_rmdup.bam
  """
}

// Mix back in libraries from skipping dedup, skipping library merging
if (!params.skip_deduplication) {
    ch_input_for_skiplibrarymerging.mix(ch_output_from_librarymerging)
        .filter { it =~/.*_rmdup.bam/ }
        .into { ch_rmdup_for_skipdamagemanipulation; ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_for_sexdeterrmine; ch_for_nuclear_contamination; ch_rmdup_for_bedtools; ch_rmdup_formtnucratio } 

} else {
    ch_input_for_skiplibrarymerging.mix(ch_output_from_librarymerging)
        .into { ch_rmdup_for_skipdamagemanipulation; ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_for_sexdeterrmine; ch_for_nuclear_contamination; ch_rmdup_for_bedtools; ch_rmdup_formtnucratio } 
}

//////////////////////////////////////////////////
/* --     POST DEDUPLICATION EVALUATION      -- */
//////////////////////////////////////////////////

// Library complexity calculation from mapped reads - could a user cost-effectively sequence deeper for more unique information?

process preseq {
    label 'sc_tiny'
    tag "${libraryid}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_preseq

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(input) from (params.skip_deduplication ? ch_rmdup_for_preseq.map{ it[0,1,2,3,4,5,6,7] } : ch_hist_for_preseq )

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${input.baseName}.ccurve") into ch_preseq_for_multiqc

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

// Optional mapping statistics for specific annotations - e.g. genes in bacterial genome

// Set up channels for annotation file
if (!params.run_bedtools_coverage){
  ch_anno_for_bedtools = Channel.empty()
} else {
  Channel
    ch_anno_for_bedtools = Channel.fromPath(params.anno_file)
}

process bedtools {
  label 'mc_small'
  tag "${libraryid}"
  publishDir "${params.outdir}/bedtools", mode: 'copy'

  when:
  params.run_bedtools_coverage

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_rmdup_for_bedtools
  file anno_file from ch_anno_for_bedtools.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*")

  script:
  """
  bedtools coverage -a ${anno_file} -b $bam | pigz -p ${task.cpus} > "${bam.baseName}".breadth.gz
  bedtools coverage -a ${anno_file} -b $bam -mean | pigz -p ${task.cpus} > "${bam.baseName}".depth.gz
  """
}

//////////////////////////////////////////////////////////////
/* --    ANCIENT DNA EVALUATION AND BAM MODIFICATION     -- */
//////////////////////////////////////////////////////////////

// Calculate typical aDNA damage frequency distribution

process damageprofiler {
    label 'sc_small'
    tag "${libraryid}"

    publishDir "${params.outdir}/damageprofiler", mode: 'copy'

    when:
    !params.skip_damage_calculation

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_rmdup_for_damageprofiler
    file fasta from ch_fasta_for_damageprofiler.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${base}/*.txt") optional true
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${base}/*.log")
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${base}/*.pdf") optional true
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("${base}/*.json") optional true into ch_damageprofiler_results

    script:
    base = "${bam.baseName}"
    """
    damageprofiler -Xmx${task.memory.toGiga()}g -i $bam -r $fasta -l ${params.damageprofiler_length} -t ${params.damageprofiler_threshold} -o . -yaxis_damageplot ${params.damageprofiler_yaxis}
    """
}

// Optionally perform further aDNA evaluation or filtering for  just reads with damage etc.

process pmdtools {
    label 'mc_small'
    tag "${libraryid}"
    publishDir "${params.outdir}/pmdtools", mode: 'copy'

    when: params.run_pmdtools

    input: 
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_rmdup_for_pmdtools
    file fasta from ch_fasta_for_pmdtools.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.bam"), file("*.{bai,csi}") into ch_output_from_pmdtools
    file "*.cpg.range*.txt"

    script:
    //Check which treatment for the libraries was used
    def treatment = udg ? (udg == 'half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(params.snpcapture){
        snpcap = (params.pmdtools_reference_mask != '') ? "--refseq ${params.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    size = "${params.large_ref}" ? '-c' : ''
    """
    #Run Filtering step 
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --threshold ${params.pmdtools_threshold} $treatment $snpcap --header | samtools view -@ ${task.cpus} -Sb - > "${libraryid}".pmd.bam
    #Run Calc Range step
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --deamination --range ${params.pmdtools_range} $treatment $snpcap -n ${params.pmdtools_max_reads} > "${libraryid}".cpg.range."${params.pmdtools_range}".txt 
    samtools index "${size}" ${libraryid}.pmd.bam
    """
}

// BAM Trimming for just non-UDG or half-UDG libraries to remove damage prior genotyping

if ( params.run_trim_bam ) {

    // You wouldn't want to make UDG treated reads even shorter, so skip trimming if UDG.
    // We assume same trim amount for both non-UDG/UDG half as could trim a bit more off half-UDG to match non-UDG if needed, with minimal effect 
    // Note: Trimming of e.g. adapters are sequencing artefacts and should be removed before mapping, so we don't account for this here.
    ch_bamutils_decision = ch_rmdup_for_bamutils.branch{
        totrim: it[6] == 'none' || it[6] == 'half' 
        notrim: it[6] == 'full'
    }

} else {

    ch_bamutils_decision = ch_rmdup_for_bamutils.branch{
        totrim: it[6] == "dummy"
        notrim: it[6] == 'full' || it[6] == 'none' || it[6] == 'half'
    }

}

process bam_trim {
    label 'mc_small'
    tag "${libraryid}" 
    publishDir "${params.outdir}/trimmed_bam", mode: 'copy'
 
    when: params.run_trim_bam

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_bamutils_decision.totrim

    output: 
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.trimmed.bam"), file("*.{bai,csi}") into ch_trimmed_from_bamutils

    script:
    softclip = "${params.bamutils_softclip}" ? '-c' : '' 
    size = "${params.large_ref}" ? '-c' : ''
    """
    bam trimBam $bam tmp.bam -L ${params.bamutils_clip_left} -R ${params.bamutils_clip_right} ${softclip}
    samtools sort -@ ${task.cpus} tmp.bam -o ${libraryid}.trimmed.bam 
    samtools index "${size}" ${libraryid}.trimmed.bam
    """
}

// Post trimming merging, because we will presume that if trimming is turned on, 'lab-removed' libraries can be combined with merged with 'in-silico damage removed' libraries to improve genotyping

ch_trimmed_formerge = ch_bamutils_decision.notrim
  .mix(ch_trimmed_from_bamutils)
  .groupTuple(by:[0,4,5])
  .map{
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]     
        def bam = it[7].flatten()
        def bai = it[8].flatten()

      [samplename, libraryid, lane, seqtype, organism, strandedness, udg, bam, bai ]
  }
  .branch{
    skip_merging: it[7].size() == 1
    merge_me: it[7].size() > 1
  }

//////////////////////////////////////////////////////////////////////////
/* --    POST aDNA BAM MODIFICATION AND FINAL MAPPING STATISTICS     -- */
//////////////////////////////////////////////////////////////////////////

process additional_library_merge {
  label 'mc_tiny'
  tag "${samplename}"
  publishDir "${params.outdir}/merged_bams/additional", mode: 'copy'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_trimmed_formerge.merge_me

  output:
  tuple samplename, val("${samplename}_libmerged"), lane, seqtype, organism, strandedness, udg, file("*_libmerged_rg_add.bam"), file("*_libmerged_rg_add.bam.{bai,csi}") into ch_output_from_trimmerge

  script:
  size = "${params.large_ref}" ? '-c' : ''
  """
  samtools merge ${samplename}_libmerged_add.bam ${bam}
  picard AddOrReplaceReadGroups I=${samplename}_libmerged_add.bam O=${samplename}_libmerged_rg_add.bam RGID=1 RGLB="${samplename}_additionalmerged" RGPL=illumina RGPU=4410 RGSM="${samplename}_additionalmerged"
  samtools index "${size}" ${samplename}_libmerged_rg_add.bam
  """
}

ch_trimmed_formerge.skip_merging
  .mix(ch_output_from_trimmerge)
  .into{ ch_output_from_bamutils; ch_addlibmerge_for_qualimap }

  // General mapping quality statistics for whole reference sequence - e.g. X and % coverage

process qualimap {
    label 'mc_small'
    tag "${samplename}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_addlibmerge_for_qualimap
    file fasta from ch_fasta_for_qualimap.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*") into ch_qualimap_results

    script:
    snpcap = ''
    if(params.snpcapture) snpcap = "-gff ${params.bedfile}"
    """
    qualimap bamqc -bam $bam -nt ${task.cpus} -outdir . -outformat "HTML" ${snpcap}
    """
}

/////////////////////////////
/* --    GENOTYPING     -- */
/////////////////////////////

// Reroute files for genotyping; we have to ensure to select lib-merged BAMs, as input channel will also contain the un-merged ones resulting in unwanted multi-sample VCFs
if ( params.run_genotyping && params.genotyping_source == 'raw' ) {
    ch_output_from_bamutils
      .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( params.run_genotyping && params.genotyping_source == "trimmed" && !params.run_trim_bam )  {
    exit 1, "[nf-core/eager] error: Cannot run genotyping with 'trimmed' source without running BAM trimming (--run_trim_bam)! Please check input parameters."

} else if ( params.run_genotyping && params.genotyping_source == "trimmed" && params.run_trim_bam )  {
    ch_output_from_bamutils
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( params.run_genotyping && params.genotyping_source == "pmd" && !params.run_run_pmdtools )  {
    exit 1, "[nf-core/eager] error: Cannot run genotyping with 'pmd' source without running pmtools (--run_pmdtools)! Please check input parameters."

} else if ( params.run_genotyping && params.genotyping_source == "pmd" && params.run_pmdtools )  {
   ch_output_from_pmdtools
     .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && !params.run_trim_bam && !params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
     .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && !params.run_trim_bam && params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
     .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

} else if ( !params.run_genotyping && params.run_trim_bam && !params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
     .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes; ch_damagemanipulation_for_genotyping_pileupcaller; ch_damagemanipulation_for_genotyping_angsd }

}

// Unified Genotyper - although not-supported, better for aDNA (because HC does de novo assembly which requires higher coverages), and needed for MultiVCFAnalyzer

if ( params.gatk_ug_jar != '' ) {
  Channel
    .fromPath( params.gatk_ug_jar )
    .set{ ch_unifiedgenotyper_jar }
} else {
  Channel
    .empty()
    .set{ ch_unifiedgenotyper_jar }
}

 process genotyping_ug {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'ug'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_ug
  file fasta from ch_fasta_for_genotyping_ug.collect()
  file jar from ch_unifiedgenotyper_jar.collect()
  file fai from ch_fai_for_ug.collect()
  file dict from ch_dict_for_ug.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*vcf.gz") into ch_ug_for_multivcfanalyzer,ch_ug_for_vcf2genome
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*realign.bam") optional true

  script:
  defaultbasequalities = params.gatk_ug_defaultbasequalities == '' ? '' : " --defaultBaseQualities ${params.gatk_ug_defaultbasequalities}" 
  keep_realign = params.gatk_ug_keep_realign_bam ? "T" : "F"
  if (params.gatk_dbsnp == '')
    """
    samtools index -b ${bam}
    java -Xmx${task.memory.toGiga()}g -jar ${jar} -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    java -Xmx${task.memory.toGiga()}g -jar ${jar} -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplename}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    java -Xmx${task.memory.toGiga()}g -jar ${jar} -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}
    
    if [[ ${keep_realign} == 'F' ]]; then
      rm ${samplename}.realign.bam
    fi
    
    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
  else if (params.gatk_dbsnp != '')
    """
    samtools index ${bam}
    java -jar ${jar} -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    java -jar ${jar} -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplenane}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    java -jar ${jar} -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --dbsnp ${params.gatk_dbsnp} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}
    
    if [[ ${keep_realign} == 'F' ]]; then
      rm ${samplename}.realign.bam
    fi
    
    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
 }

 // HaplotypeCaller as 'best practise' tool for human DNA in particular 

 process genotyping_hc {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'hc'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_hc
  file fasta from ch_fasta_for_genotyping_hc.collect()
  file fai from ch_fai_for_hc.collect()
  file dict from ch_dict_for_hc.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*vcf.gz")

  script:
  if (params.gatk_dbsnp == '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${samplename}.haplotypecaller.vcf -stand-call-conf ${params.gatk_call_conf} --sample-ploidy ${params.gatk_ploidy} --output-mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${samplename}.haplotypecaller.vcf
    """

  else if (params.gatk_dbsnp != '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${samplename}.haplotypecaller.vcf --dbsnp ${params.gatk_dbsnp} -stand-call-conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} --output_mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${samplename}.haplotypecaller.vcf
    """
 }

 // Freebayes for 'more efficient/simple' and more generic genotyping (vs HC) 

 process genotyping_freebayes {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'freebayes'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_freebayes
  file fasta from ch_fasta_for_genotyping_freebayes.collect()
  file fai from ch_fai_for_freebayes.collect()
  file dict from ch_dict_for_freebayes.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*vcf.gz")
  
  script:
  skip_coverage = "${params.freebayes_g}" == 0 ? "" : "-g ${params.freebayes_g}"
  """
  freebayes -f ${fasta} -p ${params.freebayes_p} -C ${params.freebayes_C} ${skip_coverage} ${bam} > ${samplename}.freebayes.vcf
  pigz -p ${task.cpus} ${samplename}.freebayes.vcf
  """
 }

 // pileupCaller for 'random sampling' genotyping

if (params.pileupcaller_bedfile.isEmpty()) {
  ch_bed_for_pileupcaller = 'NO_FILE_BED'
} else {
  ch_bed_for_pileupcaller = Channel.fromPath(params.pileupcaller_bedfile)
}

if (params.pileupcaller_snpfile.isEmpty ()) {
  ch_snp_for_pileupcaller = 'NO_FILE'
} else {
  ch_snp_for_pileupcaller = Channel.fromPath(params.pileupcaller_snpfile)
}

 process genotyping_pileupcaller {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'pileupcaller'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_pileupcaller
  file fasta from ch_fasta_for_genotyping_pileupcaller.collect()
  file fai from ch_fai_for_pileupcaller.collect()
  file dict from ch_dict_for_pileupcaller.collect()
  file bed from ch_bed_for_pileupcaller.collect()
  file snp from ch_snp_for_pileupcaller.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("pileupcaller.${samplename}.*")

  script:
  caller = "--${params.pileupcaller_method}"
  ssmode = strandedness == "single" ? "--singleStrandMode" : ""
  """
  samtools mpileup -B -q 30 -Q 30 -l ${bed} -f ${fasta} ${bam} | pileupCaller ${caller} ${ssmode} --sampleNames ${samplename} -f ${snp} -e pileupcaller.${samplename}
  """
 }

 process genotyping_angsd {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'angsd'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_angsd
  file fasta from ch_fasta_for_genotyping_angsd.collect()
  file fai from ch_fai_for_angsd.collect()
  file dict from ch_dict_for_angsd.collect()

  output: 
  path("${samplename}*")
  
  script:
  switch ( "${params.angsd_glmodel}" ) {
    case "samtools":
    angsd_glmodel = "1"; break
    case "gatk":
    angsd_glmodel = "2"; break
    case "soapsnp":
    angsd_glmodel = "3"; break
    case "syk":
    angsd_glmodel = "4"; break
  }

  switch ( "${params.angsd_glformat}" ) {
    case "text":
    angsd_glformat = "4"; break
    case "binary":
    angsd_glformat = "1"; break
    case "beagle":
    angsd_glformat = "2"; break
    case "binary_three":
    angsd_glformat = "3"; break
  }
  
  angsd_fasta = !params.angsd_createfasta ? '' : params.angsd_fastamethod == 'random' ? '-doFasta 1 -doCounts 1' : '-doFasta 2 -doCounts 1' 
  """
  echo ${bam} > bam.filelist
  mkdir angsd
  angsd -bam bam.filelist -nThreads ${task.cpus} -GL ${angsd_glmodel} -doGlF ${angsd_glformat} ${angsd_fasta} -out ${samplename}.angsd
  """
 }

////////////////////////////////////
/* --    CONSENSUS CALLING     -- */
////////////////////////////////////

// Generate a simple consensus-called FASTA file based on genotype VCF

process vcf2genome {
  label  'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/consensus_sequence", mode: 'copy'

  when: 
  params.run_vcf2genome

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(vcf) from ch_ug_for_vcf2genome
  file fasta from ch_fasta_for_vcf2genome.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.fasta.gz")

  script:
  out = "${params.vcf2genome_outfile}" == '' ? "${samplename}.fasta" : "${params.vcf2genome_outfile}"
  fasta_head = "${params.vcf2genome_header}" == '' ? "${samplename}" : "${params.vcf2genome_header}"
  """
  pigz -f -d -p ${task.cpus} *.vcf.gz
  vcf2genome -draft ${out}.fasta -draftname "${fasta_head}" -in ${vcf.baseName} -minc ${params.vcf2genome_minc} -minfreq ${params.vcf2genome_minfreq} -minq ${params.vcf2genome_minq} -ref ${fasta} -refMod ${out}_refmod.fasta -uncertain ${out}_uncertainy.fasta
  pigz -p ${task.cpus} *.fasta 
  pigz -p ${task.cpus} *.vcf
  """
}

// More complex consensus caller with additional filtering functionality (e.g. for heterozygous calls) to generate SNP tables and other things sometimes used in aDNA bacteria studies

// Create input channel for MultiVCFAnalyzer, possibly mixing with pre-made VCFs
if (params.additional_vcf_files == '') {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it[7] }.collect()
} else {
    ch_extravcfs_for_multivcfanalyzer = Channel.fromPath(params.additional_vcf_files)
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it [7] }.collect().mix(ch_extravcfs_for_multivcfanalyzer)
}

 process multivcfanalyzer {
  label  'mc_small'
  publishDir "${params.outdir}/MultiVCFAnalyzer", mode: 'copy'

  when:
  params.genotyping_tool == 'ug' && params.run_multivcfanalyzer && params.gatk_ploidy == '2'

  input:
  file vcf from ch_vcfs_for_multivcfanalyzer.collect()
  file fasta from ch_fasta_for_multivcfanalyzer.collect()

  output:
  file('fullAlignment.fasta.gz') into ch_output_multivcfanalyzer_fullalignment
  file('info.txt.gz') into ch_output_multivcfanalyzer_info
  file('snpAlignment.fasta.gz') into ch_output_multivcfanalyzer_snpalignment
  file('snpAlignmentIncludingRefGenome.fasta.gz') into ch_output_multivcfanalyzer_snpalignmentref
  file('snpStatistics.tsv.gz') into ch_output_multivcfanalyzer_snpstatistics
  file('snpTable.tsv.gz') into ch_output_multivcfanalyzer_snptable
  file('snpTableForSnpEff.tsv.gz') into ch_output_multivcfanalyzer_snptablesnpeff
  file('snpTableWithUncertaintyCalls.tsv.gz') into ch_output_multivcfanalyzer_snptableuncertainty
  file('structureGenotypes.tsv.gz') into ch_output_multivcfanalyzer_structuregenotypes
  file('structureGenotypes_noMissingData-Columns.tsv.gz') into ch_output_multivcfanalyzer_structuregenotypesclean
  file('MultiVCFAnalyzer.json') optional true into ch_multivcfanalyzer_for_multiqc

  script:
  write_freqs = "$params.write_allele_frequencies" ? "T" : "F"
  """
  gunzip -f *.vcf.gz
  multivcfanalyzer ${params.snp_eff_results} ${fasta} ${params.reference_gff_annotations} . ${write_freqs} ${params.min_genotype_quality} ${params.min_base_coverage} ${params.min_allele_freq_hom} ${params.min_allele_freq_het} ${params.reference_gff_exclude} *.vcf
  pigz -p ${task.cpus} *.tsv *.txt snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta fullAlignment.fasta
  """
 }

////////////////////////////////////////////////////////////
/* --    HUMAN DNA SPECIFIC ADDITIONAL INFORMATION     -- */
////////////////////////////////////////////////////////////

// Mitochondrial to nuclear ratio helps to evaluate quality of tissue sampled

 process mtnucratio {
  tag "${samplename}"
  publishDir "${params.outdir}/mtnucratio", mode: "copy"

  when: 
  params.run_mtnucratio

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_rmdup_formtnucratio

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mtnucratio")
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.json") into ch_mtnucratio_for_multiqc

  script:
  """
  mtnucratio ${bam} "${params.mtnucratio_header}"
  """
 }

// Human biological sex estimation

if (params.sexdeterrmine_bedfile == '') {
  ch_bed_for_sexdeterrmine = file('NO_FILE')
} else {
  ch_bed_for_sexdeterrmine = Channel.fromPath(params.sexdeterrmine_bedfile)
}

// As we collect all files for a single sex_deterrmine run, we DO NOT use the normal input/output tuple
process sex_deterrmine {
    label 'sc_small'
    publishDir "${params.outdir}/sex_determination", mode:"copy"
     
    input:
    file bam from ch_for_sexdeterrmine.map { it[7] }.collect()
    file bed from ch_bed_for_sexdeterrmine

    output:
    file "SexDet.txt"
    file "*.json" into ch_sexdet_for_multiqc

    when:
    params.run_sexdeterrmine
    
    script:
    def filter = bed.name != 'NO_FILE' ? "-b $bed" : ''
    """
    
    for i in *.bam; do
        echo \$i >> bamlist.txt
    done
  
    samtools depth -aa -q30 -Q30 $filter -f bamlist.txt | sexdeterrmine -f bamlist.txt > SexDet.txt
    """
}

// Human DNA nuclear contamination estimation

 process nuclear_contamination{
    label 'sc_small'
    tag "${samplename}"
    publishDir "${params.outdir}/nuclear_contamination", mode:"copy"

    // ANGSD Xcontamination will exit with status 134 when the number of SNPs is too low
    validExitStatus 0,134

    when:
    params.run_nuclear_contamination

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(input), file(bai) from ch_for_nuclear_contamination

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file('*.X.contamination.out') into ch_from_nuclear_contamination

    script:
    """
    samtools index ${input}
    angsd -i ${input} -r ${params.contamination_chrom_name}:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${input.baseName}.doCounts
    contamination -a ${input.baseName}.doCounts.icnts.gz -h ${baseDir}/assets/angsd_resources/HapMapChrX.gz 2> ${input.baseName}.X.contamination.out
    """
 }
 
// As we collect all files for a single print_nuclear_contamination run, we DO NOT use the normal input/output tuple
process print_nuclear_contamination{
    label 'sc_tiny'
    publishDir "${params.outdir}/nuclear_contamination", mode:"copy"

    when:
    params.run_nuclear_contamination

    input:
    val 'Contam' from ch_from_nuclear_contamination.map { it[7] }.collect()

    output:
    file 'nuclear_contamination.txt'

    script:
    """
    print_x_contamination.py ${Contam.join(' ')}
    """
 }

/////////////////////////////////////////////////////////
/* --    METAGENOMICS-SPECIFIC ADDITIONAL STEPS     -- */
/////////////////////////////////////////////////////////

// MALT is a super-fast BLAST replacement typically used for pathogen detection or microbiome profiling against large databases, here using off-target reads from mapping

// As we collect all files for a all metagenomic runs, we DO NOT use the normal input/output tuple!
if (params.metagenomic_tool == 'malt') {
  ch_bam_filtering_for_metagenomic
    .set {ch_bam_filtering_for_metagenomic_malt}

  ch_bam_filtering_for_metagenomic_kraken = Channel.empty()
} else if (params.metagenomic_tool == 'kraken') {
  ch_bam_filtering_for_metagenomic
    .set {ch_bam_filtering_for_metagenomic_kraken}

  ch_bam_filtering_for_metagenomic_malt = Channel.empty()
}

// As we collect all files for a single MALT run, we DO NOT use the normal input/output tuple
process malt {
  label 'mc_small'
  publishDir "${params.outdir}/metagenomic_classification/malt", mode:"copy"

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_discard_unmapped && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'malt'

  input:
  file fastqs from ch_bam_filtering_for_metagenomic_malt.map { it[7] }.collect()

  output:
  file "*.rma6" into ch_rma_for_maltExtract
  file "malt.log" into ch_malt_for_multiqc

  script:
  if ("${params.malt_min_support_mode}" == "percent") {
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  -d ${params.database} \
  -id ${params.percent_identity} \
  -m ${params.malt_mode} \
  -at ${params.malt_alignment_mode} \
  -top ${params.malt_top_percent} \
  -supp ${params.malt_min_support_percent} \
  -mq ${params.malt_max_queries} \
  --memoryMode ${params.malt_memory_mode} \
  -i ${fastqs.join(' ')} |&tee malt.log
  """
  } else if ("${params.malt_min_support_mode}" == "reads") {
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  -d ${params.database} \
  -id ${params.percent_identity} \
  -m ${params.malt_mode} \
  -at ${params.malt_alignment_mode} \
  -top ${params.malt_top_percent} \
  -sup ${params.metagenomic_min_support_reads} \
  -mq ${params.malt_max_queries} \
  --memoryMode ${params.malt_memory_mode} \
  -i ${fastqs.join(' ')} |&tee malt.log
  """
  }

}

// Create input channel for MaltExtract taxon list, to allow downloading of taxon list
if (params.maltextract_taxon_list== '') {
    ch_taxonlist_for_maltextract = Channel.empty()
} else {
    ch_taxonlist_for_maltextract = Channel.fromPath(params.maltextract_taxon_list)
}

// MaltExtract performs aDNA evaluation from the output of MALT (damage patterns, read lengths etc.)

// As we collect all files for a single MALT extract run, we DO NOT use the normal input/output tuple
process maltextract {
  label 'mc_large'
  publishDir "${params.outdir}/MaltExtract/", mode:"copy"

  when: 
  params.run_maltextract && params.metagenomic_tool == 'malt'

  input:
  file rma6 from ch_rma_for_maltExtract.collect()
  file taxon_list from ch_taxonlist_for_maltextract
  
  output:
  path "results/" type('dir')
  file "results/*_Wevid.json" optional true into ch_hops_for_multiqc 

  script:
  ncbifiles = params.maltextract_ncbifiles == '' ? "" : "-r ${params.maltextract_ncbifiles}"
  destack = params.maltextract_destackingoff ? "--destackingOff" : ""
  downsam = params.maltextract_downsamplingoff ? "--downSampOff" : ""
  dupremo = params.maltextract_duplicateremovaloff ? "--dupRemOff" : ""
  matches = params.maltextract_matches ? "--matches" : ""
  megsum = params.maltextract_megansummary ? "--meganSummary" : ""
  topaln = params.maltextract_topalignment ?  "--useTopAlignment" : ""
  ss = params.single_stranded ? "--singleStranded" : ""
  """
  MaltExtract \
  -Xmx${task.memory.toGiga()}g \
  -t ${taxon_list} \
  -i ${rma6.join(' ')} \
  -o results/ \
  ${ncbifiles} \
  -p ${task.cpus} \
  -f ${params.maltextract_filter} \
  -a ${params.maltextract_toppercent} \
  --minPI ${params.maltextract_percentidentity} \
  ${destack} \
  ${downsam} \
  ${dupremo} \
  ${matches} \
  ${megsum} \
  ${topaln} \
  ${ss}

  postprocessing.AMPS.r -r results/ -m ${params.maltextract_filter} -t ${task.cpus} -n ${taxon_list} -j
  """
}

// Kraken is offered as a replacement for MALT as MALT is _very_ resource hungry

if (params.run_metagenomic_screening && params.database.endsWith(".tar.gz") && params.metagenomic_tool == 'kraken'){
  comp_kraken = file(params.database)

  process decomp_kraken {
    input:
    file(ckdb) from comp_kraken
    
    output:
    file(dbname) into ch_krakendb
    
    script:
    dbname = params.database.tokenize("/")[-1].tokenize(".")[0]
    """
    tar xvzf $ckdb
    """
  }

} else if (! params.database.endsWith(".tar.gz") && params.run_metagenomic_screening && params.metagenomic_tool == 'kraken') {
    ch_krakendb = file(params.database)
} else {
    ch_krakendb = Channel.empty()
}

process kraken {
  tag "$prefix"
  label 'mc_huge'
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode:"copy"

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_discard_unmapped && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'kraken'

  input:
  file(fastq) from ch_bam_filtering_for_metagenomic_kraken.map { it[7] }
  file(krakendb) from ch_krakendb

  output:
  file "*.kraken.out" into ch_kraken_out
  tuple prefix, file("*.kreport") into ch_kraken_report, ch_kraken_for_multiqc

  script:
  prefix = fastq.toString().tokenize('.')[0]
  out = prefix+".kraken.out"
  kreport = prefix+".kreport"

  """
  kraken2 --db ${krakendb} --threads ${task.cpus} --output $out --report $kreport $fastq
  """
}

process kraken_parse {
  tag "$name"
  errorStrategy 'ignore'

  input:
  tuple val(name), file(kraken_r) from ch_kraken_report

  output:
  tuple val(name), file('*.kraken_parsed.csv') into ch_kraken_parsed

  script:
  out = name+".kraken_parsed.csv"
  """
  kraken_parse.py -c ${params.metagenomic_min_support_reads} -o $out $kraken_r
  """    
}

process kraken_merge {
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode:"copy"

  input:
  file csv_count from ch_kraken_parsed.map{ it[1] }.collect().dump()

  output:
  file('kraken_count_table.csv')

  script:
  out = "kraken_count_table.csv"
  """
  merge_kraken_res.py -o $out
  """    
}

//////////////////////////////////////
/* --    PIPELINE COMPLETION     -- */
//////////////////////////////////////

// Pipeline documentation for on-server guidance

process output_documentation {
    label 'sc_tiny'
    publishDir "${params.outdir}/Documentation", mode: 'copy'

    input:
    file output_docs from ch_output_docs
    file images from ch_output_docs_images

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}

// Collect all software versions for inclusion in MultiQC report

process get_software_versions {
  label 'sc_tiny'
  publishDir "${params.outdir}/SoftwareVersions", mode: 'copy'

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version &> v_pipeline.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt 2>&1 || true
    AdapterRemoval --version  &> v_adapterremoval.txt 2>&1 || true
    fastp --version &> v_fastp.txt 2>&1 || true
    bwa &> v_bwa.txt 2>&1 || true
    circulargenerator --help | head -n 1 &> v_circulargenerator.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    dedup -v &> v_dedup.txt 2>&1 || true
    ## bioconda recipe of picard is incorrectly set up and extra warning made with stderr, this ugly command ensures only version exported
    ( exec 7>&1; picard MarkDuplicates --version 2>&1 >&7 | grep -v '/' >&2 ) 2> v_markduplicates.txt || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    preseq &> v_preseq.txt 2>&1 || true
    gatk --version 2>&1 | head -n 1 > v_gatk.txt 2>&1 || true
    freebayes --version &> v_freebayes.txt 2>&1 || true
    bedtools --version &> v_bedtools.txt 2>&1 || true
    damageprofiler --version &> v_damageprofiler.txt 2>&1 || true
    bam --version &> v_bamutil.txt 2>&1 || true
    pmdtools --version &> v_pmdtools.txt 2>&1 || true
    angsd -h |& head -n 1 | cut -d ' ' -f3-4 &> v_angsd.txt 2>&1 || true 
    multivcfanalyzer --help | head -n 1 &> v_multivcfanalyzer.txt 2>&1 || true
    malt-run --help |& tail -n 3 | head -n 1 | cut -f 2 -d'(' | cut -f 1 -d ',' &> v_malt.txt 2>&1 || true
    MaltExtract --help | head -n 2 | tail -n 1 &> v_maltextract.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    vcf2genome -h |& head -n 1 &> v_vcf2genome.txt || true
    mtnucratio --help &> v_mtnucratiocalculator.txt || true
    sexdeterrmine --version &> v_sexdeterrmine.txt || true
    kraken2 --version | head -n 1 &> v_kraken.txt || true
    endorS.py --version &> v_endorSpy.txt || true
    pileupCaller --version &> v_sequencetools.txt 2>&1 || true
    bowtie2 --version | grep -a 'bowtie2-.* -fdebug' > v_bowtie2.txt || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

// MultiQC file generation for pipeline report

process multiqc {
    label 'sc_small'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file ('fastqc_raw/*') from ch_prefastqc_for_multiqc.collect().ifEmpty([])
    file('fastqc/*') from ch_fastqc_after_clipping.collect().ifEmpty([])
    file software_versions_mqc from software_versions_yaml.collect().ifEmpty([])
    file ('adapter_removal/*') from ch_adapterremoval_logs.collect().ifEmpty([])
    file ('mapping/bt2/*') from ch_bt2_for_multiqc.collect().ifEmpty([])
    file ('flagstat/*') from ch_flagstat_for_multiqc.collect().ifEmpty([])
    file ('flagstat_filtered/*') from ch_bam_filtered_flagstat_for_multiqc.collect().ifEmpty([])
    file ('preseq/*') from ch_preseq_for_multiqc.collect().ifEmpty([])
    file ('damageprofiler/dmgprof*/*') from ch_damageprofiler_results.collect().ifEmpty([])
    file ('qualimap/qualimap*/*') from ch_qualimap_results.collect().ifEmpty([])
    file ('markdup/*') from ch_markdup_results_for_multiqc.collect().ifEmpty([])
    file ('dedup*/*') from ch_dedup_results_for_multiqc.collect().ifEmpty([])
    file ('fastp/*') from ch_fastp_for_multiqc.collect().ifEmpty([])
    file ('sexdeterrmine/*') from ch_sexdet_for_multiqc.collect().ifEmpty([])
    file ('mutnucratio/*') from ch_mtnucratio_for_multiqc.collect().ifEmpty([])
    file ('endorspy/*') from ch_endorspy_for_multiqc.collect().ifEmpty([])
    file ('multivcfanalyzer/*') from ch_multivcfanalyzer_for_multiqc.collect().ifEmpty([])
    file ('malt/*') from ch_malt_for_multiqc.collect().ifEmpty([])
    file ('kraken/*') from ch_kraken_for_multiqc.collect().ifEmpty([])
    file ('hops/*') from ch_hops_for_multiqc.collect().ifEmpty([])
    file logo from ch_eager_logo

    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $multiqc_config $custom_config_file .
    """
}

// Send completion emails if requested, so user knows data is ready

workflow.onComplete {

    // On success try attach the multiqc report to a variable
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report.getClass() == ArrayList) {
                log.warn "[nf-core/eager] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
    }

    //TODO ALEX P: Not sure we need this? and for what?
    // Check if we are only sending emails on failure
    //email_address = params.email
    //if (!params.email && params.email_on_fail && !workflow.success) {
    //    email_address = params.email_on_fail
    //}

    CommonFunctions.workflowReport(summary, workflow, params, baseDir, mqc_report)
}

// Channelling the TSV file containing FASTQ or BAM 
// Header Format is: "Sample_Name  Library_ID  Lane  SeqType  Organism  Strandedness  UDG_Treatment  R1  R2  BAM  BAM_Index Group  Populations  Age"
def extract_data(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 11)
            def samplename = row.Sample_Name
            def libraryid  = row.Library_ID
            def lane = row.Lane
            def colour = row.Colour_Chemistry
            def seqtype = row.SeqType
            def organism = row.Organism
            def strandedness = row.Strandedness
            def udg = row.UDG_Treatment
            def r1 = row.R1.matches('NA') ? 'NA' : return_file(row.R1)
            def r2 = row.R2.matches('NA') ? 'NA' : return_file(row.R2)
            def bam = row.BAM.matches('NA') ? 'NA' : return_file(row.BAM)

            // Check no 'empty' rows
            if (r1.matches('NA') && r2.matches('NA') && bam.matches('NA') && bai.matches('NA')) exit 1, "[nf-core/eager] error: A row in your TSV appears to have all files defined as NA. See --help or documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Ensure BAMs aren't submitted with PE
            if (!bam.matches('NA') && seqtype.matches('PE')) exit 1, "[nf-core/eager] error: BAM input rows in TSV cannot be set as PE, only SE. See --help or documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Check valid colour chemistry
            if (!colour == 2 && !colour == 4) exit 1, "[nf-core/eager] error: Colour chemistry in TSV can either be 2 (e.g. NextSeq/NovaSeq) or 4 (e.g. HiSeq/MiSeq)"

            //  Ensure that we do not accept incompatible chemistry setup
            if (!seqtype.matches('PE') && !seqtype.matches('SE')) exit 1, "[nf-core/eager] error:  SeqType for one or more rows in TSV is neither SE nor PE! see --help or documentation under 'running the pipeline' for more information. You have: ${seqtype}"
                   
           // So we don't accept existing files that are wrong format: e.g. fasta or sam
            if ( !r1.matches('NA') && !has_extension(r1, "fastq.gz") && !has_extension(r1, "fq.gz") && !has_extension(r1, "fastq") && !has_extension(r1, "fq")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See --help or documentation under 'running the pipeline' for more information. Check: ${r1}"
            if ( !r2.matches('NA') && !has_extension(r2, "fastq.gz") && !has_extension(r2, "fq.gz") && !has_extension(r2, "fastq") && !has_extension(r2, "fq")) exit 1, "[nf-core/eager] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See --help or documentation under 'running the pipeline' for more information. Check: ${r2}"
            if ( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable BAM extension or is not NA. See --help or documentation under 'running the pipeline' for more information. Check: ${bam}"
             
            [ samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ]

         }

    }

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core/eager] error:  Invalid TSV input - malformed row (e.g. missing column) in ${row}, see --help or documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core/eager] error: Cannot find supplied FASTQ or BAM input file. If using input method TSV set to NA if no file required. See --help or documentation under 'running the pipeline' for more information. Check file: ${it}" 
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// To convert a string to an array when not an array already
// From: https://stackoverflow.com/a/55453674/11502856
def arrayify(it) {
  [] + it ?: [it]
}

// Extract FastQs from Path
// Create a channel of FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and/or _R2_ in their names.
def retrieve_input_paths(input, colour_chem, pe_se, ds_ss, udg_treat, bam_in) {

  if ( !bam_in ) {
        if( pe_se ) {
            log.info "Generating single-end FASTQ data TSV"
            Channel
                .fromFilePairs( input, size: 1 )
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty { exit 1, "[nf-core/eager] error:  Your specified FASTQ read files did not end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'. Did you forget --bam?" }
                .map { row -> [ row[0], [ row[1][0] ] ] }
                .ifEmpty { exit 1, "[nf-core/eager] error:  --input was empty - no input files supplied!" }
                .into { ch_reads_for_faketsv; ch_reads_for_validate }

                // Check we don't have any duplicated sample names due to fromFilePairs behaviour of calculating sample name from anything before R1/R2 glob
                ch_reads_for_validate
                  .groupTuple()
                  .map{
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). Did you accidentally supply paired-end data?  See --help or documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" } 
                  }

        } else if (!pe_se ){
            log.info "Generating paired-end FASTQ data TSV"

            Channel
                .fromFilePairs( input )
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty { exit 1, "[nf-core/eager] error: Your specified FASTQ read files did not end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq' " }
                .map { row -> [ row[0], [ row[1][0], row[1][1] ] ] }
                .ifEmpty { exit 1, "[nf-core/eager] error: --input was empty - no input files supplied!" }
                .into { ch_reads_for_faketsv; ch_reads_for_validate }

                // Check we don't have any duplicated sample names due to fromFilePairs behaviour of calculating sample name from anything before R1/R2 glob
                ch_reads_for_validate
                  .groupTuple()
                  .map{
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). See --help or documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" } 
                  }

        } 

    } else if ( bam_in ) {
              log.info "Generating BAM data TSV"

         Channel
            .fromFilePairs( input, size: 1 )
            .filter { it =~/.*.bam/ }
            .map { row -> [ row[0], [ row[1][0] ] ] }
            .ifEmpty { exit 1, "[nf-core/eager] error: Cannot find any bam file matching: ${input}" }
            .set { ch_reads_for_faketsv }

    }

ch_reads_for_faketsv
  .map{

      def samplename = it[0]
      def libraryid  = it[0]
      def lane = 0
      def colour = "${colour_chem}"
      def seqtype = pe_se ? 'SE' : 'PE'
      def organism = 'NA'
      def strandedness = ds_ss ? 'single' : 'double'
      def udg = udg_treat
      def r1 = !bam_in ? return_file(it[1][0]) : 'NA'
      def r2 = !bam_in && !pe_se ? return_file(it[1][1]) : 'NA'
      def bam = bam_in && pe_se ? return_file(it[1][0]) : 'NA'

      [ samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ]
  }
  .ifEmpty {exit 1, "[nf-core/eager] error: Invalid file paths with --input"}

}

// Function to check length of collection in a channel closure is as expected (e.g. with .map())
def validate_size(collection, size){
    if ( collection.size() != size ) { return false } else { return true }
}
