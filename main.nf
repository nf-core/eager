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
      --reads                       Path to input data (must be surrounded with quotes). For paired end data, the path must use '{1,2}' notation to specify read pairs
      -profile                      Institution or personal hardware config to use (e.g. standard, docker, singularity, conda, aws). Ask your system admin if unsure, or check documentation
      --singleEnd                   Specifies that the input is single end reads (required if not pairedEnd)
      --pairedEnd                   Specifies that the input is paired end reads (required if not singleEnd)
      --bam                         Specifies that the input is in BAM format
      --fasta                       Path and name of FASTA reference file (required if not iGenome reference). File suffixes can be: '.fa', '.fn', '.fna', '.fasta'
      --genome                      Name of iGenomes reference (required if not fasta reference)

    Output options:     
      --outdir                      The output directory where the results will be saved
      -w                            The directory where intermediate files will be stored. Recommended: '<outdir>/work/'

    BAM Input:
    --run_convertbam                Species to convert an input BAM file into FASTQ format before processing.

    Input Data Additional Options:
      --snpcapture                  Runs in SNPCapture mode (specify a BED file if you do this!)

    References                      Optional additional pre-made indicies, or you wish to overwrite any of the references.
      --bwa_index                   Path and name of a bwa indexed FASTA reference file with index suffixes (i.e. everything before the endings '.amb' '.ann' '.bwt'. Most likely the same value supplied with the --fasta option)
      --bedfile                     Path to BED file for SNPCapture methods
      --seq_dict                    Path to picard sequence dictionary file (typically ending in '.dict')
      --fasta_index                 Path to samtools FASTA index (typically ending in '.fai')
      --saveReference               Saves reference genome indices for later reusage

    Skipping                        Skip any of the mentioned steps
      --skip_fastqc                 Skips both pre- and post-Adapter Removal FastQC steps.
      --skip_adapterremoval         
      --skip_mapping                Note: this maybe useful when input is a BAM file
      --skip_preseq
      --skip_damage_calculation
      --skip_qualimap
      --skip_deduplication

    Complexity Filtering 
      --complexity_filter_poly_g        Run poly-G removal on FASTQ files
      --complexity_filter_poly_g_min    Specify length of poly-g min for clipping to be performed (default: 10)

    Clipping / Merging
      --clip_forward_adaptor        Specify adapter sequence to be clipped off (forward). Default: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
      --clip_reverse_adaptor        Specify adapter sequence to be clipped off (reverse). Default: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA'
      --clip_readlength             Specify read minimum length to be kept for downstream analysis. Default: 30
      --clip_min_read_quality       Specify minimum base quality for trimming off bases. Default: 20 
      --min_adap_overlap            Specify minimum adapter overlap: 1
      --skip_collapse               Skip merging forward and reverse reads together. (Only for PE samples)
      --skip_trim                   Skip adapter and quality trimming
      --preserve5p                  Skip 5p quality base trimming (n, score, window) at 5p end.
      --mergedonly                  Send downstream only merged reads (unmerged reads and singletons are discarded).

    Mapping
      --mapper                      Specify which mapper to use. Options: 'bwaaln', 'bwamem', 'circularmapper'. Default: 'bwaaln'
      --bwaalnn                     Specify the -n parameter for BWA aln. Default: 0.3
      --bwaalnk                     Specify the -k parameter for BWA aln. Default: 2
      --bwaalnl                     Specify the -l parameter for BWA aln. Default: 32
      --circularextension           Specify the number of bases to extend reference by
      --circulartarget              Specify the target chromosome for CM
      --circularfilter              Specify to filter off-target reads
 
    Stripping
      --strip_input_fastq           Create pre-Adapter Removal FASTQ files without reads that mapped to reference (e.g. for public upload of privacy sensitive non-host data)
      --strip_mode                  Stripping mode. Remove mapped reads completely from FASTQ (strip) or just mask mapped reads sequence by N (replace)
      
    BAM Filtering
      --run_bam_filtering                Turn on samtools filter for mapping quality or unmapped reads of BAM files.
      --bam_mapping_quality_threshold    Minimum mapping quality for reads filter, default 0.
      --bam_discard_unmapped             Discards unmapped reads in either FASTQ or BAM format, depending on choice (see --bam_unmapped_type).
      --bam_unmapped_type                Defines whether to discard all unmapped reads, keep only bam and/or keep only fastq format Options: 'discard', 'bam', 'fastq', 'both'.
    
    DeDuplication
      --dedupper                    Deduplication method to use. Default: dedup. Options: dedup, markduplicates
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
      
    Annotation Statistics
      --run_bedtools_coverage       Turn on ability to calculate no. reads, depth and breadth coverage of features in reference
      --anno_file                   Path to GFF or BED file containing positions of features in reference file (--fasta). Path should be enclosed in quotes

    BAM Trimming
      --run_trim_bam                Turn on BAM trimming for UDG(+ or 1/2) protocols
      --bamutils_clip_left          Specify the number of bases to clip off reads from 'left' end of read
      --bamutils_clip_right         Specify the number of bases to clip off reads from 'right' end of read
      --bamutils_softclip           Use softclip instead of hard masking

    Genotyping
      --run_genotyping                Perform genotyping on deduplicated BAMs.
      --genotyping_tool               Specify which genotyper to use either GATK UnifiedGenotyper, GATK HaplotypeCaller or Freebayes. Note: UnifiedGenotyper uses now deprecated GATK 3.5 and requires internet access. Options: 'ug', 'hc', 'freebayes'
      --genotyping_source             Specify which input BAM to use for genotyping. Options: 'raw', 'trimmed' or 'pmd' Default: 'raw'
      --gatk_call_conf                Specify GATK phred-scaled confidence threshold. Default: 30.
      --gatk_ploidy                   Specify GATK organism ploidy. Default: 2.
      --gatk_dbsnp                    Specify VCF file for output VCF SNP annotation (Optional). Gzip not accepted.
      --gatk_ug_out_mode              Specify GATK output mode. Options: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. Default: 'EMIT_VARIANTS_ONLY'. 
      --gatk_hc_out_mode              Specify GATK output mode. Options: 'EMIT_VARIANTS_ONLY', E'MIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_ACTIVE_SITES'. Default: 'EMIT_VARIANTS_ONLY'. 
      --gatk_ug_genotype_model        Specify UnifiedGenotyper likelihood model. Options: 'SNP', 'INDEL', 'BOTH', 'GENERALPLOIDYSNP', 'GENERALPLOIDYINDEL'.  Default: 'SNP'. 
      --gatk_hc_emitrefconf           Specify HaplotypeCaller mode for emitting reference confidence calls . Options: 'NONE', 'BP_RESOLUTION', 'GVCF'. Default: 'GVCF'.
      --gatk_downsample               Maximum depth coverage allowed for genotyping before downsampling is turned on. Default: 250
      --gatk_ug_defaultbasequalities  Supply a default base quality if a read is missing a base quality score. Default: -1 (turned off)
      --freebayes_C                   Specify minimum required supporting observations to consider a variant. Default: 1
      --freebayes_g                   Specify to skip over regions of high depth by discarding alignments overlapping positions where total read depth is greater than specified in --freebayes_C. Default: turned off.
      --freebayes_p                   Specify ploidy of sample in FreeBayes. Default: 2 (diploid).

    Concensus Sequence Generation
      --run_vcf2genome              Turns on ability to create a concensus sequence FASTA file based on a UnifiedGenotyper VCF file and the original reference (only considers SNPs).
      --vcf2genome_outfile          Specify name of the output FASTA file containing the concensus sequence. Do not inclvcf2 Default: '<input_vcf>'
      --vcf2genome_header           Specify the header name of the concensus sequence entry within the FASTA file. Default: '<input_vcf>'
      --vcf2genome_minc             Minimum depth coverage required for a call to be included (else N will be called). Default: 5
      --vcf2genome_minq             Minimum genotyping quality of a call to be called. Else N will be called. Default: 30
      --vcf2genome_minfreq          Minimum fraction of reads supporting a call to be included. Else N will be called. Default: 0.8

    SNP Table Generation
      --run_multivcfanalyzer        Turn on MultiVCFAnalyzer. Note: This currently only supports diploid GATK UnifiedGenotyper input. Default: false
      --write_allele_frequencies    Specify to also write allele frequencies in the SNP table. Default: turned off.
      --min_genotype_quality        Specify the minimum genotyping quality threshold for a SNP to be called. Default: 30
      --min_base_coverage           Specify the minimum number of reads a position needs to be covered to be considered for base calling. Default: 5
      --min_allele_freq_hom         Specify the minimum allele frequency that a base requires to be considered a 'homozygous' call. Default: 0.9
      --min_allele_freq_het         Specify the minimum allele frequency that a base requires to be considered a 'heterozygous' call. Default: 0.9
      --additional_vcf_files        Specify paths to additional pre-made VCF files to be included in the SNP table generation. Use wildcard(s) for multiple files. (Optional)
      --reference_gff_annotations   Specify the reference genome annotations in '.gff' format. (Optional)
      --reference_gff_exclude       Specify positions to be excluded in '.gff' format. (Optional)
      --snp_eff_results             Specify the output file from SNP effect analysis in '.txt' format. (Optional)

    Mitochondrial to Nuclear Ratio
      --run_mtnucratio              Turn on mitochondrial to nuclear ratio calculation.
      --mtnucratio_header           Specify the name of the reference FASTA entry corresponding to the mitochondrial genome (up to the first space). Default: 'MT'

    Sex Determination
      --run_sexdeterrmine           Turn on sex determination.
      --sexdeterrmine_bedfile       Specify SNP panel in bed format for error bar calculation. (Optional, see documentation)

    Nuclear Contamination for Human DNA
      --run_nuclear_contamination   Enable nuclear contamination estimation.
      --contamination_chrom_name    The name of the chromosome in your bam. 'X' for hs37d5, 'chrX' for HG19. 

    Metagenomic Screening
      --run_metagenomic_screening   Turn on metagenomic screening module for reference-unmapped reads
      --metagenomic_tool            Specify which classifier to use. Options: 'malt'. Default: 'malt'
      --database                    Specify path to classifer database directory.
      --percent_identity            Percent identity value threshold. Default: 85
      --malt_mode                   Specify which alignment method to use. Options: 'Unknown', 'BlastN', 'BlastP', 'BlastX', 'Classifier'. Default: 'BlastN'
      --malt_alignment_mode         Specify alignment method. Options: 'Local', 'SemiGlobal'. Default: 'SemiGlobal'
      --malt_top_percent            Specify the percent for LCA algorithm (see MEGAN6 CE manual). Default: 1
      --malt_min_support_mode       Specify whether to use percent or raw number of reads for minimum support required for taxon to be retained. Options: 'percent', 'reads'. Default: 'percent'
      --malt_min_support_percent    Specify the minimum percentage of reads a taxon of sample total is required to have to be retained. Default: 0.01
      --malt_min_support_reads      Specify a minimum number of reads  a taxon of sample total is required to have to be retained. Not compatible with . Default: 1
      --malt_max_queries            Specify the maximium number of queries a read can have. Default: 100
      --malt_memory_mode            Specify the memory load method. Do not use 'map' with GTFS file system. Options: 'load', 'page', 'map'. Default: 'load'

    Metagenomic Authentication
      --run_maltextract                  Turn on MaltExtract for MALT aDNA characteristics authentication
      --maltextract_taxon_list           Path to a txt file with taxa of interest (one taxon per row, NCBI taxonomy name format)
      --maltextract_ncbifiles            Path to directory containing containing NCBI resource files (ncbi.tre and ncbi.map; avaliable: https://github.com/rhuebler/HOPS/)
      --maltextract_filter               Specify which MaltExtract filter to use. Options: 'def_anc', 'ancient', 'default', 'crawl', 'scan', 'srna', 'assignment'. Default: 'def_anc' 
      --maltextract_toppercent           Specify percent of top alignments to use. Default: 0.01
      --maltextract_destackingoff        Turn off destacking.
      --maltextract_downsamplingoff      Turn off downsampling.
      --maltextract_duplicateremovaloff  Turn off duplicate removal.
      --maltextract_matches              Export alignments of hits in BLAST format. Default: off
      --maltextract_megansummary         Export MEGAN summary files. Default: off
      --maltextract_percentidentity      Minimum percent identity alignments are required to have to be reported. Recommended to set same as MALT parameter. Default: 85.0
      --maltextract_topalignment         Use top alignments per read after filtering. Default: off
      --maltextract_singlestranded       Switch damage patterns to single-stranded mode. Default: off

    Other options:     
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --max_memory                  Memory limit for each step of pipeline. Should be in form e.g. --max_memory '8.GB'
      --max_time                    Time limit for each step of the pipeline. Should be in form e.g. --max_memory '2.h'
      --max_cpus                    Maximum number of CPUs to use for each step of the pipeline. Should be in form e.g. --max_cpus 1
      --email                       Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --plaintext_email             Receive plain text emails rather than HTML
      --maxMultiqcEmailFileSize     Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)
      
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


multiqc_config = file(params.multiqc_config)
output_docs = file("$baseDir/docs/output.md")
where_are_my_files = file("$baseDir/assets/where_are_my_files.txt")

/*
* SANITY CHECKING
*/

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

// Check if genome exists in the config file. params.genomes is from igenomes.conf, params.genome specified by user
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
}

// Index files provided? Then check whether they are correct and complete
if (params.mapper != 'bwaaln' && !params.mapper == 'circularmapper' && !params.mapper == 'bwamem'){
    exit 1, "Invalid mapper option. Options are: 'bwaaln', 'bwamem', 'circularmapper'. Default: 'bwaaln'. You gave ${params.mapper}!"
}
if( params.bwa_index && (params.mapper == 'bwaaln' | params.mapper == 'bwamem')){
    lastPath = params.bwa_index.lastIndexOf(File.separator)
    bwa_dir =  params.bwa_index.substring(0,lastPath+1)
    bwa_base = params.bwa_index.substring(lastPath+1)

    Channel
        .fromPath(bwa_dir, checkIfExists: true)
        .ifEmpty { exit 1, "BWA index directory not found: ${bwa_dir}" }
        .into {bwa_index; bwa_index_bwamem}
}

// Validate not trying to run adapterremoval on a BAM file
if (params.bam && !params.run_convertbam && !params.skip_adapterremoval ) {
    exit 1, "AdapterRemoval cannot be run on BAMs. Please validate your parameters."
}

// Validate BAM is single end only
if (params.bam && !params.singleEnd){
    exit 1, "BAM input must be used with --singleEnd "
}

// Validate that you're not trying to pass FASTQs to BAM only processes
if (params.run_convertbam && params.skip_mapping) {
  exit 1, "You can't convert a BAM to FASTQ and skip mapping! Post-mapping steps require BAM input. Please validate your parameters!"
}

// Validate that you're not trying to pass FASTQs to BAM only processes
if (params.bam && !params.run_convertbam && !params.skip_mapping) {
  exit 1, "You can't directly map a BAM file! Please supply the --run_convertbam parameter!"
}

// Validate that either pairedEnd or singleEnd has been specified by the user!
if( params.singleEnd || params.pairedEnd || params.bam){
} else {
    exit 1, "Please specify either --singleEnd, --pairedEnd to execute the pipeline on FastQ files and --bam for previously processed BAM files!"
}

// Validate that skip_collapse is only set to True for pairedEnd reads!
if (params.skip_collapse  && params.singleEnd){
    exit 1, "--skip_collapse can only be set for pairedEnd samples!"
}

// Strip mode sanity checking
if (params.strip_input_fastq){
    if (!(['strip','replace'].contains(params.strip_mode))) {
        exit 1, "--strip_mode can only be set to strip or replace!"
    }

    if (params.bam && !params.run_convertbam) {
        exit 1, "--strip_input_fastq can only be used on FASTQ, but you gave BAM input and didn't specify --run_convertbam!"
    }
}

// Mapper sanity checking
if(params.mapper != "bwaaln" && params.mapper != "bwamem" && params.mapper != "circularmapper") {
    exit 1, "Please specify a valid mapper. Options: 'bwaaln', 'bwamem', 'circularmapper'. You gave: ${params.mapper}!"
}

if (params.bam_discard_unmapped && params.bam_unmapped_type == '') {
    exit 1, "Please specify valid unmapped read output format. Options: 'discard', 'bam', 'fastq', 'both'!"
}

// Bedtools sanity checking
if(params.run_bedtools_coverage && params.anno_file == ''){
  exit 1, "You have turned on bedtools coverage, but not specified a BED or GFF file with --anno_file. Please validate your parameters!"
}

// BAM filtering sanity checking - FIRST ONE CURRENTLY DOES NOT WORK!
if (params.bam_discard_unmapped && !params.run_bam_filtering) {
  "Please turn on BAM filtering before trying to indicate how to deal with unmapped reads! Give --run_bam_filtering"
}

if (params.run_bam_filtering && params.bam_discard_unmapped && params.bam_unmapped_type == '') {
  "Please specify how to deal with unmapped reads. Options: 'discard', 'bam', 'fastq', 'both'"
}

// Genotyping sanity checking
if (params.run_genotyping){
  if (params.genotyping_tool != 'ug' && params.genotyping_tool != 'hc' && params.genotyping_tool != 'freebayes') {
  exit 1, "Please specify a genotyper. Options: 'ug', 'hc', 'freebayes'. You gave: ${params.genotyping_tool}!"
  }
  
  if (params.gatk_ug_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_ug_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_ug_out_mode != 'EMIT_ALL_SITES') {
  exit 1, "Please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. You gave: ${params.gatk_out_mode}!"
  }

  if (params.gatk_hc_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_hc_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_hc_out_mode != 'EMIT_ALL_ACTIVE_SITES') {
  exit 1, "Please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. You gave: ${params.gatk_out_mode}!"
  }
  
  if (params.genotyping_tool == 'ug' && (params.gatk_ug_genotype_model != 'SNP' && params.gatk_ug_genotype_model != 'INDEL' && params.gatk_ug_genotype_model != 'BOTH' && params.gatk_ug_genotype_model != 'GENERALPLOIDYSNP' && params.gatk_ug_genotype_model != 'GENERALPLOIDYINDEL')) {
    exit 1, "Please check your UnifiedGenotyper genotype model. Options: 'SNP', 'INDEL', 'BOTH', 'GENERALPLOIDYSNP', 'GENERALPLOIDYINDEL'. You gave: ${params.gatk_ug_genotype_model}!"
  }

  if (params.genotyping_tool == 'hc' && (params.gatk_hc_emitrefconf != 'NONE' && params.gatk_hc_emitrefconf != 'GVCF' && params.gatk_hc_emitrefconf != 'BP_RESOLUTION')) {
    exit 1, "Please check your HaplotyperCaller reference confidence parameter. Options: 'NONE', 'GVCF', 'BP_RESOLUTION'. You gave: ${params.gatk_hc_emitrefconf}!"
  }
}

// Consensus sequence generation sanity checking
if (params.run_vcf2genome) {
    if (!params.run_genotyping) {
      exit 1, "Consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Please check your genotyping parameters"
    }

    if (params.genotyping_tool != 'ug') {
      exit 1, "Consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Please check your genotyping parameters"
    }
}

// MultiVCFAnalyzer sanity checking
if (params.run_multivcfanalyzer) {
  if (!params.run_genotyping) {
    exit 1, "MultiVCFAnalyzer requires genotyping on be turned on with the parameter --run_genotyping. Please check your genotyping parameters"
  }

  if (params.genotyping_tool != "ug") {
    exit 1, "MultiVCFAnalyzer only accepts VCF files from GATK UnifiedGenotyper. Please check your genotyping parameters"
  }

  if (params.gatk_ploidy != '2') {
    exit 1, "MultiVCFAnalyzer only accepts VCF files generated with a GATK ploidy set to 2. Please check your genotyping parameters"
  }
}

// MALT sanity checking
if (params.run_metagenomic_screening) {
  if ( !params.bam_discard_unmapped ) {
  exit 1, "Metagenomic classification can only run on unmapped reads. Please supply --bam_discard_unmapped and --bam_unmapped_type 'fastq'"
  }

  if (params.bam_discard_unmapped && params.bam_unmapped_type != 'fastq' ) {
  exit 1, "Metagenomic classification can only run on unmapped reads in FASTSQ format. Please supply --bam_unmapped_type 'fastq'. You gave '${params.bam_unmapped_type}'!"
  }

  if (params.metagenomic_tool != 'malt' ) {
    exit 1, "Metagenomic classification can currently only be run with 'malt'. Please check your classifer. You gave '${params.metagenomic_tool}'!"
  }

  if (params.database == '' ) {
    exit 1, "Metagenomic classification requires a path to a database directory. Please specify one with --database '/path/to/database/'."
  }

  if (params.malt_mode != 'BlastN' && params.malt_mode != 'BlastP' && params.malt_mode != 'BlastX') {
    exit 1, "Unknown MALT mode specified. Options: 'BlastN', 'BlastP', 'BlastX'. You gave '${params.malt_mode}'!"
  }

  if (params.malt_alignment_mode != 'Local' && params.malt_alignment_mode != 'SemiGlobal') {
    exit 1, "Unknown MALT alignment mode specified. Options: 'Local', 'SemiGlobal'. You gave '${params.malt_alignment_mode}'!"
  }

  if (params.malt_min_support_mode == 'percent' && params.malt_min_support_reads != 1) {
    exit 1, "Incompatible MALT min support configuration. Percent can only be used with --malt_min_support_percent. You modified --malt_min_support_reads!"
  }

  if (params.malt_min_support_mode == 'reads' && params.malt_min_support_percent != 0.01) {
    exit 1, "Incompatible MALT min support configuration. Reads can only be used with --malt_min_supportreads. You modified --malt_min_support_percent!"
  }

  if (params.malt_memory_mode != 'load' && params.malt_memory_mode != 'page' && params.malt_memory_mode != 'map') {
    exit 1, "Unknown MALT memory mode specified. Options: 'load', 'page', 'map'. You gave '${params.malt_memory_mode}'!"
  }
}

// MaltExtract Sanity checking
if (params.run_maltextract) {

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'!"
  }

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'!"
  }

  if (params.maltextract_taxon_list == '') {
    exit 1, "MaltExtract requires a taxon list specify target taxa of interest. Specify the file with --params.maltextract_taxon_list!"
  }

  if (params.maltextract_filter != 'def_anc' && params.maltextract_filter != 'default' && params.maltextract_filter != 'ancient' && params.maltextract_filter != 'scan' && params.maltextract_filter != 'crawl' && params.maltextract_filter != 'srna') {
    exit 1, "Unknown MaltExtract filter specified. Options are: 'def_anc', 'default', 'ancient', 'scan', 'crawl', 'srna'. You gave: ${params.maltextract_filter}!"
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


/*
 * Create a channel for input read files
 * Dump can be used for debugging purposes, e.g. using the -dump-channels operator on run
 */


// If read paths
//    Is single FASTQ
//    Is paired-end FASTQ
//    Is single BAM
// If NOT read paths && FASTQ
// is NOT read paths && BAM

if( params.readPaths ){
    if( params.singleEnd && !params.bam) {
        Channel
            .from( params.readPaths )
            .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
            .ifEmpty { exit 1, "Your specified FASTQ read files did not end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq' " }
            .map { row -> [ row[0], [ file( row[1][0] ) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied!" }
            .into { ch_input_for_skipconvertbam; ch_input_for_convertbam; ch_input_for_indexbam }

    } else if (!params.bam){
        Channel
            .from( params.readPaths )
            .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
            .ifEmpty { exit 1, "Your specified FASTQ read files did not end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq' " }
            .map { row -> [ row[0], [ file( row[1][0] ), file( row[1][1] ) ] ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied!" }
            .into { ch_input_for_skipconvertbam; ch_input_for_convertbam; ch_input_for_indexbam }
    } else {
        Channel
            .from( params.readPaths )
            .filter { it =~/.*.bam/ }
            .ifEmpty { exit 1, "Your specified BAM read files did not end in: '.bam' " }
            .map { row -> [ file( row )  ] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied!" }
            .dump()
            .into { ch_input_for_skipconvertbam; ch_input_for_convertbam; ch_input_for_indexbam }

    }
} else if (!params.bam){
     Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs " +
            "to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nValid input file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq'\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { ch_input_for_skipconvertbam; ch_input_for_convertbam; ch_input_for_indexbam  }

} else {
     Channel
        .fromPath( params.reads )
        .filter { it =~/.*.bam/ }
        .map { row -> [  file( row )  ] }
        .ifEmpty { exit 1, "Cannot find any bam file matching: ${params.reads}\nValid input file types: .fastq.gz', '.fq.gz', '.fastq', or '.fq'\nNB: Path needs " +
            "to be enclosed in quotes!\n" }
        .dump() //For debugging purposes
        .into { ch_input_for_skipconvertbam; ch_input_for_convertbam; ch_input_for_indexbam }

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
summary['Skipping Mapping?'] = params.skip_mapping ? 'Yes' : 'No'
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
summary['Current Home']   = "$HOME"
summary['Current User']   = "$USER"
summary['Current Path']   = "$PWD"
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

if(!params.bwa_index && !params.fasta.isEmpty() && (params.mapper == 'bwaaln' || params.mapper == 'bwamem' || params.mapper == 'circularmapper')){
process makeBWAIndex {
    label 'sc_medium'
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
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/fasta_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.fasta_index == '' && !params.fasta.isEmpty() && ( params.mapper == 'bwaaln' || params.mapper == 'bwamem' || params.mapper == 'circularmapper')

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

ch_fai_for_skipfastaindexing.mix(ch_fasta_faidx_index) 
  .into { ch_fai_for_ug; ch_fai_for_hc; ch_fai_for_freebayes }


/*
 * PREPROCESSING - Create Sequence Dictionary for FastA if not specified on CLI 
 */

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
    tag {fasta}
    publishDir path: "${params.outdir}/reference_genome/seq_dict", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.seq_dict == '' && !params.fasta.isEmpty()

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

ch_dict_for_skipdict.mix(ch_seq_dict)
  .into { ch_dict_for_ug; ch_dict_for_hc; ch_dict_for_freebayes }

/*
* PREPROCESSING - Convert BAM to FastQ if BAM input is specified instead of FastQ file(s)
*/ 

process convertBam {
    label 'mc_small'
    tag "$bam"
    
    when: 
    params.bam && params.run_convertbam

    input: 
    file bam from ch_input_for_convertbam 

    output:
    set val("${base}"), file("*.fastq.gz") into ch_output_from_convertbam

    script:
    base = "${bam.baseName}"
    """
    samtools fastq -tn ${bam} | pigz -p ${task.cpus} > ${base}.converted.fastq.gz
    """ 
}

/*
* PREPROCESSING - Index a input BAM if not being converted to FASTQ
*/

process indexinputbam {
  label 'sc_small'
  tag "$prefix"

  when: 
  params.bam && !params.run_convertbam

  input:
  file bam from ch_input_for_indexbam

  output:
  file "*.{bai,csi}" into ch_mappingindex_for_skipmapping,ch_filteringindex_for_skiprmdup

  script:
  size = "${params.large_ref}" ? '-c' : ''
  prefix = "${bam.baseName}"
  """
  samtools index "${size}" ${bam}
  """
}

// convertbam bypass
if (params.run_convertbam) {
    ch_input_for_skipconvertbam.mix(ch_output_from_convertbam)
        .filter{ it =~/.*converted.fastq.gz/}
        .into { ch_convertbam_for_fastp; ch_convertbam_for_skipfastp; ch_convertbam_for_fastqc; ch_convertbam_for_stripfastq } 
} else {
    ch_input_for_skipconvertbam
      .into { ch_convertbam_for_fastp; ch_convertbam_for_skipfastp; ch_convertbam_for_fastqc; ch_convertbam_for_stripfastq } 
}

/*
 * STEP 1a - FastQC
 */
process fastqc {
    label 'sc_tiny'
    tag "$name"
    publishDir "${params.outdir}/FastQC/input_fastq", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: 
    !params.bam && !params.skip_fastqc || params.bam && params.run_convertbam

    input:
    set val(name), file(reads) from ch_convertbam_for_fastqc

    output:
    file "*_fastqc.{zip,html}" into ch_prefastqc_for_multiqc

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
    label 'mc_small'
    tag "$name"
    publishDir "${params.outdir}/FastP", mode: 'copy'

    when: 
    !params.bam && params.complexity_filter_poly_g || params.bam && params.run_convertbam && params.complexity_filter_poly_g

    input:
    set val(name), file(reads) from ch_convertbam_for_fastp

    output:
    set val(name), file("*pG.fq.gz") into ch_output_from_fastp
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


// fastp bypass
if (params.complexity_filter_poly_g) {
    ch_convertbam_for_skipfastp.mix(ch_output_from_fastp)
        .filter { it =~/.*pG.fq.gz/ }
        .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval } 
} else {
    ch_convertbam_for_skipfastp
        .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval } 
}


/*
 * STEP 2 - Adapter Clipping / Read Merging
 */

process adapter_removal {
    label 'mc_small'
    tag "$name"
    publishDir "${params.outdir}/read_merging", mode: 'copy'

    when: 
    !params.bam && !params.skip_adapterremoval || params.bam && params.run_convertbam && !params.skip_adapterremoval

    input:
    set val(name), file(reads) from ch_fastp_for_adapterremoval

    output:
    set val(base), file("output/*.gz") into ch_output_from_adapterremoval, ch_adapterremoval_for_postfastqc
    file("output/*.settings") into ch_adapterremoval_logs

    script:
    base = reads[0].baseName
    //This checks whether we skip trimming and defines a variable respectively
    trim_me = params.skip_trim ? '' : "--trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}"
    collapse_me = params.skip_collapse ? '' : '--collapse'
    preserve5p = params.preserve5p ? '--preserve5p' : ''
    mergedonly = params.mergedonly ? "Y" : "N"
    
    //PE mode, dependent on trim_me and collapse_me the respective procedure is run or not :-) 
    if (!params.singleEnd && !params.skip_collapse && !params.skip_trim){
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} ${trim_me} --gzip --threads ${task.cpus} ${collapse_me} ${preserve5p}
    
    #Combine files
    if [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "N" ]; then 
      zcat *.collapsed.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > output/${base}.combined.fq.gz
    elif [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "Y" ] ; then
      zcat *.collapsed.gz | gzip > output/${base}.combined.fq.gz
    elif [ ${mergedonly} = "Y" ] ; then
      zcat *.collapsed.gz *.collapsed.truncated.gz | gzip > output/${base}.combined.fq.gz
    else
      zcat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > output/${base}.combined.fq.gz
    fi
   
    mv *.settings output/
    """
    //PE, don't collapse, but trim reads
    } else if (!params.singleEnd && params.skip_collapse && !params.skip_trim) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} --gzip --threads ${task.cpus} ${trim_me} ${collapse_me} ${preserve5p}
    mv *.settings ${base}.pair*.truncated.gz output/
    """
    //PE, collapse, but don't trim reads
    } else if (!params.singleEnd && !params.skip_collapse && params.skip_trim) {
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --file2 ${reads[1]} --basename ${base} --gzip --threads ${task.cpus} --basename ${base} ${collapse_me} ${trim_me}
    
    if [ ${mergedonly} = "Y" ]; then
      zcat *.collapsed.gz *.collapsed.truncated.gz | gzip > output/${base}.combined.fq.gz
    else
      zcat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz | gzip > output/${base}.combined.fq.gz
    fi

    mv *.settings output/
    """
    } else {
    //SE, collapse not possible, trim reads
    """
    mkdir -p output
    AdapterRemoval --file1 ${reads[0]} --basename ${base} --gzip --threads ${task.cpus} ${trim_me} ${preserve5p}
    
    mv *.settings *.truncated.gz output/
    """
    }
}


// Adapterremoval bypass
if (!params.skip_adapterremoval) {
    ch_output_from_adapterremoval.mix(ch_fastp_for_skipadapterremoval)
        .filter { it =~/.*combined.fq.gz|.*truncated.gz/ }
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_skipmap; ch_adapteremoval_for_bwa; ch_adapteremoval_for_cm; ch_adapteremoval_for_bwamem } 
} else {
    ch_fastp_for_skipadapterremoval
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_skipmap; ch_adapteremoval_for_bwa; ch_adapteremoval_for_cm; ch_adapteremoval_for_bwamem;  } 
}



/*
* STEP 2b - FastQC after clipping/merging (if applied!)
*/
process fastqc_after_clipping {
    label 'sc_tiny'
    tag "${name}"
    publishDir "${params.outdir}/FastQC/after_clipping", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    when: !params.bam  && !params.skip_adapterremoval && !params.skip_fastqc || params.bam && params.run_convertbam && !params.skip_adapterremoval && !params.skip_fastqc

    input:
    set val(name), file(reads) from ch_adapterremoval_for_fastqc_after_clipping

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
    label 'mc_medium'
    tag "${name}"
    publishDir "${params.outdir}/mapping/bwa", mode: 'copy'

    when: params.mapper == 'bwaaln' && !params.skip_mapping

    input:
    set val(name), file(reads) from ch_adapteremoval_for_bwa
    file index from bwa_index.collect()

    output:
    file "*.mapped.bam" into ch_output_from_bwa
    file "*.{bai,csi}" into ch_outputindex_from_bwa
    

    script:
    size = "${params.large_ref}" ? '-c' : ''
    fasta = "${index}/${bwa_base}"

    //PE data without merging, PE data without any AR applied
    if (!params.singleEnd && (params.skip_collapse || params.skip_adapterremoval)){
    prefix = "${reads[0].baseName}"
    """
    bwa aln -t ${task.cpus} $fasta ${reads[0]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r1.sai
    bwa aln -t ${task.cpus} $fasta ${reads[1]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta ${prefix}.r1.sai ${prefix}.r2.sai ${reads[0]} ${reads[1]} | samtools sort -@ ${task.cpus} -O bam - > ${prefix}.mapped.bam
    samtools index "${size}" "${prefix}".mapped.bam
    """
    } else {
    //PE collapsed, or SE data 
    prefix = "${reads.baseName}"
    """
    bwa aln -t ${task.cpus} $fasta $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $fasta ${prefix}.sai $reads | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".mapped.bam
    samtools index "${size}" "${prefix}".mapped.bam
    """
    }
    
}

process circulargenerator{
    label 'sc_tiny'
    tag "$prefix"
    publishDir "${params.outdir}/reference_genome/circularmapper_index", mode: 'copy', saveAs: { filename -> 
            if (params.saveReference) filename 
            else if(!params.saveReference && filename == "where_are_my_files.txt") filename
            else null
    }

    when: params.mapper == 'circularmapper' && !params.skip_mapping

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
    label 'mc_medium'
    tag "$prefix"
    publishDir "${params.outdir}/mapping/circularmapper", mode: 'copy'

    when: params.mapper == 'circularmapper' && !params.skip_mapping

    input:
    set val(name), file(reads) from ch_adapteremoval_for_cm
    file index from ch_circularmapper_indices.collect()
    file fasta from fasta_for_indexing

    output:
    file "*.mapped.bam" into ch_output_from_cm
    file "*.{bai,csi}" into ch_outputindex_from_cm
    
    script:
    filter = "${params.circularfilter}" ? '' : '-f true -x false'
    elongated_root = "${fasta.baseName}_${params.circularextension}.fasta"

    
    size = "${params.large_ref}" ? '-c' : ''

    if (!params.singleEnd && params.skip_collapse ){
    prefix = "${reads[0].baseName}"
    """ 
    bwa aln -t ${task.cpus} $elongated_root ${reads[0]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r1.sai
    bwa aln -t ${task.cpus} $elongated_root ${reads[1]} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $elongated_root ${prefix}.r1.sai ${prefix}.r2.sai ${reads[0]} ${reads[1]} > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > ${prefix}.mapped.bam
    samtools index "${size}" ${prefix}.mapped.bam
    """
    } else {
    prefix = reads[0].toString().tokenize('.')[0]
    """ 
    bwa aln -t ${task.cpus} $elongated_root $reads -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${prefix}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" $elongated_root ${prefix}.sai $reads > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > "${prefix}".mapped.bam
    samtools index "${size}" "${prefix}".mapped.bam
    """
    }
    
}

process bwamem {
    label 'mc_medium'
    tag "$prefix"
    publishDir "${params.outdir}/mapping/bwamem", mode: 'copy'

    when: params.mapper == 'bwamem' && !params.skip_mapping

    input:
    set val(name), file(reads) from ch_adapteremoval_for_bwamem
    file index from bwa_index_bwamem.collect()

    output:
    file "*.mapped.bam" into ch_output_from_bwamem
    file "*.{bai,csi}" into ch_outputindex_from_bwamem
    

    script:
    fasta = "${index}/${bwa_base}"
    prefix = "${reads[0].baseName}"
    size = "${params.large_ref}" ? '-c' : ''

    if (!params.singleEnd && params.skip_collapse){
    """
    bwa mem -t ${task.cpus} $fasta ${reads[0]} ${reads[1]} -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".mapped.bam
    samtools index "${size}" -@ ${task.cpus} "${prefix}".mapped.bam
    """
    } else {
    """
    bwa mem -t ${task.cpus} $fasta $reads -R "@RG\\tID:ILLUMINA-${prefix}\\tSM:${prefix}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${prefix}".mapped.bam
    samtools index "${size}" -@ ${task.cpus} "${prefix}".mapped.bam
    """
    }
    
}


// mapping bypass
if (!params.skip_mapping) {
    ch_output_from_bwa.mix(ch_output_from_bwamem, ch_output_from_cm)
        .filter { it =~/.*mapped.bam/ }
        .into { ch_mapping_for_filtering; ch_mapping_for_skipfiltering; ch_mapping_for_samtools_flagstat } 

    ch_outputindex_from_bwa.mix(ch_outputindex_from_bwamem, ch_outputindex_from_cm)
        .filter { it =~/.*mapped.bam.bai|.*mapped.bam.csi/ }
          .into {  ch_mappingindex_for_skipfiltering; ch_mappingindex_for_filtering } 

} else {
    ch_adapterremoval_for_skipmap
        .into { ch_mapping_for_skipfiltering; ch_mapping_for_filtering;  ch_mapping_for_samtools_flagstat }

     ch_mappingindex_for_skipmapping
        .into {  ch_mappingindex_for_skipfiltering; ch_mappingindex_for_filtering } 
}

/*
* Step 3b - flagstat
*/

process samtools_flagstat {
    label 'sc_tiny'
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    when:
    !params.skip_mapping

    input:
    file(bam) from ch_mapping_for_samtools_flagstat

    output:
    file "*.stats" into ch_flagstat_for_multiqc

    script:
    prefix = "$bam" - ~/(\.bam)?$/
    """
    samtools flagstat $bam > ${prefix}_flagstat.stats
    """
}


/*
* Step 4a - Keep unmapped/remove unmapped reads
*/

process samtools_filter {
    label 'mc_medium'
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
    file bam from ch_mapping_for_filtering

    output:
    file "*filtered.bam" into ch_output_from_filtering
    file "*.unmapped.fastq.gz" optional true into ch_bam_filtering_for_malt
    file "*.unmapped.bam" optional true
    file "*.{bai,csi}" into ch_outputindex_from_filtering

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
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
        samtools index "${size}" ${prefix}.filtered.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "fastq"){
        """
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam
        samtools index "${size}" ${prefix}.filtered.bam
        samtools fastq -tn ${prefix}.unmapped.bam | pigz -p ${task.cpus} > ${prefix}.unmapped.fastq.gz
        rm ${prefix}.unmapped.bam
        """
    } else if("${params.bam_discard_unmapped}" && "${params.bam_unmapped_type}" == "both"){
        """
        samtools view -h $bam | samtools view - -@ ${task.cpus} -f4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.unmapped.bam)
        samtools view -h $bam | samtools view - -@ ${task.cpus} -F4 -q ${params.bam_mapping_quality_threshold} -o ${prefix}.filtered.bam)
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


// samtools_filter bypass 
if (params.run_bam_filtering) {
    ch_mapping_for_skipfiltering.mix(ch_output_from_filtering)
        .filter { it =~/.*filtered.bam/ }
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_stripfastq; ch_filtering_for_flagstat } 

  ch_mappingindex_for_skipfiltering.mix(ch_outputindex_from_filtering)
        .filter { it =~/.*filtered.bam.bai|.*filtered.bam.csi/ }
        .into { ch_filteringindex_for_skiprmdup; ch_filteringindex_for_dedup; ch_filteringindex_for_markdup } 

} else {
    ch_mapping_for_skipfiltering
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_stripfastq; ch_filtering_for_flagstat } 

    ch_mappingindex_for_skipfiltering
        .into { ch_filteringindex_for_skiprmdup; ch_filteringindex_for_dedup; ch_filteringindex_for_markdup } 

}




process strip_input_fastq {
    label 'mc_medium'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/samtools/stripped_fastq", mode: 'copy'

    when: 
    params.strip_input_fastq

    input: 
    set val(name), file(fq) from ch_convertbam_for_stripfastq
    file bam from ch_filtering_for_stripfastq

    output:
    file "*.fq.gz" into ch_output_from_stripfastq


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
    label 'sc_tiny'
    tag "$prefix"
    publishDir "${params.outdir}/samtools/stats", mode: 'copy'

    when:
    params.run_bam_filtering

    input:
    file(bam) from ch_filtering_for_flagstat

    output:
    file "*.stats" into ch_bam_filtered_flagstat_for_multiqc

    script:
    prefix = "$bam" - ~/(\.bam)?$/
    """
    samtools flagstat $bam > ${prefix}_postfilterflagstat.stats
    """
}


/*
Step 5a: DeDup
*/ 

process dedup{
    label 'mc_small'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/deduplication/", mode: 'copy',
        saveAs: {filename -> "${prefix}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    file bam from ch_filtering_for_dedup

    output:
    file "*.hist" into ch_hist_for_preseq
    file "*.json" into ch_dedup_results_for_multiqc
    file "${prefix}_rmdup.bam" into ch_output_from_dedup
    file "*.{bai,csi}" into ch_outputindex_from_dedup

    script:
    prefix="${bam.baseName}"
    treat_merged="${params.dedup_all_merged}" ? '-m' : ''
    size = "${params.large_ref}" ? '-c' : ''
    
    if(params.singleEnd) {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix"_rmdup.bam
    samtools index "${size}" "$prefix"_rmdup.bam
    """  
    } else {
    """
    dedup -i $bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "$prefix"_rmdup.bam -o "$prefix"_rmdup.bam
    samtools index "${size}" "$prefix"_rmdup.bam
    """  
    }
}

/*
 Step 5b: MarkDuplicates
 */

process markDup{
    label 'mc_small'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/deduplication/"

    when:
    !params.skip_deduplication && params.dedupper != 'dedup'

    input:
    file bam from ch_filtering_for_markdup

    output:
    file "*.metrics" into ch_markdup_results_for_multiqc
    file "*.rmdup.bam" into ch_output_from_markdup
    file "*.{bai,csi}" into ch_outputindex_from_markdup


    script:
    prefix = "${bam.baseName}"
    size = "${params.large_ref}" ? '-c' : ''
    """
    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M MarkDuplicates INPUT=$bam OUTPUT=${prefix}._rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="${prefix}.rmdup.metrics" VALIDATION_STRINGENCY=SILENT
    samtools index "${size}" ${prefix}_rmdup.bam
    """
}


if (!params.skip_deduplication) {
    ch_filtering_for_skiprmdup.mix(ch_output_from_dedup, ch_output_from_markdup)
        .filter { it =~/.*_rmdup.bam/ }
        .into { ch_rmdup_for_skipdamagemanipulation; ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_rmdup_for_qualimap; ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_for_sexdeterrmine; ch_for_nuclear_contamination; ch_rmdup_for_bedtools; ch_rmdup_formtnucratio } 

    ch_filteringindex_for_skiprmdup.mix(ch_outputindex_from_dedup, ch_outputindex_from_markdup)
        .filter { it =~/.*_rmdup.bam.bai|.*_rmdup.bam.csi/ }
        .into { ch_rmdupindex_for_skipdamagemanipulation; ch_rmdupindex_for_damageprofiler; ch_rmdupindex_for_qualimap; ch_rmdupindex_for_pmdtools; ch_rmdupindex_for_bamutils } 

} else {
    ch_filtering_for_skiprmdup
        .into { ch_rmdup_for_skipdamagemanipulation; ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_rmdup_for_qualimap; ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_for_sexdeterrmine; ch_for_nuclear_contamination; ch_rmdup_for_bedtools; ch_rmdup_formtnucratio } 

    ch_filteringindex_for_skiprmdup
        .into { ch_rmdupindex_for_skipdamagemanipulation; ch_rmdupindex_for_damageprofiler; ch_rmdupindex_for_qualimap; ch_rmdupindex_for_pmdtools; ch_rmdupindex_for_bamutils } 
}


/*
Step 6: Preseq
*/

process preseq {
    label 'sc_tiny'
    tag "${input.baseName}"
    publishDir "${params.outdir}/preseq", mode: 'copy'

    when:
    !params.skip_preseq

    input:
    file input from (params.skip_deduplication ? ch_rmdup_for_preseq : ch_hist_for_preseq )

    output:
    file "${input.baseName}.ccurve" into ch_preseq_for_multiqc

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
Step 7a: DMG Assessment
*/ 

process damageprofiler {
    label 'sc_tiny'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/damageprofiler", mode: 'copy'

    when:
    !params.skip_damage_calculation

    input:
    file bam from ch_rmdup_for_damageprofiler
    file fasta from fasta_for_indexing
    file bai from ch_rmdupindex_for_damageprofiler
    

    output:
    file "${base}/*.txt"
    file "${base}/*.log"
    file "${base}/*.pdf"
    file "${base}/*.json" into ch_damageprofiler_results

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
    label 'mc_small'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/qualimap", mode: 'copy'

    when:
    !params.skip_qualimap

    input:
    file bam from ch_rmdup_for_qualimap
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
 Step 9: Bedtools
*/

// Set up channels for annotation file

if (!params.run_bedtools_coverage){
  ch_anno_for_bedtools = Channel.empty()
} else {
  Channel
    ch_anno_for_bedtools = Channel.fromPath(params.anno_file)
}

process bedtools {
  label 'mc_small'
  tag "${bam.baseName}"
  publishDir "${params.outdir}/bedtools", mode: 'copy'

  when:
  params.run_bedtools_coverage

  input:
  file bam from ch_rmdup_for_bedtools
  file anno_file from ch_anno_for_bedtools

  output:
  file "*"

  script:
  """
  bedtools coverage -a ${anno_file} -b $bam | pigz -p ${task.cpus} > "${bam.baseName}".breadth.gz
  bedtools coverage -a ${anno_file} -b $bam -mean | pigz -p ${task.cpus} > "${bam.baseName}".depth.gz
  """
}

/*
 Step 10: PMDtools
 */

process pmdtools {
    label 'mc_small'
    tag "${bam.baseName}"
    publishDir "${params.outdir}/pmdtools", mode: 'copy'

    when: params.run_pmdtools

    input: 
    file bam from ch_rmdup_for_pmdtools
    file fasta from fasta_for_indexing

    output:
    file "*.bam" into ch_output_from_pmdtools
    file "*.cpg.range*.txt"
    file "*.{bai,csi}" into ch_outputindex_from_pmdtools

    script:
    //Check which treatment for the libraries was used
    def treatment = params.pmd_udg_type ? (params.pmd_udg_type =='half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(params.snpcapture){
        snpcap = (params.pmdtools_reference_mask != '') ? "--refseq ${params.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    size = "${params.large_ref}" ? '-c' : ''
    prefix = "${bam.baseName}"
    """
    #Run Filtering step 
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --threshold ${params.pmdtools_threshold} $treatment $snpcap --header | samtools view -@ ${task.cpus} -Sb - > "${bam.baseName}".pmd.bam
    #Run Calc Range step
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --deamination --range ${params.pmdtools_range} $treatment $snpcap -n ${params.pmdtools_max_reads} > "${bam.baseName}".cpg.range."${params.pmdtools_range}".txt 
    samtools index "${size}" ${prefix}.pmd.bam
    """
}

/*
* Step 11 - BAM Trimming step using bamUtils 
* Can be used for UDGhalf protocols to clip off -n bases of each read
*/

process bam_trim {
    label 'mc_small'
    tag "${prefix}" 
    publishDir "${params.outdir}/trimmed_bam", mode: 'copy'
 
    when: params.run_trim_bam

    input:
    file bam from ch_rmdup_for_bamutils

    output: 
    file "*.trimmed.bam" into ch_output_from_bamutils
    file "*.{bai,csi}" into ch_outputindex_from_bamutils

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


if ( params.run_genotyping && params.genotyping_source == 'raw' ) {
    ch_rmdup_for_skipdamagemanipulation.mix(ch_output_from_pmdtools,ch_output_from_bamutils)
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes }

    ch_rmdupindex_for_skipdamagemanipulation.mix(ch_outputindex_from_pmdtools,ch_outputindex_from_bamutils)
        .into { ch_damagemanipulationindex_for_skipgenotyping; ch_damagemanipulationindex_for_genotyping_hc; ch_damagemanipulationindex_for_genotyping_freebayes }

} else if ( params.run_genotyping && params.genotyping_source == "trimmed" )  {
    ch_rmdup_for_skipdamagemanipulation.mix(ch_output_from_pmdtools,ch_output_from_bamutils)
        .filter { it =~/.*trimmed.bam/ }
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes } 

    ch_rmdupindex_for_skipdamagemanipulation.mix(ch_outputindex_from_pmdtools,ch_outputindex_from_bamutils)
        .filter { it =~/.*trimmed.bam.bai|.*.trimmed.bam.csi/ }
        .into { ch_damagemanipulationindex_for_skipgenotyping; ch_damagemanipulationindex_for_genotyping_hc; ch_damagemanipulationindex_for_genotyping_freebayes }

} else if ( params.run_genotyping && params.genotyping_source == "pmd" )  {
    ch_rmdup_for_skipdamagemanipulation.mix(ch_output_from_pmdtools,ch_output_from_bamutils)
        .filter { it =~/.*pmd.bam/ }
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes } 

    ch_rmdupindex_for_skipdamagemanipulation.mix(ch_outputindex_from_pmdtools,ch_outputindex_from_bamutils)
        .filter { it =~/.*pmd.bam.bai|.*.pmd.bam.csi/ }
        .into { ch_damagemanipulationindex_for_skipgenotyping; ch_damagemanipulationindex_for_genotyping_hc; ch_damagemanipulationindex_for_genotyping_freebayes }

} else if ( !params.run_genotyping && !params.run_trim_bam && !params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes } 

    ch_rmdupindex_for_skipdamagemanipulation
        .into { ch_damagemanipulationindex_for_skipgenotyping; ch_damagemanipulationindex_for_genotyping_hc; ch_damagemanipulationindex_for_genotyping_freebayes } 
} else if ( !params.run_genotyping && !params.run_trim_bam && params.run_pmdtools )  {
    ch_rmdup_for_skipdamagemanipulation
        .into { ch_damagemanipulation_for_skipgenotyping; ch_damagemanipulation_for_genotyping_ug; ch_damagemanipulation_for_genotyping_hc; ch_damagemanipulation_for_genotyping_freebayes } 

    ch_rmdupindex_for_skipdamagemanipulation
        .into { ch_damagemanipulationindex_for_skipgenotyping; ch_damagemanipulationindex_for_genotyping_hc; ch_damagemanipulationindex_for_genotyping_freebayes } 
}


/*
 Step 12a: Genotyping - UnifiedGenotyper Downloading
 NB: GATK 3.5 is the last release with VCF output in "old" VCF format, not breaking downstream tools. Therefore we need it (for now at least until downstream tools can read proper 4.2 VCFs... )
    
 */

ch_gatk_download = Channel.value("download")

 process download_gatk_v3_5 {
    label 'sc_tiny'
    when: params.run_genotyping && params.genotyping_tool == 'ug'

    input: 
    val "download" from ch_gatk_download

    output:
    file "*.jar" into ch_unifiedgenotyper_jar,ch_unifiedgenotyper_versions_jar

    """
    wget -O GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 --referer https://software.broadinstitute.org/ 'https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=3.5-0-g36282e4' 
    tar xjf GenomeAnalysisTK-3.5-0-g36282e4.tar.bz2 
    """

 }

/*
 Step 12b: Genotyping - UG
*/

 process genotyping_ug {
  'mc_small'
  tag "${prefix}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'ug'

  input:
  file fasta from fasta_for_indexing
  file jar from ch_unifiedgenotyper_jar
  file bam from ch_damagemanipulation_for_genotyping_ug
  file fai from ch_fai_for_ug
  file dict from ch_dict_for_ug

  output: 
  file "*vcf.gz" into ch_ug_for_multivcfanalyzer,ch_ug_for_vcf2genome

  script:
  prefix="${bam.baseName}"
  defaultbasequalities = params.gatk_ug_defaultbasequalities == '' ? '' : " --defaultBaseQualities ${params.gatk_ug_defaultbasequalities}" 
  if (params.gatk_dbsnp == '')
    """
    samtools index -b ${bam}
    java -jar ${jar} -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${bam}.intervals 
    java -jar ${jar} -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${bam}.intervals -o ${bam}.realign.bam
    java -jar ${jar} -T UnifiedGenotyper -R ${fasta} -I ${bam}.realign.bam -o ${bam}.unifiedgenotyper.vcf -nt ${task.cpus} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}
    pigz -p ${task.cpus} ${bam}.unifiedgenotyper.vcf
    """
  else if (params.gatk_dbsnp != '')
    """
    samtools index ${bam}
    java -jar ${jar} -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${bam}.intervals 
    java -jar ${jar} -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${bam}.intervals -o ${bam}.realign.bam
    java -jar ${jar} -T UnifiedGenotyper -R ${fasta} -I ${bam}.realign.bam -o ${bam}.unifiedgenotyper.vcf -nt ${task.cpus} --dbsnp ${params.gatk_dbsnp} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}

    pigz -p ${task.cpus} ${bam}.unifiedgenotyper.vcf
    """
 }

  process genotyping_hc {
  label  'mc_small'
  tag "${prefix}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'hc'

  input:
  file fasta from fasta_for_indexing
  file bam from ch_damagemanipulation_for_genotyping_hc
  file fai from ch_fai_for_hc
  file dict from ch_dict_for_hc
  file bai from ch_damagemanipulationindex_for_genotyping_hc

  output: 
  file "*vcf.gz" into ch_vcf_hc

  script:
  prefix="${bam.baseName}"
  if (params.gatk_dbsnp == '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${bam}.haplotypecaller.vcf -stand-call-conf ${params.gatk_call_conf} --sample-ploidy ${params.gatk_ploidy} --output-mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${bam}.haplotypecaller.vcf
    """

  else if (params.gatk_dbsnp != '')
    """
    gatk HaplotypeCaller -R ${fasta} -I ${bam} -O ${bam}.haplotypecaller.vcf --dbsnp ${params.gatk_dbsnp} -stand-call-conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} --output_mode ${params.gatk_hc_out_mode} --emit-ref-confidence ${params.gatk_hc_emitrefconf}
    pigz -p ${task.cpus} ${bam}.haplotypecaller.vcf
    """
 }

 /*
 *  Step 12c: FreeBayes genotyping, should probably add in some options for users to set 
 */ 
 process genotyping_freebayes {
  tag "${prefix}"
  publishDir "${params.outdir}/genotyping", mode: 'copy'

  when:
  params.run_genotyping && params.genotyping_tool == 'freebayes'

  input:
  file fasta from fasta_for_indexing
  file bam from ch_damagemanipulation_for_genotyping_freebayes
  file fai from ch_fai_for_freebayes
  file dict from ch_dict_for_freebayes
  file bai from ch_damagemanipulationindex_for_genotyping_freebayes

  output: 
  file "*vcf.gz" into ch_vcf_freebayes
  
  script:
  prefix="${bam.baseName}"
  skip_coverage = "${params.freebayes_g}" == 0 ? "" : "-g ${params.freebayes_g}"
  """
  freebayes -f ${fasta} -p ${params.freebayes_p} -C ${params.freebayes_C} ${skip_coverage} ${bam} > ${bam.baseName}.vcf 
  pigz -p ${task.cpus} ${bam.baseName}.vcf
  """
 }


/*
 * Step 13: VCF2Genome
*/


process vcf2genome {
  label  'mc_small'
  tag "${prefix}"
  publishDir "${params.outdir}/consensus_sequence", mode: 'copy'

  when: 
  params.run_vcf2genome

  input:
  file vcf from ch_ug_for_vcf2genome
  file fasta from fasta_for_indexing

  output:
  file "*.fasta.gz"

  script:
  prefix = "${vcf.baseName}"
  out = "${params.vcf2genome_outfile}" == '' ? "${prefix}.fasta" : "${params.vcf2genome_outfile}"
  fasta_head = "${params.vcf2genome_header}" == '' ? "${prefix}" : "${params.vcf2genome_header}"
  """
  pigz -f -d -p ${task.cpus} *.vcf.gz
  vcf2genome -draft ${out}.fasta -draftname "${fasta_head}" -in ${vcf.baseName} -minc ${params.vcf2genome_minc} -minfreq ${params.vcf2genome_minfreq} -minq ${params.vcf2genome_minq} -ref ${fasta} -refMod ${out}_refmod.fasta -uncertain ${out}_uncertainy.fasta
  pigz -p ${task.cpus} *.fasta 
  pigz -p ${task.cpus} *.vcf
  """
}


/*
 * Step 13: SNP Table Generation
 */

// Create input channel for MultiVCFAnalyzer, possibly mixing with pre-made VCFs
if (params.additional_vcf_files == '') {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.collect()
} else {
    ch_extravcfs_for_multivcfanalyzer = Channel.fromPath(params.additional_vcf_files)
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.mix(ch_extravcfs_for_multivcfanalyzer).collect()
}

 process multivcfanalyzer {
  label  'mc_small'
  publishDir "${params.outdir}/MultiVCFAnalyzer", mode: 'copy'

  when:
  params.genotyping_tool == 'ug' && params.run_multivcfanalyzer && params.gatk_ploidy == '2'

  input:
  file fasta from fasta_for_indexing
  file vcf from ch_vcfs_for_multivcfanalyzer

  output:
  file 'fullAlignment.fasta.gz' into ch_output_multivcfanalyzer_fullalignment
  file 'info.txt.gz' into ch_output_multivcfanalyzer_info
  file 'snpAlignment.fasta.gz' into ch_output_multivcfanalyzer_snpalignment
  file 'snpAlignmentIncludingRefGenome.fasta.gz' into ch_output_multivcfanalyzer_snpalignmentref
  file 'snpStatistics.tsv.gz' into ch_output_multivcfanalyzer_snpstatistics
  file 'snpTable.tsv.gz' into ch_output_multivcfanalyzer_snptable
  file 'snpTableForSnpEff.tsv.gz' into ch_output_multivcfanalyzer_snptablesnpeff
  file 'snpTableWithUncertaintyCalls.tsv.gz' into ch_output_multivcfanalyzer_snptableuncertainty
  file 'structureGenotypes.tsv.gz' into ch_output_multivcfanalyzer_structuregenotypes
  file 'structureGenotypes_noMissingData-Columns.tsv.gz' into ch_output_multivcfanalyzer_structuregenotypesclean

  script:
  write_freqs = "$params.write_allele_frequencies" ? "T" : "F"
  """
  gunzip -f *.vcf.gz
  multivcfanalyzer ${params.snp_eff_results} ${fasta} ${params.reference_gff_annotations} . ${write_freqs} ${params.min_genotype_quality} ${params.min_base_coverage} ${params.min_allele_freq_hom} ${params.min_allele_freq_het} ${params.reference_gff_exclude} *.vcf
  pigz -p ${task.cpus} *.tsv *.txt snpAlignment.fasta snpAlignmentIncludingRefGenome.fasta fullAlignment.fasta
  rm *.vcf
  """
 }

 /*
  * Step 14 Mitochondrial to Nuclear Ratio
 */

 process mtnucratio {
  tag "${prefix}"
  publishDir "${params.outdir}/mtnucratio", mode: "copy"

  when: 
  params.run_mtnucratio

  input:
  file bam from ch_rmdup_formtnucratio

  output:
  file '*.mtnucratio'
  file '*.json' into ch_mtnucratio_for_multiqc

  script:
  prefix="${bam.baseName}"
  """
  mtnucratio ${bam} "${params.mtnucratio_header}"
  """
 }

 /*
  * Step 15 Sex determintion with error bar calculation.
  */

if (params.sexdeterrmine_bedfile == '') {
  ch_bed_for_sexdeterrmine = file('NO_FILE')
} else {
  ch_bed_for_sexdeterrmine = Channel.fromPath(params.sexdeterrmine_bedfile)
}



 process sex_deterrmine {
    label 'sc_small'
    publishDir "${params.outdir}/sex_determination", mode:"copy"
    
     when:
     params.run_sexdeterrmine
    
     input:
     file bam from ch_for_sexdeterrmine.collect()
     file bed from ch_bed_for_sexdeterrmine

     output:
     file 'SexDet.txt'
     file '*.json' into ch_sexdet_for_multiqc
     
     script:
     if (params.sexdeterrmine_bedfile == '') {
         """
         for i in *.bam; do
             echo \$i >> bamlist.txt
         done
        
         samtools depth -aa -q30 -Q30 -f bamlist.txt| sexdeterrmine -f bamlist.txt >SexDet.txt
         """
         } else {
         """
         for i in *.bam; do
             echo \$i >> bamlist.txt
         done
        
         samtools depth -aa -q30 -Q30 -b ${bed} -f bamlist.txt | sexdeterrmine -f bamlist.txt >SexDet.txt
         """
     }
 }

 /* 
  * Step 16 Nuclear contamination for Human DNA based on chromosome X heterozygosity.
  */
 process nuclear_contamination{
    label 'sc_small'
    publishDir "${params.outdir}/nuclear_contamination", mode:"copy"
    validExitStatus 0,134
    /*
     * ANGSD Xcontamination will exit with status 134 when the number of SNPs
     *     is not large enough for estimation.
     */

    when:
    params.run_nuclear_contamination

    input:
    file input from ch_for_nuclear_contamination

    output:
    file '*.X.contamination.out' into ch_from_nuclear_contamination

    script:
    """
    samtools index ${input}
    angsd -i ${input} -r ${params.contamination_chrom_name}:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${input.baseName}.doCounts
    contamination -a ${input.baseName}.doCounts.icnts.gz -h ${baseDir}/assets/angsd_resources/HapMapChrX.gz 2> ${input.baseName}.X.contamination.out
    """
 }
 
process print_nuclear_contamination{
    label 'sc_tiny'
    publishDir "${params.outdir}/nuclear_contamination", mode:"copy"

    when:
    params.run_nuclear_contamination

    input:
    val 'Contam' from ch_from_nuclear_contamination.collect()

    output:
    file 'nuclear_contamination.txt'

    script:
    """
    print_x_contamination.py ${Contam.join(' ')}
    """
 }

/*
 * Step 17: Metagenomic screening of unmapped reads
*/

process malt {
  label 'mc_huge'
  publishDir "${params.outdir}/metagenomic_classification", mode:"copy"

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_discard_unmapped && params.bam_unmapped_type == 'fastq'

  input:
  file fastqs from ch_bam_filtering_for_malt.collect()

  output:
  file "*.rma6" into ch_rma_for_maltExtract
  file "malt.log"

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
  -sup ${params.malt_min_support_reads} \
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

process maltextract {
  label 'mc_large'
  publishDir "${params.outdir}/MaltExtract/", mode:"copy"

  when: 
  params.run_maltextract

  input:
  file rma6 from ch_rma_for_maltExtract.collect()
  file taxon_list from ch_taxonlist_for_maltextract
  
  output:
  path "results/" type('dir')

  script:
  ncbifiles = params.maltextract_ncbifiles == '' ? "" : "-r ${params.maltextract_ncbifiles}"
  destack = params.maltextract_destackingoff ? "--destackingOff" : ""
  downsam = params.maltextract_downsamplingoff ? "--downSampOff" : ""
  dupremo = params.maltextract_duplicateremovaloff ? "--dupRemOff" : ""
  matches = params.maltextract_matches ? "--matches" : ""
  megsum = params.maltextract_megansummary ? "--meganSummary" : ""
  topaln = params.maltextract_topalignment ?  "--useTopAlignment" : ""
  ss = params.maltextract_singlestranded ? "--singleStranded" : ""
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
  """
}


/*
Genotyping tools:
- snpAD
- sequenceTools

Downstream VCF tools:
- gencons?
- READ/mcMLKin?
- popGen output? PLINK? 
*/

/*
 * Step 18a - Output Description HTML
 */
process output_documentation {
    label 'sc_tiny'
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
 * Step 18b - Parse software version numbers
 */
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
    picard MarkDuplicates --version &> v_markduplicates.txt  2>&1 || true
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

    ## Hardcoded as no --version flag or equivalent
    echo 'version 3.5-0-g36282e4' > v_gatk3_5.txt

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

/*
 * Step 18c - MultiQC
 */
process multiqc {
    label 'sc_tiny'

    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc_raw/*') from ch_prefastqc_for_multiqc.collect().ifEmpty([])
    file('fastqc/*') from ch_fastqc_after_clipping.collect().ifEmpty([])
    file software_versions_mqc from software_versions_yaml.collect().ifEmpty([])
    file ('adapter_removal/*') from ch_adapterremoval_logs.collect().ifEmpty([])
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
 * Step 18d - Completion e-mail notification
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

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
      log.info "${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}"
      log.info "${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}"
      log.info "${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}"
    }

    if (workflow.success) {
        log.info "${c_purple}nf-core/eager${c_green} Pipeline completed successfully${c_reset}"
    } else {
        checkHostname()
        log.info "${c_purple}nf-core/eager${c_red} Pipeline completed with errors${c_reset}"
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
