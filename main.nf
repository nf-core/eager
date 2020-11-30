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


def helpMessage() {
    log.info nfcoreHeader()
    log.info"""
    =========================================
    eager v${workflow.manifest.version}
    =========================================
    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/eager -profile <docker/singularity/conda> --reads'*_R{1,2}.fastq.gz' --fasta '<your_reference>.fasta'

    Mandatory arguments:
      -profile [str]                  Configuration profile to use. Can use multiple (comma separated). Ask system administrator if unsure.
                                      Available: conda, docker, singularity, test, awsbatch, <institute> and more
    Input
      --input [file]            Either paths or URLs to FASTQ/BAM data (must be surrounded with quotes). Indicate multiple files with wildcards (*). For paired end data, the path must use '{1,2}' notation to specify read pairs.
                                OR 
                                A path to a TSV file (ending .tsv) containing file paths and sequencing/sample metadata. Allows for merging of multiple lanes/libraries/samples. Please see documentation for template.

      --udg_type [str]          Specify here if you have UDG treated libraries, Set to 'half' for partial treatment, or 'full' for UDG. If not set, libraries are assumed to have no UDG treatment ('none'). Not required for TSV input. Default: ${params.udg_type}
      --single_stranded [bool]  Specifies that libraries are single stranded. Only effects MaltExtract and genotyping pileupCaller. Not required for TSV input.
      --single_end [bool]       Specifies that the input is single end reads. Not required for TSV input.
      --colour_chemistry [num]  Specifies what Illumina sequencing chemistry was used. Used to inform whether to poly-G trim if turned on (see below). Not required for TSV input. Options: 2, 4. Default: ${params.colour_chemistry}
      --bam [bool]              Specifies that the input is in BAM format. Not required for TSV input.


    Additional Options:
      --snpcapture_bed [file]       If library result of SNP capture, path to BED file containing SNPs positions on reference genome.
      --run_convertinputbam [bool]  Turns on conversion of an input BAM file into FASTQ format to allow re-preprocessing (e.g. AdapterRemoval etc.).

    References
      --fasta [file]          Path or URL to a FASTA reference file (required if not iGenome reference). File suffixes can be: '.fa', '.fn', '.fna', '.fasta'
      --genome [str]          Name of iGenomes reference (required if not FASTA reference).
      --bwa_index [dir]       Path to directory containing pre-made BWA indices (i.e. everything before the endings '.amb' '.ann' '.bwt'. Most likely the same path as --fasta). If not supplied will be made for you.
      --bt2_index [dir]       Path to directory containing pre-made Bowtie2 indices (i.e. everything before the endings e.g. '.1.bt2', '.2.bt2', '.rev.1.bt2'. Most likely the same value as --fasta). If not supplied will be made for you.
      --fasta_index [file]    Path to samtools FASTA index (typically ending in '.fai').
      --seq_dict [file]       Path to picard sequence dictionary file (typically ending in '.dict').
      --large_ref [bool]      Specify to generate more recent '.csi' BAM indices. If your reference genome is larger than 3.5GB, this is recommended due to more efficient data handling with the '.csi' format over the older '.bai'.
      --save_reference [bool] Turns on saving reference genome indices for later re-usage.

    Output options:     
      --outdir [dir]            The output directory where the results will be saved. Default: ${params.outdir}
      -w [dir]                  The directory where intermediate files will be stored. Recommended: '<outdir>/work/'

    Skipping                            Skip any of the mentioned steps.
      --skip_fastqc [bool]              Skips both pre- and post-Adapter Removal FastQC steps.
      --skip_adapterremoval [bool]  
      --skip_preseq [bool]      
      --skip_deduplication [bool] 
      --skip_damage_calculation [bool]
      --skip_qualimap [bool]

    Complexity Filtering 
      --complexity_filter_poly_g [bool]     Turn on running poly-G removal on FASTQ files. Will only be performed on 2 colour chemistry machine sequenced libraries.
      --complexity_filter_poly_g_min [num]  Specify length of poly-g min for clipping to be performed. Default: ${params.complexity_filter_poly_g_min}

    Clipping / Merging
      --clip_forward_adaptor [str]  Specify adapter sequence to be clipped off (forward strand). Default: '${params.clip_forward_adaptor}'
      --clip_reverse_adaptor [str]  Specify adapter sequence to be clipped off (reverse strand). Default: '${params.clip_reverse_adaptor}'
      --clip_readlength [num]       Specify read minimum length to be kept for downstream analysis. Default: ${params.clip_readlength}
      --clip_min_read_quality [num] Specify minimum base quality for trimming off bases. Default: ${params.clip_min_read_quality}
      --min_adap_overlap [num]      Specify minimum adapter overlap: Default: ${params.min_adap_overlap}
      --skip_collapse [bool]        Skip merging forward and reverse reads together. Only applicable for paired-end libraries.
      --skip_trim [bool]            Skip adapter and quality trimming
      --preserve5p [bool]           Skip 5p quality base trimming (n, score, window) of 5 prime end.
      --mergedonly [bool]           Only use merged reads downstream (un-merged reads and singletons are discarded).

    Mapping
      --mapper [str]             Specify which mapper to use. Options: 'bwaaln', 'bwamem', 'circularmapper', 'bowtie2'. Default: '${params.mapper}'
      --bwaalnn [num]            Specify the -n parameter for BWA aln, i.e. amount of allowed mismatches in alignments. Default: ${params.bwaalnn}
      --bwaalnk [num]            Specify the -k parameter for BWA aln, i.e. maximum edit distance allowed in a seed. Default: ${params.bwaalnk}
      --bwaalnl [num]            Specify the -l parameter for BWA aln, i.e. length of seeds to be used. Set to 1024 for whole read. Default: ${params.bwaalnl}
      --circularextension [num]  Specify the number of bases to extend reference by (circularmapper only). Default: ${params.circularextension}
      --circulartarget [chr]     Specify the FASTA header of the target chromosome to extend(circularmapper only). Default: '${params.circulartarget}'
      --circularfilter [bool]    Turn on to filter off-target reads (circularmapper only).
      --bt2_alignmode [str]      Specify the bowtie2 alignment mode. Options:  'local', 'end-to-end'. Default: '${params.bt2_alignmode}'
      --bt2_sensitivity [str]    Specify the level of sensitivity for the bowtie2 alignment mode. Options: 'no-preset', 'very-fast', 'fast', 'sensitive', 'very-sensitive'. Default: '${params.bt2_sensitivity}'
      --bt2n [num]               Specify the -N parameter for bowtie2 (mismatches in seed). This will override defaults from alignmode/sensitivity. Default: ${params.bt2n}
      --bt2l [num]               Specify the -L parameter for bowtie2 (length of seed substrings). This will override defaults from alignmode/sensitivity. Default: ${params.bt2l}
      --bt2_trim5 [num]          Specify number of bases to trim off from 5' (left) end of read before alignment. Default: ${params.bt2_trim5}
      --bt2_trim3 [num]          Specify number of bases to trim off from 3' (right) end of read before alignment. Default: ${params.bt2_trim3}

    Host removal
      --hostremoval_input_fastq [bool]    Turn on creating pre-Adapter Removal FASTQ files without reads that mapped to reference (e.g. for public upload of privacy sensitive non-host data)
      --hostremoval_mode [str]            Host DNA Removal mode. Remove mapped reads completely from FASTQ (remove) or just mask mapped reads sequence by N (replace). Default: '${params.hostremoval_mode}'
      
    BAM Filtering
      --run_bam_filtering [bool]               Turn on filtering of mapping quality, read lengths, or unmapped reads of BAM files.
      --bam_mapping_quality_threshold [num]    Minimum mapping quality for reads filter. Default: ${params.bam_mapping_quality_threshold}
      --bam_filter_minreadlength [num]         Specify minimum read length to be kept after mapping.
      --bam_unmapped_type [str]                Defines whether to discard all unmapped reads, keep both mapped and unmapped together, or save as bam and/or only fastq format Options: 'discard', 'bam', 'keep', 'fastq', 'both'. Default: '${params.bam_unmapped_type}'

    DeDuplication
      --dedupper [str]                    Deduplication method to use. Options: 'markduplicates', 'dedup'. Default: '${params.dedupper}'
      --dedup_all_merged [bool]           Turn on treating all reads as merged reads.

    Library Complexity Estimation
      --preseq_step_size [num]            Specify the step size of Preseq. Default: ${params.preseq_step_size}

    (aDNA) Damage Analysis
      --damageprofiler_length [num]       Specify length filter for DamageProfiler. Default: ${params.damageprofiler_length}
      --damageprofiler_threshold [num]    Specify number of bases of each read to consider for DamageProfiler calculations. Default: ${params.damageprofiler_threshold}
      --damageprofiler_yaxis [float]      Specify the maximum misincorporation frequency that should be displayed on damage plot. Set to 0 to 'autoscale'. Default: ${params.damageprofiler_yaxis} 
      --run_pmdtools [bool]               Turn on PMDtools
      --pmdtools_range [num]              Specify range of bases for PMDTools. Default: ${params.pmdtools_range} 
      --pmdtools_threshold [num]          Specify PMDScore threshold for PMDTools. Default: ${params.pmdtools_threshold} 
      --pmdtools_reference_mask [file]    Specify a path to reference mask for PMDTools.
      --pmdtools_max_reads [num]          Specify the maximum number of reads to consider for metrics generation. Default: ${params.pmdtools_max_reads}
      
    Annotation Statistics
      --run_bedtools_coverage [bool]       Turn on ability to calculate no. reads, depth and breadth coverage of features in reference.
      --anno_file [file]                   Path to GFF or BED file containing positions of features in reference file (--fasta). Path should be enclosed in quotes.

    BAM Trimming
      --run_trim_bam [bool]                  Turn on BAM trimming. Will only run on full-UDG or half-UDG libraries.
      --bamutils_clip_half_udg_left [num]    Specify the number of bases to clip off reads from 'left' end of read for half-UDG libraries. Default: ${params.bamutils_clip_half_udg_left}
      --bamutils_clip_half_udg_right [num]   Specify the number of bases to clip off reads from 'right' end of read for half-UDG libraries. Default: ${params.bamutils_clip_half_udg_right}
      --bamutils_clip_none_udg_left [num]    Specify the number of bases to clip off reads from 'left' end of read for non-UDG libraries. Default: ${params.bamutils_clip_none_udg_left}
      --bamutils_clip_none_udg_right [num]   Specify the number of bases to clip off reads from 'right' end of read for non-UDG libraries. Default: ${params.bamutils_clip_none_udg_right}
      --bamutils_softclip [bool]             Turn on using softclip instead of hard masking.

    Genotyping
      --run_genotyping [bool]                Turn on genotyping of BAM files.
      --genotyping_tool [str]                Specify which genotyper to use either GATK UnifiedGenotyper, GATK HaplotypeCaller, Freebayes, or pileupCaller. Options: 'ug', 'hc', 'freebayes', 'pileupcaller', 'angsd'.
      --genotyping_source [str]              Specify which input BAM to use for genotyping. Options: 'raw', 'trimmed' or 'pmd'. Default: '${params.genotyping_source}'
      --gatk_call_conf [num]                 Specify GATK phred-scaled confidence threshold. Default: ${params.gatk_call_conf}
      --gatk_ploidy [num]                    Specify GATK organism ploidy. Default: ${params.gatk_ploidy}
      --gatk_downsample [num]                Maximum depth coverage allowed for genotyping before down-sampling is turned on. Default: ${params.gatk_downsample}
      --gatk_dbsnp [file]                    Specify VCF file for output VCF SNP annotation. Optional. Gzip not accepted.
      --gatk_hc_out_mode [str]               Specify GATK output mode. Options: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_ACTIVE_SITES'. Default: '${params.gatk_hc_out_mode}'
      --gatk_hc_emitrefconf [str]            Specify HaplotypeCaller mode for emitting reference confidence calls . Options: 'NONE', 'BP_RESOLUTION', 'GVCF'. Default: '${params.gatk_hc_emitrefconf}'
      --gatk_ug_out_mode [str]               Specify GATK output mode. Options: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. Default: '${params.gatk_ug_out_mode}'
      --gatk_ug_genotype_model [str]         Specify UnifiedGenotyper likelihood model. Options: 'SNP', 'INDEL', 'BOTH', 'GENERALPLOIDYSNP', 'GENERALPLOIDYINDEL'. Default: '${params.gatk_ug_genotype_model}'
      --gatk_ug_keep_realign_bam [bool]      Specify to keep the BAM output of re-alignment around variants from GATK UnifiedGenotyper.
      --gatk_ug_defaultbasequalities [num]   Supply a default base quality if a read is missing a base quality score. Setting to -1 turns this off.
      --freebayes_C [num]                    Specify minimum required supporting observations to consider a variant. Default: ${params.freebayes_C}
      --freebayes_g [num]                    Specify to skip over regions of high depth by discarding alignments overlapping positions where total read depth is greater than specified in --freebayes_C. Default: ${params.freebayes_g}
      --freebayes_p [num]                    Specify ploidy of sample in FreeBayes. Default: ${params.freebayes_p}
      --pileupcaller_bedfile [file]          Specify path to SNP panel in bed format for pileupCaller.
      --pileupcaller_snpfile [file]          Specify path to SNP panel in EIGENSTRAT format for pileupCaller.
      --pileupcaller_method [str]            Specify calling method to use. Options: 'randomHaploid', 'randomDiploid', 'majorityCall'. Default: '${params.pileupcaller_method}'
      --pileupcaller_transitions_mode [str]  Specify the calling mode for transitions. Options: 'AllSites', 'TransitionsMissing', 'SkipTransitions'. Default: '${params.pileupcaller_transitions_mode}'
      --angsd_glmodel [str]                  Specify which ANGSD genotyping likelihood model to use. Options: 'samtools', 'gatk', 'soapsnp', 'syk'. Default: '${params.angsd_glmodel}'
      --angsd_glformat [str]                 Specify which output type to output ANGSD genotyping likelihood results: Options: 'text', 'binary', 'binary_three', 'beagle'. Default: '${params.angsd_glformat}' 
      --angsd_createfasta [bool]             Turn on creation of FASTA from ANGSD genotyping likelhoood.
      --angsd_fastamethod [str]              Specify which genotype type of 'base calling' to use for ANGSD FASTA generation. Options: 'random', 'common'. Default: '${params.angsd_fastamethod}'

    Consensus Sequence Generation
      --run_vcf2genome [bool]       Turns on ability to create a consensus sequence FASTA file based on a UnifiedGenotyper VCF file and the original reference (only considers SNPs).
      --vcf2genome_outfile [str]    Specify name of the output FASTA file containing the consensus sequence. Do not include `.vcf` in the file name. Default: '<input_vcf>'
      --vcf2genome_header [str]     Specify the header name of the consensus sequence entry within the FASTA file. Default: '<input_vcf>'
      --vcf2genome_minc [num]       Minimum depth coverage required for a call to be included (else N will be called). Default: ${params.vcf2genome_minc}
      --vcf2genome_minq [num]       Minimum genotyping quality of a call to be called. Else N will be called. Default: ${params.vcf2genome_minq}
      --vcf2genome_minfreq [float]  Minimum fraction of reads supporting a call to be included. Else N will be called. Default: ${params.vcf2genome_minfreq}

    SNP Table Generation
      --run_multivcfanalyzer [bool]       Turn on MultiVCFAnalyzer. Note: This currently only supports diploid GATK UnifiedGenotyper input.
      --write_allele_frequencies [bool]   Turn on writing write allele frequencies in the SNP table.
      --min_genotype_quality [num]        Specify the minimum genotyping quality threshold for a SNP to be called. Default: ${params.min_genotype_quality}
      --min_base_coverage [num]           Specify the minimum number of reads a position needs to be covered to be considered for base calling. Default: ${params.min_base_coverage}
      --min_allele_freq_hom [float]       Specify the minimum allele frequency that a base requires to be considered a 'homozygous' call. Default: ${params.min_allele_freq_hom}
      --min_allele_freq_het [float]       Specify the minimum allele frequency that a base requires to be considered a 'heterozygous' call. Default: ${params.min_allele_freq_het}
      --additional_vcf_files [file]       Specify paths to additional pre-made VCF files to be included in the SNP table generation. Use wildcard(s) for multiple files. Optional.
      --reference_gff_annotations [file]  Specify path to the reference genome annotations in '.gff' format. Optional.
      --reference_gff_exclude [file]      Specify path to the positions to be excluded in '.gff' format. Optional.
      --snp_eff_results [file]            Specify path to the output file from SNP effect analysis in '.txt' format. Optional.

    Mitochondrial to Nuclear Ratio
      --run_mtnucratio [bool]    Turn on mitochondrial to nuclear ratio calculation.
      --mtnucratio_header [str]  Specify the name of the reference FASTA entry corresponding to the mitochondrial genome (up to the first space). Default: '${params.mtnucratio_header}'

    Sex Determination
      --run_sexdeterrmine [bool]      Turn on sex determination for human reference genomes.
      --sexdeterrmine_bedfile [file]  Specify path to SNP panel in bed format for error bar calculation. Optional (see documentation).

    Nuclear Contamination for Human DNA
      --run_nuclear_contamination [bool]  Turn on nuclear contamination estimation for human reference genomes.
      --contamination_chrom_name [str]    The name of the X chromosome in your bam or FASTA header. 'X' for hs37d5, 'chrX' for HG19. Default: '${params.contamination_chrom_name}'

    Metagenomic Screening
      --run_metagenomic_screening [bool]     Turn on metagenomic screening module for reference-unmapped reads
      --metagenomic_tool [str]               Specify which classifier to use. Options: 'malt', 'kraken'. Default: '${params.contamination_chrom_name}'
      --database [dir]                       Specify path to classifer database directory. For Kraken2 this can also be a `.tar.gz` of the directory.
      --metagenomic_min_support_reads [num]  Specify a minimum number of reads a taxon of sample total is required to have to be retained. Not compatible with . Default: ${params.metagenomic_min_support_reads}
      --percent_identity [num]               Percent identity value threshold for MALT. Default: ${params.percent_identity}
      --malt_mode [str]                      Specify which alignment method to use for MALT. Options: 'Unknown', 'BlastN', 'BlastP', 'BlastX', 'Classifier'. Default: '${params.malt_mode}'
      --malt_alignment_mode [str]            Specify alignment method for MALT. Options: 'Local', 'SemiGlobal'. Default: '${params.malt_alignment_mode}'
      --malt_top_percent [num]               Specify the percent for LCA algorithm for MALT (see MEGAN6 CE manual). Default: ${params.malt_top_percent}
      --malt_min_support_mode [str]          Specify whether to use percent or raw number of reads for minimum support required for taxon to be retained for MALT. Options: 'percent', 'reads'. Default: '${params.malt_min_support_mode}'
      --malt_min_support_percent [num]       Specify the minimum percentage of reads a taxon of sample total is required to have to be retained for MALT. Default: Default: ${params.malt_min_support_percent}
      --malt_max_queries [num]               Specify the maximium number of queries a read can have for MALT. Default: ${params.malt_max_queries}
      --malt_memory_mode [str]               Specify the memory load method. Do not use 'map' with GPFS file systems for MALT as can be very slow. Options: 'load', 'page', 'map'. Default: '${params.malt_memory_mode}'
      --malt_sam_output [bool]               Specify to also produce SAM alignment files. Note this includes both aligned and unaligned reads, and are gzipped. Note this will result in very large file sizes.

    Metagenomic Authentication
      --run_maltextract [bool]                  Turn on MaltExtract for MALT aDNA characteristics authentication
      --maltextract_taxon_list [file]           Path to a txt file with taxa of interest (one taxon per row, NCBI taxonomy name format)
      --maltextract_ncbifiles [dir]             Path to directory containing containing NCBI resource files (ncbi.tre and ncbi.map; avaliable: https://github.com/rhuebler/HOPS/)
      --maltextract_filter [str]                Specify which MaltExtract filter to use. Options: 'def_anc', 'ancient', 'default', 'crawl', 'scan', 'srna', 'assignment'. Default: '${params.maltextract_filter}'
      --maltextract_toppercent [num]            Specify percent of top alignments to use. Default: ${params.maltextract_toppercent}
      --maltextract_destackingoff [bool]        Turn off destacking.
      --maltextract_downsamplingoff [bool]      Turn off downsampling.
      --maltextract_duplicateremovaloff [bool]  Turn off duplicate removal.
      --maltextract_matches [bool]              Turn on exporting alignments of hits in BLAST format.
      --maltextract_megansummary [bool]         Turn on export of MEGAN summary files.
      --maltextract_percentidentity [num]       Minimum percent identity alignments are required to have to be reported. Recommended to set same as MALT parameter. Default: ${params.maltextract_percentidentity}
      --maltextract_topalignment [int]          Turn on using top alignments per read after filtering.

    Other options:     
      -name [str]                     Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
      --max_memory [str]              Memory limit for each step of pipeline. Should be in form e.g. --max_memory '8.GB'. Default: '${params.max_memory}'
      --max_time [str]                Time limit for each step of the pipeline. Should be in form e.g. --max_time '2.h'. Default: '${params.max_time}'
      --max_cpus [str]                Maximum number of CPUs to use for each step of the pipeline. Should be in form e.g. Default: '${params.max_cpus}'
      --publish_dir_mode [str]        Mode for publishing results in the output directory. Available: symlink, rellink, link, copy, copyNoFollow, move. Default: '${params.publish_dir_mode}'
      --email [email]                 Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail [email]         Same as --email, except only send mail if the workflow is not successful
      --plaintext_email [email]       Receive plain text emails rather than HTML
      --max_multiqc_email_size [str]  Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached (Default: 25MB)

    AWSBatch options:
      --awsqueue [str]   The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion [str]  The AWS Region for your AWS Batch job to run on
      --awscli [str]     Path to the AWS CLI tool
      
    For a full list and more information of available parameters, consider the documentation (https://github.com/nf-core/eager/).
    """.stripIndent()
}

///////////////////////////////////////////////////////////////////////////////
/* --                SET UP CONFIGURATION VARIABLES                       -- */
///////////////////////////////////////////////////////////////////////////////

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Small console separator to make it easier to read errors after launch
println ""

////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////
/**FASTA input handling
**/ 

if (params.fasta) {
    file(params.fasta, checkIfExists: true)
    lastPath = params.fasta.lastIndexOf(File.separator)
    lastExt = params.fasta.lastIndexOf(".")
    fasta_base = params.fasta.substring(lastPath+1)
    index_base = params.fasta.substring(lastPath+1,lastExt)
    if (params.fasta.endsWith('.gz')) {
        fasta_base = params.fasta.substring(lastPath+1,lastExt)
        index_base = fasta_base.substring(0,fasta_base.lastIndexOf("."))

    }
} else {
    exit 1, "[nf-core/eager] error: please specify --fasta with the path to your reference"
}

// Validate reference inputs
if("${params.fasta}".endsWith(".gz")){
    process unzip_reference{
        tag "${zipped_fasta}"

        input:
        path zipped_fasta from file(params.fasta) // path doesn't like it if a string of an object is not prefaced with a root dir (/), so use file() to resolve string before parsing to `path` 

        output:
        path "$unzip" into ch_fasta into ch_fasta_for_bwaindex,ch_fasta_for_bt2index,ch_fasta_for_faidx,ch_fasta_for_seqdict,ch_fasta_for_circulargenerator,ch_fasta_for_circularmapper,ch_fasta_for_damageprofiler,ch_fasta_for_qualimap,ch_fasta_for_pmdtools,ch_fasta_for_genotyping_ug,ch_fasta_for_genotyping_hc,ch_fasta_for_genotyping_freebayes,ch_fasta_for_genotyping_pileupcaller,ch_fasta_for_vcf2genome,ch_fasta_for_multivcfanalyzer,ch_fasta_for_genotyping_angsd

        script:
        unzip = zipped_fasta.toString() - '.gz'
        """
        pigz -f -d -p ${task.cpus} $zipped_fasta
        """
        }
    } else {
    fasta_for_indexing = Channel
    .fromPath("${params.fasta}", checkIfExists: true)
    .into{ ch_fasta_for_bwaindex; ch_fasta_for_bt2index; ch_fasta_for_faidx; ch_fasta_for_seqdict; ch_fasta_for_circulargenerator; ch_fasta_for_circularmapper; ch_fasta_for_damageprofiler; ch_fasta_for_qualimap; ch_fasta_for_pmdtools; ch_fasta_for_genotyping_ug; ch_fasta__for_genotyping_hc; ch_fasta_for_genotyping_hc; ch_fasta_for_genotyping_freebayes; ch_fasta_for_genotyping_pileupcaller; ch_fasta_for_vcf2genome; ch_fasta_for_multivcfanalyzer;ch_fasta_for_genotyping_angsd }
}

// Check that fasta index file path ends in '.fai'
if (params.fasta_index && !params.fasta_index.endsWith(".fai")) {
    exit 1, "The specified fasta index file (${params.fasta_index}) is not valid. Fasta index files should end in '.fai'."
}

// Check if genome exists in the config file. params.genomes is from igenomes.conf, params.genome specified by user
if ( params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "[nf-core/eager] error: the provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}."
}

// Mapper validation
if (params.mapper != 'bwaaln' && !params.mapper == 'circularmapper' && !params.mapper == 'bwamem' && !params.mapper == "bowtie2"){
    exit 1, "[nf-core/eager] error: invalid mapper option. Options are: 'bwaaln', 'bwamem', 'circularmapper', 'bowtie2'. Default: 'bwaaln'. Found parameter: --mapper '${params.mapper}'."
}

if (params.mapper == 'bowtie2' && params.bt2_alignmode != 'local' && params.bt2_alignmode != 'end-to-end' ) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 alignment mode. Options: 'local', 'end-to-end'. Found parameter: --bt2_alignmode '${params.bt2_alignmode}'"
}

if (params.mapper == 'bowtie2' && params.bt2_sensitivity != 'no-preset' && params.bt2_sensitivity != 'very-fast' && params.bt2_sensitivity != 'fast' && params.bt2_sensitivity != 'sensitive' && params.bt2_sensitivity != 'very-sensitive' ) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 sensitivity mode. Options: 'no-preset', 'very-fast', 'fast', 'sensitive', 'very-sensitive'. Options are for both alignmodes Found parameter: --bt2_sensitivity '${params.bt2_sensitivity}'."
}

if (params.bt2n != 0 && params.bt2n != 1) {
    exit 1, "[nf-core/eager] error: invalid bowtie2 --bt2n (-N) parameter. Options: 0, 1. Found parameter: --bt2n ${params.bt2n}."

}

// Index files provided? Then check whether they are correct and complete
if( params.bwa_index != '' && (params.mapper == 'bwaaln' | params.mapper == 'bwamem' | params.mapper == 'circularmapper')){
    Channel
        .fromPath(params.bwa_index, checkIfExists: true)
        .ifEmpty { exit 1, "[nf-core/eager] error: bwa indices not found in: ${index_base}." }
        .into {bwa_index; bwa_index_bwamem}

    bt2_index = Channel.empty()
}

if( params.bt2_index != '' && params.mapper == 'bowtie2' ){
    lastPath = params.bt2_index.lastIndexOf(File.separator)
    bt2_dir =  params.bt2_index.substring(0,lastPath+1)
    bt2_base = params.bt2_index.substring(lastPath+1)

    Channel
        .fromPath(params.bt2_index, checkIfExists: true)
        .ifEmpty { exit 1, "[nf-core/eager] error: bowtie2 indices not found in: ${bt2_dir}." }
        .into {bt2_index; bt2_index_bwamem}

    bwa_index = Channel.empty()
    bwa_index_bwamem = Channel.empty()
}

// Validate BAM input isn't set to paired_end
if ( params.bam && !params.single_end ) {
  exit 1, "[nf-core/eager] error: bams can only be specified with --single_end. Please check input command."
}

// Validate that skip_collapse is only set to True for paired_end reads!
if (!has_extension(params.input, "tsv") && params.skip_collapse  && params.single_end){
    exit 1, "[nf-core/eager] error: --skip_collapse can only be set for paired_end samples."
}

// Host removal mode validation
if (params.hostremoval_input_fastq){
    if (!(['remove','replace'].contains(params.hostremoval_mode))) {
        exit 1, "[nf-core/eager] error: --hostremoval_mode can only be set to 'remove' or 'replace'."
    }
}

if (params.bam_unmapped_type == '') {
    exit 1, "[nf-core/eager] error: please specify valid unmapped read output format. Options: 'discard', 'keep', 'bam', 'fastq', 'both'. Found parameter: --bam_unmapped_type '${params.bam_unmapped_type}'."
}

// Bedtools validation
if(params.run_bedtools_coverage && params.anno_file == ''){
  exit 1, "[nf-core/eager] error: you have turned on bedtools coverage, but not specified a BED or GFF file with --anno_file. Please validate your parameters."
}

// Set up channels for annotation file
if (!params.run_bedtools_coverage){
  ch_anno_for_bedtools = Channel.empty()
} else {
  Channel
    ch_anno_for_bedtools = Channel.fromPath(params.anno_file, checkIfExists: true)
    .ifEmpty { exit 1, "[nf-core/eager] error: bedtools annotation file not found. Supplied parameter: --anno_file ${params.anno_file}."}
}

// BAM filtering validation
if (!params.run_bam_filtering && params.bam_mapping_quality_threshold != 0) {
  exit 1, "[nf-core/eager] error: please turn on BAM filtering if you want to perform mapping quality filtering! Provide: --run_bam_filtering."
}

if (params.run_bam_filtering && params.bam_unmapped_type != 'discard' && params.bam_unmapped_type != 'keep' && params.bam_unmapped_type != 'bam' && params.bam_unmapped_type != 'fastq' && params.bam_unmapped_type != 'both' ) {
  exit 1, "[nf-core/eager] error: please specify how to deal with unmapped reads. Options: 'discard', 'keep', 'bam', 'fastq', 'both'."
}

// Deduplication validation
if (params.dedupper != 'dedup' && params.dedupper != 'markduplicates') {
  exit 1, "[nf-core/eager] error: Selected deduplication tool is not recognised. Options: 'dedup' or 'markduplicates'. Found parameter: --dedupper '${params.dedupper}'."
}

if (params.dedupper == 'dedup' && !params.mergedonly) {
    log.warn "[nf-core/eager] Warning: you are using DeDup but without specifying --mergedonly for AdapterRemoval, dedup will likely fail! See documentation for more information."
}

// SexDetermination channel set up and bedfile validation
if (params.sexdeterrmine_bedfile == '') {
  ch_bed_for_sexdeterrmine = Channel.fromPath("$projectDir/assets/nf-core_eager_dummy.txt")
} else {
  ch_bed_for_sexdeterrmine = Channel.fromPath(params.sexdeterrmine_bedfile, checkIfExists: true)
}

// Genotyping validation
if (params.run_genotyping){
  if (params.genotyping_tool != 'ug' && params.genotyping_tool != 'hc' && params.genotyping_tool != 'freebayes' && params.genotyping_tool != 'pileupcaller' && params.genotyping_tool != 'angsd' ) {
  exit 1, "[nf-core/eager] error: please specify a genotyper. Options: 'ug', 'hc', 'freebayes', 'pileupcaller'. Found parameter: --genotyping_tool '${params.genotyping_tool}'."
  }
  
  if (params.gatk_ug_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_ug_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_ug_out_mode != 'EMIT_ALL_SITES') {
  exit 1, "[nf-core/eager] error: please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. Found parameter: --gatk_ug_out_mode '${params.gatk_out_mode}'."
  }

  if (params.gatk_hc_out_mode != 'EMIT_VARIANTS_ONLY' && params.gatk_hc_out_mode != 'EMIT_ALL_CONFIDENT_SITES' && params.gatk_hc_out_mode != 'EMIT_ALL_ACTIVE_SITES') {
  exit 1, "[nf-core/eager] error: please check your GATK output mode. Options are: 'EMIT_VARIANTS_ONLY', 'EMIT_ALL_CONFIDENT_SITES', 'EMIT_ALL_SITES'. Found parameter: --gatk_out_mode '${params.gatk_out_mode}'."
  }
  
  if (params.genotyping_tool == 'ug' && (params.gatk_ug_genotype_model != 'SNP' && params.gatk_ug_genotype_model != 'INDEL' && params.gatk_ug_genotype_model != 'BOTH' && params.gatk_ug_genotype_model != 'GENERALPLOIDYSNP' && params.gatk_ug_genotype_model != 'GENERALPLOIDYINDEL')) {
    exit 1, "[nf-core/eager] error: please check your UnifiedGenotyper genotype model. Options: 'SNP', 'INDEL', 'BOTH', 'GENERALPLOIDYSNP', 'GENERALPLOIDYINDEL'. Found parameter: --gatk_ug_genotype_model '${params.gatk_ug_genotype_model}'."
  }

  if (params.genotyping_tool == 'hc' && (params.gatk_hc_emitrefconf != 'NONE' && params.gatk_hc_emitrefconf != 'GVCF' && params.gatk_hc_emitrefconf != 'BP_RESOLUTION')) {
    exit 1, "[nf-core/eager] error: please check your HaplotyperCaller reference confidence parameter. Options: 'NONE', 'GVCF', 'BP_RESOLUTION'. Found parameter: --gatk_hc_emitrefconf '${params.gatk_hc_emitrefconf}'."
  }

  if (params.genotyping_tool == 'pileupcaller' && ! ( params.pileupcaller_method == 'randomHaploid' || params.pileupcaller_method == 'randomDiploid' || params.pileupcaller_method == 'majorityCall' ) ) {
    exit 1, "[nf-core/eager] error: please check your pileupCaller method parameter. Options: 'randomHaploid', 'randomDiploid', 'majorityCall'. Found parameter: --pileupcaller_method '${params.pileupcaller_method}'."
  }

  if (params.genotyping_tool == 'pileupcaller' && ( params.pileupcaller_bedfile == '' || params.pileupcaller_snpfile == '' ) ) {
    exit 1, "[nf-core/eager] error: please check your pileupCaller bed file and snp file parameters. You must supply a bed file and a snp file."
  }

  if (params.genotyping_tool == 'angsd' && ! ( params.angsd_glmodel == 'samtools' || params.angsd_glmodel == 'gatk' || params.angsd_glmodel == 'soapsnp' || params.angsd_glmodel == 'syk' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD genotyping model! Options: 'samtools', 'gatk', 'soapsnp', 'syk'. Found parameter: --angsd_glmodel' ${params.angsd_glmodel}'."
  }

  if (params.genotyping_tool == 'angsd' && ! ( params.angsd_glformat == 'text' || params.angsd_glformat == 'binary' || params.angsd_glformat == 'binary_three' || params.angsd_glformat == 'beagle' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD genotyping model! Options: 'text', 'binary', 'binary_three', 'beagle'. Found parameter: --angsd_glmodel '${params.angsd_glmodel}'."
  }

  if ( !params.angsd_createfasta && params.angsd_fastamethod != 'random' ) {
    exit 1, "[nf-core/eager] error: to output a ANGSD FASTA file, please turn on FASTA creation with --angsd_createfasta."
  }

  if ( params.angsd_createfasta && !( params.angsd_fastamethod == 'random' || params.angsd_fastamethod == 'common' ) ) {
    exit 1, "[nf-core/eager] error: please check your ANGSD FASTA file creation method. Options: 'random', 'common'. Found parameter: --angsd_fastamethod '${params.angsd_fastamethod}'."
  }

  if (params.genotyping_tool == 'pileupcaller' && ! ( params.pileupcaller_transitions_mode == 'AllSites' || params.pileupcaller_transitions_mode == 'TransitionsMissing' || params.pileupcaller_transitions_mode == 'SkipTransitions') ) {
    exit 1, "[nf-core/eager] error: please check your pileupCaller transitions mode parameter. Options: 'AllSites', 'TransitionsMissing', 'SkipTransitions'. Found parameter: --pileupcaller_transitions_mode '${params.pileupcaller_transitions_mode}'"
  }
}

 // pileupCaller channel generation and input checks for 'random sampling' genotyping
if (params.pileupcaller_bedfile.isEmpty()) {
  ch_bed_for_pileupcaller = Channel.fromPath("$projectDir/assets/nf-core_eager_dummy.txt")
} else {
  ch_bed_for_pileupcaller = Channel.fromPath(params.pileupcaller_bedfile, checkIfExists: true)
}

if (params.pileupcaller_snpfile.isEmpty ()) {
  ch_snp_for_pileupcaller = Channel.fromPath("$projectDir/assets/nf-core_eager_dummy2.txt")
} else {
  ch_snp_for_pileupcaller = Channel.fromPath(params.pileupcaller_snpfile, checkIfExists: true)
}

// Consensus sequence generation validation
if (params.run_vcf2genome) {
    if (!params.run_genotyping) {
      exit 1, "[nf-core/eager] error: consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Please check your genotyping parameters."
    }

    if (params.genotyping_tool != 'ug') {
      exit 1, "[nf-core/eager] error: consensus sequence generation requires genotyping via UnifiedGenotyper on be turned on with the parameter --run_genotyping and --genotyping_tool 'ug'. Found parameter: --genotyping_tool '${params.genotyping_tool}'."
    }
}

// MultiVCFAnalyzer validation
if (params.run_multivcfanalyzer) {
  if (!params.run_genotyping) {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer requires genotyping to be turned on with the parameter --run_genotyping. Please check your genotyping parameters."
  }

  if (params.genotyping_tool != "ug") {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer only accepts VCF files from GATK UnifiedGenotyper. Found parameter: --genotyping_tool '${params.genotyping_tool}'."
  }

  if (params.gatk_ploidy != 2) {
    exit 1, "[nf-core/eager] error: MultiVCFAnalyzer only accepts VCF files generated with a GATK ploidy set to 2. Found parameter: --gatk_ploidy ${params.gatk_ploidy}."
  }

  if (params.additional_vcf_files != '') {
      ch_extravcfs_for_multivcfanalyzer = Channel.fromPath(params.additional_vcf_files, checkIfExists: true)
  }
}

// Metagenomic validation
if (params.run_metagenomic_screening) {
  if ( params.bam_unmapped_type == "discard" ) {
  exit 1, "[nf-core/eager] error: metagenomic classification can only run on unmapped reads. Please supply --bam_unmapped_type 'fastq'. Supplied: --bam_unmapped_type '${params.bam_unmapped_type}'."
  }

  if (params.bam_unmapped_type != 'fastq' ) {
  exit 1, "[nf-core/eager] error: metagenomic classification can only run on unmapped reads in FASTQ format. Please supply --bam_unmapped_type 'fastq'. Found parameter: --bam_unmapped_type '${params.bam_unmapped_type}'."
  }

  if (params.metagenomic_tool != 'malt' &&  params.metagenomic_tool != 'kraken') {
    exit 1, "[nf-core/eager] error: metagenomic classification can currently only be run with 'malt' or 'kraken' (kraken2). Please check your classifier. Found parameter: --metagenomic_tool '${params.metagenomic_tool}'."
  }

  if (params.database == '' ) {
    exit 1, "[nf-core/eager] error: metagenomic classification requires a path to a database directory. Please specify one with --database '/path/to/database/'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_mode != 'BlastN' && params.malt_mode != 'BlastP' && params.malt_mode != 'BlastX') {
    exit 1, "[nf-core/eager] error: unknown MALT mode specified. Options: 'BlastN', 'BlastP', 'BlastX'. Found parameter: --malt_mode '${params.malt_mode}'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_alignment_mode != 'Local' && params.malt_alignment_mode != 'SemiGlobal') {
    exit 1, "[nf-core/eager] error: unknown MALT alignment mode specified. Options: 'Local', 'SemiGlobal'. Found parameter: --malt_alignment_mode '${params.malt_alignment_mode}'."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_min_support_mode == 'percent' && params.metagenomic_min_support_reads != 1) {
    exit 1, "[nf-core/eager] error: incompatible MALT min support configuration. Percent can only be used with --malt_min_support_percent. You modified: --metagenomic_min_support_reads."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_min_support_mode == 'reads' && params.malt_min_support_percent != 0.01) {
    exit 1, "[nf-core/eager] error: incompatible MALT min support configuration. Reads can only be used with --malt_min_supportreads. You modified: --malt_min_support_percent."
  }

  if (params.metagenomic_tool == 'malt' && params.malt_memory_mode != 'load' && params.malt_memory_mode != 'page' && params.malt_memory_mode != 'map') {
    exit 1, "[nf-core/eager] error: unknown MALT memory mode specified. Options: 'load', 'page', 'map'. Found parameter: --malt_memory_mode '${params.malt_memory_mode}'."
  }

  if (!params.metagenomic_min_support_reads.toString().isInteger()){
    exit 1, "[nf-core/eager] error: incompatible min_support_reads configuration. min_support_reads can only be used with integers. --metagenomic_min_support_reads Found parameter: ${params.metagenomic_min_support_reads}."
  }
}

// Create input channel for MALT database directory, checking directory exists
if ( params.database == '') {
    ch_db_for_malt = Channel.empty()
} else {
    ch_db_for_malt = Channel.fromPath(params.database, checkIfExists: true)
}

// MaltExtract validation
if (params.run_maltextract) {

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "[nf-core/eager] error: MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'. Found parameter: --metagenomic_tool '${params.metagenomic_tool}'"
  }

  if (params.run_metagenomic_screening && params.metagenomic_tool != 'malt') {
    exit 1, "[nf-core/eager] error: MaltExtract can only accept MALT output. Please supply --metagenomic_tool 'malt'. Found parameter: --metagenomic_tool '${params.metagenomic_tool}'"
  }

  if (params.maltextract_taxon_list == '') {
    exit 1, "[nf-core/eager] error: MaltExtract requires a taxon list specifying the target taxa of interest. Specify the file with --params.maltextract_taxon_list."
  }

  if (params.maltextract_filter != 'def_anc' && params.maltextract_filter != 'default' && params.maltextract_filter != 'ancient' && params.maltextract_filter != 'scan' && params.maltextract_filter != 'crawl' && params.maltextract_filter != 'srna') {
    exit 1, "[nf-core/eager] error: unknown MaltExtract filter specified. Options are: 'def_anc', 'default', 'ancient', 'scan', 'crawl', 'srna'. Found parameter: --maltextract_filter '${params.maltextract_filter}'."
  }

}

// Create input channel for MaltExtract taxon list, to allow downloading of taxon list, checking file exists.
if ( params.maltextract_taxon_list== '' ) {
    ch_taxonlist_for_maltextract = Channel.empty()
} else {
    ch_taxonlist_for_maltextract = Channel.fromPath(params.maltextract_taxon_list, checkIfExists: true)
}

// Create input channel for MaltExtract NCBI files, checking files exists.
if ( params.maltextract_ncbifiles == '' ) {
    ch_ncbifiles_for_maltextract = Channel.empty()
} else {
    ch_ncbifiles_for_maltextract = Channel.fromPath(params.maltextract_ncbifiles, checkIfExists: true)
}


// Has the run name been specified by the user?
// this has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) {
    custom_runName = workflow.runName
}

// Check AWS batch settings
if (workflow.profile.contains('awsbatch')) {
    // AWSBatch sanity checking
    if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
    // Check outdir paths to be S3 buckets if running on AWSBatch
    // related: https://github.com/nextflow-io/nextflow/issues/813
    if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
    // Prevent trace files to be stored on S3 since S3 does not support rolling files.
    if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
}

////////////////////////////////////////////////////
/* --          CONFIG FILES                    -- */
////////////////////////////////////////////////////

ch_multiqc_config = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_eager_logo = file("$projectDir/docs/images/nf-core_eager_logo.png")
ch_output_docs = file("$projectDir/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("$projectDir/docs/images/", checkIfExists: true)
where_are_my_files = file("$projectDir/assets/where_are_my_files.txt")

///////////////////////////////////////////////////
/* --    INPUT FILE LOADING AND VALIDATING    -- */
///////////////////////////////////////////////////

// check if we have valid --reads or --input
if (params.input == null) {
  exit 1, "[nf-core/eager] error: --input was not supplied! Please check '--help' or documentation under 'running the pipeline' for details"
}

// Read in files properly from TSV file
tsv_path = null
if (params.input && (has_extension(params.input, "tsv"))) tsv_path = params.input

ch_input_sample = Channel.empty()
if (tsv_path) {

    tsv_file = file(tsv_path)
    
    if (tsv_file instanceof List) exit 1, "[nf-core/eager] error: can only accept one TSV file per run."
    if (!tsv_file.exists()) exit 1, "[nf-core/eager] error: input TSV file could not be found. Does the file exist and is it in the right place? You gave the path: ${params.input}"

    ch_input_sample = extract_data(tsv_path)

} else if (params.input && !has_extension(params.input, "tsv")) {

    log.info ""
    log.info "No TSV file provided - creating TSV from supplied directory."
    log.info "Reading path(s): ${params.input}\n"
    inputSample = retrieve_input_paths(params.input, params.colour_chemistry, params.single_end, params.single_stranded, params.udg_type, params.bam)
    ch_input_sample = inputSample

} else exit 1, "[nf-core/eager] error: --input file(s) not correctly not supplied or improperly defined, see '--help' flag and documentation under 'running the pipeline' for details."

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

// Also need to send raw files for lane merging, if we want to host removed fastq
ch_fastq_channel
  .into { ch_input_for_skipconvertbam; ch_input_for_lanemerge_hostremovalfastq }

///////////////////////////////////////////////////
/* --             HEADER LOG INFO             -- */
///////////////////////////////////////////////////

log.info nfcoreHeader()
def summary = [:]
summary['Pipeline Name']  = 'nf-core/eager'
summary['Pipeline Version'] = workflow.manifest.version
if (workflow.revision) summary['Pipeline Release'] = workflow.revision
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
  summary['Skip Read Merging'] = params.bam_unmapped_type
}
summary['Run Fastq Host Removal'] = params.hostremoval_input_fastq ? 'Yes' : 'No'
if (params.hostremoval_input_fastq){
    summary['Host removal mode'] = params.hostremoval_mode
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
summary['Max Resources']    = "$params.max_memory memory, $params.max_cpus cpus, $params.max_time time per job"
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
summary['User']             = workflow.userName
if (workflow.profile.contains('awsbatch')) {
    summary['AWS Region']   = params.awsregion
    summary['AWS Queue']    = params.awsqueue
    summary['AWS CLI']      = params.awscli
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
checkHostname()

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
    publishDir path: "${params.outdir}/reference_genome/bwa_index", mode: params.publish_dir_mode, saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    path fasta from ch_fasta_for_bwaindex
    path where_are_my_files

    output:
    path "BWAIndex" into (bwa_index, bwa_index_bwamem)
    path "where_are_my_files.txt"

    script:
    """
    bwa index $fasta
    mkdir BWAIndex && mv ${fasta}* BWAIndex
    """
    }
    
    bt2_index = Channel.empty()
}

// bowtie2 Index
if(params.bt2_index == '' && !params.fasta.isEmpty() && params.mapper == "bowtie2"){
  process makeBT2Index {
    label 'sc_medium'
    tag "${fasta}"
    publishDir path: "${params.outdir}/reference_genome/bt2_index", mode: params.publish_dir_mode, saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    path fasta from ch_fasta_for_bt2index
    path where_are_my_files

    output:
    path "BT2Index" into (bt2_index)
    path "where_are_my_files.txt"

    script:
    """
    bowtie2-build $fasta $fasta
    mkdir BT2Index && mv ${fasta}* BT2Index
    """
    }

  bwa_index = Channel.empty()
  bwa_index_bwamem = Channel.empty()

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
    publishDir path: "${params.outdir}/reference_genome/fasta_index", mode: params.publish_dir_mode, saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.fasta_index == '' && !params.fasta.isEmpty() && ( params.mapper == 'bwaaln' || params.mapper == 'bwamem' || params.mapper == 'circularmapper')

    input:
    path fasta from ch_fasta_for_faidx
    path where_are_my_files

    output:
    path "*.fai" into ch_fasta_faidx_index
    path "where_are_my_files.txt"

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
    publishDir path: "${params.outdir}/reference_genome/seq_dict", mode: params.publish_dir_mode, saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }
    
    when: params.seq_dict == '' && !params.fasta.isEmpty()

    input:
    path fasta from ch_fasta_for_seqdict
    path where_are_my_files

    output:
    path "*.dict" into ch_seq_dict
    path "where_are_my_files.txt"

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
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, path(bam) from ch_input_for_convertbam 

    output:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, path("*fastq.gz"), val('NA') into ch_output_from_convertbam

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
  tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, path(bam) from ch_input_for_indexbam 

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), file("*.{bai,csi}")  into ch_indexbam_for_filtering

  script:
  def size = params.large_ref ? '-c' : ''
  """
  samtools index ${bam} ${size}
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
    publishDir "${params.outdir}/fastqc/input_fastq", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }


    input:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_convertbam_for_fastqc

    output:
    path "*_fastqc.{zip,html}" into ch_prefastqc_for_multiqc

    when: 
    !params.skip_fastqc

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
    publishDir "${params.outdir}/FastP", mode: params.publish_dir_mode

    when: 
    params.complexity_filter_poly_g

    input:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_input_for_fastp.twocol

    output:
    tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, path("*.pG.fq.gz") into ch_output_from_fastp
    path("*.json") into ch_fastp_for_multiqc

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
      def r2 = seqtype == "PE" ? it[9] : file("$projectDir/assets/nf-core_eager_dummy.txt")
      
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
    def r1 = it[8] instanceof ArrayList ? it[8].sort()[0] : it[8]
    def r2 = seqtype == "PE" ? it[8].sort()[1] : file("$projectDir/assets/nf-core_eager_dummy.txt")

    [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

  }
  .set{ ch_fastp_for_merge }

ch_skipfastp_for_merge.mix(ch_fastp_for_merge)
  .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval } 

// Sequencing adapter clipping and optional paired-end merging in preparation for mapping

process adapter_removal {
    label 'mc_small'
    tag "${libraryid}_L${lane}"
    publishDir "${params.outdir}/adapterremoval", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_fastp_for_adapterremoval

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("output/*{combined.fq,.se.truncated,pair1.truncated}.gz") into ch_output_from_adapterremoval_r1
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("output/*pair2.truncated.gz") optional true into ch_output_from_adapterremoval_r2
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("output/*.settings") into ch_adapterremoval_logs

    when: 
    !params.skip_adapterremoval

    script:
    base = "${r1.baseName}_L${lane}"
    //This checks whether we skip trimming and defines a variable respectively
    def trim_me = params.skip_trim ? '' : "--trimns --trimqualities --adapter1 ${params.clip_forward_adaptor} --adapter2 ${params.clip_reverse_adaptor} --minlength ${params.clip_readlength} --minquality ${params.clip_min_read_quality} --minadapteroverlap ${params.min_adap_overlap}"
    def collapse_me = params.skip_collapse ? '' : '--collapse'
    def preserve5p = params.preserve5p ? '--preserve5p' : ''
    def mergedonly = params.mergedonly ? "Y" : "N"
    
    //PE mode, dependent on trim_me and collapse_me the respective procedure is run or not :-) 
    if ( seqtype == 'PE'  && !params.skip_collapse && !params.skip_trim ){
    """
    mkdir -p output
    AdapterRemoval --file1 ${r1} --file2 ${r2} --basename ${base}.pe ${trim_me} --gzip --threads ${task.cpus} ${collapse_me} ${preserve5p}
    
    #Combine files
    if [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "N" ]; then 
      cat *.collapsed.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    elif [ ${preserve5p}  = "--preserve5p" ] && [ ${mergedonly} = "Y" ] ; then
      cat *.collapsed.gz > output/${base}.pe.combined.tmp.fq.gz
    elif [ ${mergedonly} = "Y" ] ; then
      cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    else
      cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    fi

    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus} > output/${base}.pe.combined.fq.gz
   
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
      cat *.collapsed.gz *.collapsed.truncated.gz > output/${base}.pe.combined.tmp.fq.gz
    else
      cat *.collapsed.gz *.collapsed.truncated.gz *.singleton.truncated.gz *.pair1.truncated.gz *.pair2.truncated.gz  > output/${base}.pe.combined.tmp.fq.gz
    fi

    ## Add R_ and L_ for unmerged reads for DeDup compatibility
    AdapterRemovalFixPrefix output/${base}.pe.combined.tmp.fq.gz | pigz -p ${task.cpus} > output/${base}.pe.combined.fq.gz

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
        def r2 = seqtype == "PE" ? file(it[7].sort()[1]) : file("$projectDir/assets/nf-core_eager_dummy.txt")

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
        def r2 = file("$projectDir/assets/nf-core_eager_dummy.txt")

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]
    }
    .into { ch_output_from_adapterremoval; ch_adapterremoval_for_postfastqc }
}

// AdapterRemoval bypass when not running it
if (!params.skip_adapterremoval) {
    ch_output_from_adapterremoval.mix(ch_fastp_for_skipadapterremoval)
        .filter { it =~/.*combined.fq.gz|.*truncated.gz/ }
        .dump(tag: "AR Bypass")
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_lanemerge; } 
} else {
    ch_fastp_for_skipadapterremoval
        .into { ch_adapterremoval_for_fastqc_after_clipping; ch_adapterremoval_for_lanemerge; } 
}

// Lane merging for libraries sequenced over multiple lanes (e.g. NextSeq)
ch_branched_for_lanemerge = ch_adapterremoval_for_lanemerge
  .groupTuple(by: [0,1,3,4,5,6])
  .map {
    it ->
      def samplename = it[0]
      def libraryid  = it[1]
      def lane = it[2]
      def seqtype = it[3]
      def organism = it[4]
      def strandedness = it[5]
      def udg = it[6]
      def r1 = it[7]
      def r2 = it[8]

      [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

  }
  .dump(tag: "LaneMerge Bypass")
  .branch {
    skip_merge: it[7].size() == 1 // Can skip merging if only single lanes
    merge_me: it[7].size() > 1
  }

ch_branched_for_lanemerge_skipme = ch_branched_for_lanemerge.skip_merge
  .map{
    it -> 
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = it[7][0]
        def r2 = it[8][0]

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]
  }
  .dump(tag: "LaneMerge Reconfigure")


ch_branched_for_lanemerge_ready = ch_branched_for_lanemerge.merge_me
  .map{
      it -> 
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = it[7]

        // find and remove duplicate dummies to prevent file collision error
        def r2 = it[8]*.toString()
        r2.removeAll{ it == "$projectDir/assets/nf-core_eager_dummy.txt" }

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]
  }

process lanemerge {
  label 'sc_tiny'
  tag "${libraryid}"
  publishDir "${params.outdir}/lanemerging", mode: params.publish_dir_mode

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(r1), path(r2) from ch_branched_for_lanemerge_ready

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*_R1_lanemerged.fq.gz") into ch_lanemerge_for_mapping_r1
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*_R2_lanemerged.fq.gz") optional true into ch_lanemerge_for_mapping_r2

  script:
  if ( seqtype == 'PE' && ( params.skip_collapse || params.skip_adapterremoval ) ){
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

// Ensuring always valid R2 file even if doesn't exist for AWS
if ( ( params.skip_collapse || params.skip_adapterremoval ) ) {
  ch_lanemerge_for_mapping_r1
    .dump(tag: "Post LaneMerge Reconfigure")
    .mix(ch_lanemerge_for_mapping_r2)
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
        def r2 = seqtype == "PE" ? file(it[7].sort()[1]) : file("$projectDir/assets/nf-core_eager_dummy.txt")

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

    }
    .mix(ch_branched_for_lanemerge_skipme)
    .into { ch_lanemerge_for_skipmap; ch_lanemerge_for_bwa; ch_lanemerge_for_cm; ch_lanemerge_for_bwamem; ch_lanemerge_for_bt2 }
} else {
  ch_lanemerge_for_mapping_r1
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
        def r2 = file("$projectDir/assets/nf-core_eager_dummy.txt")

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]
    }
    .mix(ch_branched_for_lanemerge_skipme)
    .into { ch_lanemerge_for_skipmap; ch_lanemerge_for_bwa; ch_lanemerge_for_cm; ch_lanemerge_for_bwamem; ch_lanemerge_for_bt2 }
}

// ENA upload doesn't do separate lanes, so merge raw FASTQs for mapped-reads removal 

// Per-library lane grouping done within process
process lanemerge_hostremoval_fastq {
  label 'sc_tiny'
  tag "${libraryid}"

  when: 
  params.hostremoval_input_fastq

  input:
  tuple samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_input_for_lanemerge_hostremovalfastq.groupTuple(by: [0,1,3,4,5,6,7])

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.fq.gz") into ch_fastqlanemerge_for_hostremovalfastq

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
    publishDir "${params.outdir}/fastqc/after_clipping", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"
                }


    when: !params.skip_adapterremoval && !params.skip_fastqc

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_adapterremoval_for_fastqc_after_clipping

    output:
    path("*_fastqc.{zip,html}") into ch_fastqc_after_clipping

    script:
    if ( params.skip_collapse && seqtype == 'PE' ) {
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
    publishDir "${params.outdir}/mapping/bwa", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(r1), path(r2) from ch_lanemerge_for_bwa.dump(tag: "input_tuple")
    path index from bwa_index.collect().dump(tag: "input_index")

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.mapped.bam"), path("*.{bai,csi}") into ch_output_from_bwa   

    when: 
    params.mapper == 'bwaaln'

    script:
    def size = params.large_ref ? '-c' : ''
    def fasta = "${index}/${fasta_base}"

    //PE data without merging, PE data without any AR applied
    if ( seqtype == 'PE' && ( params.skip_collapse || params.skip_adapterremoval ) ){
    """
    bwa aln -t ${task.cpus} $fasta ${r1} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r1.sai
    bwa aln -t ${task.cpus} $fasta ${r2} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $fasta ${libraryid}.r1.sai ${libraryid}.r2.sai ${r1} ${r2} | samtools sort -@ ${task.cpus} -O bam - > ${libraryid}_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size}
    """
    } else {
    //PE collapsed, or SE data 
    """
    bwa aln -t ${task.cpus} ${fasta} ${r1} -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $fasta ${libraryid}.sai $r1 | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size}
    """
    }
    
}

// bwa mem for more complex or for modern data mapping

process bwamem {
    label 'mc_medium'
    tag "$libraryid"
    publishDir "${params.outdir}/mapping/bwamem", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_bwamem
    path index from bwa_index_bwamem.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.mapped.bam"), path("*.{bai,csi}") into ch_output_from_bwamem

    when: 
    params.mapper == 'bwamem'

    script:
    def fasta = "${index}/${fasta_base}"
    def size = params.large_ref ? '-c' : ''

    if (!params.single_end && params.skip_collapse){
    """
    bwa mem -t ${task.cpus} $fasta $r1 $r2 -R "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index ${size} -@ ${task.cpus} "${libraryid}".mapped.bam
    """
    } else {
    """
    bwa mem -t ${task.cpus} $fasta $r1 -R "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" | samtools sort -@ ${task.cpus} -O bam - > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index -@ ${task.cpus} "${libraryid}"_"${seqtype}".mapped.bam ${size} 
    """
    }
    
}

// CircularMapper reference preparation and mapping for circular genomes e.g. mtDNA

process circulargenerator{
    label 'sc_tiny'
    tag "$prefix"
    publishDir "${params.outdir}/reference_genome/circularmapper_index", mode: params.publish_dir_mode, saveAs: { filename -> 
            if (params.save_reference) filename 
            else if(!params.save_reference && filename == "where_are_my_files.txt") filename
            else null
    }

    input:
    file fasta from ch_fasta_for_circulargenerator

    output:
    file "${prefix}.{amb,ann,bwt,sa,pac}" into ch_circularmapper_indices
    file "*_elongated" into ch_circularmapper_elongatedfasta

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
    publishDir "${params.outdir}/mapping/circularmapper", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_cm
    file index from ch_circularmapper_indices.collect()
    file fasta from ch_fasta_for_circularmapper.collect()
    file elongated from ch_circularmapper_elongatedfasta.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.mapped.bam"), file("*.{bai,csi}") into ch_output_from_cm

    when: 
    params.mapper == 'circularmapper'

    script:
    def filter = params.circularfilter ? '' : '-f true -x false'
    def elongated_root = "${fasta.baseName}_${params.circularextension}.fasta"
    def size = params.large_ref ? '-c' : ''

    if (!params.single_end && params.skip_collapse ){
    """
    bwa aln -t ${task.cpus} $elongated_root $r1 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r1.sai
    bwa aln -t ${task.cpus} $elongated_root $r2 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.r2.sai
    bwa sampe -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $elongated_root ${libraryid}.r1.sai ${libraryid}.r2.sai $r1 $r2 > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > ${libraryid}_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size} 
    """
    } else {
    """ 
    bwa aln -t ${task.cpus} $elongated_root $r1 -n ${params.bwaalnn} -l ${params.bwaalnl} -k ${params.bwaalnk} -f ${libraryid}.sai
    bwa samse -r "@RG\\tID:ILLUMINA-${libraryid}\\tSM:${libraryid}\\tPL:illumina" $elongated_root ${libraryid}.sai $r1 > tmp.out
    realignsamfile -e ${params.circularextension} -i tmp.out -r $fasta $filter 
    samtools sort -@ ${task.cpus} -O bam tmp_realigned.bam > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size}
    """
    }
    
}

process bowtie2 {
    label 'mc_medium'
    tag "${libraryid}"
    publishDir "${params.outdir}/mapping/bt2", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(r1), file(r2) from ch_lanemerge_for_bt2
    path index from bt2_index.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.mapped.bam"), path("*.{bai,csi}") into ch_output_from_bt2
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*_bt2.log") into ch_bt2_for_multiqc

    when: 
    params.mapper == 'bowtie2'

    script:
    def size = params.large_ref ? '-c' : ''
    def fasta = "${index}/${fasta_base}"
    def trim5 = params.bt2_trim5 != 0 ? "--trim5 ${params.bt2_trim5}" : ""
    def trim3 = params.bt2_trim3 != 0 ? "--trim3 ${params.bt2_trim3}" : ""
    def bt2n = params.bt2n != 0 ? "-N ${params.bt2n}" : ""
    def bt2l = params.bt2l != 0 ? "-L ${params.bt2l}" : ""

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
    bowtie2 -x ${fasta} -1 ${r1} -2 ${r2} -p ${task.cpus} ${sensitivity} ${bt2n} ${bt2l} ${trim5} ${trim3} 2> "${libraryid}"_bt2.log | samtools sort -@ ${task.cpus} -O bam > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size}
    """
    } else {
    //PE collapsed, or SE data 
    """
    bowtie2 -x ${fasta} -U ${r1} -p ${task.cpus} ${sensitivity} ${bt2n} ${bt2l} ${trim5} ${trim3} 2> "${libraryid}"_bt2.log | samtools sort -@ ${task.cpus} -O bam > "${libraryid}"_"${seqtype}".mapped.bam
    samtools index "${libraryid}"_"${seqtype}".mapped.bam ${size}
    """
    }
    
}

// Gather all mapped BAMs from all possible mappers into common channels to send downstream
ch_output_from_bwa.mix(ch_output_from_bwamem, ch_output_from_cm, ch_indexbam_for_filtering, ch_output_from_bt2)
  .into { ch_mapping_for_hostremovalfastq; ch_mapping_for_seqtype_merging }

// Synchronise the mapped input FASTQ and input non-remapped BAM channels
ch_fastqlanemerge_for_hostremovalfastq
    .map {
        def samplename = it[0]
        def libraryid  = it[1]
        def lane = it[2]
        def seqtype = it[3]
        def organism = it[4]
        def strandedness = it[5]
        def udg = it[6]
        def r1 = seqtype == "PE" ? file(it[7].sort()[0]) : file(it[7])
        def r2 = seqtype == "PE" ? file(it[7].sort()[1]) : file("$projectDir/assets/nf-core_eager_dummy.txt")

        [ samplename, libraryid, lane, seqtype, organism, strandedness, udg, r1, r2 ]

    }
    .mix(ch_mapping_for_hostremovalfastq)
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
    .set { ch_synced_for_hostremovalfastq }

// Remove mapped reads from original (lane merged) input FASTQ e.g. for sensitive host data when running metagenomic data

process hostremoval_input_fastq {
    label 'mc_medium'
    tag "${libraryid}"
    publishDir "${params.outdir}/hostremoved_fastq", mode: params.publish_dir_mode

    when: 
    params.hostremoval_input_fastq

    input: 
    tuple samplename, libraryid, seqtype, organism, strandedness, udg, file(r1), file(r2), file(bam), file(bai) from ch_synced_for_hostremovalfastq

    output:
    tuple samplename, libraryid, seqtype, organism, strandedness, udg, file("*.fq.gz") into ch_output_from_hostremovalfastq

    script:
    if ( seqtype == 'SE' ) {
        out_fwd = bam.baseName+'.hostremoved.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${r1} -m ${params.hostremoval_mode} -of $out_fwd -p ${task.cpus}
        """
    } else {
        out_fwd = bam.baseName+'.hostremoved.fwd.fq.gz'
        out_rev = bam.baseName+'.hostremoved.rev.fq.gz'
        """
        samtools index $bam
        extract_map_reads.py $bam ${r1} -rev ${r2} -m  ${params.hostremoval_mode} -of $out_fwd -or $out_rev -p ${task.cpus}
        """ 
    }
    
}

// Seqtype merging to combine paired end with single end  sequenceing data of the same libraries
// goes here, goes into flagstat, filter etc. Important: This type of merge of this isn't technically valid for DeDup!
// and should only be used with markduplicates!
ch_branched_for_seqtypemerge = ch_mapping_for_seqtype_merging
  .groupTuple(by: [0,1,4,5,6])
  .map {
    it ->
      def samplename = it[0]
      def libraryid  = it[1]
      def lane = it[2]
      def seqtype = it[3].unique() // How to deal with this?
      def organism = it[4]
      def strandedness = it[5]
      def udg = it[6]
      def r1 = it[7]
      def r2 = it[8]

      // We will assume if mixing it is better to set as PE as this is informative
      // for DeDup (and markduplicates doesn't care), but will throw a warning!
      def seqtype_new = seqtype.flatten().size() > 1 ? 'PE' : seqtype 
                      
      if ( seqtype.flatten().size() > 1 &&  params.dedupper == 'dedup' ) {
        log.warn "[nf-core/eager] Warning: you are running DeDup on BAMs with a mixture of PE/SE data for library: ${libraryid}. DeDup is designed for PE data only, deduplication maybe suboptimal!"
      }
      
      [ samplename, libraryid, lane, seqtype_new, organism, strandedness, udg, r1, r2 ]

  }
  .dump(tag: "Seqtype")
  .branch {
    skip_merge: it[7].size() == 1 // Can skip merging if only single lanes
    merge_me: it[7].size() > 1
  }

  process seqtype_merge {

    label 'sc_tiny'
    tag "$libraryid"

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_branched_for_seqtypemerge.merge_me

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*_seqtypemerged_rg.bam"), file("*_seqtypemerged_rg*.{bai,csi}")  into ch_seqtypemerge_for_filtering

    script:
    def size = params.large_ref ? '-c' : ''
    """
    samtools merge ${libraryid}_seqtypemerged.bam ${bam}
    ## Have to set validation as lenient because of BWA issue: "I see a read stands out the end of a chromosome and is flagged as unmapped (flag 0x4). [...]" http://bio-bwa.sourceforge.net/
    picard AddOrReplaceReadGroups I=${libraryid}_seqtypemerged.bam O=${libraryid}_seqtypemerged_rg.bam RGID=1 RGLB="${libraryid}_seqtypemerged" RGPL=illumina RGPU=4410 RGSM="${libraryid}_seqtypemerged" VALIDATION_STRINGENCY=LENIENT
    samtools index ${libraryid}_seqtypemerged_rg.bam ${size}
    """
    
  }

ch_seqtypemerge_for_filtering
  .mix(ch_branched_for_seqtypemerge.skip_merge)
  .into { ch_seqtypemerged_for_skipfiltering; ch_seqtypemerged_for_samtools_filter; ch_seqtypemerged_for_samtools_flagstat } 

// Post-mapping QC

process samtools_flagstat {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/samtools/stats", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_seqtypemerged_for_samtools_flagstat


    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*stats") into ch_flagstat_for_multiqc,ch_flagstat_for_endorspy

    script:
    """
    samtools flagstat $bam > ${libraryid}_flagstat.stats
    """
}


// BAM filtering e.g. to extract unmapped reads for downstream or stricter mapping quality

process samtools_filter {
    label 'mc_medium'
    tag "$libraryid"
    publishDir "${params.outdir}/samtools/filter", mode: params.publish_dir_mode,
    saveAs: {filename ->
            if (filename.indexOf(".fq.gz") > 0) "unmapped/$filename"
            else if (filename.indexOf(".unmapped.bam") > 0) "unmapped/$filename"
            else if (filename.indexOf(".filtered.bam")) filename
            else null
    }

    when: 
    params.run_bam_filtering

    input: 
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_seqtypemerged_for_samtools_filter

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*filtered.bam"), file("*.{bai,csi}") into ch_output_from_filtering
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.unmapped.fastq.gz") optional true into ch_bam_filtering_for_metagenomic
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.unmapped.bam") optional true

    // Using shell block rather than script because we are playing with awk
    shell:
    size = !{params.large_ref} ? '-c' : ''
    
    if ( "${params.bam_unmapped_type}" == "keep" ) {
        '''
        ## Unmapped and MAPQ filtering
        samtools view -h -b !{bam} -@ !{task.cpus} -q !{params.bam_mapping_quality_threshold} -o tmp_mapped.bam

        ## Mapped LEN filtering
        if [[ !{params.bam_filter_minreadlength} -eq 0 ]]; then
            mv tmp_mapped.bam !{libraryid}.filtered.bam
        else
            filter_bam_fragment_length.py -a -l !{params.bam_filter_minreadlength} -o !{libraryid} tmp_mapped.bam
        fi

        samtools index !{libraryid}.filtered.bam !{size}
        '''
    } else if("${params.bam_unmapped_type}" == "discard"){
        '''
        ## Unmapped and MAPQ filtering
        samtools view -h -b !{bam} -@ !{task.cpus} -F4 -q !{params.bam_mapping_quality_threshold} -o tmp_mapped.bam

        ## Mapped LEN filtering
        if [[ !{params.bam_filter_minreadlength} -eq 0 ]]; then
            mv tmp_mapped.bam !{libraryid}.filtered.bam
        else
            filter_bam_fragment_length.py -a -l !{params.bam_filter_minreadlength} -o !{libraryid} tmp_mapped.bam
        fi

        samtools index !{libraryid}.filtered.bam !{size}
        '''
    } else if("${params.bam_unmapped_type}" == "bam"){
        '''
        ## Unmapped and MAPQ filtering
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -f4 -o !{libraryid}.unmapped.bam
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -F4 -q !{params.bam_mapping_quality_threshold} -o tmp_mapped.bam

        ## Mapped LEN filtering
        if [[ !{params.bam_filter_minreadlength} -eq 0 ]]; then
            mv tmp_mapped.bam !{libraryid}.filtered.bam
        else
            filter_bam_fragment_length.py -a -l !{params.bam_filter_minreadlength} -o !{libraryid} tmp_mapped.bam
        fi

        samtools index !{libraryid}.filtered.bam !{size}
        '''
    } else if("${params.bam_unmapped_type}" == "fastq"){
        '''
        ## Unmapped and MAPQ filtering
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -f4 -o !{libraryid}.unmapped.bam
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -F4 -q !{params.bam_mapping_quality_threshold} -o tmp_mapped.bam

        ## Mapped LEN filtering
        if [[ !{params.bam_filter_minreadlength} -eq 0 ]]; then
            mv tmp_mapped.bam !{libraryid}.filtered.bam
        else
            filter_bam_fragment_length.py -a -l !{params.bam_filter_minreadlength} -o !{libraryid} tmp_mapped.bam
        fi

        samtools index !{libraryid}.filtered.bam !{size}

        ## FASTQ
        samtools fastq -tn !{libraryid}.unmapped.bam | pigz -p !{task.cpus} > !{libraryid}.unmapped.fastq.gz
        rm !{libraryid}.unmapped.bam
        '''
    } else if("${params.bam_unmapped_type}" == "both"){
        '''
        ## Unmapped and MAPQ filtering
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -f4 -o !{libraryid}.unmapped.bam
        samtools view -h !{bam} | samtools view - -@ !{task.cpus} -F4 -q !{params.bam_mapping_quality_threshold} -o tmp_mapped.bam

        ## Mapped LEN filtering
        if [[ !{params.bam_filter_minreadlength} -eq 0 ]]; then
            mv tmp_mapped.bam !{libraryid}.filtered.bam
        else
            filter_bam_fragment_length.py -a -l !{params.bam_filter_minreadlength} -o !{libraryid} tmp_mapped.bam
        fi

        samtools index !{libraryid}.filtered.bam !{size}
        samtools fastq -tn !{libraryid}.unmapped.bam | pigz -p !{task.cpus} > !{libraryid}.unmapped.fastq.gz
        '''
    }
}

// samtools_filter bypass in case not run
if (params.run_bam_filtering) {
    ch_seqtypemerged_for_skipfiltering.mix(ch_output_from_filtering)
        .filter { it =~/.*filtered.bam/ }
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_flagstat; ch_skiprmdup_for_libeval; ch_mapped_for_preseq } 

} else {
    ch_seqtypemerged_for_skipfiltering
        .into { ch_filtering_for_skiprmdup; ch_filtering_for_dedup; ch_filtering_for_markdup; ch_filtering_for_flagstat; ch_skiprmdup_for_libeval; ch_mapped_for_preseq } 

}

// Post filtering mapping QC - particularly to help see how much was removed from mapping quality filtering

process samtools_flagstat_after_filter {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/samtools/filtered_stats", mode: params.publish_dir_mode

    when:
    params.run_bam_filtering

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_filtering_for_flagstat

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.stats") into ch_bam_filtered_flagstat_for_multiqc, ch_bam_filtered_flagstat_for_endorspy

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
        def poststats = file("$projectDir/assets/nf-core_eager_dummy.txt")

      [samplename, libraryid, lane, seqtype, organism, strandedness, udg, stats, poststats ] }
    .set{ ch_allflagstats_for_endorspy }
}

// Endogenous DNA calculator to say how much of a library contained 'on-target' DNA

process endorSpy {
    label 'sc_tiny'
    tag "$libraryid"
    publishDir "${params.outdir}/endorspy", mode: params.publish_dir_mode

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(stats), path(poststats) from ch_allflagstats_for_endorspy

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_endorspy_for_multiqc

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
    publishDir "${params.outdir}/deduplication/", mode: params.publish_dir_mode,
        saveAs: {filename -> "${libraryid}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'dedup'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_filtering_for_dedup

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.hist") into ch_hist_for_preseq
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_dedup_results_for_multiqc
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${libraryid}_rmdup.bam"), path("*.{bai,csi}") into ch_output_from_dedup, ch_dedup_for_libeval

    script:
    def outname = "${bam.baseName}"
    def treat_merged = params.dedup_all_merged ? '-m' : ''
    def size = params.large_ref ? '-c' : ''
    
    """
    ## To make sure direct BAMs have a clean name
    if [[ "${bam}" != "${libraryid}.bam" ]]; then
      mv ${bam} ${libraryid}.bam
    fi
    
    dedup -Xmx${task.memory.toGiga()}g -i ${libraryid}.bam $treat_merged -o . -u 
    mv *.log dedup.log
    samtools sort -@ ${task.cpus} "${libraryid}"_rmdup.bam -o "${libraryid}"_rmdup.bam
    samtools index "${libraryid}"_rmdup.bam ${size}
    """
}

process markduplicates{
    label 'mc_small'
    tag "${libraryid}"
    publishDir "${params.outdir}/deduplication/", mode: params.publish_dir_mode,
        saveAs: {filename -> "${libraryid}/$filename"}

    when:
    !params.skip_deduplication && params.dedupper == 'markduplicates'

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_filtering_for_markdup

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.metrics") into ch_markdup_results_for_multiqc
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${libraryid}_rmdup.bam"), path("*.{bai,csi}") into ch_output_from_markdup, ch_markdup_for_libeval

    script:
    def outname = "${bam.baseName}"
    def size = params.large_ref ? '-c' : ''
    """
    ## To make sure direct BAMs have a clean name
    if [[ "${bam}" != "${libraryid}.bam" ]]; then
      mv ${bam} ${libraryid}.bam
    fi

    picard -Xmx${task.memory.toMega()}M -Xms${task.memory.toMega()}M MarkDuplicates INPUT=${libraryid}.bam OUTPUT=${libraryid}_rmdup.bam REMOVE_DUPLICATES=TRUE AS=TRUE METRICS_FILE="${libraryid}_rmdup.metrics" VALIDATION_STRINGENCY=SILENT
    samtools index ${libraryid}_rmdup.bam ${size}
    """
}

// This is for post-deduplcation per-library evaluation steps _without_ any 
// form of library merging. 
if ( params.skip_deduplication ) {
  ch_skiprmdup_for_libeval.mix(ch_dedup_for_libeval, ch_markdup_for_libeval)
    .into{ ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_for_nuclear_contamination; ch_rmdup_formtnucratio }
} else {
  ch_dedup_for_libeval.mix(ch_markdup_for_libeval)
    .into{ ch_rmdup_for_preseq; ch_rmdup_for_damageprofiler; ch_for_nuclear_contamination; ch_rmdup_formtnucratio }
}

// Merge independent libraries sequenced but with same treatment (often done to 
// improve complexity) with the same _sample_ name. Different strand/UDG libs 
// not merged because bamtrim/pmdtools/genotyping needs that info.

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
      def bam = it[7].flatten()
      def bai = it[8].flatten()

      [it[0], libraryid, it[2], it[3], it[4], it[5], it[6], bam, bai ]
    }
  .set { ch_input_for_skiplibrarymerging }

ch_input_for_librarymerging.merge_me
  .map{
    it ->
      def libraryid = it[1][0]
      def seqtype = "merged"
      def bam = it[7].flatten()
      def bai = it[8].flatten()

      [it[0], libraryid, it[2], seqtype, it[4], it[5], it[6], bam, bai ]
    }
  .set { ch_fixedinput_for_librarymerging }

process library_merge {
  label 'sc_tiny'
  tag "${samplename}"
  publishDir "${params.outdir}/merged_bams/initial", mode: params.publish_dir_mode

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_fixedinput_for_librarymerging.dump(tag: "Input Tuple Library Merge")

  output:
  tuple samplename, val("${samplename}_libmerged"), lane, seqtype, organism, strandedness, udg, path("*_libmerged_rg_rmdup.bam"), path("*_libmerged_rg_rmdup.bam.{bai,csi}") into ch_output_from_librarymerging

  script:
  def size = params.large_ref ? '-c' : ''
  """
  samtools merge ${samplename}_libmerged_rmdup.bam ${bam}
  ## Have to set validation as lenient because of BWA issue: "I see a read stands out the end of a chromosome and is flagged as unmapped (flag 0x4). [...]" http://bio-bwa.sourceforge.net/
  picard AddOrReplaceReadGroups I=${samplename}_libmerged_rmdup.bam O=${samplename}_libmerged_rg_rmdup.bam RGID=1 RGLB="${samplename}_merged" RGPL=illumina RGPU=4410 RGSM="${samplename}_merged" VALIDATION_STRINGENCY=LENIENT
  samtools index ${samplename}_libmerged_rg_rmdup.bam ${size}
  """
}

// Mix back in libraries from skipping dedup, skipping library merging
if (!params.skip_deduplication) {
    ch_input_for_skiplibrarymerging.mix(ch_output_from_librarymerging)
        .filter { it =~/.*_rmdup.bam/ }
        .into { ch_rmdup_for_skipdamagemanipulation;  ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_rmdup_for_bedtools } 

} else {
    ch_input_for_skiplibrarymerging.mix(ch_output_from_librarymerging)
        .into { ch_rmdup_for_skipdamagemanipulation; ch_rmdup_for_pmdtools; ch_rmdup_for_bamutils; ch_rmdup_for_bedtools } 
}

//////////////////////////////////////////////////
/* --     POST DEDUPLICATION EVALUATION      -- */
//////////////////////////////////////////////////

// Library complexity calculation from mapped reads - could a user cost-effectively sequence deeper for more unique information?
if ( params.skip_deduplication ) {
  ch_input_for_preseq = ch_rmdup_for_preseq.map{ it[0,1,2,3,4,5,6,7] }

} else if ( !params.skip_deduplication && params.dedupper == "markduplicates" ) {
  ch_input_for_preseq = ch_mapped_for_preseq.map{ it[0,1,2,3,4,5,6,7] }

} else if ( !params.skip_deduplication && params.dedupper == "dedup" ) {
  ch_input_for_preseq = ch_hist_for_preseq

}

process preseq {
    label 'sc_tiny'
    tag "${libraryid}"
    publishDir "${params.outdir}/preseq", mode: params.publish_dir_mode

    when:
    !params.skip_preseq

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(input) from ch_input_for_preseq

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${input.baseName}.ccurve") into ch_preseq_for_multiqc

    script:
    pe_mode = params.skip_collapse && seqtype == "PE" ? '-P' : ''
    if(!params.skip_deduplication && params.dedupper == "dedup"){
    """
    preseq c_curve -s ${params.preseq_step_size} -o ${input.baseName}.ccurve -H ${input}
    """
    } else if( !params.skip_deduplication && params.dedupper == "markduplicates"){
    """
    preseq c_curve -s ${params.preseq_step_size} -o ${input.baseName}.ccurve -B ${input} ${pe_mode}
    """
    } else if ( params.skip_deduplication ) {
    """
    preseq c_curve -s ${params.preseq_step_size} -o ${input.baseName}.ccurve -B ${input} ${pe_mode}
    """
    }
}

// Optional mapping statistics for specific annotations - e.g. genes in bacterial genome

process bedtools {
  label 'mc_small'
  tag "${libraryid}"
  publishDir "${params.outdir}/bedtools", mode: params.publish_dir_mode

  when:
  params.run_bedtools_coverage

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_rmdup_for_bedtools
  file anno_file from ch_anno_for_bedtools.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*")

  script:
  """
  bedtools coverage -nonamecheck -a ${anno_file} -b $bam | pigz -p ${task.cpus} > "${bam.baseName}".breadth.gz
  bedtools coverage -nonamecheck -a ${anno_file} -b $bam -mean | pigz -p ${task.cpus} > "${bam.baseName}".depth.gz
  """
}

//////////////////////////////////////////////////////////////
/* --    ANCIENT DNA EVALUATION AND BAM MODIFICATION     -- */
//////////////////////////////////////////////////////////////

// Calculate typical aDNA damage frequency distribution

process damageprofiler {
    label 'sc_small'
    tag "${libraryid}"

    publishDir "${params.outdir}/damageprofiler", mode: params.publish_dir_mode

    when:
    !params.skip_damage_calculation

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_rmdup_for_damageprofiler
    file fasta from ch_fasta_for_damageprofiler.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${base}/*.txt") optional true
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${base}/*.log")
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${base}/*.pdf") optional true
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("${base}/*.json") optional true into ch_damageprofiler_results

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
    publishDir "${params.outdir}/pmdtools", mode: params.publish_dir_mode

    when: params.run_pmdtools

    input: 
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_rmdup_for_pmdtools
    file fasta from ch_fasta_for_pmdtools.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.pmd.bam"), path("*.pmd.bam.{bai,csi}") into ch_output_from_pmdtools
    file "*.cpg.range*.txt"

    script:
    //Check which treatment for the libraries was used
    def treatment = udg ? (udg == 'half' ? '--UDGhalf' : '--CpG') : '--UDGminus'
    if(params.snpcapture_bed != ''){
        snpcap = (params.pmdtools_reference_mask != '') ? "--refseq ${params.pmdtools_reference_mask}" : ''
        log.info"######No reference mask specified for PMDtools, therefore ignoring that for downstream analysis!"
    } else {
        snpcap = ''
    }
    def size = params.large_ref ? '-c' : ''
    """
    #Run Filtering step 
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --threshold ${params.pmdtools_threshold} $treatment $snpcap --header | samtools view -@ ${task.cpus} -Sb - > "${libraryid}".pmd.bam
    
    #Run Calc Range step
    ## To allow early shut off of pipe: https://github.com/nextflow-io/nextflow/issues/1564
    trap 'if [[ \$? == 141 ]]; then echo "Shutting samtools early due to -n parameter" && samtools index ${libraryid}.pmd.bam ${size}; exit 0; fi' EXIT
    samtools calmd -b $bam $fasta | samtools view -h - | pmdtools --deamination --range ${params.pmdtools_range} $treatment $snpcap -n ${params.pmdtools_max_reads} > "${libraryid}".cpg.range."${params.pmdtools_range}".txt
    
    echo "Running indexing"
    samtools index ${libraryid}.pmd.bam ${size}
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
    publishDir "${params.outdir}/trimmed_bam", mode: params.publish_dir_mode

    when: params.run_trim_bam

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_bamutils_decision.totrim

    output: 
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.trimmed.bam"), file("*.trimmed.bam.{bai,csi}") into ch_trimmed_from_bamutils

    script:
    def softclip = params.bamutils_softclip ? '-c' : '' 
    def size = params.large_ref ? '-c' : ''
    def left_clipping = udg == "half" ? "${params.bamutils_clip_half_udg_left}" : "${params.bamutils_clip_none_udg_left}"
    def right_clipping = udg == "half" ? "${params.bamutils_clip_half_udg_right}" : "${params.bamutils_clip_none_udg_right}"
    """
    bam trimBam $bam tmp.bam -L ${left_clipping} -R ${right_clipping} ${softclip}
    samtools sort -@ ${task.cpus} tmp.bam -o ${libraryid}.trimmed.bam 
    samtools index ${libraryid}.trimmed.bam ${size}
    """
}

// Post trimming merging of libraries to single samples, except for SS/DS 
// libraries as they should be genotyped separately, because we will assume 
// that if trimming is turned on, 'lab-removed' libraries can be combined with 
// merged with 'in-silico damage removed' libraries to improve genotyping

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
  label 'sc_tiny'
  tag "${samplename}"
  publishDir "${params.outdir}/merged_bams/additional", mode: params.publish_dir_mode

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_trimmed_formerge.merge_me

  output:
  tuple samplename, val("${samplename}_libmerged"), lane, seqtype, organism, strandedness, udg, path("*_libmerged_rg_add.bam"), path("*_libmerged_rg_add.bam.{bai,csi}") into ch_output_from_trimmerge

  script:
  def size = params.large_ref ? '-c' : ''
  """
  samtools merge ${samplename}_libmerged_add.bam ${bam}
  picard AddOrReplaceReadGroups I=${samplename}_libmerged_add.bam O=${samplename}_libmerged_rg_add.bam RGID=1 RGLB="${samplename}_additionalmerged" RGPL=illumina RGPU=4410 RGSM="${samplename}_additionalmerged" VALIDATION_STRINGENCY=LENIENT
  samtools index ${samplename}_libmerged_rg_add.bam ${size}
  """
}

ch_trimmed_formerge.skip_merging
  .mix(ch_output_from_trimmerge)
  .into{ ch_output_from_bamutils; ch_addlibmerge_for_qualimap; ch_for_sexdeterrmine }

  // General mapping quality statistics for whole reference sequence - e.g. X and % coverage

process qualimap {
    label 'mc_small'
    tag "${samplename}"
    publishDir "${params.outdir}/qualimap", mode: params.publish_dir_mode

    when:
    !params.skip_qualimap

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_addlibmerge_for_qualimap
    file fasta from ch_fasta_for_qualimap.collect()

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*") into ch_qualimap_results

    script:
    def snpcap = params.snpcapture_bed != '' ? "-gff ${params.snpcapture_bed}" : ''
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

} else if ( params.run_genotyping && params.genotyping_source == "pmd" && !params.run_pmdtools )  {
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

process genotyping_ug {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode, pattern: '*{.vcf.gz,.realign.bam,realign.bai}'

  when:
  params.run_genotyping && params.genotyping_tool == 'ug'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_ug
  file fasta from ch_fasta_for_genotyping_ug.collect()
  file fai from ch_fai_for_ug.collect()
  file dict from ch_dict_for_ug.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*vcf.gz") into ch_ug_for_multivcfanalyzer,ch_ug_for_vcf2genome
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file("*.realign.{bam,bai}") optional true

  script:
  def defaultbasequalities = params.gatk_ug_defaultbasequalities == '' ? '' : " --defaultBaseQualities ${params.gatk_ug_defaultbasequalities}" 
  def keep_realign = params.gatk_ug_keep_realign_bam ? "samtools index ${samplename}.realign.bam" : "rm ${samplename}.realign.{bam,bai}"
  if (params.gatk_dbsnp == '')
    """
    samtools index -b ${bam}
    gatk3 -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    gatk3 -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplename}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    gatk3 -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}
    
    $keep_realign
    
    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
  else if (params.gatk_dbsnp != '')
    """
    samtools index ${bam}
    gatk3 -T RealignerTargetCreator -R ${fasta} -I ${bam} -nt ${task.cpus} -o ${samplename}.intervals ${defaultbasequalities}
    gatk3 -T IndelRealigner -R ${fasta} -I ${bam} -targetIntervals ${samplenane}.intervals -o ${samplename}.realign.bam ${defaultbasequalities}
    gatk3 -T UnifiedGenotyper -R ${fasta} -I ${samplename}.realign.bam -o ${samplename}.unifiedgenotyper.vcf -nt ${task.cpus} --dbsnp ${params.gatk_dbsnp} --genotype_likelihoods_model ${params.gatk_ug_genotype_model} -stand_call_conf ${params.gatk_call_conf} --sample_ploidy ${params.gatk_ploidy} -dcov ${params.gatk_downsample} --output_mode ${params.gatk_ug_out_mode} ${defaultbasequalities}
    
    $keep_realign
    
    pigz -p ${task.cpus} ${samplename}.unifiedgenotyper.vcf
    """
}

 // HaplotypeCaller as 'best practise' tool for human DNA in particular 

process genotyping_hc {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'hc'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_hc
  file fasta from ch_fasta_for_genotyping_hc.collect()
  file fai from ch_fai_for_hc.collect()
  file dict from ch_dict_for_hc.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*vcf.gz")

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
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'freebayes'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, file(bam), file(bai) from ch_damagemanipulation_for_genotyping_freebayes
  file fasta from ch_fasta_for_genotyping_freebayes.collect()
  file fai from ch_fai_for_freebayes.collect()
  file dict from ch_dict_for_freebayes.collect()

  output: 
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*vcf.gz")
  
  script:
  def skip_coverage = "${params.freebayes_g}" == 0 ? "" : "-g ${params.freebayes_g}"
  """
  freebayes -f ${fasta} -p ${params.freebayes_p} -C ${params.freebayes_C} ${skip_coverage} ${bam} > ${samplename}.freebayes.vcf
  pigz -p ${task.cpus} ${samplename}.freebayes.vcf
  """
}


 // Branch channel by strandedness
ch_damagemanipulation_for_genotyping_pileupcaller
  .branch{
      singleStranded: it[5] == "single"
      doubleStranded: it[5] == "double"
  }
  .set{ch_input_for_genotyping_pileupcaller}

 // Create pileupcaller input tuples
ch_input_for_genotyping_pileupcaller.singleStranded
  .groupTuple(by:[5])
  .set {ch_prepped_for_pileupcaller_single}

ch_input_for_genotyping_pileupcaller.doubleStranded
  .groupTuple(by:[5])
  .set {ch_prepped_for_pileupcaller_double}

process genotyping_pileupcaller {
  label 'mc_small'
  tag "${strandedness}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

  when:
  params.run_genotyping && params.genotyping_tool == 'pileupcaller'

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, bam, bai from ch_prepped_for_pileupcaller_double.mix(ch_prepped_for_pileupcaller_single)
  file fasta from ch_fasta_for_genotyping_pileupcaller.collect()
  file fai from ch_fai_for_pileupcaller.collect()
  file dict from ch_dict_for_pileupcaller.collect()
  path(bed) from ch_bed_for_pileupcaller.collect()
  path(snp) from ch_snp_for_pileupcaller.collect().dump(tag: "Pileupcaller SNP file")

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("pileupcaller.${strandedness}.*") into ch_for_eigenstrat_snp_coverage

  script:
  def use_bed = bed.getName() != 'nf-core_eager_dummy.txt' ? "-l ${bed}" : ''
  def use_snp = snp.getName() != 'nf-core_eager_dummy2.txt' ? "-f ${snp}" : ''

  def transitions_mode = strandedness == "single" ? "" : "${params.pileupcaller_transitions_mode}" == 'SkipTransitions' ? "--skipTransitions" : "${params.pileupcaller_transitions_mode}" == 'TransitionsMissing' ? "--transitionsMissing" : ""
  def caller = "--${params.pileupcaller_method}"
  def ssmode = strandedness == "single" ? "--singleStrandMode" : ""
  def bam_list = bam.flatten().join(" ")
  def sample_names = samplename.flatten().join(",")
  """
  samtools mpileup -B -q 30 -Q 30 ${use_bed} -f ${fasta} ${bam_list} | pileupCaller ${caller} ${ssmode} ${transitions_mode} --sampleNames ${sample_names} ${use_snp} -e pileupcaller.${strandedness}
  """
 }
 
process eigenstrat_snp_coverage {
   label 'mc_tiny'
   tag "${strandedness}"
   publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode
   
   when:
   params.run_genotyping && params.genotyping_tool == 'pileupcaller'
   
   input:
   tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*") from ch_for_eigenstrat_snp_coverage.dump()
   
   output:
   tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_eigenstrat_snp_cov_for_multiqc
   path("*_eigenstrat_coverage.txt")
   
   script:
   /* 
   The following code block can be swapped in once the eigenstratdatabasetools MultiQC module becomes available.
   """
   eigenstrat_snp_coverage -i pileupcaller.${strandedness} -s ".txt" >${strandedness}_eigenstrat_coverage.txt -j ${strandedness}_eigenstrat_coverage_mqc.json
   """
   */
   """
   eigenstrat_snp_coverage -i pileupcaller.${strandedness} -s ".txt" >${strandedness}_eigenstrat_coverage.txt
   parse_snp_cov.py ${strandedness}_eigenstrat_coverage.txt
   """
 }
 
process genotyping_angsd {
  label 'mc_small'
  tag "${samplename}"
  publishDir "${params.outdir}/genotyping", mode: params.publish_dir_mode

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
  
  def angsd_fasta = !params.angsd_createfasta ? '' : params.angsd_fastamethod == 'random' ? '-doFasta 1 -doCounts 1' : '-doFasta 2 -doCounts 1' 
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
  publishDir "${params.outdir}/consensus_sequence", mode: params.publish_dir_mode

  when: 
  params.run_vcf2genome

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(vcf) from ch_ug_for_vcf2genome
  file fasta from ch_fasta_for_vcf2genome.collect()

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.fasta.gz")

  script:
  def out = "${params.vcf2genome_outfile}" == '' ? "${samplename}.fasta" : "${params.vcf2genome_outfile}"
  def fasta_head = "${params.vcf2genome_header}" == '' ? "${samplename}" : "${params.vcf2genome_header}"
  """
  pigz -f -d -p ${task.cpus} *.vcf.gz
  vcf2genome -draft ${out}.fasta -draftname "${fasta_head}" -in ${vcf.baseName} -minc ${params.vcf2genome_minc} -minfreq ${params.vcf2genome_minfreq} -minq ${params.vcf2genome_minq} -ref ${fasta} -refMod ${out}_refmod.fasta -uncertain ${out}_uncertainy.fasta
  pigz -p ${task.cpus} *.fasta 
  pigz -p ${task.cpus} *.vcf
  """
}

// More complex consensus caller with additional filtering functionality (e.g. for heterozygous calls) to generate SNP tables and other things sometimes used in aDNA bacteria studies

// Create input channel for MultiVCFAnalyzer, possibly mixing with pre-made VCFs.
if (params.additional_vcf_files == '') {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it[7] }.collect()
} else {
    ch_vcfs_for_multivcfanalyzer = ch_ug_for_multivcfanalyzer.map{ it [7] }.collect().mix(ch_extravcfs_for_multivcfanalyzer)
}

process multivcfanalyzer {
  label  'mc_small'
  publishDir "${params.outdir}/multivcfanalyzer", mode: params.publish_dir_mode

  when:
  params.genotyping_tool == 'ug' && params.run_multivcfanalyzer && params.gatk_ploidy == '2'

  input:
  file vcf from ch_vcfs_for_multivcfanalyzer.collect()
  file fasta from ch_fasta_for_multivcfanalyzer.collect()

  output:
  file('fullAlignment.fasta.gz')
  file('info.txt.gz')
  file('snpAlignment.fasta.gz')
  file('snpAlignmentIncludingRefGenome.fasta.gz')
  file('snpStatistics.tsv.gz')
  file('snpTable.tsv.gz')
  file('snpTableForSnpEff.tsv.gz')
  file('snpTableWithUncertaintyCalls.tsv.gz')
  file('structureGenotypes.tsv.gz')
  file('structureGenotypes_noMissingData-Columns.tsv.gz')
  file('MultiVCFAnalyzer.json') optional true into ch_multivcfanalyzer_for_multiqc

  script:
  def write_freqs = params.write_allele_frequencies ? "T" : "F"
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
  publishDir "${params.outdir}/mtnucratio", mode: params.publish_dir_mode

  when: 
  params.run_mtnucratio

  input:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(bam), path(bai) from ch_rmdup_formtnucratio

  output:
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.mtnucratio")
  tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path("*.json") into ch_mtnucratio_for_multiqc

  script:
  """
  mtnucratio ${bam} "${params.mtnucratio_header}"
  """
 }

// Human biological sex estimation

// As we collect all files for a single sex_deterrmine run, we DO NOT use the normal input/output tuple
process sex_deterrmine {
    label 'sc_small'
    publishDir "${params.outdir}/sex_determination", mode: params.publish_dir_mode
     
    input:
    path bam from ch_for_sexdeterrmine.map { it[7] }.collect()
    path(bed) from ch_bed_for_sexdeterrmine

    output:
    file "SexDet.txt"
    file "*.json" into ch_sexdet_for_multiqc

    when:
    params.run_sexdeterrmine
    
    script:
    def filter = bed.getName() != 'nf-core_eager_dummy.txt' ? "-b $bed" : ''
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
    publishDir "${params.outdir}/nuclear_contamination", mode: params.publish_dir_mode

    when:
    params.run_nuclear_contamination

    input:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path(input), path(bai) from ch_for_nuclear_contamination

    output:
    tuple samplename, libraryid, lane, seqtype, organism, strandedness, udg, path('*.X.contamination.out') into ch_from_nuclear_contamination

    script:
    """
    samtools index ${input}
    angsd -i ${input} -r ${params.contamination_chrom_name}:5000000-154900000 -doCounts 1 -iCounts 1 -minMapQ 30 -minQ 30 -out ${libraryid}.doCounts
    contamination -a ${libraryid}.doCounts.icnts.gz -h ${projectDir}/assets/angsd_resources/HapMapChrX.gz 2> ${libraryid}.X.contamination.out
    """
 }
 
// As we collect all files for a single print_nuclear_contamination run, we DO NOT use the normal input/output tuple
process print_nuclear_contamination{
    label 'sc_tiny'
    publishDir "${params.outdir}/nuclear_contamination", mode: params.publish_dir_mode

    when:
    params.run_nuclear_contamination

    input:
    val 'Contam' from ch_from_nuclear_contamination.map { it[7] }.collect()

    output:
    file 'nuclear_contamination.txt'
    file 'nuclear_contamination_mqc.json' into ch_nuclear_contamination_for_multiqc

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
} else if ( params.metagenomic_tool == '' ) {
  ch_bam_filtering_for_metagenomic_malt = Channel.empty()
  ch_bam_filtering_for_metagenomic_kraken = Channel.empty()

}

// As we collect all files for a single MALT run, we DO NOT use the normal input/output tuple
process malt {
  label 'mc_small'
  publishDir "${params.outdir}/metagenomic_classification/malt", mode: params.publish_dir_mode

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'malt'

  input:
  file fastqs from ch_bam_filtering_for_metagenomic_malt.map { it[7] }.collect()
  file db from ch_db_for_malt

  output:
  path("*.rma6") into ch_rma_for_maltExtract
  path("*.sam.gz") optional true
  path("malt.log") into ch_malt_for_multiqc

  script:
  if ( "${params.malt_min_support_mode}" == "percent" ) {
    min_supp = "-supp ${params.malt_min_support_percent}" 
  } else if ( "${params.malt_min_support_mode}" == "reads" ) {
    min_supp = "-sup ${params.metagenomic_min_support_reads}"
  }
  def sam_out = params.malt_sam_output ? "-a . -f SAM" : ""
  """
  malt-run \
  -J-Xmx${task.memory.toGiga()}g \
  -t ${task.cpus} \
  -v \
  -o . \
  -d ${db} \
  ${sam_out} \
  -id ${params.percent_identity} \
  -m ${params.malt_mode} \
  -at ${params.malt_alignment_mode} \
  -top ${params.malt_top_percent} \
  ${min_supp} \
  -mq ${params.malt_max_queries} \
  --memoryMode ${params.malt_memory_mode} \
  -i ${fastqs.join(' ')} |&tee malt.log
  """
}

// MaltExtract performs aDNA evaluation from the output of MALT (damage patterns, read lengths etc.)

// As we collect all files for a single MALT extract run, we DO NOT use the normal input/output tuple
process maltextract {
  label 'mc_medium'
  publishDir "${params.outdir}/maltextract/", mode: params.publish_dir_mode

  when: 
  params.run_maltextract && params.metagenomic_tool == 'malt'

  input:
  file rma6 from ch_rma_for_maltExtract.collect()
  file taxon_list from ch_taxonlist_for_maltextract
  file ncbifiles from ch_ncbifiles_for_maltextract
  
  output:
  path "results/" type('dir')
  file "results/*_Wevid.json" optional true into ch_hops_for_multiqc 

  script:
  def destack = params.maltextract_destackingoff ? "--destackingOff" : ""
  def downsam = params.maltextract_downsamplingoff ? "--downSampOff" : ""
  def dupremo = params.maltextract_duplicateremovaloff ? "--dupRemOff" : ""
  def matches = params.maltextract_matches ? "--matches" : ""
  def megsum = params.maltextract_megansummary ? "--meganSummary" : ""
  def topaln = params.maltextract_topalignment ?  "--useTopAlignment" : ""
  def ss = params.single_stranded ? "--singleStranded" : ""
  """
  MaltExtract \
  -Xmx${task.memory.toGiga()}g \
  -t ${taxon_list} \
  -i ${rma6.join(' ')} \
  -o results/ \
  -r ${ncbifiles} \
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
    path(ckdb) from comp_kraken
    
    output:
    path(dbname) into ch_krakendb
    
    script:
    dbname = params.database.tokenize("/")[-1].tokenize(".")[0]
    """
    tar xvzf $ckdb
    """
  }

} else if (! params.database.endsWith(".tar.gz") && params.run_metagenomic_screening && params.metagenomic_tool == 'kraken') {
    ch_krakendb = path(params.database)
} else {
    ch_krakendb = Channel.empty()
}

process kraken {
  tag "$prefix"
  label 'mc_huge'
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode: params.publish_dir_mode

  when:
  params.run_metagenomic_screening && params.run_bam_filtering && params.bam_unmapped_type == 'fastq' && params.metagenomic_tool == 'kraken'

  input:
  path(fastq) from ch_bam_filtering_for_metagenomic_kraken.map { it[7] }
  path(krakendb) from ch_krakendb

  output:
  file "*.kraken.out" into ch_kraken_out
  tuple prefix, path("*.kreport") into ch_kraken_report, ch_kraken_for_multiqc

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
  tuple val(name), path(kraken_r) from ch_kraken_report

  output:
  tuple val(name), path('*.kraken_parsed.csv') into ch_kraken_parsed

  script:
  out = name+".kraken_parsed.csv"
  """
  kraken_parse.py -c ${params.metagenomic_min_support_reads} -o $out $kraken_r
  """    
}

process kraken_merge {
  publishDir "${params.outdir}/metagenomic_classification/kraken", mode: params.publish_dir_mode

  input:
  file csv_count from ch_kraken_parsed.map{ it[1] }.collect()

  output:
  path('kraken_count_table.csv')

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
    publishDir "${params.outdir}/documentation", mode: params.publish_dir_mode

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
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      if (filename.indexOf(".csv") > 0) filename
                      else null
                }

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml
    file "software_versions.csv"

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
    angsd -h |& head -n 1 | cut -d ' ' -f3-4 &> v_angsd.txt 2>&1 || true 
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
    eigenstrat_snp_coverage --version | cut -d ' ' -f2 >v_eigenstrat_snp_coverage.txt || true

    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

// MultiQC file generation for pipeline report

process multiqc {
    label 'sc_medium'

    publishDir "${params.outdir}/multiqc", mode: params.publish_dir_mode

    input:
    file multiqc_config from ch_multiqc_config
    file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
    file software_versions_mqc from software_versions_yaml.collect().ifEmpty([])
    file logo from ch_eager_logo
    file ('fastqc_raw/*') from ch_prefastqc_for_multiqc.collect().ifEmpty([])
    path('fastqc/*') from ch_fastqc_after_clipping.collect().ifEmpty([])
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
    file ('nuclear_contamination/*') from ch_nuclear_contamination_for_multiqc.collect().ifEmpty([])
    file ('genotyping/*') from ch_eigenstrat_snp_cov_for_multiqc.collect().ifEmpty([])

    file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")

    output:
    file "*multiqc_report.html" into ch_multiqc_report
    file "*_data"

    script:
    def rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    def rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    def custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f $rtitle $rfilename $multiqc_config $custom_config_file .
    """
}

// Send completion emails if requested, so user knows data is ready

workflow.onComplete {

    // Set up the e-mail variables
    def subject = "[nf-core/eager] Successful: $workflow.runName"
    if (!workflow.success) {
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
    if (workflow.repository) email_fields['summary']['Pipeline repository Git URL'] = workflow.repository
    if (workflow.commitId) email_fields['summary']['Pipeline repository Git Commit'] = workflow.commitId
    if (workflow.revision) email_fields['summary']['Pipeline Git branch/tag'] = workflow.revision
    email_fields['summary']['Nextflow Version'] = workflow.nextflow.version
    email_fields['summary']['Nextflow Build'] = workflow.nextflow.build
    email_fields['summary']['Nextflow Compile Timestamp'] = workflow.nextflow.timestamp

    // On success try attach the multiqc report
    def mqc_report = null
    try {
        if (workflow.success) {
            mqc_report = ch_multiqc_report.getVal()
            if (mqc_report instanceof ArrayList) {
                log.warn "[nf-core/eager] Found multiple reports from process 'multiqc', will use only one"
                mqc_report = mqc_report[0]
            }
        }
    } catch (all) {
        log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
    }

    // Check if we are only sending emails on failure
    email_address = params.email
    if (!params.email && params.email_on_fail && !workflow.success) {
        email_address = params.email_on_fail
    }

    // Render the TXT template
    def engine = new groovy.text.GStringTemplateEngine()
    def tf = new File("$projectDir/assets/email_template.txt")
    def txt_template = engine.createTemplate(tf).make(email_fields)
    def email_txt = txt_template.toString()

    // Render the HTML template
    def hf = new File("$projectDir/assets/email_template.html")
    def html_template = engine.createTemplate(hf).make(email_fields)
    def email_html = html_template.toString()

    // Render the sendmail template
    def smail_fields = [ email: email_address, subject: subject, email_txt: email_txt, email_html: email_html, projectDir: "$projectDir", mqcFile: mqc_report, mqcMaxSize: params.max_multiqc_email_size.toBytes() ]
    def sf = new File("$projectDir/assets/sendmail_template.txt")
    def sendmail_template = engine.createTemplate(sf).make(smail_fields)
    def sendmail_html = sendmail_template.toString()

    // Send the HTML e-mail
    if (email_address) {
        try {
            if (params.plaintext_email) { throw GroovyException('Send plaintext e-mail, not HTML') }
            // Try to send HTML e-mail using sendmail
            [ 'sendmail', '-t' ].execute() << sendmail_html
            log.info "[nf-core/eager] Sent summary e-mail to $email_address (sendmail)"
        } catch (all) {
            // Catch failures and try with plaintext
            def mail_cmd = [ 'mail', '-s', subject, '--content-type=text/html', email_address ]
            if (mqc_report == NULL) {
                log.warn "[nf-core/eager] Could not attach MultiQC report to summary email"
            } else if ( mqc_report.size() <= params.max_multiqc_email_size.toBytes() ) {
                mail_cmd += [ '-A', mqc_report ]
            }
            mail_cmd.execute() << email_html 
            log.info "[nf-core/eager] Sent summary e-mail to $email_address (mail)"
        }
    }

    // Write summary e-mail HTML to a file
    def output_d = new File("${params.outdir}/pipeline_info/")
    if (!output_d.exists()) {
        output_d.mkdirs()
    }
    def output_hf = new File(output_d, "pipeline_report.html")
    output_hf.withWriter { w -> w << email_html }
    def output_tf = new File(output_d, "pipeline_report.txt")
    output_tf.withWriter { w -> w << email_txt }

    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_red = params.monochrome_logs ? '' : "\033[0;31m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";

    if (workflow.stats.ignoredCount > 0 && workflow.success) {
        log.info "-${c_purple}Warning, pipeline completed, but with errored process(es) ${c_reset}-"
        log.info "-${c_red}Number of ignored errored process(es) : ${workflow.stats.ignoredCount} ${c_reset}-"
        log.info "-${c_green}Number of successfully ran process(es) : ${workflow.stats.succeedCount} ${c_reset}-"
    }

    if (workflow.success) {
        log.info "-${c_purple}[nf-core/eager]${c_green} Pipeline completed successfully${c_reset}-"
    } else {
        checkHostname()
        log.info "-${c_purple}[nf-core/eager]${c_red} Pipeline completed with errors${c_reset}-"
    }
    
}

/////////////////////////////////////
/* --    AUXILARY FUNCTIONS     -- */
/////////////////////////////////////

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue = params.monochrome_logs ? '' : "\033[0;34m";
    c_cyan = params.monochrome_logs ? '' : "\033[0;36m";
    c_dim = params.monochrome_logs ? '' : "\033[2m";
    c_green = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset = params.monochrome_logs ? '' : "\033[0m";
    c_white = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
    ${c_purple}  nf-core/eager v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}


def checkHostname() {
    def c_reset = params.monochrome_logs ? '' : "\033[0m"
    def c_white = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
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

// Channelling the TSV file containing FASTQ or BAM 
def extract_data(tsvFile) {
    Channel.fromPath(tsvFile)
        .splitCsv(header: true, sep: '\t')
        .dump()
        .map { row ->

            def expected_keys = ['Sample_Name', 'Library_ID', 'Lane', 'Colour_Chemistry', 'SeqType', 'Organism', 'Strandedness', 'UDG_Treatment', 'R1', 'R2', 'BAM']
            if ( !row.keySet().containsAll(expected_keys) ) exit 1, "[nf-core/eager] error: Invalid TSV input - malformed column names. Please check input TSV. Column names should be: ${expected_keys.join(", ")}"

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

            // check no empty metadata fields
            if (samplename == '' || libraryid == '' || lane == '' || colour == '' || seqtype == '' || seqtype == '' || udg == '' || r1 == '' || r2 == '') exit 1, "[nf-core/eager] error: a field does not contain any information. Ensure all cells are filled or contain 'NA' for optional fields. Check row:\n ${row}"

            // Check no 'empty' rows
            if (r1.matches('NA') && r2.matches('NA') && bam.matches('NA') && bai.matches('NA')) exit 1, "[nf-core/eager] error: A row in your TSV appears to have all files defined as NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Ensure BAMs aren't submitted with PE
            if (!bam.matches('NA') && seqtype.matches('PE')) exit 1, "[nf-core/eager] error: BAM input rows in TSV cannot be set as PE, only SE. See '--help' flag and documentation under 'running the pipeline' for more information. Check row for: ${samplename}"

            // Check valid UDG treatment
            if (!udg.matches('none') && !udg.matches('half') && !udg.matches('full')) exit 1, "[nf-core/eager] error: UDG treatment can only be 'none', 'half' or 'full'. See '--help' flag and documentation under 'running the pipeline' for more information. You have '${udg}'"

            // Check valid colour chemistry
            if (!colour == 2 && !colour == 4) exit 1, "[nf-core/eager] error: Colour chemistry in TSV can either be 2 (e.g. NextSeq/NovaSeq) or 4 (e.g. HiSeq/MiSeq)"

            //  Ensure that we do not accept incompatible chemistry setup
            if (!seqtype.matches('PE') && !seqtype.matches('SE')) exit 1, "[nf-core/eager] error:  SeqType for one or more rows in TSV is neither SE nor PE! see '--help' flag and documentation under 'running the pipeline' for more information. You have: '${seqtype}'"
            
           // So we don't accept existing files that are wrong format: e.g. fasta or sam
            if ( !r1.matches('NA') && !has_extension(r1, "fastq.gz") && !has_extension(r1, "fq.gz") && !has_extension(r1, "fastq") && !has_extension(r1, "fq")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r1}"
            if ( !r2.matches('NA') && !has_extension(r2, "fastq.gz") && !has_extension(r2, "fq.gz") && !has_extension(r2, "fastq") && !has_extension(r2, "fq")) exit 1, "[nf-core/eager] error: A specified R2 file either has a non-recognizable FASTQ extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${r2}"
            if ( !bam.matches('NA') && !has_extension(bam, "bam")) exit 1, "[nf-core/eager] error: A specified R1 file either has a non-recognizable BAM extension or is not NA. See '--help' flag and documentation under 'running the pipeline' for more information. Check: ${bam}"
            
            [ samplename, libraryid, lane, colour, seqtype, organism, strandedness, udg, r1, r2, bam ]

        }

    }

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "[nf-core/eager] error:  Invalid TSV input - malformed row (e.g. missing column) in ${row}, see '--help' flag and documentation under 'running the pipeline' for more information"
    return true
}

// Return file if it exists
def return_file(it) {
    if (!file(it).exists()) exit 1, "[nf-core/eager] error: Cannot find supplied FASTQ or BAM input file. If using input method TSV set to NA if no file required. See '--help' flag and documentation under 'running the pipeline' for more information. Check file: ${it}" 
    return file(it)
}

// Check file extension
def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
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
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). Did you accidentally supply paired-end data?  see '--help' flag and documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" } 
                  }

        } else if (!pe_se ){
            log.info "Generating paired-end FASTQ data TSV"

            Channel
                .fromFilePairs( input )
                .filter { it =~/.*.fastq.gz|.*.fq.gz|.*.fastq|.*.fq/ }
                .ifEmpty { exit 1, "[nf-core/eager] error: Files could not be found. Do the specified FASTQ read files end in: '.fastq.gz', '.fq.gz', '.fastq', or '.fq'? Did you forget --single_end?" }
                .map { row -> [ row[0], [ row[1][0], row[1][1] ] ] }
                .ifEmpty { exit 1, "[nf-core/eager] error: --input was empty - no input files supplied!" }
                .into { ch_reads_for_faketsv; ch_reads_for_validate }

                // Check we don't have any duplicated sample names due to fromFilePairs behaviour of calculating sample name from anything before R1/R2 glob
                ch_reads_for_validate
                  .groupTuple()
                  .map{
                    if ( validate_size(it[1], 1) ) { null } else { exit 1, "[nf-core/eager] error: You have supplied non-unique sample names (text before R1/R2 indication). See '--help' flag and documentation under 'running the pipeline' for more information. Check duplicates of: ${it[0]}" } 
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
