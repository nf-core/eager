// Prepare various reference FASTA index files for downstream steps

params.gunzip_options                           = [:]
params.bwa_index_options                        = [:]
params.bowtie2_build_options                    = [:]
params.samtools_faidx_options                   = [:]
params.gatk4_createsequencedictionary_options   = [:]

include { GUNZIP as GUNZIP_FASTA }          from '../../modules/nf-core/modules/gunzip/main'                            addParams( options: params.gunzip_options )
include { BWA_INDEX }                       from '../../modules/nf-core/modules/bwa/index/main'                         addParams( options: params.bwa_index_options )
include { BOWTIE2_BUILD }                   from '../../modules/nf-core/modules/bowtie2/build/main'                     addParams( options: params.bowtie2_build_options )
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/modules/samtools/faidx/main'                    addParams( options: params.samtools_faidx_options )
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'    addParams( options: params.gatk4_createsequencedictionary_options )

workflow INDEX_FASTA {

    // Requires dummy to uncompress, as FASTA needs to be given as string from params (From ViralRecon)
    take:
    dummy_file

    main:
    // Uncompress genome fasta file if required
    if ( params.fasta.endsWith('.gz') ) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
    } else {
        ch_fasta = file( params.fasta )
    }

    // Run BWA indexing
    ch_bwa_index = Channel.empty()
    if ( !params.bwa_index_dir ) {
        ch_bwa_index = BWA_INDEX( ch_fasta ).index
    } else {
        ch_bwa_index = file( params.bwa_index_dir )
    }

    // Run Bowtie2 Indexing
    ch_bt2_build = Channel.empty()
    if ( !params.bt2_index_dir ) {
        ch_bt2_build = BOWTIE2_BUILD( ch_fasta ).index
    } else {
        ch_bt2_build = file( params.bt2_index_dir )
    }

    // Run samtools faidx
    ch_samtools_faidx = Channel.empty()
    if ( !params.fasta_fai ) {
        ch_samtools_faidx = SAMTOOLS_FAIDX( ch_fasta ).fai
    } else {
        ch_samtools_faidx = file( params.fasta_fai )
    }

    // Run GATK4/Picard CreateqSequenceDictionary
    ch_gatk_seqdict = Channel.empty()
    if ( !params.fasta_dict ) {
        ch_gatk_seqdict = GATK4_CREATESEQUENCEDICTIONARY( ch_fasta ).dict
    } else {
        ch_gatk_seqdict = file( params.fasta_dict )
    }

    emit:
    fasta           = ch_fasta
    bwa_index       = ch_bwa_index
    bt2_index       = ch_bt2_build
    samtools_faidx  = ch_samtools_faidx
    gatk_seqdict    = ch_gatk_seqdict

}



