// Prepare various reference FASTA index files for downstream steps
include { GUNZIP as GUNZIP_FASTA }          from '../../modules/nf-core/modules/gunzip/main'
include { BWA_INDEX }                       from '../../modules/nf-core/modules/bwa/index/main'
include { BOWTIE2_BUILD }                   from '../../modules/nf-core/modules/bowtie2/build/main'
include { SAMTOOLS_FAIDX }                  from '../../modules/nf-core/modules/samtools/faidx/main'
include { GATK4_CREATESEQUENCEDICTIONARY }  from '../../modules/nf-core/modules/gatk4/createsequencedictionary/main'

workflow INDEX_FASTA {

    // Requires dummy to uncompress, as FASTA needs to be given as string from params (From ViralRecon)
    take:
    dummy_file

    main:
    // Set versions channel
    ch_versions = Channel.empty()

    // Uncompress genome fasta file if required
    if ( params.fasta.endsWith('.gz') ) {
        ch_fasta = GUNZIP_FASTA ( params.fasta ).gunzip
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file( params.fasta, checkIfExists: true )
    }

    // Run BWA indexing
    // TODO add check that expected files are present and filename matches reference?
    ch_bwa_index = Channel.empty()
    if (  ( params.mapping_mapper == "bwaaln" | params.mapping_mapper == "bwamem" | params.mapping_mapper == 'circularmapper' ) && params.index_bwa_index_dir ) {
        ch_bwa_index = file( params.index_bwa_index_dir, checkIfExists: true )
    } else if ( ( params.mapping_mapper == "bwaaln" | params.mapping_mapper == "bwamem" | params.mapping_mapper == 'circularmapper' ) && !params.index_bwa_index_dir )  {
        ch_bwa_index = BWA_INDEX( ch_fasta ).index
        ch_versions = ch_versions.mix(BWA_INDEX.out.versions)
    }

    // Run Bowtie2 Indexing
    ch_bt2_build = Channel.empty()
    if ( params.mapping_mapper == "bowtie2" && params.mapping_mapper && params.index_bt2_index_dir ) {
        ch_bt2_build = file( params.index_bt2_index_dir, checkIfExists: true )
    } else if ( params.mapping_mapper == "bowtie2" && params.mapping_mapper && !params.index_bt2_index_dir ) {
        ch_bt2_build = BOWTIE2_BUILD( ch_fasta ).index
        ch_versions = ch_versions.mix(BOWTIE2_BUILD.out.versions)
    }

    // Run samtools faidx
    ch_samtools_faidx = Channel.empty()
    if ( !params.index_fasta_fai ) {
        ch_samtools_faidx = SAMTOOLS_FAIDX( [[], ch_fasta] ).fai
        ch_versions = ch_versions.mix( SAMTOOLS_FAIDX.out.versions )
    } else {
        ch_samtools_faidx = file( params.index_fasta_fai, checkIfExists: true )
    }

    // Run GATK4/Picard CreateqSequenceDictionary
    ch_gatk_seqdict = Channel.empty()
    if ( !params.index_fasta_dict ) {
        ch_gatk_seqdict = GATK4_CREATESEQUENCEDICTIONARY( ch_fasta ).dict
        ch_versions = ch_versions.mix(GATK4_CREATESEQUENCEDICTIONARY.out.versions)
    } else {
        ch_gatk_seqdict = file( params.index_fasta_dict, checkIfExists: true )
    }

    emit:
    fasta           = ch_fasta
    bwa_index       = ch_bwa_index
    bt2_index       = ch_bt2_build
    //samtools_faidx  = ch_samtools_faidx
    gatk_seqdict    = ch_gatk_seqdict
    versions        = ch_versions
}
