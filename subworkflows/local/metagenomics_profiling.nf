//
// Complexity filtering and metagenomics screening of sequencing reads
//

// Much taken from nf-core/taxprofile subworkflows/local/profiling.nf

include { MALT_RUN                       } from '../../modules/nf-core/malt/run/main'
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { METAPHLAN_METAPHLAN            } from '../../modules/nf-core/metaphlan/metaphlan/main'

workflow METAGENOMICS_PROFILING {

    take:
    reads
    database

    main:

    ch_versions                  = Channel.empty()
    ch_raw_classifications       = Channel.empty()
    ch_raw_profiles              = Channel.empty()
    ch_multiqc_files             = Channel.empty()
    ch_postprocessing_input      = Channel.empty()

    /*
        PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool has a slightly different input structure and generally separate
    // input channels for reads vs database. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/database
    // channel element order in sync with each other

    if ( params.metagenomics_profiling_tool == 'malt' ) {

        // Reset entire input meta for MALT to just database name,
        // as we don't run run on a per-sample basis due to huge databases
        // so all samples are in one run and so sample-specific metadata
        // unnecessary. Set as database name to prevent `null` job ID and prefix.

        if ( params.metagenomics_malt_group_size > 0 ) {
            ch_input_for_malt_tmp =  reads
                .map { meta, reads -> reads }
                .collate( params.metagenomics_malt_group_size ) //collate into bins of defined lengths
                .map{
                    reads ->
                    // add new meta with db-name as id
                    [[id: file(params.metagenomics_profiling_database).getBaseName() ], reads]
                }
                .combine(database) //combine with database
                .multiMap{
                    // and split apart again
                    meta, reads, database ->
                        reads: [meta, reads]
                        database: database
                }
            ch_input_for_malt = ch_input_for_malt_tmp.reads
            database = ch_input_for_malt_tmp.database
        }

        else {
            ch_input_for_malt =  reads
                .map { meta, reads -> reads }
                .collect()
                .map{
                    // make sure id is db_name for publishing purposes.
                    reads ->
                    [[id: file(params.metagenomics_profiling_database).getBaseName() ], reads]
                }
        }

        // Run MALT
        MALT_RUN ( ch_input_for_malt, database )

        ch_maltrun_for_megan = MALT_RUN.out.rma6
            .transpose()
            .map {
                meta, rma ->
                    // re-extract meta from file names, use filename without rma to
                    // ensure we keep paired-end information in downstream filenames
                    // when no pair-merging
                    def meta_new = meta.clone()
                    meta_new['db_name'] = meta.id
                    meta_new['id'] = rma.baseName
                [ meta_new, rma ]
            }

        ch_versions                 = ch_versions.mix( MALT_RUN.out.versions.first() )
        ch_raw_classifications      = ch_raw_classifications.mix( ch_maltrun_for_megan )
        ch_multiqc_files            = ch_multiqc_files.mix( MALT_RUN.out.log )
        ch_postprocessing_input     = ch_postprocessing_input.mix( ch_maltrun_for_megan )
    }

    else if ( params.metagenomics_profiling_tool == 'metaphlan' ) {

        reads = reads.combine(database)
        metaphlan_reads = reads.map{ meta, reads, database -> [meta, reads] }
        metaphlan_db = reads.map{ meta, reads, database -> [database] }

        METAPHLAN_METAPHLAN ( metaphlan_reads , metaphlan_db )
        ch_versions             = ch_versions.mix( METAPHLAN_METAPHLAN.out.versions.first() )
        ch_raw_profiles         = ch_raw_profiles.mix( METAPHLAN_METAPHLAN.out.profile )
        ch_postprocessing_input = ch_postprocessing_input.mix( METAPHLAN_METAPHLAN.out.profile )

    }

    if ( params.metagenomics_profiling_tool == 'krakenuniq' ) {
        // run kraken uniq per sample, to preserve the meta-data

        reads = reads.combine(database)
        krakenuniq_reads = reads.map{ meta, reads, database -> [meta, reads] }
        krakenuniq_db = reads.map{ meta, reads, database -> [database] }

        KRAKENUNIQ_PRELOADEDKRAKENUNIQ (
            krakenuniq_reads,
            krakenuniq_db,
            params.metagenomics_krakenuniq_ram_chunk_size,
            params.metagenomics_kraken_save_reads,
            true,
            params.metagenomics_kraken_save_read_classifications
        )

        ch_multiqc_files           = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_versions                = ch_versions.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions.first() )
        ch_raw_classifications     = ch_raw_classifications.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.classified_assignment )
        ch_raw_profiles            = ch_raw_profiles.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_multiqc_files           = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_postprocessing_input    = ch_postprocessing_input.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
    }

    else if ( params.metagenomics_profiling_tool == 'kraken2' ) {
        // run kraken2 per sample

        reads = reads.combine(database)
        kraken2_reads = reads.map{meta, reads, database -> [meta, reads]}
        kraken2_db = reads.map{meta, reads, database -> [database]}

        KRAKEN2_KRAKEN2 (
            kraken2_reads,
            kraken2_db,
            params.metagenomics_kraken_save_reads,
            params.metagenomics_kraken_save_read_classifications
        )

        ch_multiqc_files            = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions                 = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications      = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles             = ch_raw_profiles.mix( KRAKEN2_KRAKEN2.out.report )
        ch_postprocessing_input     = ch_postprocessing_input.mix( KRAKEN2_KRAKEN2.out.report )
    }

    emit:
    versions             = ch_versions          // channel: [ versions.yml ]
    classifications      = ch_raw_classifications
    profiles             = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    postprocessing_input = ch_postprocessing_input // channel: [ val(meta), [ inputs_for_postprocessing_tools ] ] // see info at metagenomics_postprocessing
    mqc                  = ch_multiqc_files

}
