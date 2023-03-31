//
// Complexity filtering and metagenomics screening of sequencing reads
//

// Much taken from nf-core/taxprofile subworkflows/local/profiling.nf

include { MALT_RUN                                      } from '../../modules/nf-core/malt/run/main'
include { KRAKEN2_KRAKEN2                               } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ                } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { METAPHLAN3_METAPHLAN3                         } from '../../modules/nf-core/metaphlan3/metaphlan3/main'

workflow METAGENOMICS_PROFILING {

    take:
    reads // channel: [ [ meta ] , [ reads ] ]
    database // channel: [ [ meta ] , path ]

    main:

    ch_versions                  = Channel.empty()
    ch_raw_classifications       = Channel.empty()
    ch_raw_profiles              = Channel.empty()
    ch_multiqc_files             = Channel.empty()
    // TODO: malt, metaphylan, kraken2, krakenuniq
    // TODO: maltextract, krakenparse

    /*
        PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool as a slightly different input structure and generally separate
    // input channels for reads vs database. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/database
    // channel element order in sync with each other

    if ( params.metagenomics_profiling_tool == 'malt' ) {

        if ( params.metagenomics_profiling_malt_group_size > 0 ) {
            ch_input_for_malt =  reads
                .map {
                    meta, reads ->

                        // Reset entire input meta for MALT to just database name,
                        // as we don't run run on a per-sample basis due to huge datbaases
                        // so all samples are in one run and so sample-specific metadata
                        // unnecessary. Set as database name to prevent `null` job ID and prefix.

                        def temp_meta = [ id: database ]

                        // Combine reduced sample metadata with updated database parameters metadata,
                        // make sure id is db_name for publishing purposes.

                        [ temp_meta, reads, database ]

                }
                .groupTuple(by: [0,2], size: params.metagenomics_profiling_malt_group_size, remainder: true)
                .multiMap {
                    meta, reads, database ->
                        reads: [ meta, reads ]
                        database: database
                }
        }

        else {
            ch_input_for_malt =  reads
                .map {
                    meta, reads ->

                        // Reset entire input meta for MALT to just database name,
                        // as we don't run run on a per-sample basis due to huge datbaases
                        // so all samples are in one run and so sample-specific metadata
                        // unnecessary. Set as database name to prevent `null` job ID and prefix.

                        def temp_meta = [ id: database ]

                        // Combine reduced sample metadata with updated database parameters metadata,
                        // make sure id is db_name for publishing purposes.

                        [ temp_meta, reads, database ]

                }
                .groupTuple(by: [0,2])
                .multiMap {
                    meta, reads, database ->
                        reads: [ meta, reads ]
                        database: database
                }
        }

        // MALT: We groupTuple to have all samples in one channel for MALT as database
        // loading takes a long time, so we only want to run it once per database, unless otherwise specified (eg grouping samples)

        MALT_RUN ( ch_input_for_malt.reads, ch_input_for_malt.database )

        ch_maltrun_for_megan = MALT_RUN.out.rma6
                                .transpose()
                                .map{
                                    meta, rma ->
                                            // re-extract meta from file names, use filename without rma to
                                            // ensure we keep paired-end information in downstream filenames
                                            // when no pair-merging
                                            def meta_new = meta.clone()
                                            meta_new['db_name'] = meta.id
                                            meta_new['id'] = rma.baseName
                                        [ meta_new, rma ]
                                }

        ch_versions            = ch_versions.mix( MALT_RUN.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( ch_maltrun_for_megan )
        ch_multiqc_files       = ch_multiqc_files.mix( MALT_RUN.out.log )
    }

    if ( params.metagenomics_profiling_tool == 'metaphlan3' ) {

        METAPHLAN3_METAPHLAN3 ( reads , database )
        ch_versions        = ch_versions.mix( METAPHLAN3_METAPHLAN3.out.versions.first() )
        ch_raw_profiles    = ch_raw_profiles.mix( METAPHLAN3_METAPHLAN3.out.profile )

    }

    if ( params.metagenomics_profiling_tool == 'krakenuniq' ) {
        ch_input_for_krakenuniq =  ch_input_for_profiling.krakenuniq
                                    .map {
                                        meta, reads, db_meta, db ->
                                            [[id: db_meta.db_name, single_end: meta.single_end], reads, db_meta, db]
                                    }
                                    .groupTuple(by: [0,2,3])
                                    .multiMap {
                                        single_meta, reads, db_meta, db ->
                                            reads: [ single_meta + db_meta, reads.flatten() ]
                                            db: db
                                }
        // Hardcode to _always_ produce the report file (which is our basic output, and goes into)
        KRAKENUNIQ_PRELOADEDKRAKENUNIQ ( reads , database , params.metagenomics_profiling_krakenuniq_ram_chunk_size, params.metagenomics_profiling_krakenuniq_save_reads, true, params.metagenomics_profiling_krakenuniq_save_read_classifications )
        ch_versions            = ch_versions.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.classified_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report )

    }

    if ( params.metagenomics_profiling_tool == 'kraken2' ) {

        KRAKEN2_KRAKEN2 ( reads, database, params.metagenomics_profiling_kraken2_save_reads, params.metagenomics_profiling_kraken2_save_readclassification )
        ch_multiqc_files       = ch_multiqc_files.mix( KRAKEN2_KRAKEN2.out.report )
        ch_versions            = ch_versions.mix( KRAKEN2_KRAKEN2.out.versions.first() )
        ch_raw_classifications = ch_raw_classifications.mix( KRAKEN2_KRAKEN2.out.classified_reads_assignment )
        ch_raw_profiles        = ch_raw_profiles.mix(
            KRAKEN2_KRAKEN2.out.report
                // Set the tool to be strictly 'kraken2' instead of potentially 'bracken' for downstream use.
                // Will remain distinct from 'pure' Kraken2 results due to distinct database names in file names.
                .map { meta, report -> [meta + [tool: 'kraken2'], report]}
        )

    }

    emit:
    versions          = ch_versions          // channel: [ versions.yml ]
    classifications   = ch_raw_classifications
    profiles          = ch_raw_profiles    // channel: [ val(meta), [ reads ] ] - should be text files or biom
    mqc               = ch_multiqc_files
}
