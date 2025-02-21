//
// Complexity filtering and metagenomics screening of sequencing reads
//

// Much taken from nf-core/taxprofile subworkflows/local/profiling.nf

include { MALT_RUN                       } from '../../modules/nf-core/malt/run/main'
include { KRAKEN2_KRAKEN2                } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENUNIQ_PRELOADEDKRAKENUNIQ } from '../../modules/nf-core/krakenuniq/preloadedkrakenuniq/main'
include { METAPHLAN_METAPHLAN            } from '../../modules/nf-core/metaphlan/metaphlan/main'
include { CAT_CAT as CAT_CAT_MALT        } from '../../modules/nf-core/cat/cat/main'
include { UNTAR                          } from '../../modules/nf-core/untar/main'

workflow METAGENOMICS_PROFILING {

    take:
    ch_reads
    ch_database

    main:
    ch_versions             = Channel.empty()
    ch_multiqc_files        = Channel.empty()
    ch_postprocessing_input = Channel.empty()

    /*
        UNTAR THE DATABASE IF NECESSARY
    */

    ch_database = ch_database
    .branch{
        untar: it ==~ /.*.tar.gz/
        base:true
    }

    // untar the database
    ch_untar_input = ch_database.untar.map{ [[], it] }

    UNTAR( ch_untar_input )
    ch_untar_output = UNTAR.out.untar.map{ it[1] }

    // back to the original database channel...
    ch_database = ch_database.base.mix(ch_untar_output)

    /*
        PREPARE PROFILER INPUT CHANNELS & RUN PROFILING
    */

    // Each tool has a slightly different input structure and generally separate
    // input channels for reads vs database. We restructure the channel tuple
    // for each tool and make liberal use of multiMap to keep reads/database
    // channel element order in sync with each other

    ch_reads.view()

    if ( params.metagenomics_profiling_tool == 'malt' ) {

        // Optional parallel run of malt available:
        // If parallel execution, split into groups with meta id of the first library id of group
        // Merging of maltlog will be done by concatenation

        // If no parallel execution (default):
        // Reset entire input meta for MALT to just database name,
        // as we don't run run on a per-sample basis due to huge databases
        // so all samples are in one run and so sample-specific metadata
        // unnecessary. Set as database name to prevent `null` job ID and prefix.

        def label = file(params.metagenomics_profiling_database).getBaseName()

        // For the next step we need the number of analysis-groups for the spezified number of input files
        // since we work with channels, we need a channel that stores that information
        ch_tmp_groups =  params.metagenomics_malt_group_size > 0 ? ch_reads.collate(params.metagenomics_malt_group_size).count() : Channel.of(1)
        // this is for enumerating the channel-entries in the ch_reads channel
        def groups_counter = 0

        //replace the meta in a way that groupTuple splits the entries
        //by strandedness and metagenomics_malt_group_size
        //NOTE: known limitations
        // this method splits the entries into groups of size malt_group_size, but if there is a mix of ds and ss entries
        // the groups are split again by ds and ss with groupTuple. So they might end up smaller than malt_group_size
        // could be prevented by branching early and running the lower part twice for ss and ds individually
        // but this is an edge-case and might never be relevant...
        ch_input_for_malt = ch_reads.combine(ch_tmp_groups).map{ meta, reads, n_groups ->
            [
                [
                    label: label,
                    strandedness:meta.strandedness,
                    id:"${meta.strandedness}stranded_${groups_counter++%n_groups}"
                ],
                reads
            ]
        }
        .groupTuple(by:0)

        // We might have multiple chunks in the reads_channel
        // each of which requires the database
        ch_input_for_malt = ch_input_for_malt.combine(ch_database)

        // Split Channels into reads and database
        ch_input_for_malt = ch_input_for_malt
            .multiMap{ meta, reads, database ->
                reads: [meta, reads]
                database: database
            }

        // Run MALT
        MALT_RUN ( ch_input_for_malt.reads, ch_input_for_malt.database )

        // Recombine log files for outputting
        ch_log_for_cat =
            MALT_RUN.out.log
                .map {
                    meta,log -> log
                }
                .collect()
                .map {
                    log -> [['id': label], log]
                }

        CAT_CAT_MALT ( ch_log_for_cat )

        ch_versions             = MALT_RUN.out.versions.first()
        ch_multiqc_files        = MALT_RUN.out.log
        ch_postprocessing_input = MALT_RUN.out.rma6
    }

    else if ( params.metagenomics_profiling_tool == 'metaphlan' ) {

        ch_metaphlan_input = ch_reads
            .combine(ch_database)
            .multiMap{ meta, reads, db ->
                reads: [meta, reads]
                database: db
            }

        METAPHLAN_METAPHLAN ( ch_metaphlan_input.reads , ch_metaphlan_input.database )
        ch_versions             = METAPHLAN_METAPHLAN.out.versions.first()
        ch_postprocessing_input = METAPHLAN_METAPHLAN.out.profile

    }

    else if ( params.metagenomics_profiling_tool == 'krakenuniq' ) {
        // run krakenuniq once for all samples, unless non-merged PE vs SE data

        ch_krakenuniq_input = ch_reads
            .map{ meta, file ->
                [
                    ['single_end':meta['single_end']], // retain single_end vs paired_end bools for input splitting
                    file
                ]
            }
            .groupTuple(by:0)
                .map { meta, files ->
                [
                    meta, files.flatten()
                ]}

        ch_krakenuniq_input.view()

        ch_krakenuniq_input = ch_krakenuniq_input.combine(ch_database)
            .multiMap{
                meta, files, database ->
                    meta_files_input: [meta, files]
                    database: database
            }

        KRAKENUNIQ_PRELOADEDKRAKENUNIQ (
            ch_krakenuniq_input.meta_files_input,
            'fastq', // only fastq files can get to the input channel
            ch_krakenuniq_input.database,
            params.metagenomics_krakenuniq_ramchunksize,
            params.metagenomics_kraken2_savereads,
            true, // save read assignments
            params.metagenomics_kraken2_savereadclassifications
        )

        ch_versions             = KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.versions
        ch_multiqc_files        = KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report
        ch_postprocessing_input = KRAKENUNIQ_PRELOADEDKRAKENUNIQ.out.report

    }

    else if ( params.metagenomics_profiling_tool == 'kraken2' ) {
        // run kraken2 per sample
        // it lacks the option of krakenuniq

        ch_kraken2_input = ch_reads
            .combine(ch_database)
            .multiMap{ meta, reads, db ->
                reads: [meta, reads]
                database: db
            }

        KRAKEN2_KRAKEN2 (
            ch_kraken2_input.reads,
            ch_kraken2_input.database,
            params.metagenomics_kraken2_savereads,
            params.metagenomics_kraken2_savereadclassifications
        )

        ch_multiqc_files        = KRAKEN2_KRAKEN2.out.report
        ch_versions             = KRAKEN2_KRAKEN2.out.versions.first()
        ch_postprocessing_input = KRAKEN2_KRAKEN2.out.report
    }

    emit:
    versions             = ch_versions             // channel: [ versions.yml ]
    postprocessing_input = ch_postprocessing_input // channel: [ val(meta), [ inputs_for_postprocessing_tools ] ] // see info at metagenomics_postprocessing
    mqc                  = ch_multiqc_files

}
