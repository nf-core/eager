include { MALTEXTRACT          } from '../../modules/nf-core/maltextract/main'
include { AMPS                 } from '../../modules/nf-core/amps/main'
include { TAXPASTA_MERGE       } from '../../modules/nf-core/taxpasta/merge/main'
include { TAXPASTA_STANDARDISE } from '../../modules/nf-core/taxpasta/standardise/main'
include { MEGAN_RMA2INFO       } from '../../modules/nf-core/megan/rma2info/main'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    ch_postprocessing_input // different between each profiling --> postprocessing tool,
                            // defined in metagenomics profiling subworkflow

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // For MALT we have an additional step that includes maltextract+amps
    if ( params.metagenomics_run_postprocessing && params.metagenomics_profiling_tool == 'malt' ) {

        tax_list = Channel.fromPath(params.metagenomics_maltextract_taxonlist)
        ncbi_dir = Channel.fromPath(params.metagenomics_maltextract_ncbidir)

        // Malt could have been executed multiple times (group_size paramter and strandedness)
        // We want to combine the chunks, but run MaltExtract on double and singlestranded individually
        ch_strandedness = ch_postprocessing_input
            .transpose()
            .map{ meta, reads ->
                [
                    meta + [ 'id': "${meta.strandedness}stranded" ],
                    reads
                ]
            }
            .groupTuple(by:0)

        // could no be two entries in the channel, so combine with the tax_list and ncbi
        ch_maltextract_input = ch_strandedness
            .combine(tax_list)
            .combine(ncbi_dir)
            .multiMap{
                rma6:[it[0],it[1]]
                tax_list:it[2]
                ncbi_dir:it[3]
            }

        //RUN MaltExtract
        MALTEXTRACT ( ch_maltextract_input.rma6, ch_maltextract_input.tax_list, ch_maltextract_input.ncbi_dir )

        // now we need to run AMPS for each MALTEXTRACT output
        ch_amps_input = MALTEXTRACT.out.results.map{ it[1] }

        AMPS ( ch_amps_input, ch_maltextract_input.tax_list, params.metagenomics_maltextract_filter )

        //Now, prepare Malt rma6 output for taxpasta by running rma2info
        ch_rma2info_input = ch_postprocessing_input
            .transpose()
            .map {
                meta, rma ->
                    // re-extract meta from file names, use filename without rma to
                    // ensure we keep paired-end information in downstream filenames
                    // when no pair-merging
                    [
                        meta + ['db_name': meta.id, 'id': rma.baseName ],
                        rma
                    ]
            }

        MEGAN_RMA2INFO( ch_rma2info_input, true )
        ch_postprocessing_input = MEGAN_RMA2INFO.out.txt

        ch_versions      = ch_versions.mix( MALTEXTRACT.out.versions.first(), AMPS.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( AMPS.out.json )
    }

    // Run taxpasta for everything!
    // We need to know how many reports we have, so that we can run either taxpasta standardise or taxpasta merge
    ch_report_count = ch_postprocessing_input.transpose().count()

    ch_postprocessing_input = ch_postprocessing_input
        .transpose()
        .map{
            meta, report ->
            [
                [
                    "id":"${params.metagenomics_profiling_tool}_profiles_all_samples_merged_taxpasta",
                    "profiler":params.metagenomics_profiling_tool == 'malt' ? 'megan6' : params.metagenomics_profiling_tool
                ],
                report
            ]
        }
        .groupTuple(by:0)
        .combine(ch_report_count)
        .branch{
            standardise: it[2] == 1
            merge: true
        }

    ch_standardise_input = ch_postprocessing_input.standardise.map{ meta, reports, count ->
        [meta, reports]
    }

    TAXPASTA_STANDARDISE( ch_standardise_input, [] )
    ch_versions = ch_versions.mix(TAXPASTA_STANDARDISE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TAXPASTA_STANDARDISE.out.standardised_profile)

    ch_merge_input = ch_postprocessing_input.merge.map{ meta, reports, count ->
        [meta, reports]
    }

    TAXPASTA_MERGE( ch_merge_input, [], [] )
    ch_versions = ch_versions.mix(TAXPASTA_MERGE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TAXPASTA_MERGE.out.merged_profiles)


    emit:
    versions = ch_versions
    mqc      = ch_multiqc_files

}
