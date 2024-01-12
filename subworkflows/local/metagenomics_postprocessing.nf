include { MALTEXTRACT     } from '../../modules/nf-core/maltextract/main'
include { AMPS            } from '../../modules/nf-core/amps/main'
include { TAXPASTA_MERGE  } from '../../modules/nf-core/taxpasta/merge/main'
include { MEGAN_RMA2INFO  } from '../../modules/nf-core/megan/rma2info/main'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    ch_postprocessing_input // different between each profiling --> postprocessing tool,
                            // defined in metagenomics profiling subworkflow

    main:
    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // For MALT we have an additional step that includes maltextract+amps
    if ( params.metagenomics_run_postprocessing && params.metagenomics_profiling_tool == 'malt' ) {

        //maltextract doesnt accepts a meta param in the first input channel, so remove it
        ch_maltextract_input = ch_postprocessing_input.map{it[1]}

        tax_list = Channel.fromPath(params.metagenomics_maltextract_taxonlist)
        ncbi_dir = Channel.fromPath(params.metagenomics_maltextract_ncbidir)

        MALTEXTRACT ( ch_maltextract_input, tax_list, ncbi_dir)

        AMPS ( MALTEXTRACT.out.results, tax_list, params.metagenomics_maltextract_filter )


        //Also, prepare Malt for taxpasta by running rma2info

        MEGAN_RMA2INFO( ch_postprocessing_input, true )
        ch_postprocessing_input = MEGAN_RMA2INFO.out.txt

        ch_versions      = ch_versions.mix( MALTEXTRACT.out.versions.first(), AMPS.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( AMPS.out.json )
    }

   // Run taxpasta for everything!

    ch_postprocessing_input = ch_postprocessing_input
    .map{
        meta, report ->
        [report]
    }
    .collect()
    .map{
        reports ->
        [
            [
                "id":"${params.metagenomics_profiling_tool}_profiles_all_samples_merged_taxpasta",
                "profiler":params.metagenomics_profiling_tool == 'malt' ? 'megan6' : params.metagenomics_profiling_tool
            ],
            reports
        ]
    }

    TAXPASTA_MERGE( ch_postprocessing_input, [], [] )
    ch_versions = ch_versions.mix(TAXPASTA_MERGE.out.versions)
    ch_multiqc_files = ch_multiqc_files.mix(TAXPASTA_MERGE.out.merged_profiles)

    emit:
    versions = ch_versions
    mqc      = ch_multiqc_files

}

