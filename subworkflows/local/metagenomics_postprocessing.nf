// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { MALTEXTRACT } from '../../../modules/nf-core/maltextract/main'
include { KRAKENPARSE } from '../../../modules/local/krakenparse'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    // TODO nf-core: edit input (take) channels
    ch_postprocessing_input // different between kraken and malt

    main:

    ch_versions      = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( params.metagenomics_postprocessing_tool == 'maltextract') {
        MALTEXTRACT ( ch_postprocessing_input, params.taxon_list, params.ncbi_dir )
        ch_versions = ch_versions.mix( MALTEXTRACT.out.versions.first() )
        ch_results  = ch_results.mix( MALTEXTRACT.out.results )
    }
    else if ( params.metagenomics_postprocessing_tool == 'krakenparse' ) {
        KRAKENPARSE ( ch_postprocessing_input )
        ch_versions = ch_versions.mix( KRAKENPARSE.out.versions.first() )
        ch_results  = ch_results.mix( KRAKENPARSE.out.results )
    }
// TODO check how to actually emit krakenparse output channels into one directory
// TODO check if necessary to have merge_kraken_parsed
// TODO add paths for multiqc files (maltextract and maybe kraken parse)
    emit:
    // TODO nf-core: edit emitted channels
    versions          = ch_versions
    results_directory = ch_results
    mqc               = ch_multiqc_files

}

