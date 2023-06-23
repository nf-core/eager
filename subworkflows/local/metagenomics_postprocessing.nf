include { MALTEXTRACT                    } from '../../modules/nf-core/maltextract/main'
include { AMPS                           } from '../../modules/nf-core/amps/main'
include { KRAKENPARSE                    } from '../../modules/local/krakenparse'
include { KRAKENMERGE                    } from '../../modules/local/krakenmerge'
include { METAPHLAN_MERGEMETAPHLANTABLES } from '../modules/nf-core/metaphlan/mergemetaphlantables/main'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    ch_postprocessing_input // different between each profiling --> postprocessing tool,
                            // defined in metagenomics profiling subworkflow

    main:
    ch_versions      = Channel.empty()
    ch_results       = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( params.metagenomics_postprocessing_tool == 'maltextract' ) {

        MALTEXTRACT ( ch_postprocessing_input, params.metagenomics_maltextract_taxon_list, params.metagenomics_maltextract_ncbi_dir )

        AMPS ( MALTEXTRACT.out.results, params.taxon_list, params.metagenomics_maltextract_filter )

        ch_versions      = ch_versions.mix( MALTEXTRACT.out.versions.first(), AMPS.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( AMPS.out.results.json )

    }

    elif ( params.metagenomics_postprocessing_tool == 'krakenmerge' || ['kraken2', 'krakenuniq'].contains(params.metagenomics_profiling_tool) ) {

        KRAKENPARSE ( ch_postprocessing_input )

        ch_list_of_kraken_parse_reads  = KRAKENPARSE.out.read_kraken_parsed.map {
            meta, read_out -> [ read_out ]
        }
        ch_list_of_kraken_parse_kmer = KRAKENPARSE.out.kmer_kraken_parsed.map {
            meta, kmer_out -> [ kmer_out ]
        }

        KRAKENMERGE ( ch_list_of_kraken_parse_reads.collect() , ch_list_of_kraken_parse_kmer.collect() )

        ch_versions      = ch_versions.mix( KRAKENPARSE.out.versions.first(), KRAKENMERGE.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( KRAKENMERGE.out.read_count_table, KRAKENMERGE.out.kmer_duplication_table )

    }

    elif ( params.metagenomics_postprocessing_tool == 'mergemetaphlantables' ) {
        METAPHLAN_MERGEMETAPHLANTABLES ( ch_postprocessing_input , params.metagenomics_profiling_database )

        ch_versions      = ch_versions.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.txt )

    }

    emit:
    versions = ch_versions
    mqc      = ch_multiqc_files

}

