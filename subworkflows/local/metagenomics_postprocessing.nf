include { MALTEXTRACT                    } from '../../modules/nf-core/maltextract/main'
include { AMPS                           } from '../../modules/nf-core/amps/main'
include { KRAKENPARSE                    } from '../../modules/local/krakenparse'
include { KRAKENMERGE                    } from '../../modules/local/krakenmerge'
include { METAPHLAN_MERGEMETAPHLANTABLES } from '../../modules/nf-core/metaphlan/mergemetaphlantables/main'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    ch_postprocessing_input // different between each profiling --> postprocessing tool,
                            // defined in metagenomics profiling subworkflow

    main:
    ch_versions      = Channel.empty()
    ch_results       = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( params.metagenomics_postprocessing_tool == 'maltextract' ) {

        //maltextract doesnt accepts a meta param in the first input channel, so remove it
        ch_maltextract_input = ch_postprocessing_input.map{it[1]}

        tax_list = Channel.fromPath(params.metagenomics_maltextract_taxon_list)
        ncbi_dir = Channel.fromPath(params.metagenomics_maltextract_ncbi_dir)
        maltex_filter = Channel.fromPath(params.metagenomics_maltextract_filter)

        MALTEXTRACT ( ch_maltextract_input, tax_list, ncbi_dir)

        AMPS ( MALTEXTRACT.out.results, tax_list, maltex_filter )

        ch_versions      = ch_versions.mix( MALTEXTRACT.out.versions.first(), AMPS.out.versions.first() )
        ch_results       = ch_results.mix( AMPS.out.candidate_pdfs, AMPS.out.tsv, AMPS.out.summary_pdf )
        ch_multiqc_files = ch_multiqc_files.mix( AMPS.out.json )

    }

    else if ( ['kraken2', 'krakenuniq'].contains(params.metagenomics_profiling_tool) ) {

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

    else if ( params.metagenomics_postprocessing_tool == 'mergemetaphlantables' ) {

        METAPHLAN_MERGEMETAPHLANTABLES ( ch_postprocessing_input.map{ [[id:"metaphlan_profiles_all_samples_merged"], it[1]] }.groupTuple() )

        ch_versions      = ch_versions.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( METAPHLAN_MERGEMETAPHLANTABLES.out.txt )

    }

    emit:
    versions = ch_versions
    mqc      = ch_multiqc_files

}

