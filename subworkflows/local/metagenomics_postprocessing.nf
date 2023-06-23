// TODO: publish the files in ch_results directly?

include { MALTEXTRACT } from '../../modules/nf-core/maltextract/main'
include { AMPS        } from '../../modules/nf-core/amps/main'
include { KRAKENPARSE } from '../../modules/local/krakenparse'
include { KRAKENMERGE } from '../../modules/local/krakenmerge'

workflow METAGENOMICS_POSTPROCESSING {

    take:
    ch_postprocessing_input // different between kraken and malt

    main:
    ch_versions      = Channel.empty()
    ch_results       = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if ( params.metagenomics_postprocessing_tool == 'maltextract') {

        MALTEXTRACT ( ch_postprocessing_input, params.metagenomics_maltextract_taxon_list, params.metagenomics_maltextract_ncbi_dir )

        AMPS ( MALTEXTRACT.out.results, params.metagenomics_maltextract_taxon_list, params.metagenomics_maltextract_filter )

        ch_versions      = ch_versions.mix( MALTEXTRACT.out.versions.first(), AMPS.out.versions.first() )
        ch_results       = ch_results.mix( AMPS.out.candidate_pdfs, AMPS.out.tsv, AMPS.out.summary_pdf )
        ch_multiqc_files = ch_multiqc_files.mix( AMPS.out.json )

    }

    else {

        KRAKENPARSE ( ch_postprocessing_input )

        ch_list_of_kraken_parse_reads  = KRAKENPARSE.out.read_kraken_parsed.map {
            meta, read_out -> [ read_out ]
        }
        ch_list_of_kraken_parse_kmer = KRAKENPARSE.out.kmer_kraken_parsed.map {
            meta, kmer_out -> [ kmer_out ]
        }

        KRAKENMERGE ( ch_list_of_kraken_parse_reads.collect() , ch_list_of_kraken_parse_kmer.collect() )

        ch_versions      = ch_versions.mix( KRAKENPARSE.out.versions.first(), KRAKENMERGE.out.versions.first() )
        ch_results       = ch_results.mix( KRAKENMERGE.out.read_count_table, KRAKENMERGE.out.kmer_duplication_table )
        ch_multiqc_files = ch_multiqc_files.mix( KRAKENMERGE.out.read_count_table, KRAKENMERGE.out.kmer_duplication_table )

    }

    emit:
    versions = ch_versions
    results  = ch_results
    mqc      = ch_multiqc_files

}

