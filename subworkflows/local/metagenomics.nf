include { METAGENOMICS_COMPLEXITYFILTER } from './metagenomics_complexityfilter'
include { METAGENOMICS_PROFILING        } from './metagenomics_profiling'
include { METAGENOMICS_POSTPROCESSING   } from './metagenomics_postprocessing'

workflow METAGENOMICS {
    take: ch_bamfiltered_for_metagenomics

    main:
    // Define channels
    ch_multiqc_files                = Channel.empty()
    ch_versions                     = Channel.empty()
    ch_bamfiltered_for_metagenomics = ch_bamfiltered_for_metagenomics
        .map{ meta, fastq ->
            [meta+['single_end':true], fastq]
        }

    //
    // Run the complexity filter subworkflow
    //

    if ( params.run_metagenomics_complexityfiltering ) {
        METAGENOMICS_COMPLEXITYFILTER( ch_bamfiltered_for_metagenomics )
        ch_reads_for_metagenomics = METAGENOMICS_COMPLEXITYFILTER.out.fastq
        ch_versions = ch_versions.mix(METAGENOMICS_COMPLEXITYFILTER.out.versions.first())
        ch_multiqc_files = ch_multiqc_files.mix(METAGENOMICS_COMPLEXITYFILTER.out.fastq.collect{it[1]}.ifEmpty([]))
    } else {
        ch_reads_for_metagenomics = ch_bamfiltered_for_metagenomics
    }

    //
    // Run the profiling subworkflow
    //

    database = Channel.fromPath(params.metagenomics_profiling_database)

    METAGENOMICS_PROFILING( ch_reads_for_metagenomics, database )

    ch_versions      = ch_versions.mix( METAGENOMICS_PROFILING.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( METAGENOMICS_PROFILING.out.mqc.collect{it[1]}.ifEmpty([]) )

    //
    // Run the post profiling subworkflow (optionally run for malt, mandatory for kraken2/krakenuniq)
    //

    if ( params.metagenomics_postprocessing_tool || ['kraken2', 'krakenuniq'].contains(params.metagenomics_profiling_tool) ) {

        METAGENOMICS_POSTPROCESSING ( METAGENOMICS_PROFILING.out.postprocessing_input )

        ch_versions      = ch_versions.mix( METAGENOMICS_POSTPROCESSING.out.versions.first() )
        ch_multiqc_files = ch_multiqc_files.mix( METAGENOMICS_POSTPROCESSING.out.mqc.collect{it[1]}.ifEmpty([]) )
    }


    emit:
    versions      = ch_versions
    ch_multiqc_files = ch_multiqc_files

}
