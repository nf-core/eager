//
// Prepare reference indexing for downstream
//

include { SEQKIT_SPLIT2                                 } from '../../modules/nf-core/seqkit/split2/main'
include { FASTQ_ALIGN_BWAALN                            } from '../../subworkflows/nf-core/fastq_align_bwaaln/main'
include { BWA_MEM                                       } from '../../modules/nf-core/bwa/mem/main'
include { BOWTIE2_ALIGN                                 } from '../../modules/nf-core/bowtie2/align/main'
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LANES        } from '../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_SORT  as SAMTOOLS_SORT_MERGED_LANES  } from '../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MEM          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_BT2          } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_INDEX as SAMTOOLS_INDEX_MERGED_LANES } from '../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_FLAGSTAT as SAMTOOLS_FLAGSTAT_MAPPED } from '../../modules/nf-core/samtools/flagstat/main'
include { CIRCULARMAPPER                                } from '../../subworkflows/local/circularmapper'

workflow MAP {
    take:
    reads          // [ [meta], [read1, reads2] ] or [ [meta], [read1] ]
    index          // [ [meta], [ index ], [ fasta ] ]
    elogated_index // [ [meta], [ index ], [ fasta ], [ circular_target ] ]

    main:
    ch_versions       = Channel.empty()
    ch_multiqc_files  = Channel.empty()

    if ( params.run_fastq_sharding ) {

        ch_input_for_sharding = reads

        SEQKIT_SPLIT2( ch_input_for_sharding )
        ch_versions        = ch_versions.mix ( SEQKIT_SPLIT2.out.versions.first() )

        sharded_reads = SEQKIT_SPLIT2.out.reads
            .transpose()
            .map {
                meta, reads ->
                    new_meta = meta.clone()
                    new_meta.shard_number = reads.getName().replaceAll(/.*(part_\d+).(?:fastq|fq).gz/, '$1')
                    [ new_meta, reads ]
            }
            .groupTuple()

        ch_input_for_mapping = sharded_reads
            .combine(index.map{ meta, index, fasta -> [ meta, index ] })
            .multiMap {
                meta, reads, meta2, index ->
                    new_meta = meta.clone()
                    new_meta.reference = meta2.id
                    reads: [ new_meta, reads ]
                    index: [ meta2, index ]
            }

    } else {
        ch_input_for_mapping = reads
            .combine(index.map{ meta, index, fasta -> [ meta, index ] })
            .multiMap {
                meta, reads, meta2, index ->
                    new_meta = meta.clone()
                    new_meta.reference = meta2.id
                    reads: [ new_meta, reads ]
                    index: [ meta2, index ]
            }
    }

    if ( params.mapping_tool == 'bwaaln' ) {
        ch_index_for_mapping = index.map{ meta, index, fasta -> [ meta, index ] }
        ch_reads_for_mapping = reads

        FASTQ_ALIGN_BWAALN ( ch_reads_for_mapping, ch_index_for_mapping )
        ch_versions        = ch_versions.mix ( FASTQ_ALIGN_BWAALN.out.versions.first() )
        ch_mapped_lane_bam = FASTQ_ALIGN_BWAALN.out.bam
                                .map{
                                    // create meta consistent with rest of workflow
                                    meta, bam ->
                                    new_meta = meta + [ reference: meta.id_index ]
                                [ new_meta, bam ]
                                }

        ch_mapped_lane_bai = params.fasta_largeref ? FASTQ_ALIGN_BWAALN.out.csi : FASTQ_ALIGN_BWAALN.out.bai

    } else if ( params.mapping_tool == 'bwamem' ) {
        ch_input_for_mapping = reads
                            .combine( index )
                            .multiMap {
                                meta, reads, meta2, index, fasta ->
                                    new_meta = meta + [ reference: meta2.id ]
                                    reads: [ new_meta, reads ]
                                    index: [ meta2, index ]
                                    fasta: [ meta2, fasta ]
                            }

        BWA_MEM ( ch_input_for_mapping.reads, ch_input_for_mapping.index, ch_input_for_mapping.fasta, true )
        ch_versions        = ch_versions.mix ( BWA_MEM.out.versions.first() )
        ch_mapped_lane_bam = BWA_MEM.out.bam

        SAMTOOLS_INDEX_MEM ( ch_mapped_lane_bam )
        ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_MEM.out.versions.first())
        ch_mapped_lane_bai = params.fasta_largeref ? SAMTOOLS_INDEX_MEM.out.csi : SAMTOOLS_INDEX_MEM.out.bai

    } else if ( params.mapping_tool == 'bowtie2' ) {
        ch_input_for_mapping = reads
                            .combine( index.map{ meta, index, fasta -> [ meta, index ] } )
                            .multiMap {
                                meta, reads, meta2, index ->
                                    new_meta = meta + [ reference: meta2.id ]
                                    reads: [ new_meta, reads ]
                                    index: [ meta2, index ]
                            }

        BOWTIE2_ALIGN ( ch_input_for_mapping.reads, ch_input_for_mapping.index, false, true )
        ch_versions        = ch_versions.mix ( BOWTIE2_ALIGN.out.versions.first() )
        ch_mapped_lane_bam = BOWTIE2_ALIGN.out.aligned

        SAMTOOLS_INDEX_BT2 ( ch_mapped_lane_bam )
        ch_versions        = ch_versions.mix(SAMTOOLS_INDEX_BT2.out.versions.first())
        ch_mapped_lane_bai = params.fasta_largeref ? SAMTOOLS_INDEX_BT2.out.csi : SAMTOOLS_INDEX_BT2.out.bai

    } else if ( params.mapping_tool == 'circularmapper' ) {
        // Reference elongation and indexing takes place in the reference_indexing swf.
        // Circularmapper takes non-elongated AND elongated references and reads as input (i think. wait for Alex's reply).

        // ch_input_for_circularmapper = reads
        //                             .combine(index.map{ meta, index, fasta -> [ meta, fasta ] })
        //                             .dump(tag:"CM Inputs", pretty:true)
        //                             .multiMap {
        //                                 meta, reads, meta2, fasta ->
        //                                     reads: [ meta, reads ]
        //                                     reference: [ meta2, fasta ]
        //                             }
        // CIRCULARMAPPER( ch_input_for_circularmapper.reads, params.elongation_factor, ch_input_for_circularmapper.reference )
        // ch_versions        = ch_versions.mix ( CIRCULARMAPPER.out.versions )
        // // TODO - Update SWF outputs
        // ch_mapped_lane_bam      = CIRCULARMAPPER.out.bam
        // ch_mapped_lane_bai      = Channel.empty() // Circularmapper doesn't give a bai


    }


    // Only run merge lanes if we have more than one BAM to merge!
    ch_input_for_lane_merge = ch_mapped_lane_bam
                                .map {
                                    meta, bam ->
                                    new_meta = meta.clone().findAll{ it.key !in ['lane', 'colour_chemistry', 'shard_number'] }
                                    [ new_meta, bam ]
                                }
                                .groupTuple()
                                .branch {
                                    meta, bam ->
                                        merge: bam.size() > 1
                                        skip: true
                                }

    SAMTOOLS_MERGE_LANES ( ch_input_for_lane_merge.merge, [[], []], [[], []] )
    ch_versions.mix( SAMTOOLS_MERGE_LANES.out.versions )

    // Then mix back merged and single lane libraries for everything downstream
    ch_mergedlanes_for_sorting = ch_input_for_lane_merge.skip
                                    .mix( SAMTOOLS_MERGE_LANES.out.bam )

    SAMTOOLS_SORT_MERGED_LANES ( ch_mergedlanes_for_sorting )
    ch_mapped_bam = SAMTOOLS_SORT_MERGED_LANES.out.bam
    ch_versions.mix( SAMTOOLS_SORT_MERGED_LANES.out.versions )

    SAMTOOLS_INDEX_MERGED_LANES( ch_mapped_bam )
    ch_mapped_bai =  params.fasta_largeref ? SAMTOOLS_INDEX_MERGED_LANES.out.csi : SAMTOOLS_INDEX_MERGED_LANES.out.bai
    ch_versions.mix( SAMTOOLS_INDEX_MERGED_LANES.out.versions )

    ch_input_for_flagstat = SAMTOOLS_SORT_MERGED_LANES.out.bam.join( SAMTOOLS_INDEX_MERGED_LANES.out.bai, failOnMismatch: true )

    SAMTOOLS_FLAGSTAT_MAPPED ( ch_input_for_flagstat )
    ch_versions.mix( SAMTOOLS_FLAGSTAT_MAPPED.out.versions.first() )
    ch_multiqc_files = ch_multiqc_files.mix( SAMTOOLS_FLAGSTAT_MAPPED.out.flagstat )

    emit:
    bam        = ch_mapped_bam                            // [ [ meta ], bam ]
    bai        = ch_mapped_bai                            // [ [ meta ], bai ]
    flagstat   = SAMTOOLS_FLAGSTAT_MAPPED.out.flagstat    // [ [ meta ], stats ]
    mqc        = ch_multiqc_files
    versions   = ch_versions

}
