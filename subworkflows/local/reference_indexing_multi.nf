//
// Index input reference as required
//

include { GUNZIP as GUNZIP_FASTA          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
include { MAPAD_INDEX                     } from '../../modules/nf-core/mapad/index/main'
include { samplesheetToList               } from 'plugin/nf-schema'

workflow REFERENCE_INDEXING_MULTI {
    take:
    referencesheet // file: /path/to/name.{csv,tsv}

    main:
    ch_versions = Channel.empty()

    // Import reference sheet and change empty arrays to empty strings for compatibility with single reference input
    ch_splitreferencesheet_for_branch = Channel
        .fromList(samplesheetToList(referencesheet, "${projectDir}/assets/schema_fasta.json"))
        .map { meta, fasta, fai, dict, mapper_index, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex, mitochondrion, capture_bed, pileupcaller_bed, pileupcaller_snp, hapmap, pmd_masked_fasta, pmd_bed_for_masking, sexdet_bed, bedtools_feature, genotyping_gatk_dbsnp ->
            meta.ploidy = meta.genotyping_ploidy != null ? meta.genotyping_ploidy : params.genotyping_reference_ploidy
            fai = fai != [] ? fai : ""
            dict = dict != [] ? dict : ""
            mapper_index = mapper_index != [] ? mapper_index : ""
            circular_target = circular_target != [] ? circular_target : ""
            circularmapper_elongatedfasta = circularmapper_elongatedfasta != [] ? circularmapper_elongatedfasta : ""
            circularmapper_elongatedindex = circularmapper_elongatedindex != [] ? circularmapper_elongatedindex : ""
            mitochondrion = mitochondrion != [] ? mitochondrion : ""
            capture_bed = capture_bed != [] ? capture_bed : ""
            pileupcaller_bed = pileupcaller_bed != [] ? pileupcaller_bed : ""
            pileupcaller_snp = pileupcaller_snp != [] ? pileupcaller_snp : ""
            hapmap = hapmap != [] ? hapmap : ""
            pmd_masked_fasta = pmd_masked_fasta != [] ? pmd_masked_fasta : ""
            pmd_bed_for_masking = pmd_bed_for_masking != [] ? pmd_bed_for_masking : ""
            sexdet_bed = sexdet_bed != [] ? sexdet_bed : ""
            bedtools_feature = bedtools_feature != [] ? bedtools_feature : ""
            genotyping_gatk_dbsnp = genotyping_gatk_dbsnp != [] ? genotyping_gatk_dbsnp : ""
            [meta - meta.subMap('genotyping_ploidy'), fasta, fai, dict, mapper_index, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex, mitochondrion, capture_bed, pileupcaller_bed, pileupcaller_snp, hapmap, pmd_masked_fasta, pmd_bed_for_masking, sexdet_bed, bedtools_feature, genotyping_gatk_dbsnp]
        }

    // GENERAL DESCRIPTION FOR NEXT SECTIONS
    // This will be the same scheme for all other generation steps, i.e.
    // for those that need to be processed send them to a  multiMap (others skip via branch)
    // take only those files needed for the generation step in one sub-channel, and all other
    // channel elements go into 'remainder'. After generation, join these back, and then mix back the
    // ones that don't skip

    //
    // DECOMPRESSION
    //

    ch_input_from_referencesheet = ch_splitreferencesheet_for_branch.multiMap { meta, fasta, fai, dict, mapper_index, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex, mitochondrion, capture_bed, pileupcaller_bed, pileupcaller_snp, hapmap, pmd_masked_fasta, pmd_bed_for_masking, sexdet_bed, bedtools_feature, genotyping_gatk_dbsnp ->
        generated: [meta, fasta, fai, dict, mapper_index]
        circularmapper: [meta, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex]
        mitochondrion_header: [meta, mitochondrion]
        angsd_hapmap: [meta, hapmap]
        pmd_masked_fasta: [meta, pmd_masked_fasta]
        pmd_bed_for_masking: [meta, pmd_bed_for_masking]
        snp_bed: [meta, capture_bed]
        pileupcaller_bed_snp: [meta, pileupcaller_bed, pileupcaller_snp]
        sexdeterrmine_bed: [meta, sexdet_bed]
        bedtools_feature: [meta, bedtools_feature]
        dbsnp: [meta, genotyping_gatk_dbsnp]
    }

    // Detect if fasta is gzipped or not
    ch_fasta_for_gunzip = ch_input_from_referencesheet.generated.branch { meta, fasta, fai, dict, mapper_index ->
        forgunzip: fasta.extension == "gz"
        skip: true
    }

    // Pull out name/file to match cardinality for GUNZIP module
    ch_gunzip_input = ch_fasta_for_gunzip.forgunzip.multiMap { meta, fasta, fai, dict, mapper_index ->
        gunzip: [meta, fasta]
        remainder: [meta, fai, dict, mapper_index]
    }


    GUNZIP_FASTA(ch_gunzip_input.gunzip)
    ch_version = ch_versions.mix(GUNZIP_FASTA.out.versions)

    // Mix back gunzipped fasta with remaining files, and then mix back with pre-gunzipped references
    ch_gunzippedfasta_formix = GUNZIP_FASTA.out.gunzip.join(ch_gunzip_input.remainder, failOnMismatch: true)
    ch_fasta_for_faiindexing = ch_fasta_for_gunzip.skip.mix(ch_gunzippedfasta_formix)

    //
    // INDEXING: fai
    //

    // Separate out non-faidxed references
    ch_fasta_for_faidx = ch_fasta_for_faiindexing.branch { meta, fasta, fai, dict, mapper_index ->
        forfaidx: fai == ""
        skip: true
    }

    // Split channel to ensure cardinality matching
    ch_faidx_input = ch_fasta_for_faidx.forfaidx.multiMap { meta, fasta, fai, dict, mapper_index ->
        faidx: [meta, fasta]
        remainder: [meta, fasta, dict, mapper_index]
    }

    SAMTOOLS_FAIDX(ch_faidx_input.faidx, [[], []])
    ch_version = ch_versions.mix(SAMTOOLS_FAIDX.out.versions)

    // Rejoin output channel with main reference indicies channel elements
    ch_faidxed_formix = SAMTOOLS_FAIDX.out.fai
        .join(ch_faidx_input.remainder, failOnMismatch: true)
        .map { meta, fai, fasta, dict, mapper_index ->

            [meta, fasta, fai, dict, mapper_index]
        }

    // Mix back newly faidx'd references with the pre-indexed ones
    ch_fasta_for_dictindexing = ch_fasta_for_faidx.skip.mix(ch_faidxed_formix)

    //
    // INDEXING: dict
    //

    ch_fasta_for_dict = ch_fasta_for_dictindexing.branch { meta, fasta, fai, dict, mapper_index ->
        fordict: dict == ""
        skip: true
    }

    ch_dict_input = ch_fasta_for_dict.fordict.multiMap { meta, fasta, fai, dict, mapper_index ->
        dict: [meta, fasta]
        remainder: [meta, fasta, fai, mapper_index]
    }

    PICARD_CREATESEQUENCEDICTIONARY(ch_dict_input.dict)
    ch_version = ch_versions.mix(PICARD_CREATESEQUENCEDICTIONARY.out.versions)

    ch_dicted_formix = PICARD_CREATESEQUENCEDICTIONARY.out.reference_dict
        .join(ch_dict_input.remainder, failOnMismatch: true)
        .map { meta, dict, fasta, fai, mapper_index ->

            [meta, fasta, fai, dict, mapper_index]
        }

    ch_dict_formapperindexing = ch_fasta_for_dict.skip.mix(ch_dicted_formix)

    //
    // INDEXING: Mapping indicies
    //

    // Generate mapper indicies if not supplied, and if supplied generate meta

    ch_fasta_for_mapperindex = ch_dict_formapperindexing.branch { meta, fasta, fai, dict, mapper_index ->
        forindex: mapper_index == ""
        skip: true
    }

    ch_mapindex_input = ch_fasta_for_mapperindex.forindex.multiMap { meta, fasta, fai, dict, mapper_index ->
        index: [meta, fasta]
        remainder: [meta, fasta, fai, dict]
    }

    if (params.mapping_tool == "bwaaln" || params.mapping_tool == "bwamem" || params.mapping_tool == "circularmapper") {
        BWA_INDEX(ch_mapindex_input.index)
        ch_version = ch_versions.mix(BWA_INDEX.out.versions)
        ch_indexed_forremap = BWA_INDEX.out.index
    }
    else if (params.mapping_tool == "bowtie2") {
        BOWTIE2_BUILD(ch_mapindex_input.index)
        ch_version = ch_versions.mix(BOWTIE2_BUILD.out.versions)
        ch_indexed_forremap = BOWTIE2_BUILD.out.index
    } else if (params.mapping_tool == "mapad") {
        MAPAD_INDEX (ch_mapindex_input.index)
        ch_version = ch_versions.mix( MAPAD_INDEX.out.versions )
        ch_indexed_forremap = MAPAD_INDEX.out.index
    }

    ch_indexed_formix = ch_indexed_forremap
        .join(ch_mapindex_input.remainder, failOnMismatch: true)
        .map { meta, mapper_index, fasta, fai, dict ->

            [meta, fasta, fai, dict, mapper_index]
        }

    ch_indexmapper_for_reference = ch_fasta_for_mapperindex.skip.mix(ch_indexed_formix)

    emit:
    reference            = ch_indexmapper_for_reference // [ meta, fasta, fai, dict, mapindex ]
    elongated_reference  = ch_input_from_referencesheet.circularmapper // [ meta, circular_target, circularmapper_elongatedfasta, circularmapper_elongatedindex ]
    mitochondrion_header = ch_input_from_referencesheet.mitochondrion_header // [ meta, mitochondrion ]
    hapmap               = ch_input_from_referencesheet.angsd_hapmap // [ meta, hapmap ]
    pmd_masked_fasta     = ch_input_from_referencesheet.pmd_masked_fasta // [ meta, pmd_masked_fasta ]
    pmd_bed_for_masking  = ch_input_from_referencesheet.pmd_bed_for_masking // [ meta, pmd_bed_for_masking ]
    snp_capture_bed      = ch_input_from_referencesheet.snp_bed // [ meta, capture_bed ]
    pileupcaller_bed_snp = ch_input_from_referencesheet.pileupcaller_bed_snp // [ meta, pileupcaller_bed, pileupcaller_snp ]
    sexdeterrmine_bed    = ch_input_from_referencesheet.sexdeterrmine_bed // [ meta, sexdet_bed ]
    bedtools_feature     = ch_input_from_referencesheet.bedtools_feature // [ meta, bedtools_feature ]
    dbsnp                = ch_input_from_referencesheet.dbsnp // [ meta, genotyping_gatk_dbsnp ]
    versions             = ch_versions
}
