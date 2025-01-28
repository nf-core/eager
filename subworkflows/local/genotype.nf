//
// Genotype the input data using the requested genotyper.
//

include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_PILEUPCALLER } from '../../modules/nf-core/samtools/mpileup/main'
include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE     } from '../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main'
include { SEQUENCETOOLS_PILEUPCALLER                        } from '../../modules/nf-core/sequencetools/pileupcaller/main'
include { COLLECT_GENOTYPES                                 } from '../../modules/local/collect_genotypes'
include { GATK_REALIGNERTARGETCREATOR                       } from '../../modules/nf-core/gatk/realignertargetcreator/main'
include { GATK_INDELREALIGNER                               } from '../../modules/nf-core/gatk/indelrealigner/main'
include { GATK_UNIFIEDGENOTYPER                             } from '../../modules/nf-core/gatk/unifiedgenotyper/main'
include { GATK4_HAPLOTYPECALLER                             } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { FREEBAYES                                         } from '../../modules/nf-core/freebayes/main'
include { ANGSD_GL                                          } from '../../modules/nf-core/angsd/gl/main'
include { BCFTOOLS_STATS as BCFTOOLS_STATS_GENOTYPING       } from '../../modules/nf-core/bcftools/stats/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_UG               } from '../../modules/nf-core/bcftools/index/main'
include { BCFTOOLS_INDEX as BCFTOOLS_INDEX_FREEBAYES        } from '../../modules/nf-core/bcftools/index/main'
include { addNewMetaFromAttributes                          } from '../../subworkflows/local/utils_nfcore_eager_pipeline/main'

workflow GENOTYPE {
    take:
    ch_bam_bai                // [ [ meta ], bam , bai ]
    ch_fasta_plus             // [ [ meta ], fasta, fai, dict ]
    ch_pileupcaller_aux_files // [ [ meta ], bed, snp ]
    ch_dbsnp                  // [ [ meta ], dbsnp ]

    main:
    ch_versions                        = Channel.empty()
    ch_multiqc_files                   = Channel.empty()
    ch_pileupcaller_genotypes          = Channel.empty()
    ch_eigenstrat_coverage_stats       = Channel.empty()
    ch_angsd_genotype_likelihoods      = Channel.empty()
    ch_genotypes_vcf                   = Channel.empty()
    ch_bcftools_stats                  = Channel.empty()

    if ( params.genotyping_tool == 'pileupcaller' ) {

        // Compile together all reference based files
        ch_refs_prep = ch_fasta_plus
            // Because aux files are optional, the channel can be [[],[],[]]. remainder:true will output both the empty list and the fasta_plus channel with an added 'null'.
            .join( ch_pileupcaller_aux_files, remainder: true ) // [ [ref_meta], fasta, fai, dict, bed, snp ]
            // Also filter out the empty list aux_files (meta == [])
            .filter { it[0] != [] }
            // .branch to separate succesfully joined from unsuccesfully joined elements
            .branch {
                it ->
                no_aux: it[4] == null
                has_aux: true
            }

        // mix the two branches back together after fixing cardinality
        ch_refs_for_mpileup_pileupcaller = ch_refs_prep.no_aux
            .map {
                ref_meta, fasta, fai, dict, empty ->
                [ ref_meta, fasta, fai, dict, [], [] ]
            }
            .mix( ch_refs_prep.has_aux )
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "id" , "reference" , false )
            } // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, bed, snp ]

        // Prepare collect bams for mpileup
        ch_mpileup_inputs_bams = ch_bam_bai
            .map {
                addNewMetaFromAttributes( it, ["reference", "strandedness"] , ["reference", "strandedness"] , false )
            }
            .groupTuple()
            .map {
                combo_meta, metas, bams, bais ->
                def ids = metas.collect { meta -> meta.sample_id }
                [ combo_meta + [sample_id: ids], bams ] // Drop bais
            } // Collect all IDs into a list in meta.sample_id. Useful when running pileupCaller later

            // Combine prepped bams and references
            ch_mpileup_inputs = ch_mpileup_inputs_bams
                .map {
                    addNewMetaFromAttributes( it, "reference", "reference" , false )
                }
                .combine( ch_refs_for_mpileup_pileupcaller , by:0 )
                // do not run if no bed file is provided
                .filter { it[7] != []}
                .multiMap {
                    ignore_me, combo_meta, bams, ref_meta, fasta, fai, dict, bed, snp ->
                        def bedfile = bed != "" ? bed : []
                        bams:  [ combo_meta, bams, bedfile ]
                        fasta: [ fasta ]
                }

            SAMTOOLS_MPILEUP_PILEUPCALLER(
                ch_mpileup_inputs.bams,
                ch_mpileup_inputs.fasta,
            )
            ch_versions = ch_versions.mix( SAMTOOLS_MPILEUP_PILEUPCALLER.out.versions.first() )

            ch_pileupcaller_input = SAMTOOLS_MPILEUP_PILEUPCALLER.out.mpileup
                .map {
                    addNewMetaFromAttributes( it, "reference", "reference" , false )
                }
                .combine( ch_refs_for_mpileup_pileupcaller, by:0 )
                .multiMap {
                    ignore_me, meta, mpileup, ref_meta, fasta, fai, dict, bed, snp ->
                        // def snpfile = snp != "" ? snp : []
                        mpileup: [ meta, mpileup ]
                        snpfile: snp
                }

            // Run PileupCaller
            SEQUENCETOOLS_PILEUPCALLER(
                ch_pileupcaller_input.mpileup,
                ch_pileupcaller_input.snpfile,
                []
            )
            ch_versions = ch_versions.mix( SEQUENCETOOLS_PILEUPCALLER.out.versions.first() )

            // Merge/rename genotyping datasets
            ch_final_genotypes = SEQUENCETOOLS_PILEUPCALLER.out.eigenstrat
                .map {
                    addNewMetaFromAttributes( it, "reference" , "reference" , false )
                }
                .groupTuple()
                .map {
                    combo_meta, metas, geno, snp, ind ->
                    [ combo_meta, geno, snp, ind ]
                }

            COLLECT_GENOTYPES( ch_final_genotypes )
            // Add genotyper info to the meta
            ch_pileupcaller_genotypes = COLLECT_GENOTYPES.out.collected
            .map {
                meta, geno, snp, ind ->
                [ meta + [ genotyper: "pileupcaller" ], geno , snp, ind ]
            }
            ch_versions               = ch_versions.mix( COLLECT_GENOTYPES.out.versions.first() )

            // Calculate coverage stats for collected eigenstrat dataset
            EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE(
                ch_pileupcaller_genotypes
            )
            ch_eigenstrat_coverage_stats = EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE.out.tsv
            ch_versions                  = ch_versions.mix( EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE.out.versions.first() )
            ch_multiqc_files             = ch_multiqc_files.mix( EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE.out.json )
    }

    if ( params.genotyping_tool == 'ug' ) {
        // Use correct reference for each input bam/bai pair.
        ch_bams_for_multimap = ch_bam_bai
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }

        ch_fasta_for_multimap = ch_fasta_plus
            // Because dbsnp is optional, the channel can be [[],[]]. remainder:true will output both the empty list and the fasta_plus channel with an added 'null'.
            .join( ch_dbsnp, remainder:true ) // [ [ref_meta], fasta, fai, dict, dbsnp ]
            // Also filter out the empty list dbsnp (meta == [])
            .filter { it[0] != [] }
            // convert added null dbsnp into an empty list
            .map {
                ref_meta, fasta, fai, dict, dbsnp ->
                def final_dbsnp = dbsnp != null ? dbsnp : []
                [ ref_meta, fasta, fai, dict, final_dbsnp ]
            }
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "id" , "reference" , false )
            } // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, dbsnp ]

        ch_input_for_targetcreator = ch_bams_for_multimap
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam , bai ]
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
                    dict:  [ ref_meta, dict ]
            }

        GATK_REALIGNERTARGETCREATOR(
            ch_input_for_targetcreator.bam,
            ch_input_for_targetcreator.fasta,
            ch_input_for_targetcreator.fai,
            ch_input_for_targetcreator.dict,
            [[], []] // No known_vcf
        )
        ch_versions = ch_versions.mix( GATK_REALIGNERTARGETCREATOR.out.versions.first() )

        // Join the bam/bai pairs to the intervals file, then redo multiMap to get the correct ordering for each bam/reference/intervals set.
        ch_input_for_indelrealigner = ch_bam_bai
            .join( GATK_REALIGNERTARGETCREATOR.out.intervals )
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, intervals, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam, bai, intervals ]
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
                    dict:  [ ref_meta, dict ]
            }

        GATK_INDELREALIGNER(
            ch_input_for_indelrealigner.bam,
            ch_input_for_indelrealigner.fasta,
            ch_input_for_indelrealigner.fai,
            ch_input_for_indelrealigner.dict,
            [[], []] // No known_vcf
        )
        ch_versions = ch_versions.mix( GATK_INDELREALIGNER.out.versions.first() )

        // Use realigned bams as input for UG. combine with reference info to get correct ordering.
        ch_bams_for_ug = GATK_INDELREALIGNER.out.bam
            .map {
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam, bai ]
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
                    dict:  [ ref_meta, dict ]
                    dbsnp: [ ref_meta, dbsnp ]
            }

        GATK_UNIFIEDGENOTYPER(
            ch_bams_for_ug.bam,
            ch_bams_for_ug.fasta,
            ch_bams_for_ug.fai,
            ch_bams_for_ug.dict,
            [[], []], // No intervals
            [[], []], // No contamination
            ch_bams_for_ug.dbsnp,
            [[], []]  // No comp
        )
        ch_gatk_ug_vcf = GATK_UNIFIEDGENOTYPER.out.vcf
        ch_versions                        = ch_versions.mix( GATK_UNIFIEDGENOTYPER.out.versions.first() )

        // Index the VCFs
        BCFTOOLS_INDEX_UG( ch_gatk_ug_vcf )
        ch_versions = ch_versions.mix( BCFTOOLS_INDEX_UG.out.versions.first() )

        // Add genotyper info to the meta
        ch_genotypes_vcf = ch_gatk_ug_vcf.join(BCFTOOLS_INDEX_UG.out.tbi)
            .map {
                meta, vcf, tbi ->
                [ meta + [ genotyper: "ug" ], vcf , tbi ]
            }
    }

    if ( params.genotyping_tool == 'hc' ) {
        ch_bams_for_multimap = ch_bam_bai
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }

        ch_fasta_for_multimap = ch_fasta_plus
            // Because dbsnp is optional, the channel can be [[],[]]. remainder:true will output both the empty list and the fasta_plus channel with an added 'null'.
            .join( ch_dbsnp, remainder:true ) // [ [ref_meta], fasta, fai, dict, dbsnp ]
            // Also filter out the empty list dbsnp (meta == [])
            .filter { it[0] != [] }
            // convert added null dbsnp into an empty list
            .map {
                ref_meta, fasta, fai, dict, dbsnp ->
                def final_dbsnp = dbsnp != null ? dbsnp : []
                [ ref_meta, fasta, fai, dict, final_dbsnp ]
            }
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "id" , "reference" , false )
            } // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, dbsnp ]

        ch_input_for_hc = ch_bams_for_multimap
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam , bai, [], [] ] // No intervals or dragSTR model inputs to HC module
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
                    dict:  [ ref_meta, dict ]
                    dbsnp: [ ref_meta, dbsnp ]
            }

        GATK4_HAPLOTYPECALLER(
            ch_input_for_hc.bam,
            ch_input_for_hc.fasta,
            ch_input_for_hc.fai,
            ch_input_for_hc.dict,
            ch_input_for_hc.dbsnp,
            [[], []] // No dbsnp_tbi
        )
        // Add genotyper info to the meta
        ch_genotypes_vcf = GATK4_HAPLOTYPECALLER.out.vcf.join( GATK4_HAPLOTYPECALLER.out.tbi )
            .map {
                meta, vcf, tbi ->
                [ meta + [ genotyper: "hc" ], vcf , tbi ]
            }
        ch_versions                       = ch_versions.mix( GATK4_HAPLOTYPECALLER.out.versions.first() )
    }

    if ( params.genotyping_tool == 'freebayes' ) {
        ch_bams_for_multimap = ch_bam_bai
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }

        // NOTE: dbsnp is not used by Freebayes, but we need to provide it to the module anyway, to ensure correct cardinality of the fasta channel within the BCFTOOLS_STATS channel operations.
        //     i.e. to keep the definition of the ch_fasta_for_multimap channel consistent regardless of genotyper, so the `combine -> multiMap` in lines 327-328 work.
        ch_fasta_for_multimap = ch_fasta_plus
            // Because dbsnp is optional, the channel can be [[],[]]. remainder:true will output both the empty list and the fasta_plus channel with an added 'null'.
            .join( ch_dbsnp, remainder:true ) // [ [ref_meta], fasta, fai, dict, dbsnp ]
            // Also filter out the empty list dbsnp (meta == [])
            .filter { it[0] != [] }
            // convert added null dbsnp into an empty list
            .map {
                ref_meta, fasta, fai, dict, dbsnp ->
                def final_dbsnp = dbsnp != null ? dbsnp : []
                [ ref_meta, fasta, fai, dict, final_dbsnp ]
            }
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "id" , "reference" , false )
            } // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, dbsnp ]

        ch_input_for_freebayes = ch_bams_for_multimap
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam , bai, [], [], [] ] // No second bam, second bai, or regions-BED file
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
            }

        FREEBAYES(
            ch_input_for_freebayes.bam,
            ch_input_for_freebayes.fasta,
            ch_input_for_freebayes.fai,
            [ [], [] ], // No samples file
            [ [], [] ], // No populations file
            [ [], [] ]  // No CNV file
            )
        ch_freebayes_genotypes = FREEBAYES.out.vcf
        ch_versions            = ch_versions.mix( FREEBAYES.out.versions.first() )

        // Index the VCFs
        BCFTOOLS_INDEX_FREEBAYES( ch_freebayes_genotypes )
        ch_versions = ch_versions.mix( BCFTOOLS_INDEX_FREEBAYES.out.versions.first() )

        // Add genotyper info to the meta
        ch_genotypes_vcf = ch_freebayes_genotypes.join(BCFTOOLS_INDEX_FREEBAYES.out.tbi)
            .map {
                meta, vcf, tbi ->
                [ meta + [ genotyper: "freebayes" ], vcf , tbi ]
            }
    }

    if ( params.genotyping_tool == 'angsd' ) {
        ch_bams_for_multimap = ch_bam_bai
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .groupTuple()
                .map {
                    combo_meta, metas, bams, bais ->
                    def new_map = [:]
                    def ids = metas.collect { meta -> meta.sample_id }
                    def strandedness = metas.collect { meta -> meta.strandedness }
                    def single_ends = metas.collect { meta -> meta.single_end }
                    def reference = combo_meta.reference
                    new_meta = [ sample_id: ids, strandedness: strandedness, single_end: single_ends, reference: reference ]

                    [ combo_meta, new_meta, bams, bais ] // Drop bais
                } // Collect all IDs into a list in meta.sample_id.

        ch_fasta_for_multimap = ch_fasta_plus
            // Because dbsnp is optional, the channel can be [[],[]]. remainder:true will output both the empty list and the fasta_plus channel with an added 'null'.
            .join( ch_dbsnp, remainder:true ) // [ [ref_meta], fasta, fai, dict, dbsnp ]
            // Also filter out the empty list dbsnp (meta == [])
            .filter { it[0] != [] }
            // convert added null dbsnp into an empty list
            .map {
                ref_meta, fasta, fai, dict, dbsnp ->
                def final_dbsnp = dbsnp != null ? dbsnp : []
                [ ref_meta, fasta, fai, dict, final_dbsnp ]
            }
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                addNewMetaFromAttributes( it, "id" , "reference" , false )
            } // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, dbsnp ]

        ch_input_for_angsd = ch_bams_for_multimap
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict, dbsnp ->
                    bam:   [ meta, bam ]
                    fasta: [ ref_meta, fasta ]
            }

        ANGSD_GL(
            ch_input_for_angsd.bam,
            ch_input_for_angsd.fasta,
            [[], []], // No errors file
        )
        ch_angsd_genotype_likelihoods = ANGSD_GL.out.genotype_likelihood
        ch_versions                   = ch_versions.mix( ANGSD_GL.out.versions.first() )

        // Add genotyper info to the meta
        ch_angsd_genotype_likelihoods = ch_angsd_genotype_likelihoods
            .map {
                meta, glf ->
                [ meta + [ genotyper: "angsd" ], glf ]
            }
    }

    // Run BCFTOOLS_STATS on output from GATK UG, HC and Freebayes
    if ( !params.skip_bcftools_stats && ( params.genotyping_tool == 'hc' || params.genotyping_tool == 'ug' || params.genotyping_tool == 'freebayes' ) ) {
        ch_bcftools_input= ch_genotypes_vcf
            .map {
                addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, vcf, tbi, ref_meta, fasta, fai, dict, dbsnp ->
                    vcf:   [ meta, vcf, tbi ] // bcftools stats module expects a tbi file with the vcf.
                    fasta: [ ref_meta, fasta ]
            }

        BCFTOOLS_STATS_GENOTYPING(
            ch_bcftools_input.vcf,  // vcf
            [ [], [] ],             // regions
            [ [], [] ],             // targets
            [ [], [] ],             // samples
            [ [], [] ],             // exons
            ch_bcftools_input.fasta // fasta
        )
        ch_bcftools_stats = BCFTOOLS_STATS_GENOTYPING.out.stats
        ch_multiqc_files  = ch_multiqc_files.mix(BCFTOOLS_STATS_GENOTYPING.out.stats)
        ch_versions       = ch_versions.mix( BCFTOOLS_STATS_GENOTYPING.out.versions.first() )
    }

    emit:
    eigenstrat          = ch_pileupcaller_genotypes     // [ [ meta ], geno, snp, ind ]
    vcf                 = ch_genotypes_vcf              // [ [ meta ], vcf ] ]
    vcf_stats           = ch_bcftools_stats             // [ [ meta ], stats ]
    glf                 = ch_angsd_genotype_likelihoods // [ [ meta ], glf ]
    eigenstrat_coverage = ch_eigenstrat_coverage_stats  // [ [ meta ], stats ]
    versions            = ch_versions
    mqc                 = ch_multiqc_files
}
