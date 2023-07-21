//
// Genotype the input data using the requested genotyper.
//

include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_PILEUPCALLER } from '../../modules/nf-core/samtools/mpileup/main'
include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE     } from '../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main'
include { SEQUENCETOOLS_PILEUPCALLER                        } from '../../modules/nf-core/sequencetools/pileupcaller/main'
include { GATK_REALIGNERTARGETCREATOR                       } from '../../modules/nf-core/gatk/realignertargetcreator/main'
include { GATK_INDELREALIGNER                               } from '../../modules/nf-core/gatk/indelrealigner/main'
include { GATK_UNIFIEDGENOTYPER                             } from '../../modules/nf-core/gatk/unifiedgenotyper/main'
include { GATK4_HAPLOTYPECALLER                             } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { FREEBAYES                                         } from '../../modules/nf-core/freebayes/main'
// TODO Add ANGSD GTL module. The current module does not pick up the .glf.gz output files.

workflow GENOTYPE {
    take:
    ch_bam_bai              // [ [ meta ], bam , bai ]
    ch_fasta_plus           // [ [ meta ], fasta, fai, dict ] // TODO add dbSNP
    ch_snpcapture_bed       // [ [ meta ], bed ]
    ch_pileupcaller_bedfile // [ [ meta ], bed ]
    ch_pileupcaller_snpfile // [ [ meta ], snp ]

    main:
    ch_versions                        = Channel.empty()
    ch_multiqc_files                   = Channel.empty()
    ch_pileupcaller_genotypes          = Channel.empty()
    ch_gatk_haplotypecaller_genotypes  = Channel.empty()
    ch_gatk_unifiedgenotyper_genotypes = Channel.empty()
    ch_freebayes_genotypes             = Channel.empty()
    ch_angsd_genotypes                 = Channel.empty()


    if ( params.genotyping_tool == 'pileupcaller' ) {
        // SAMTOOLS_MPILEUP_PILEUPCALLER( ch_bam_bai, ch_fasta_plus )

    /*
    // TODO - this is not working yet. Need snpcapture Bed and pileupcaller snp file to add here.
    SEQUENCETOOLS_PILEUPCALLER( ch_bam_bai, ch_fasta_plus, ch_versions, ch_multiqc_files )
    */
    }

    if ( params.genotyping_tool == 'ug' ) {
        // Use correct reference for each input bam/bai pair.
        ch_bams_for_multimap = ch_bam_bai
            .map {
            // Prepend a new meta that contains the meta.reference value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }

        ch_fasta_for_multimap = ch_fasta_plus
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }

        ch_input_for_targetcreator = ch_bams_for_multimap
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, ref_meta, fasta, fai, dict -> // TODO add dbSNP
                    bam: [ meta, bam , bai ]
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
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, intervals, ref_meta, fasta, fai, dict -> // TODO add dbSNP
                    bam: [ meta, bam, bai, intervals ]
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
        ch_versions = ch_versions.mix( GATK_INDELREALIGNER.out.versions.first() ) // TODO is this actually needed, since all GATK modules have the same version?

        // Use realigned bams as input for UG. combine with reference info to get correct ordering.
        ch_bams_for_ug = GATK_INDELREALIGNER.out.bam
            .map {
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, bam, bai, intervals, ref_meta, fasta, fai, dict -> // TODO add dbSNP
                    bam: [ meta, bam, bai, intervals ]
                    fasta: [ ref_meta, fasta ]
                    fai:   [ ref_meta, fai ]
                    dict:  [ ref_meta, dict ]
                    // dbsnp: [ ref_meta, dbsnp ] // TODO add dbsnp
            }

        GATK_UNIFIEDGENOTYPER(
            ch_bams_for_ug.bam,
            ch_bams_for_ug.fasta,
            ch_bams_for_ug.fai,
            ch_bams_for_ug.dict,
            [[], []], // No intervals
            [[], []], // No contamination
            [[], []], // No dbsnp //TODO add dbsnp // ch_bams_for_ug.dbsnp
            [[], []]  // No comp
        )
    }

    if ( params.genotyping_tool == 'hc' ) {
        // TODO
    }

    if ( params.genotyping_tool == 'freebayes' ) {
        // TODO
    }

    if ( params.genotyping_tool == 'angsd' ) {
        // TODO
    }

    emit:
    pileupcaller_genotypes          = ch_pileupcaller_genotypes          // [ [ meta ], geno, snp, ind ]
    gatk_haplotypecaller_genotypes  = ch_gatk_haplotypecaller_genotypes  // [ [ meta ], vcf ] ]
    gatk_unifiedgenotyper_genotypes = ch_gatk_unifiedgenotyper_genotypes // [ [ meta ], vcf ] ]
    freebayes_genotypes             = ch_freebayes_genotypes             // [ [ meta ], vcf ] ]
    angsd_genotypes                 = ch_angsd_genotypes                 // [ [ meta ], glf ] ]
    versions                        = ch_versions
    mqc                             = ch_multiqc_files

}
