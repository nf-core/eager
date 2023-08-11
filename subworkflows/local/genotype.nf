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
include { BCFTOOLS_STATS as BCFTOOLS_STATS_GENOTYPING       } from '../../modules/nf-core/bcftools/stats/main'
// TODO Add ANGSD GTL module. The current module does not pick up the .glf.gz output files.

workflow GENOTYPE {
    take:
    ch_bam_bai              // [ [ meta ], bam , bai ]
    ch_fasta_plus           // [ [ meta ], fasta, fai, dict ] // TODO add dbSNP
    ch_snpcapture_bed       // [ [ meta ], bed ]
    ch_pileupcaller_bedfile // [ [ meta ], bed ]
    ch_pileupcaller_snpfile // [ [ meta ], snp ]
    ch_dbsnp                // [ dbsnp ]

    main:
    ch_versions                        = Channel.empty()
    ch_multiqc_files                   = Channel.empty()
    ch_pileupcaller_genotypes          = Channel.empty()
    ch_gatk_haplotypecaller_genotypes  = Channel.empty()
    ch_gatk_unifiedgenotyper_genotypes = Channel.empty()
    ch_freebayes_genotypes             = Channel.empty()
    ch_angsd_genotypes                 = Channel.empty()
    ch_bcftools_stats                  = Channel.empty()

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

        // RESULT: [ [combination_meta], [ref_meta], fasta, fai, dict, dbsnp ]
        ch_fasta_for_multimap = ch_fasta_plus
            .map {
            // Prepend a new meta that contains the meta.id value as the new_meta.reference attribute
                WorkflowEager.addNewMetaFromAttributes( it, "id" , "reference" , false )
            }
            .combine( ch_dbsnp ) // TODO This will need tweaking ('by:0'?) once the dbsnp channel gets a meta.
            .dump(tag:"dbsnp")

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
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
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
        ch_versions = ch_versions.mix( GATK_INDELREALIGNER.out.versions.first() ) // TODO is this actually needed, since all GATK modules have the same version?

        // Use realigned bams as input for UG. combine with reference info to get correct ordering.
        ch_bams_for_ug = GATK_INDELREALIGNER.out.bam
            .map {
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
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
        ch_versions = ch_versions.mix( GATK_UNIFIEDGENOTYPER.out.versions.first() )
        ch_gatk_unifiedgenotyper_genotypes = GATK_UNIFIEDGENOTYPER.out.vcf

    if ( ! params.skip_bcftools_stats ) {
        // TODO this section could be moved outside the UG specific section into its own if clause and take input from HC and FB as well.
        ch_bcftools_input= ch_gatk_unifiedgenotyper_genotypes
            .map {
                WorkflowEager.addNewMetaFromAttributes( it, "reference" , "reference" , false )
            }
            .combine( ch_fasta_for_multimap , by:0 )
            .multiMap {
                ignore_me, meta, vcf, ref_meta, fasta, fai, dict, dbsnp ->
                    vcf:   [ meta, vcf, [] ] // bcftools stats module expects a tbi file with the vcf.
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
        ch_versions = ch_versions.mix( BCFTOOLS_STATS_GENOTYPING.out.versions.first() )
    }
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
    geno_pileupcaller = ch_pileupcaller_genotypes          // [ [ meta ], geno, snp, ind ]
    geno_gatk_hc      = ch_gatk_haplotypecaller_genotypes  // [ [ meta ], vcf ] ]
    geno_gatk_ug      = ch_gatk_unifiedgenotyper_genotypes // [ [ meta ], vcf ] ]
    geno_freebayes    = ch_freebayes_genotypes             // [ [ meta ], vcf ] ]
    geno_angsd        = ch_angsd_genotypes                 // [ [ meta ], glf ] ]
    versions          = ch_versions
    mqc               = ch_multiqc_files

}
