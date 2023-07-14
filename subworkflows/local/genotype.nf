//
// Genotype the input data using the requested genotyper.
//

include { SAMTOOLS_MPILEUP as SAMTOOLS_MPILEUP_PILEUPCALLER } from '../../modules/nf-core/samtools/mpileup/main'
include { EIGENSTRATDATABASETOOLS_EIGENSTRATSNPCOVERAGE     } from '../../modules/nf-core/eigenstratdatabasetools/eigenstratsnpcoverage/main'
include { SEQUENCETOOLS_PILEUPCALLER                        } from '../../modules/nf-core/sequencetools/pileupcaller/main'
include { GATK_INDELREALIGNER                               } from '../../modules/nf-core/gatk/indelrealigner/main'
include { GATK_REALIGNERTARGETCREATOR                       } from '../../modules/nf-core/gatk/realignertargetcreator/main'
include { GATK_UNIFIEDGENOTYPER                             } from '../../modules/nf-core/gatk/unifiedgenotyper/main'
include { GATK4_HAPLOTYPECALLER                             } from '../../modules/nf-core/gatk4/haplotypecaller/main'
include { FREEBAYES                                         } from '../../modules/nf-core/freebayes/main'
// TODO Add ANGSD GTL module. The current module does not pick up the .glf.gz output files.

workflow GENOTYPE {
    take:
    ch_bam_bai              // [ [ meta ], bam , bai ]
    ch_fasta                // [ [ meta ], fasta ]
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
        // SAMTOOLS_MPILEUP_PILEUPCALLER( ch_bam_bai, ch_fasta )

    /*
    // TODO - this is not working yet. Need snpcapture Bed and pileupcaller snp file to add here.
    SEQUENCETOOLS_PILEUPCALLER( ch_bam_bai, ch_fasta, ch_versions, ch_multiqc_files )
    */
    }

    if ( params.genotyping_tool == 'unifiedgenotyper' ) {
        // TODO
    }

    if ( params.genotyping_tool == 'haplotypecaller' ) {
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
