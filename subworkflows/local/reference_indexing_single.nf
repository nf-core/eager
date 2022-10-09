//
// Index input reference as required
//

include { GUNZIP                          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
// missing: circulargeneraotr?

workflow REFERENCE_INDEXING_SINGLE {

    take:
    fasta                // file: /path/to/name.{fasta,fa,fna,fas,fasta.gz,fa.gz,fna.gz,fas.gz}
    fasta_fai            // file: /path/to/name.{fasta,fa,fna,fas}.fai
    fasta_dict           // file: /path/to/fasta.dict
    fasta_mapperindexdir // file: /path/to/

    main:

    // Detect if fasta is gzipped or not, unzip if necessary, and generate meta ID by sanitizing file
    if ( fasta.extension == 'gz' ) {
        ch_gz_ref = Channel.fromPath(fasta).map{[ [], it ] }
        GUNZIP ( ch_gz_ref )
        ch_ungz_ref = GUNZIP.out.gunzip
                    .map{
                        meta, fasta ->
                            def new_meta = [:]
                            def clean_name = fasta.name.toString() - ('.' + fasta.extension)
                            new_meta['id'] = clean_name.replaceAll("\\.", "_")

                        [ new_meta, fasta ]
                    }
    } else {
        ch_ungz_ref = Channel.fromPath(fasta)
            .map{
                fasta ->
                    def meta = [:]
                    def clean_name = fasta.name.toString() - ('.' + fasta.extension)
                    meta['id'] = clean_name.replaceAll("\\.", "_")

                [ meta, fasta ]
            }
    }

    // Generate FAI if not supplied, and if supplied generate meta ID
    if ( !fasta_fai ) {
        ch_fasta_fai = SAMTOOLS_FAIDX ( ch_ungz_ref ).fai
    } else {
        ch_fasta_fai = Channel.fromPath(fasta_fai)
            .map {
                fasta_fai ->
                    def meta = [:]
                    def clean_name = fasta_fai.name.toString() - ('.' + 'fai') // TODO problem: how to remove when gzipped extension?
                    meta['id'] = clean_name.replaceAll("\\.", "_")

                [ meta, fasta_fai ]
            }
    }

    // TODO Q: What if user only supplies one, how to ensure meta generated properly?
    // Just specify naming requirements must be same as FASTA file?

    // Join s
    ch_reference_for_mapping = ch_ungz_ref.join(ch_fasta_fai).dump(tag: "refout")

    emit:
    reference = ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex ]

}
