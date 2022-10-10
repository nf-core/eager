//
// Index input reference as required
//

include { GUNZIP                          } from '../../modules/nf-core/gunzip/main'
include { BWA_INDEX                       } from '../../modules/nf-core/bwa/index/main'
include { BOWTIE2_BUILD                   } from '../../modules/nf-core/bowtie2/build/main'
include { SAMTOOLS_FAIDX                  } from '../../modules/nf-core/samtools/faidx/main'
include { PICARD_CREATESEQUENCEDICTIONARY } from '../../modules/nf-core/picard/createsequencedictionary/main'
// TODO missing: circulargeneraotr?

workflow REFERENCE_INDEXING_SINGLE {

    take:
    fasta                // file: /path/to/name.{fasta,fa,fna,fas,fasta.gz,fa.gz,fna.gz,fas.gz}
    fasta_fai            // file: /path/to/name.{fasta,fa,fna,fas}.fai
    fasta_dict           // file: /path/to/fasta.dict
    fasta_mapperindexdir // file: /path/to/

    main:

    def fasta_ext = grab_ungzipped_extension(fasta)
    def clean_name = fasta.name.toString() - fasta_ext

    // TODO for all sections check if saving based on --save_reference

    // Detect if fasta is gzipped or not, unzip if necessary, and generate meta ID by sanitizing file
    if ( fasta.extension == 'gz' ) {
        ch_gz_ref = Channel.fromPath(fasta).map{[[], it]}
        GUNZIP ( ch_gz_ref )
        ch_ungz_ref = GUNZIP.out.gunzip.map{[[id: clean_name], it[1] ]}
    } else {
        ch_ungz_ref = Channel.fromPath(fasta).map{[[id: clean_name], it ]}
    }

    // Generate FAI if not supplied, and if supplied generate meta ID
    if ( !fasta_fai ) {
        ch_fasta_fai = SAMTOOLS_FAIDX ( ch_ungz_ref ).fai.map{[ [id: clean_name - '.fai'], it[1] ] }
    } else {
        ch_fasta_fai = Channel.fromPath(fasta_fai).map{[[id: clean_name], it ]}
    }

    // Generate DICT if not supplied, and if supplied generate meta
    if ( !fasta_dict ) {
        ch_fasta_dict = PICARD_CREATESEQUENCEDICTIONARY ( ch_ungz_ref ).reference_dict.map{[ [id: clean_name - '.fai'], it[1] ] }
    } else {
        ch_fasta_dict = Channel.fromPath(ch_fasta_dict).map{[[id: clean_name], it ]}
    }

    // Generate mapper indicies if not supplied, and if supplied generate meta
    if ( params.mapping_tool == 'bwaaln' ){

        if ( !fasta_mapperindexdir ) {
            ch_fasta_mapperindexdir = BWA_INDEX ( ch_ungz_ref ).index
        } else {
            ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
        }

    } else if ( params.mapping_tool == "bowtie2" ) {

        if ( !fasta_mapperindexdir ) {
            ch_fasta_mapperindexdir = BOWTIE2_BUILD ( ch_ungz_ref ).index
        } else {
            ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}
        }

    }


    // TODO: document that the 'base name' of all indicies must be the same, i.e. correspond to the FASTA

    // Join all together into a single map. Include failOnMismatch as a check if
    // a user supplies indicies with different 'base' names.
    ch_reference_for_mapping = ch_ungz_ref.join(ch_fasta_fai, failOnMismatch: true).join(ch_fasta_dict, failOnMismatch: true).join(ch_fasta_mapperindexdir, failOnMismatch: true)

    emit:
    reference = ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex ]

}

def grab_ungzipped_extension (Path infile) {

    def split_name = infile.toString().tokenize('.')
    def output = split_name.reverse().first() == 'gz' ? split_name.reverse()[1,0].join('.') : split_name.reverse()[0]

    return '.' + output

}
