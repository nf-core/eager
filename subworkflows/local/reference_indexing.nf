//
// Prepare reference indexing for downstream
//

include { REFERENCE_INDEXING_SINGLE } from '../../subworkflows/local/reference_indexing_single.nf'
include { REFERENCE_INDEXING_MULTI  } from '../../subworkflows/local/reference_indexing_multi.nf'

workflow REFERENCE_INDEXING {
    take:
    fasta // file: /path/to/samplesheet.csv
    fasta_fai
    fasta_dict
    fasta_mapperindexdir

    main:
    ch_versions = Channel.empty()
    // TODO add WARN: if fasta.ext == csv && fai/dict/mapperindexdir supplied, then latter will be ignored with preference for info in csv!


    if ( fasta_fai && fasta_dict && fasta_mapperindexdir ) {
        def fasta_ext = grabUngzippedExtension(fasta)
        def clean_name = fasta.name.toString() - fasta_ext

        ch_fasta = Channel.fromPath(fasta).map{[[id: clean_name], it ]}
        ch_fasta_fai = Channel.fromPath(fasta_fai).map{[[id: clean_name], it ]}
        ch_fasta_dict = Channel.fromPath(fasta_dict).map{[[id: clean_name], it ]}
        ch_fasta_mapperindexdir = Channel.fromPath(fasta_mapperindexdir).map{[[id: clean_name], it ]}

        ch_reference_for_mapping = ch_fasta
                                    .join(ch_fasta_fai, failOnMismatch: true)
                                    .join(ch_fasta_dict, failOnMismatch: true)
                                    .join(ch_fasta_mapperindexdir, failOnMismatch: true)
                                    .map{
                                        meta, fasta, fai, dict, mapper_index ->
                                        [ meta, fasta, fai, dict, mapper_index, params.fasta_circular_target, params.fasta_mitochondrion_header ]
                                    }

    } else if ( fasta.extension == 'csv' | fasta.extension == 'tsv' ) {
        REFERENCE_INDEXING_MULTI ( fasta )
        ch_reference_for_mapping = REFERENCE_INDEXING_MULTI.out.reference
        ch_versions = ch_versions.mix( REFERENCE_INDEXING_MULTI.out.versions )
    } else {
        ch_reference_for_mapping = REFERENCE_INDEXING_SINGLE ( fasta, fasta_fai, fasta_dict, fasta_mapperindexdir ).reference
        ch_versions = ch_versions.mix( REFERENCE_INDEXING_SINGLE.out.versions )
    }

    emit:
    reference = ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex ]
    versions  = ch_versions

}

// TODO Move to lib?
def grabUngzippedExtension (Path infile) {

    def split_name = infile.toString().tokenize('.')
    def output = split_name.reverse().first() == 'gz' ? split_name.reverse()[1,0].join('.') : split_name.reverse()[0]

    return '.' + output

}
