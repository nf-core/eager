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

    // Warn user if they've given a reference sheet that already includes fai/dict/mapper index etc.
    if ( ( fasta.extension == 'csv' || fasta.extension == 'tsv' ) && (fasta_fai || fasta_dict || fasta_mapperindexdir)) log.warn("A TSV or CSV has been supplied to `--fasta` as well as e.g. `--fasta_fai`. --fasta CSV/TSV takes priority and --fasta_* parameters will be ignored.")

    if ( fasta.extension == 'csv' || fasta.extension == 'tsv' ) {
        // If input (multi-)reference sheet supplied
        REFERENCE_INDEXING_MULTI ( fasta )
        ch_reference_for_mapping = REFERENCE_INDEXING_MULTI.out.reference
        ch_versions = ch_versions.mix( REFERENCE_INDEXING_MULTI.out.versions )
    } else {
        // If input FASTA and/or indicies supplied
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
