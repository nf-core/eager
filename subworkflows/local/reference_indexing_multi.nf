//
// Parse fasta CSV and index where necessary
//

workflow REFERENCE_INDEXING_MULTI {
    take:
    fasta // file: /path/to/samplesheet.csv
    fasta_fai
    fasta_dict
    fasta_mapperindexdir

    main:

    // TODO add WARN: if fasta.ext == csv && fai/dict/mapperindexdir supplied, then latter will be ignored with preference for info in csv!

    if ( reference.fasta.ext == csv ) {
        ch_reference_for_mapping = REFERENCE_INDEXING_MULTI ( fasta ).reference
    } else {
        ch_reference_for_mapping = REFERENCE_INDEXING_SINGLE ( fasta, fasta_fai, fasta_dict, fasta_mapperindexdir ).reference
    }

    emit:
    reference = ch_reference_for_mapping // [ meta, fasta, fai, dict, mapindex ]

}
