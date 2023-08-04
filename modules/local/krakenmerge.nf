process KRAKENMERGE {
    label 'process_single'

    conda "conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:1.5.2' :
        'quay.io/biocontainers/pandas:1.5.2' }"
    input:

    path kraken_parse_reads
    path kraken_parse_kmers

    output:
    path "kraken_read_count_table.csv" , emit: read_count_table
    path "kraken_kmer_duplication.csv" , emit: kmer_duplication_table
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def read_out = "kraken_read_count_table.csv"
    def kmer_out = "kraken_kmer_duplication.csv"
    """
    merge_kraken_res.py \\
        -or $read_out \\
        -ok $kmer_out \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
