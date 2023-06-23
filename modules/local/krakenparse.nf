process KRAKENPARSE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
    input:
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*read_kraken_parsed.csv"), emit: read_kraken_parsed
    tuple val(meta), path("*kmer_kraken_parsed.csv"), emit: kmer_kraken_parsed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def read_out = "${meta.id}.read_kraken_parsed.csv"
    def kmer_out = "${meta.id}.kmer_kraken_parsed.csv"
    """
    kraken_parse.py \\
        -c ${params.metagenomics_min_support_reads} \\
        -or $read_out \\
        -ok $kmer_out \\
        $report

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
