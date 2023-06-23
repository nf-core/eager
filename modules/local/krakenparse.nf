process KRAKENPARSE {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'quay.io/biocontainers/python:3.8.3' }"
    input:
    // TODO nf-core: Where applicable all sample-specific information e.g. "id", "single_end", "read_group"
    //               MUST be provided as an input via a Groovy Map called "meta".
    //               This information may not be required in some instances e.g. indexing reference genome files:
    //               https://github.com/nf-core/modules/blob/master/modules/nf-core/bwa/index/main.nf
    // TODO nf-core: Where applicable please provide/convert compressed files as input/output
    //               e.g. "*.fastq.gz" and NOT "*.fastq", "*.bam" and NOT "*.sam" etc.
    tuple val(meta), path(report)

    output:
    tuple val(meta), path("*read_kraken_parsed.csv"), emit: read_kraken_parsed
    tuple val(meta), path("*kmer_kraken_parsed.csv"), emit: kmer_kraken_parsed
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO default value for KRAKEN_PARSE min reads shared with MALTEXTRACT, but recommended defaults in tools is 50 vs 1, respectively: add check and warning?
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
