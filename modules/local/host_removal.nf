process HOST_REMOVAL {
    tag "$meta.id"
    label 'process_medium'

    //TODO:check if correct conda and containers used
    conda "bioconda::pysam=0.20.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.20.0--py39h9abd093_0' :
        'quay.io/biocontainers/pysam:0.20.0--py39h9abd093_0' }"

    input:
    tupple val(meta), path(bam), path(forwardfastq), path(reversefastq)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fastqs
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def VERSION_PYSAM = '0.16.0'
    def VERSION_XOPEN = '1.1.0'
    def inputR2 = reversefastq ? "-rev ${reversefastq}" : ''
    def outputR2 = reversefastq ? "-or ${prefix}_hostremoved" : ''


    """
    extract_map_reads.py \\
        $args \\
        $inputR2 \\
        $outputR2 \\
        -of ${prefix}_hostremoved \\
        -t $task.cpus \\
        $bam \\
        $forwardfastq


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: $VERSION_PYSAM
        xopen: $VERSION_XOPEN
    END_VERSIONS
    """

}
