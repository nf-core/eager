process FILTER_BAM_FRAGMENT_LENGTH {
    label 'process_single'

    conda "bioconda::pysam=0.20.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pysam:0.20.0--py39h9abd093_0' :
        'quay.io/biocontainers/pysam:0.20.0--py39h9abd093_0' }"

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*_lengthfiltered.bam"), emit: bam
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    filter_bam_fragment_length.py \\
        -a \\
        $args \\
        -o ${prefix}_lengthfiltered.bam \\
        $bam

    ## TODO GET ALL PACKAGES
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

}
