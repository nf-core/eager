process HOST_REMOVAL {
    tag "$meta.id"
    label 'process_medium'

    //TODO:check if correct conda and containers used
    conda "bioconda::xopen=1.1.0 bioconda::pysam=0.16.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-adbddd275b72bafc166f5c77e6b1a8c3d200cfe7:2c601eb41a331051f5d90b147919a69ac4e17e19-0' :
        'quay.io/biocontainers/mulled-v2-adbddd275b72bafc166f5c77e6b1a8c3d200cfe7:2c601eb41a331051f5d90b147919a69ac4e17e19-0' }"

    input:
    // [ meta_bam, bam, meta_fastqs, fastqs ]
    tuple val(meta), path(bam), path(bai), val(meta_fastqs), path(fastqs)

    output:
    tuple val(meta), path("*.fq.gz"), emit: fastqs
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}_${meta.reference}"
    def VERSION_PYSAM = '0.16.0'
    def VERSION_XOPEN = '1.1.0'
    def fastqs_input = [fastqs].flatten().sort().size() > 1 ? "${fastqs[0]} -rev ${fastqs[1]}" : "${fastqs[0]}"
    def outputR2 = [fastqs].flatten().sort().size() > 1 ? "-or ${prefix}_hostremoved.r2.fq.gz" : ''

    """
    extract_map_reads.py \\
        $args \\
        $outputR2 \\
        -of ${prefix}_hostremoved.r1.fq.gz \\
        -t $task.cpus \\
        $bam \\
        $fastqs_input


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pysam: $VERSION_PYSAM
        xopen: $VERSION_XOPEN
    END_VERSIONS
    """

}
