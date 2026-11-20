process PRINT_CONTAMINATION_ANGSD {
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    path angsd_output

    output:
    path "nuclear_contamination.txt"     , emit: txt
    path "nuclear_contamination_mqc.json", emit: json
    path "versions.yml"                  , emit: versions

    when:
    params.run_contamination_estimation_angsd

    script:

    """
    print_x_contamination.py ${angsd_output.join(' ')}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
