process COLLECT_GENOTYPES {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::python=3.8.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.8.3' :
        'biocontainers/python:3.8.3' }"

    input:
    tuple val(meta), path(geno), path(snp), path(ind)

    output:
    tuple val(meta), path("*.geno"), path("*.snp"), path("*.ind")  , emit: collected
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix   = task.ext.prefix ?: "${meta.id}"
    // If there are multiple genotype datasets, then merge them, else just rename the output for consistency.
    println "geno = ${geno.toList().size()}"
    println "${geno.toList()}"
    if ( geno.toList().size() == 1 ) {
        """
            mv ${geno[0]} ${prefix}.geno
            mv ${snp[0]}  ${prefix}.snp
            mv ${ind[0]}  ${prefix}.ind

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            collect_genotypes.py: \$(collect_genotypes.py -v)
        END_VERSIONS
        """
    } else {
        """
            collect_genotypes.py \
                --genoFn1 ${geno[0]} \
                --snpFn1 ${snp[0]} \
                --indFn1 ${ind[0]} \
                --genoFn2 ${geno[1]} \
                --snpFn2 ${snp[1]} \
                --indFn2 ${ind[1]} \
                --output ${prefix}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            collect_genotypes.py: \$(collect_genotypes.py -v)
        END_VERSIONS
        """
    }
}
