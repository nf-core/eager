process CAT_ADAPTERREMOVAL {
    //tag "${meta.id}" TODO WORK OUT WHY CAN'T FIND META?
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::sed=4.7" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'biocontainers/biocontainers:v1.2.0_cv1' }"

    input:
    // TODO Doesn't like a list going into reads for some reason? Should do tho?
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${prefix}.combined.fq.gz"), emit: reads
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    println(meta)
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    /*

    T if SE                                                            ONLY truncated
    T if PE skip_trim merged_only preserve 5p                          ONLY collpased (combined with below)
    T if PE           merged_only preserve 5p                          ONLY collapsed (combined with above)
    T if PE           merged_only,                                     ONLY collapsd, collapsed_tuncated
    T if PE skip_trim,                                                 ONLY collpased, pair1 trunc, pair2_trun (pair files, generated, so including but expect empty)
    T if PE                                   skip_collapse,           ONLY pair1, pair2 !!SEPARATE (no merging! paired-end output/mapping)
    T if PE                       preserve5p,                          all EXCEPT collpased_truncated
    T if PE merge trim                                                 all!

    */

    if ( meta.single_end  ) {
            // single
            """
            cat ${prefix}.truncated.gz > ${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && !params.clipmerge_mergedonly && !params.clipmerge_skipcollapse && params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, all, merge, trim, preserve5p
            """
            cat ${prefix}.collapsed.gz > ${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && params.clipmerge_mergedonly && !params.clipmerge_skipcollapse && !params.clipmerge_skiptrim && !params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, mergedonly, merge, trim, clip5p
            """
            cat ${prefix}.collapsed.gz ${prefix}.collapsed.truncated.gz > ${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && !params.clipmerge_mergedonly && !params.clipmerge_skipcollapse && params.clipmerge_skiptrim && !params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, all, merge, skiptrim, clip5p
            """
            cat ${prefix}.collapsed.gz ${prefix}.pair1.truncated.gz ${prefix}.pair2.truncated.gz > ${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && !params.clipmerge_mergedonly && params.clipmerge_skipcollapse && !params.clipmerge_skiptrim && !params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, all, nomerge, trim, clip5p
            """
            cat ${prefix}.pair1.truncated.gz > ${prefix}_1.combined.fq.gz
            cat ${prefix}.pair2.truncated.gz > ${prefix}_2.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && !params.clipmerge_mergedonly && !params.clipmerge_skipcollapse && !params.clipmerge_skiptrim && params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, all, merge, trim, preserve5p
            """
            cat ${prefix}.collapsed.gz ${prefix}.singleton.truncated.gz ${prefix}.pair1.truncated.gz ${prefix}.pair2.truncated.gz > output/${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    } else if ( !meta.single_end && !params.clipmerge_mergedonly && !params.clipmerge_skipcollapse && !params.clipmerge_skiptrim && !params.clipmerge_adapterremoval_preserve5p  ) {
            // paired, all, merge, trim, clip5p
            """
            cat ${prefix}.collapsed.gz ${prefix}.collapsed.truncated.gz ${prefix}.singleton.truncated.gz ${prefix}.pair1.truncated.gz ${prefix}.pair2.truncated.gz > output/${prefix}.combined.fq.gz

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                cat: \$(echo \$(cat --version 2>&1) | sed 's/^.*coreutils) //; s/ .*\$//')
            END_VERSIONS
            """
    }
}
