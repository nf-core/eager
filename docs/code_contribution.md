# Code Contribution Guidelines

To make the EAGER2 code and processing logic more understandable for new contributers, and to ensure quality. We are making an attempt to somewhat-standardise the way the code is written.

If you wish to contribute a new module, please use the following coding standards.

## Default Values

Default values should go in `nextflow.config` under the `params` scope.

## Default resource processes

Defining recommended 'minimum' resource requiements (CPUs/Memory) for a process should be defined in `conf/base.config`. This can be utilised within the process using `${task.cpu}` or `${task.memory}` variables in the `script:` block.

## Process Concept

We are providing a highly configurable pipeline, with many options to turn on and off different processes in different combinations. This can make a very complex graph structure that can cause a large amount of duplicated channels coming out of every process to account for each possible combination.

The EAGER pipeline can currently be broken down into the following 'stages', where a stage is a collection of  non-terminal mututally exclusive processes, which is the output of which is used for another file reporting module (but not reporting!) .

* Input
* Convert BAM
* PolyG Clipping
* AdapterRemoval
* Mapping (either: `bwa`, `bwamem`, or `circularmapper`)
* BAM Filtering
* Deduplication (either `dedup` or `markduplicates`
* BAM Trimming
* PMDtools
* Genotyping

Every step can potentially be skipped, therefore the output of a previous stage must be able to be passed to the next stage, if the given stage is not run.

To somewhat simplify this logic, we have implemented the following structure.

The concept is as follows:

* Every 'stage' of the pipeline (i.e. collection of mutually exclusive processes) must always have a if else statement following it.
* This if else 'bypass' statement collects and standardises all possible input files into single channel(s) for the next stage.
* Importantly - within the bypass statement, a channel from the previous stage's bypass mixes into these output channels. This additional channel is named `ch_previousstage_for_skipcurrentstage`. This contains the output from the previous stage, i.e. not the modified version from the current stage.
* The bypass statement works as follows:
  * If the current stage is turned on: will mix the previous stage and current stage output and filter for file suffixes unique to the current stage output
  * If the current stage is turned off or skipped: will mix the previous stage and current stage output. However as there there is no files in the output channel from the current stage, no filtering is required and the files in the 'ch_XXX_for_skipXXX' stage will be used.
  
 This ensures the same channel inputs to the next stage is 'homogenous' - i.e. all comes from the same source (the bypass statement)
  
 An example schematic can be given as follows

```nextflow
 // PREVIOUS STAGE OUTPUT
if (params.run_bam_filtering) {
    ch_input_for_skipconvertbam.mix(ch_output_ch_convertbam)
        .filter{ it =~/.*converted.fq/}
        .into { ch_convertbam_for_fastp; ch_convertbam_for_skipfastp }
} else {
    ch_input_for_skipconvertbam
        .into { ch_convertbam_for_fastp; ch_convertbam_for_skipfastp }
}

// SKIPPABLE CURRENT STAGE PROCESS
process fastp {
    publishDir "${params.outdir}/fastp", mode: 'copy'

    when:
    params.run_fastp

    input:
    file fq from ch_convertbam_for_fastp

    output:
    file "*pG.fq" into ch_output_from_fastp

    script:
    """
    echo "I have been fastp'd" > ${fq}  
    mv ${fq} ${fq}.pG.fq
    """
}

// NEXT STAGE INPUT PREPARATION
if (params.run_fastp) {
    ch_convertbam_for_skipfastp.mix(ch_output_from_fastp)
        .filter { it =~/.*pG.fq/ }
        .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval }
} else {
    ch_convertbam_for_skipfastp
        .into { ch_fastp_for_adapterremoval; ch_fastp_for_skipadapterremoval }
}

 ```

## Naming Schemes

Please use the following naming schemes, to make it easy to understand what is going where.

* process output: `ch_output_from_<process>`(this should always go into the bypass statement described above).
* skipped process output: `ch_<previousstage>_for_<skipprocess>`(this goes out of the bypass statement described above)
* process inputs: `ch_<previousstage>_for_<process>` (this goes into a process)

## Nextflow Version Bumping

If you have agreement from reviewers, you may bump the 'default' minimum version of nextflow (e.g. for testing).

For this, you need to update the in the `manifest{}` scope of `nextflow.config`, and also in `.travis.yml` and `.github/workflows/nf-core_eager.yml` 