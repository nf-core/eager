# nf-core/eager: Contributing Guidelines

Hi there!
Many thanks for taking an interest in improving nf-core/eager.

We try to manage the required tasks for nf-core/eager using GitHub issues, you probably came to this page when creating one.
Please use the pre-filled template to save time.

However, don't be put off by this template - other more general issues and suggestions are welcome!
Contributions to the code are even more welcome ;)

> If you need help using or modifying nf-core/eager then the best place to ask is on the nf-core Slack [#eager](https://nfcore.slack.com/channels/eager) channel ([join our Slack here](https://nf-co.re/join/slack)).

## Contribution workflow

If you'd like to write some code for nf-core/eager, the standard workflow is as follows:

1. Check that there isn't already an issue about your idea in the [nf-core/eager issues](https://github.com/nf-core/eager/issues) to avoid duplicating work
    * If there isn't one already, please create one so that others know you're working on this
2. [Fork](https://help.github.com/en/github/getting-started-with-github/fork-a-repo) the [nf-core/eager repository](https://github.com/nf-core/eager) to your GitHub account
3. Make the necessary changes / additions within your forked repository (following [code contribution guidelines](https://github.com/nf-core/eager/blob/dev/.github/CONTRIBUTING.md))
4. Use `nf-core schema build .` and add any new parameters to the pipeline JSON schema (requires nf-core tools >= 1.10).
5. Submit a Pull Request against the `dev` branch and wait for the code to be reviewed and merged

If you're not used to this workflow with git, you can start with some [docs from GitHub](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests) or even their [excellent `git` resources](https://try.github.io/).

## Tests

When you create a pull request with changes, [GitHub Actions](https://github.com/features/actions) will run automatic tests.
Typically, pull-requests are only fully reviewed when these tests are passing, though of course we can help out before then.

There are typically two types of tests that run:

### Lint Tests

`nf-core` has a [set of guidelines](https://nf-co.re/developers/guidelines) which all pipelines must adhere to.
To enforce these and ensure that all pipelines stay in sync, we have developed a helper tool which runs checks on the pipeline code. This is in the [nf-core/tools repository](https://github.com/nf-core/tools) and once installed can be run locally with the `nf-core lint <pipeline-directory>` command.

If any failures or warnings are encountered, please follow the listed URL for more documentation.

### Pipeline Tests

Each `nf-core` pipeline should be set up with a minimal set of test-data.
`GitHub Actions` then runs the pipeline on this data to ensure that it exits successfully.
If there are any failures then the automated tests fail.
These tests are run both with the latest available version of `Nextflow` and also the minimum required version that is stated in the pipeline code.

## Patch

:warning: Only in the unlikely and regretful event of a release happening with a bug.

* On your own fork, make a new branch `patch` based on `upstream/master`.
* Fix the bug, and bump version (X.Y.Z+1).
* A PR should be made on `master` from patch to directly this particular bug.

## Getting help

For further information/help, please consult the [nf-core/eager documentation](https://nf-co.re/eager/docs) and don't hesitate to get in touch on the nf-core Slack [#eager](https://nfcore.slack.com/channels/eager) channel ([join our Slack here](https://nf-co.re/join/slack)).

# Code Contribution Guidelines

To make the EAGER2 code and processing logic more understandable for new contributors, and to ensure quality. We are making an attempt to somewhat-standardise the way the code is written.

If you wish to contribute a new module, please use the following coding standards.

The typical workflow for adding a new module is as follows:

1. Define the corresponding input channel into your new process from the expected previous process channel (or re-routing block, see below).
2. Write the process block (see below).
3. Define the output channel if needed (see below).
4. Add any new flags/options to `nextflow.config` with a default (see below).
5. Add any new flags/options to the help message (for integer/text parameters, print to help the corresponding `nextflow.config` parameter).
6. Add sanity checks for all relevant parameters.
7. Add any new software to the `scrape_software_versions.py` script in `bin/` and the version command to the `scrape_software_versions` process in `main.nf`.
8. Do local tests that the new code works properly and as expected.
9. Add a new test command in `.github/workflow/ci.yaml`.
10. If applicable add a [MultiQC](https://https://multiqc.info/) module.
11. Update MultiQC config `assests/multiqc_config.yaml` so relevant suffixes, name clean up, General Statistics Table column order, and module figures are in the right order.
12. Add new flags/options to 'usage' documentation under `docs/usage.md`.
13. Add any descriptions of MultiQC report sections and output files to `docs/output.md`.

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

## Software Version Reporting

If you add a new tool to the pipeline, please ensure you add the information of the tool to the `get_software_version` process.

Add to the script block of the process, something like the following:

```bash
<YOUR_TOOL> --version &> v_<YOUR_TOOL>.txt 2>&1 || true
```

or

```bash
<YOUR_TOOL> --help | head -n 1 &> v_<YOUR_TOOL>.txt 2>&1 || true
```

You then need to edit the script `bin/scrape_software_versions.py` to

1. add a (python) regex for your tools --version output (as in stored in the `v_<YOUR_TOOL>.txt` file), to ensure the version is reported as a `v` and the version number e.g. `v2.1.1`
2. add a HTML block entry to the `OrderedDict` for formatting in MultiQC.

> If a tool does not unfortunately offer any printing of version data, you may add this 'manually' e.g. with `echo "v1.1" > v_<YOUR_TOOL>.txt`

## Images and Figures

For all internal nf-core/eager documentation images we are using the 'Kalam' font by the Indian Type Foundry and licensed under the Open Font License.

For the overview image we follow the nf-core [style and guidelines](https://nf-co.re/developers/design_guidelines).

These can be found for download [here](https://fonts.google.com/specimen/Kalam)
