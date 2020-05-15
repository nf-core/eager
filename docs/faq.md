# Frequently Asked Questions

## Warning about sticked on revision

If you get a warning like the following:

```bash
Project nf-core/eager currently is sticked on revision: dev -- you need to specify explicitly a revision with the option -r to use it
```

This is a nextflow error. This means that you have multiple versions of nf-core/eager downloaded (e.g. 2.0.0, 2.1.0, 2.1.1, dev etc.). With every `nextflow run nf-core/eager` command you must always specify which one you want to use with `-r` - there is no 'default' version it will use.

For example you must run :

```bash
nextflow run nf-core/eager -r 2.1.0 --reads '/<path>/<to>/data/*_{R1,R2}_*.fq.gz' <...>
```

Specifying the version of the run you are using is highly recommended, as it helps in full reproducibility. Whereby if you record the whole command _with_ the version for your publication or internal reports, then anyone who wants to check your work can use the exact version you also used (including all internal tools).

You can see more information on the nextflow documentation [here](https://www.nextflow.io/docs/latest/sharing.html?highlight=revision#handling-revisions).

## I specified a module and it didn't produce the expected output

Possible options:

1. Check there if you have a typo in the parameter name. Nextflow _does not_
check for this
2. Check that an upstream module was turned on (if a module requires the output
of a previous module, it will not be activated unless it receives the output)

## The pipeline crashes almost immediately with an early pipeline step

Sometimes a newly downloaded and set up nf-core/eager pipeline will encounter
an issue where a run almost immediately crashes
(e.g. at `fastqc`, `outputdocumentation` etc.) saying the tool could not be
found or similar.

If you're running singularity, it could be that nextflow cannot access your
singularity image properly - often due to missing bind paths.

See [here](https://nf-co.re/usage/troubleshooting#cannot-find-input-files-when-using-singularity)
for more information.

## The pipeline has crashed with an error but nextflow is still running

If this happens, you can either wait until all other already running jobs to
safely finish, or if nextflow _still_ does not stop press `ctrl + c` on your
keyboard (or equivalent) to stop the nextflow run.

> :warning: if you do this, and do not plan to fix the run make sure to delete
the output folder. Otherwise you may end up a lot of large intermediate files
being left! You can clean a nextflow run of all intermediate files with
`nextflow clean -f -k` or delete the `work/` directory.

## I get a file name collision error during library merging

When using TSV input, nf-core EAGER will attempt to merge all files with the
same `Sample ID`. However, if you have specified the same `Library_ID` for two
sets of FASTQ files you will likely receive an error such as

```bash
Error executing process > 'library_merge (JK2782)'
Caused by:
  Process `library_merge` input file name collision -- There are multiple input files for each of the following file names: JK2782.mapped_rmdup.bam.csi, JK2782.mapped_rmdup.bam
Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`
Execution cancelled -- Finishing pending tasks before exit
```

This could happen when you specify the same `Library_ID` but with different
sequencing configurations (PE vs SE).

In this case you must modify your `Library_ID` accordingly to make them unique.

## How to investigate a failed run

As with most pipelines, nf-core/eager can sometimes fail either through a
problem with the pipeline itself, but also sometimes through an issue of the
program being run at the given step.

To help try and identify what has caused the error, you can perform the
following steps before reporting the issue:

### 1a Nextflow reports an 'error executing process' with command error

Firstly, take a moment to read the terminal output that is printed by an
nf-core/eager command.

When reading the following, you can see that the actual _command_ failed. When
you get this error, this would suggest that an actual program used by the
pipeline has failed. This is identifiable that you get an `exit status` and
a `Command error:`, the latter of which is what is reported by the failed
program itself.

```bash
ERROR ~ Error executing process > 'circulargenerator (hg19_complete_500.fasta)'

Caused by:
  Process `circulargenerator (hg19_complete_500.fasta)` terminated with an error exit status (1)

Command executed:

  circulargenerator -e 500 -i hg19_complete.fasta -s MT
  bwa index hg19_complete_500.fasta

Command exit status:
  1

Command output:
  (empty)

Command error:
  Exception in thread "main" java.lang.OutOfMemoryError: Java heap space
        at java.util.Arrays.copyOf(Arrays.java:3332)
        at java.lang.AbstractStringBuilder.ensureCapacityInternal(AbstractStringBuilder.java:124)
        at java.lang.AbstractStringBuilder.append(AbstractStringBuilder.java:448)
        at java.lang.StringBuffer.append(StringBuffer.java:270)
  rk dirat CircularGenerator.extendFastA(CircularGenerator.java:155)
        at CircularGenerator.main(CircularGenerator.java:119)

Work dir:
  /projects1/microbiome_calculus/RIII/03-preprocessing/mtCap_preprocessing/work/7f/52f33fdd50ed2593d3d62e7c74e408

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`

 -- Check '.nextflow.log' file for details
```

If you find it is a common error try and fix it yourself by changing your
options in your nf-core/eager run - it could be a
configuration error on your part - however in some cases it could be an error
in the way we've set up the process in nf-core/eager.

To further investigate, go to step 2.

### 1b Nextflow reports an 'error executing process' with no command error

Alternatively, you may get an error with nextflow itself. The most common one
would be a 'process fails' and it looks like the following.

```bash
Error executing process > 'library_merge (JK2782)'
Caused by:
  Process `library_merge` input file name collision -- There are multiple input files for each of the following file names: JK2782.mapped_rmdup.bam.csi, JK2782.mapped_rmdup.bam
Tip: you can try to figure out what's wrong by changing to the process work dir and showing the script file named `.command.sh`
Execution cancelled -- Finishing pending tasks before exit
```

However in this case, there is no `exit status` or `Command error:` message. In
this case this is a nextflow issue.

The example above is because a user has specified multiple sequencing runs of
different libraries but with the same library name. In this case nextflow could
not identify which is the correct file to merge because they have the same name.

This again can also be a user or nextflow error, but the errors are often more
abstract and less clear how to solve (unless you are familiar with nextflow).

Try to investigate a bit further and see if you can understand what the
error refers to, but if you cannot - please ask on the #eager channel on the
[nf-core slack](https://nf-co.re/join/slack) or leave a
[github issue](https://github.com/nf-core/eager/issues).

### 2 Investigating an failed process's `work/` directory

If you haven't found a clear solution to the failed process from the reported
error's, you can next go into the directory where the process was working in,
and investigate the log and error messages that are produced by each command
of the process.

For example, in the error in
[1a](#1a-Nextflow-reports-an-error-executing-process-with-command-error) you can
see the following line

```bash
Work dir:
  /projects1/microbiome_calculus/RIII/03-preprocessing/mtCap_preprocessing/work/7f/52f33fdd50ed2593d3d62e7c74e408
```

If you change into this with `cd` and run `ls -la` you should see a collection
normal files, symbolic links (symlinks) and hidden files (indicated with `.`
at the beginning of the file name).

- Symbolic links: are typically input files from previous processes.
- Normal files: are typically successfully completed output files from some of
some of the commands in the process
- Hidden files are nextflow generated files and include the submission commands
as well as log files

When you have an error run, you can firstly check the contents of the output
files to see if they are empty or not (e.g. with `cat` or `zcat`),
interpretation of which will depend on the program thus dependent on the user
knowledge.

Next, you can investigate `.command.err` and `.command.out`, or `.command.log`.
These represent the standard out or err (in the case of `.log`, both combined)
of all the commands/programs in the process - i.e. what would be printed
to screen if you were running the command/program yourself. Again, view
these with e.g. `cat` and see if you can identify the error of the program
itself.

Finally, you can also try running the commands _yourself_. You can firstly
try to do this by loading your given nf-core/eager environment (e.g.
`singularity shell /<path>/<to>/nf-core-eager-X-X-X.img` or
`conda activate nf-core-eager-X.X.X`), then running `bash .command.sh`.

If this doesn't work this suggests either there is something wrong with the
nf-core/eager environment confirugration, _or_ there still a problem with the
program itself. To confirm the former, try running the command within the
`.command.sh` file (viewable with `cat`) but with locally installed versions
of programs you may already have on your system. If the command still doesn't
work - it is a problem with the program or your specified configuration. If
it does work locally, please report as a [github issue](https://github.com/nf-core/eager/issues).

If it does, please ask the developer of the tool (although we will endevour to
help as much as we can via the [nf-core slack](https://nf-co.re/join/slack) in
the #eager channel).
