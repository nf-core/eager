# nf-core/eager: Troubleshooting

## My pipeline update doesn't seem to do anything

To download a new version of a pipeline, you can use the following, replacing `<VERSION>` the the corresponding verison.

```bash
nextflow pull nf-core/eager -r <VERSION>
```

However, in very rare cases, minor fixes to a version will be pushed out without a version number bump. This can confuse nextflow slightly, as it thinks you already have the 'broken' version from your original pipeline download.

If when running the pipeline you don't see any changes in the fixed version when running it, you can try removing your nextflow EAGER cache typically stored in your home directory with

```bash
rm -r ~/.nextflow/assets/nf-core/eager
```

And re-pull the pipeline with the command above. This will install a fresh version of the version with the fixes.

## Input files not found

If no file, only one input file, or only read one and not read two is picked up then something is wrong with your input file declaration

1. The path must be enclosed in quotes (`'` or `"`)
2. The path must have at least one `*` wildcard character. This is even if you are only running one paired end sample.
3. When using the pipeline with paired end data, the path must use `{1,2}` or `{R1,R2}` notation to specify read pairs.
4. If you are running Single end data make sure to specify `--single_end`

If the pipeline can't find your files then you will get the following error

```bash
ERROR ~ Cannot find any reads matching: *{1,2}.fastq.gz
```

Note that if your sample name is "messy" then you have to be very particular with your glob specification. A file name like `L1-1-D-2h_S1_L002_R1_001.fastq.gz` can be difficult enough for a human to read. Specifying `*{1,2}*.gz` won't work give you what you want whilst `*{R1,R2}*.gz` will.

## Data organization

The pipeline can't take a list of multiple input files - it takes a glob expression. If your input files are scattered in different paths then we recommend that you generate a directory with symlinked files. If running in paired end mode please make sure that your files are sensibly named so that they can be properly paired. See the previous point.

## Extra resources and getting help

If you still have an issue with running the pipeline then feel free to contact us.
Have a look at the [pipeline website](https://github.com/nf-core/eager) to find out how.

If you have problems that are related to Nextflow and not our pipeline then check out the [Nextflow gitter channel](https://gitter.im/nextflow-io/nextflow) or the [google group](https://groups.google.com/forum/#!forum/nextflow).
