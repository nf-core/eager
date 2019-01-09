# Frequently Asked Questions

## My pipeline update doesn't seem to do anything!

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
