# Naming and formatting conventions for the pipeline code

## Module names

EVERY non-specific module being imported into a local subworkflow (i.e. anywhere not in the main eager.nf script) HAS to be assigned an alias that clearly identifies it. This is to:

- Avoid any confusion between modules with the same name in different subworkflows.
- Avoid annoying warnings popping up in the console when running the pipeline, saying that `No process matching configuration .DEDUPLICATION:SAMTOOLS_VIEW was found.` when users do not run parts of the pipeline.

### Example

Within the map.nf we used to have the following:

```nextflow
include { SAMTOOLS_MERGE } from '../../modules/nf-core/samtools/merge/main'
```

This was renamed to:

```nextflow
include { SAMTOOLS_MERGE as SAMTOOLS_MERGE_LANES } from '../../modules/nf-core/samtools/merge/main'
```

The alias should ideally make it intuitive to understand which subworkflow the module is being imported into. In this case, the module is being imported into the map.nf subworkflow, and since we merge lanes after mapping, the name `SAMTOOLS_MERGE_LANES` should be clear enough.

## Module configuration

- The unique module names specified above should make it possible to always configure modules without the need for a regex/glob when using `withName`. Exception to this is modules named within nf-core subworkflows, which should be configured with a regex/glob.
- The order of attributes within configuration blocks should always be the following:
  1.  tag (mandatory)
  2.  ext.args\* (optional. Followed by ext.args{2,3,...} in ascending order)
  3.  ext.prefix (optional)
  4.  publishDir (optional)
  5.  any other attributes go to the end.
- NEVER use `meta.id` in module configuration (`tag`,`ext.*`), but instead the full explicit combination of unique attributes expected. `meta.sample_id` is fine to use and is equivalent to `meta.id`, but should be supplemented by `meta.library_id` and `meta.lane` etc, as required.
- Every process that is reference-specific MUST include `${meta.reference}` in its `tag` and `ext.prefix` attributes. This is to avoid confusion when running the pipeline with multiple references.
  - Tags that include reference and sample information should be formatted as `${meta.reference}|${meta.sample_id}_*`. Reference specific attributes go on the left-hand-side of the tag, data-specific attributes on the right-hand-side.

### Example

NOT LIKE THIS:

```nextflow
withName ".*:MAP:SAMTOOLS_MERGE_LANES" {
  tag = { "${meta.id}" }
  ext.prefix = { "${meta.id}_${meta.library_id}" }
}
```

BUT LIKE THIS:

```nextflow
withName SAMTOOLS_MERGE_LANES {
        tag = { "${meta.reference}|${meta.sample_id}_${meta.library_id}" }
        ext.prefix = { "${meta.sample_id}_${meta.library_id}_${meta.reference}" }
}
```
