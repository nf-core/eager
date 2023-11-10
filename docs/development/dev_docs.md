# Documentation and how tos for developing eager

## How to add new input files and options to the reference sheet

To add new input files or option to the reference sheet, you have to complete all the following steps.

### Single-reference input workflow

1. Add your new parameter to nextflow.config.
2. Add parameter description to schema (nf-core schema build).
3. Read in new parameter (params.NEW) as input within the reference_indexing_single local subworkflow.
   1. Add new line to the large .map{} operation starting on line 80 and add check if the file exists.
      `def PARAM_NAME = params.NEW != null ? file(params.NEW, checkIfExists: true ) : ""`
   2. Add PARAM_NAME to the result of the map operation. Double-check the order!
   3. With the ch_ref_index_single.multiMap{} below you add the reference name as a meta. You can also combine your new parameter with others if useful for the workflow step.
      `NEW_SUBCHANNEL: [ meta, PARAM_NAME ]`
   4. Add your ch_ref_index_single.NEW_SUBCHANNEL to the final emit.
      `NEW_EMIT = ch_ref_index_single.NEW_SUBCHANNEL`

### Multi-reference input workflow

1. Add new column named SOFTWARE_FILETYPE and test data to the test reference sheet (https://github.com/nf-core/test-datasets/blob/eager/reference/reference_sheet_multiref.csv).
2. Read in new input within the reference_indexing_multi local subworkflow.
   1. Add new line to the large .map{} operation starting on line 30. Add check if the file exists if appropriate.
      `def PARAM_NAME = row["SOFTWARE_FILETYPE"] != "" ? file(row["SOFTWARE_FILETYPE"], checkIfExists: true) : ""`
   2. Add PARAM_NAME to the result of the map operation. Double-check the order!
   3. With the ch_input_from_referencesheet.multiMap{} below you add the reference name as a meta. You can also combine your new parameter with others if useful for the workflow step.
      `NEW_SUBCHANNEL: [ meta, PARAM_NAME ]`
   4. Add ch_input_from_referencesheet.NEW_SUBCHANNEL to the final emit.
      `NEW_EMIT = ch_input_from_referencesheet.NEW_SUBCHANNEL`

### Combining in the Reference Indexing workflow

1. Add you new parameter channel to the if condition selecting between the direct parameter input or the reference sheet input.
   1. below line 25 for reference sheet input
      `ch_NEW_CHANNEL = REFERENCE_INDEXING_MULTI.out.NEW_EMIT`
   2. below "REFERENCE_INDEXING_SINGLE"
      `ch_NEW_CHANNEL = REFERENCE_INDEXING_SINGLE.out.NEW_EMIT`
   3. Filter out options that have not been provided.
      `ch_NEW_CHANNEL = ch_NEW_CHANNEL.filter{ it[1] != "" }`
   4. Add unzipping of zipped input files with GUNZIP.
   5. Add ch_NEW_CHANNEL to the final emit.
      `NEW_EMIT = ch_NEW_CHANNEL`
   6. Call new inputs within the main eager.nf with `REFERENCE_INDEXING.out.NEW_EMIT`.
