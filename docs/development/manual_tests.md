# Manual Tests

Here is a list of manual tests we can run with the expect output commands

## Preprocessing

All possible params

```
skip_preprocessing                               = false
preprocessing_tool                               = 'fastp'
preprocessing_skippairmerging                    = false
preprocessing_excludeunmerged                    = false
preprocessing_skipadaptertrim                    = false
preprocessing_adapter1                           = null
preprocessing_adapter2                           = null
preprocessing_adapterlist                        = null
preprocessing_minlength                          = 25
preprocessing_savepreprocessedreads              = false

preprocessing_fastp_complexityfilter             = false
preprocessing_trim5p                             = 0      // WARN: slightly different behaviour between fastp & adapterremoval
preprocessing_trim3p                             = 0      // WARN: slightly different behaviour between fastp & adapterremoval
preprocessing_fastp_complexityfilter_threshold   = 30
preprocessing_adapterremoval_preserve5p          = false
preprocessing_adapterremoval_adapteroverlap      = 1
preprocessing_adapterremoval_skipqualitytrimming = false
preprocessing_adapterremoval_trimbasequalitymin  = 20
preprocessing_adapterremoval_skipntrimming       = false
preprocessing_adapterremoval_qualitymax          = 41

skip_deduplication                               = false
deduplication_tool                               = 'markduplicates'
```

General Combinations:

- All numeric values change ✅✅

- with/without preprocessing ✅✅
- with/without skipmerging ✅✅
- with/without merged only ✅✅

- with/without skip adapter trimming ✅✅
- with adapterlist ✅✅
- with custom adapterseq (1&2) ✅✅
- overriding order of skipping+trim/adapterlist/custom adapterseq ✅✅
- with default ✅✅

- with/without preprocessed reads ✅✅

Tool Specific combinations

- fastp
  - with/without complexity filtering ✅
- AdapterRemoval

  - with/without skipqualitytim ✅
  - with/without skipntrimming ✅

- Markduplicates

  - With FastP
    - SE&PE data ✅
    - SE&PE data + preprocessing_excludeunmerged ✅
    - PE_only + preprocessing_excludeunmerged ✅
  - With AdapterRemoval
    - With FastP
    - SE&PE data ✅
    - SE&PE data + preprocessing_excludeunmerged ✅
    - PE_only + preprocessing_excludeunmerged ✅

- Dedup
  - With FastP
    - SE&PE data ✅ (expected failure)
    - SE&PE data + preprocessing_excludeunmerged ✅ (expected failure)
    - PE_only + preprocessing_excludeunmerged ✅
  - With AdapterRemoval
    - SE&PE data ✅ (expected failure)
    - SE&PE data + preprocessing_excludeunmerged ✅ (expected failure)
    - PE_only + preprocessing_excludeunmerged ✅

### AdapterRemoval

```bash
## IMPORTANT: CHECK COMMANDS IN BOTH SE/PE EXAMPLES!

## numeric tests: check all numeric parameters are present and match those in the command:
##                 should also show `--collapse`, and single `_merged` file in dumped channel
##                 .settings should indicate two adapter sequences on lines 7/8
##                check no other results files (i.e., FASTQs) than settings
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_minlength 30 --preprocessing_trim5p 5 --preprocessing_trim3p 3 --preprocessing_adapterremoval_adapteroverlap 3 --preprocessing_adapterremoval_trimbasequalitymin 25 --preprocessing_adapterremoval_qualitymax 47

## Skip preprocessing: no adapterremoval executed
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --skip_preprocessing

## Skip merging: command for paired does not have `--collapse`, and  a.dump() on preprocessing out shows two FASTQs
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_skippairmerging

## Merged only:  CAT_FASTQ command should only have collapsed and collapsed truncated, no pair1/pair2 truncated
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_excludeunmerged

## Skip adapter trimming: .command.sh should have adapter1/2 specifies with empty string in quotes adapter sequences in .settings log file should be empty,
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_skipadaptertrim

## Adapter List: adapterlist file should be present in the command, settings should have multiple sets of 1/2 sadapter sequences
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_adapterlist 'data/adapterlist_ar.txt'

## Custom adapter trimming sequences: .command.sh should display the input sequences as does the `.settings` file
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check skipping adapter trimming overrides adapterlist/sequences
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_skipadaptertrim --preprocessing_adapterlist 'data/adapterlist_ar.txt' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check adapterlist overrides custom sequences
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_adapterlist 'data/adapterlist_ar.txt' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check default trimming works: command only contains file and defaults numeric and flag parameters as specified in nextflow.config. Settings file should contain tool-default adapter sequences
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval'

## Check saving preprocessed reads: check output directory contains preprocessed FASTQ reads for both files
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_savepreprocessedreads

## Check quality trimming and trimm n is is turned off: flags should not be in .command.sh and settings should saythey are set as 'No'
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'adapterremoval' --preprocessing_adapterremoval_skipqualitytrimming --preprocessing_adapterremoval_skipntrimming
```

### fastp

```bash
## numeric tests: check all numeric parameters are present and match those in the command:
##                 should also show `-m`, and single `_merged` file in dumped channel for PqE data
##                 should see warning that --include_unmerged is specified so not dumping to separate file
##                 .log should indicate two default Illumrina read1/read2 adapter sequences
##                 should see 'Merged and filtered' section
##                check no other results files (i.e., FASTQs) than settings
##                  Should see command right at the end of log has correct numberic values
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_minlength 30 --preprocessing_trim5p 5 --preprocessing_trim3p 3

## Skip preprocessing: no adapterremoval executed
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --skip_preprocessing

## Skip merging: command for paired does not have `-m`, and  a.dump() on preprocessing out shows two FASTQs, single end unmodified
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_skippairmerging

## Merged only:  no `--include_unmerged` in .command.sh or log file (and can check merged.fastq.gz is small than default run)
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_excludeunmerged

## Skip adapter trimming: .command.sh should have --disable_adapter_trimming, and no references to any adapters at the top of the log
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_skipadaptertrim

## Adapter List: adapterlist file should be present in the command and in command in log... (no other specific stats on this :( ))
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_adapterlist 'data/adapterlist_fastp.fa'

## Custom adapter trimming sequences: .command.sh should display the input sequences and no 'detected' reads in log; HTML should show AAAA adapters of all lengths or ends of 'detected' adapters
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check saving preprocessed reads: check output directory contains preprocessed FASTQ reads for both files
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_minlength 30 --preprocessing_trim5p 5 --preprocessing_trim3p 3 --preprocessing_savepreprocessedreads
7d/dd2a00

## Check skippting adapter trimming overrides adapterlist/sequences, i.e. no adapter trimming references in log and --disabled in command!
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_skipadaptertrim --preprocessing_adapterlist 'data/adapterlist_fastp.fa' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check adapterlist overrides custom sequences: no custom sequences in log/command, only adapter list (auto-detection allowed)
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_adapterlist 'data/adapterlist_fastp.fa' --preprocessing_adapter1 AAAA --preprocessing_adapter2 AAAAAAA

## Check default trimming works: command only contains --detect_adpater_for_pe (no other adapter related stuff). Log file should contain tool-default auto-detection (PE only)
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp'

## Check complexity filter: in command and in log (should say: reads failed due to low complexity)
nextflow run ../main.nf -profile test,singularity --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --outdir ./results -ansi-log false -dump-channels --preprocessing_tool 'fastp' --preprocessing_fastp_complexityfilter --preprocessing_fastp_complexityfilter_threshold 90
```

## Mapping

### General

```bash
## Large Ref: if --fasta_largeref specified, `csi` file generated and published alongisde BAM file; .command.sh contains -c`
nextflow run ../main.nf -profile test,singularity --outdir ./results -resume -dump-channels -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --fasta_largeref
```

### BWA ALN

```bash
## BWA ALN (with merging): SAMSE run twice, BWA ALN .command.sh read gorup header is not `null` & contains -n -k -l -o values, `.bai` file present in `mapping/` directory
## Samtools flagstat figures &  all versions preent in MultiQC report
nextflow run ../main.nf -profile test,singularity --outdir ./results -resume -dump-channels -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta

## BWA ALN (with merging): Both SAMSE and SAMPE executed
nextflow run ../main.nf -profile test,singularity --outdir ./results -resume -dump-channels -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --preprocessing_skippairmerging

```

## Deduplication

In cases where deduplication is expected, check that both `mapping` and `dedup` bams exist, and that the latter contains one read less overlapping position `NC_007596.2:187-187`, where a PCR duplicate exists.

### MARKDUPLICATES

#### With FastP (implicit)

```bash
## MARKDUPLICATES with default parameters
## Expect: deduplication directory with bam,bai,flagstat for each library (6 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --outdir ./results/markduplicates -dump-channels -ansi-log false --deduplication_tool 'markduplicates' -resume

## MARKDUPLICATES run only on merged reads
## Expect: deduplication directory with bam,bai,flagstat for each library (6 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --outdir ./results/markduplicates_merged -dump-channels -ansi-log false --deduplication_tool 'markduplicates' --preprocessing_excludeunmerged -resume

#### MARKDUPLICATES run only merged only PE
## Use as input a version of the TSV in the test profile that only includes the PE data row.
## Expect: deduplication directory with a single bam,bai,flagstat for the library (3 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --input ~/eager_dsl2_testing/input/only_PE/pe_only.tsv --outdir ./results/markduplicates_merged_PE_only  -dump-channels -ansi-log false --deduplication_tool 'markduplicates' --preprocessing_excludeunmerged -resume
```

#### With AdapterRemoval

```bash
#### MARKDUPLICATES with AR on all reads
## Expect: deduplication directory with bam,bai,flagstat for each library (6 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --outdir ./results/AR_markduplicates -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool 'markduplicates' -resume

#### MARKDUPLICATES with AR on merged reads
## Expect: deduplication directory with bam,bai,flagstat for each library (6 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --outdir ./results/AR_markduplicates_merged -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool 'markduplicates' --preprocessing_excludeunmerged -resume

#### MARKDUPLICATES with AR on merged reads from PE run only.
## Use as input a version of the TSV in the test profile that only includes the PE data row.
## Expect: deduplication directory with a single bam,bai,flagstat for the library (3 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --input ~/eager_dsl2_testing/input/only_PE/pe_only.tsv --outdir ./results/AR_markduplicates_merged_PE_only -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool 'markduplicates' --preprocessing_excludeunmerged -resume

```

### DEDUP

#### With FastP (implicit)

```bash
#### DEDUP all reads from FastP
## Expect: Run fails with error "Dedup can only be used on collapsed (i.e. merged) PE reads. For all other cases, please set --deduplication_tool to 'markduplicates'."
nextflow run main.nf -profile docker,test --outdir ./results/dedup -dump-channels -ansi-log false --deduplication_tool 'dedup' -resume

#### DEDUP only on merged from FastP
## Expect: Run fails with error "[nf-core] Error: Invalid input/parameter combination: '--deduplication_tool' cannot be 'dedup' on runs that include SE data. Use  'markduplicates' for all or separate SE and PE data into separate runs."
nextflow run main.nf -profile docker,test --outdir ./results/dedup_merged -dump-channels -ansi-log false --deduplication_tool 'dedup' --preprocessing_excludeunmerged -resume

#### DEDUP only on merged from PE runs from  FastP
## Use as input a version of the TSV in the test profile that only includes the PE data row.
## Expect: deduplication directory with a single bam,bai,flagstat for the library (3 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --input ~/eager_dsl2_testing/input/only_PE/pe_only.tsv --outdir ./results/dedup_merged_PE_only -dump-channels -ansi-log false --deduplication_tool 'dedup' --preprocessing_excludeunmerged -resume

```

#### With AdapterRemoval

```bash
#### DEDUP on all reads from Adapter Removal
## Expect: Run fails with error "Dedup can only be used on collapsed (i.e. merged) PE reads. For all other cases, please set --deduplication_tool to 'markduplicates'."
nextflow run main.nf -profile docker,test --outdir ./results/AR_dedup -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool dedup -resume

#### DEDUP only on merged from AdapterRemoval
## Expect: Run fails with error "[nf-core] Error: Invalid input/parameter combination: '--deduplication_tool' cannot be 'dedup' on runs that include SE data. Use  'markduplicates' for all or separate SE and PE data into separate runs."
nextflow run main.nf -profile docker,test --outdir ./results/AR_dedup_merged -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool 'dedup' --preprocessing_excludeunmerged -resume

#### DEDUP only on merged from PE runs from  AdapterRemoval
## Use as input a version of the TSV in the test profile that only includes the PE data row.
## Expect: deduplication directory with a single bam,bai,flagstat for the library (3 files total). Flagstat for each library should include fewer mapped reads than the mapped bam version. Check that duplicate at NC_007596.2:187-187 is removed.
nextflow run main.nf -profile docker,test --input ~/eager_dsl2_testing/input/only_PE/pe_only.tsv --outdir ./results/AR_dedup_merged_PE_only -dump-channels -ansi-log false --preprocessing_tool 'adapterremoval' --deduplication_tool 'dedup' --preprocessing_excludeunmerged -resume
```
