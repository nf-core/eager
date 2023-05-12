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

General Combinations (test logging):

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

### Multi-reference tests

```bash
## Test: (1) Two references, only FASTAs ✅
## Expect: Expect all of fai file (x2 SAMTOOLS_FAIDX processes), dict file (x2 PICARD_CREATESEQUENCEDICTIONARY), bwa index directory (x2 BWA_INDEX) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test01.csv -ansi-log false -dump-channels --save_reference

## Test: (2) Two reference FASTAs, one also has fai ✅
## Expect: Expect one fai file (x1 SAMTOOLS_FAIDX processes), 2 dict (x2 PICARD_CREATESEQUENCEDICTIONARY), 2 bwa  index directory (x2 BWA_INDEX) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test02.csv -ansi-log false -dump-channels --save_reference

## Test: (3) Two reference FASTAs, one also has dict ✅
## Expect: Expect two fai file (x2 SAMTOOLS_FAIDX processes), 1 dict file (x1 PICARD_CREATESEQUENCEDICTIONARY), 2 bwa index directory (x2 BWA_INDEX) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test03.csv -ansi-log false -dump-channels --save_reference

## Test: (4) Two reference FASTAs, one also has bwa_index ref ✅
## Expect: Expect two fai (x2 SAMTOOLS_FAIDX processes), two dict (x2 PICARD_CREATESEQUENCEDICTIONARY), 1 bwa  index directory (x1 BWA_INDEX) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test04.csv -ansi-log false -dump-channels --save_reference

## Test: (5) Two reference FASTAs, one also has bowtie2_index ref ✅
## Expect: Expect two fai (x2 SAMTOOLS_FAIDX processes), two dict (x2 PICARD_CREATESEQUENCEDICTIONARY), 1 bowtie2  index directory (x1 BOWTIE2_BUILD) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test05.csv -ansi-log false -dump-channels --save_reference --mapping_tool bowtie2

## Test: (6) Mapper index mismatch with `--mapping_tool` ✅
## Expect: Expect FAIL at mapping step for Mammoth two fai (x2 SAMTOOLS_FAIDX processes), two dict (x2 PICARD_CREATESEQUENCEDICTIONARY), 1 bowtie2  index directory (BOWTIE2_BUILD) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test06.csv -ansi-log false -dump-channels --save_reference --mapping_tool bowtie2

## Test: (7) Mammoth has all pre-supplied ✅
## Expect: Expect one fai (x1 SAMTOOLS_FAIDX processes), one dict (x1 PICARD_CREATESEQUENCEDICTIONARY), 1 bowtie2  index directory (BOWTIE2_BUILD) etc. to be generated and present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test07.csv -ansi-log false -dump-channels --save_reference

## Test: (8) No indexing necessary, all already supplied ✅
## Expect: Expect no files to be generated/processes executed, nor results present in per reference results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test08.csv -ansi-log false -dump-channels --save_reference

## Test: (9) All but Human FAI provided ❓️ - fails at step not involving FAI?
## Expect: Expect one fai (x1 SAMTOOLS_FAIDX processes), and nothing else results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test09.csv -ansi-log false -dump-channels --save_reference

## Test: (10) All but Human dict provided ✅
## Expect: Expect one dict (x1 PICARD_CREATESEQUENCEDICTIONARY processes), and nothing else results/reference folder
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test10.csv -ansi-log false -dump-channels --save_reference

## Test: (11) Broken path correctly fails pipeline ✅
## Expect: Expect fail
nextflow run ../main.nf -profile singularity,test --outdir ./results --input samplesheet.tsv --fasta reference_sheet_multiref_test11.csv -ansi-log false -dump-channels --save_reference
```

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

## BAM filtering

All possible parameters

```
    // BAM Filtering
    run_bamfiltering                      = false
    bamfiltering_minreadlength            = 0
    bamfiltering_mappingquality           = 0
    bamfiltering_generateunmappedfastq    = false
    bamfiltering_generatemappedfastq      = false
    bamfiltering_retainunmappedgenomicbam = false // downstream genomics only
    bamfiltering_savefilteredbam          = false // can include unmapped reads if --bamfiltering_retainunmappedgenomicbam specified

    // Metagenomic Screening
    run_metagenomicscreening   = false
    metagenomicscreening_input = 'unmapped' // mapped, all, unmapped -> mapped vs all specified in SAMTOOLS_FASTQ_MAPPED in modules.conf, unmapped hardcoded SAMTOOLS_FASTQ_UMAPPED
```

Tests

```bash
## Check no BAM filtering
## Expect: full completion of pipeline without any bam filtering execution
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta

## Check BAM filtering, mapped reads only in downstream BAM
## Expect to see FILTER_BAM workflow with VIEW and FLAGSTAT, no results directory
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering

## Check BAM filtering, mapped reads only in downstream BAM
## Expect to see FILTER_BAM workflow with VIEW and FLAGSTAT, only quality filtered BAMs in results directory, `samtools view` should show only mapped reads
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams

## Check BAM filtering, mapped reads only in downstream BAM with length filtering
## Expect: to see filtered and length filtered BAM in results dir, `samtools stats` on lengthonly, mapped and unmapped , and grep RL shortest read at 50 and only mapped reads, straight BAM should be both mapped only and RL shortest at 50, multiQC should have flagstats
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50

## Check BAM filtering mapped reads only in downstream BAM, no length filtering, but with quality filtering
## Expect: *filtered.bam, samtools stats has < 50 bp, mapped only and samtools view -q 0 produces reads
## TODO: BEFORE RELEASE UPDATE MULTIQC TO PICK UP BOTH BEFORE AND AFTER FILTERING
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams

## Check BAM filtering mapped reads only in downstream BAM, length filtering and quality filtering
## Expect: *filtered.bam, samtools stats has >= 50 bp, mapped only, and samtools stats MQ0: 0
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37

## Check BAM filtering without length/quality filtering, with retained unmapped reads in genomic BAM (why would one need this!? But w/e)
## Expect: 'filtered' BAM still contains unmapped reads in samtools stats
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_retainunmappedgenomicbam

## Check BAM filtering with length/quality filtering, with retained unmapped reads in genomic BAM
## Expect: filtered and length filtered BAMs, final BAM (not length only) should have min RL of 50, MQO: 0, but retain unmapped reads
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_retainunmappedgenomicbam --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37

## Check BAM filtering with length filtering but no quality filtering, with retained unmapped reads in genomic BAM
## Expect: filtered and length filtered BAMs, final BAM (not length only) should have min RL of 50, MQO: >=1, but and include retain unmapped reads
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_retainunmappedgenomicbam --bamfiltering_minreadlength 50

## Check BAM filtering with unmapped FASTQ generation with length/quality filtering of genomic BAM
## Expect: BAM and unmapped FASTQ files in results, samtools stats should have RL >= 50 and MPQ:0
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --bamfiltering_generateunmappedfastq

## Check BAM filtering with mapped FASTQ generation with length/quality filtering of genomic BAM
## Expect: BAM and FASTQ files in results, FASTQ should have more number of reads than the BAM (as fastq should not have length nor quality filtering)
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --bamfiltering_generatemappedfastq

## Check BAM filtering with mapped FASTQ generation WITHOUT length/quality filtering of genomic BAM
## Expect: BAM and FASTQ files in results, FASTQ should have save number of reads as in the BAM
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams  --bamfiltering_generatemappedfastq

## Check BAM filtering with unmmaped AND mapped FASTQ generation with length/quality filtering of genomic BAM
## Expect: filtered BAM and mapped AND unmapped FASTQ files, and larger FASTQ unmapped files, and BAM files are length/quality filtered
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --bamfiltering_generateunmappedfastq --bamfiltering_generatemappedfastq

## Check BAM filtering (mapped only/length/quality on genomic bam) with metagenomics screening, with unmapped reads to metagenomics
# Expect: filtered BAM (samtools stats | grep SN total/mapped same), and a dump() on the ch_bam_for_metagenomics channel should report unmapped_other. Nr. of reads in dumped FASTQ should match approx unmmaped reads in results/mapping/*.flagstat
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --run_metagenomicscreening -dump-channels

## Check BAM filtering (mapped only/length/quality on genomic bam) with metagenomics screening, with mapped only reads going to metagenomics
# Expect: filtered BAM (samtools stats | grep SN total/mapped same), and a dump() on the ch_bam_for_metagenomics channel should report mapped_other. Nr. of reads in dumped FASTQ should match approx mmaped reads in results/mapping/*.flagstat
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --run_metagenomicscreening --metagenomicscreening_input 'mapped' -dump-channels

## Check BAM filtering (mapped only/length/quality on genomic bam) with metagenomics screening, with all reads going to metagenomics
# Expect: filtered BAM (samtools stats | grep SN total/mapped same), and a dump() on the ch_bam_for_metagenomics channel should report mapped_other. Nr. of reads in dumped FASTQ should match total reads in results/mapping/*.flagstat
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --bamfiltering_minreadlength 50 --bamfiltering_mappingquality 37 --run_metagenomicscreening --metagenomicscreening_input 'all' -dump-channels

## Check BAM filtering NO LENGTH/QAULITY with metagenomics screening, with unmapped reads to metagenomics
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is different  and  RL reads min <50), and a dump() on the ch_bam_for_metagenomics channel should report mapped_other. Nr. of reads in dumped FASTQ should match unmapped reads as calculated from results/mapping/*.flagstat. Note: No filtered flagstat expected!
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'unmapped' -dump-channels

## Check BAM filtering NO LENGTH/QAULITY with metagenomics screening, with unmapped reads to metagenomics and save unmapped FASTQ
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is different  and  RL reads min <50), and a dump() on the ch_bam_for_metagenomics channel should report unmapped_other. Nr. of reads in dumped FASTQ should match unmapped reads as calculated from results/mapping/*.flagstat; and unmapped other fASTQ in bam_filtering directoryt. Note: No filtered flagstat expected!
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'unmapped' -dump-channels

## Check BAM filtering NO LENGTH/QAULITY with metagenomics screening, with mapped only reads going to metagenomics
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is different  and  RL reads min <50), and a dump() on the ch_bam_for_metagenomics channel should report mapped_other. Nr. of reads in dumped FASTQ should be roughly matching mappd reads as calculated from results/mapping/*.flagstatt. Note: No filtered flagstat expected!
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'mapped' -dump-channels

## Check BAM filtering NO LENGTH/QAULITY with metagenomics screening, with all reads going to metagenomics
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is different  and  RL reads min <50), and a dump() on the ch_bam_for_metagenomics channel should report mapped_other. Nr. of reads in dumped FASTQ should be roughly matching total reads as calculated from results/mapping/*.flagstatt. Note: No filtered flagstat expected!
## Some reads lost, not 100% why command looks OK... but not just unmapped as more than that
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'all' -dump-channels

## Check BAM filtering ONLY length filtering, with metagenomics screening, with unmapped reads to metagenomics and save unmapped FASTQ
## Metagenomics with length only
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is different  and  RL reads min >= 50), and a dump() on the ch_bam_for_metagenomics channel should report unmapped_other. Nr. of reads in dumped FASTQ should match unmapped reads as calculated from results/mapping/*.flagstat; and unmapped other fASTQ in bam_filtering directoryt.
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'unmapped' -dump-channels --bamfiltering_minreadlength 50

## Check BAM filtering ONLY length filtering, with metagenomics screening, with unmapped reads to metagenomics and save unmapped FASTQ
## Metagenomics with length only
# Expect: filtered BAM (samtools stats SN quality average < 36.7 or view -q 0 vs. -q 37 is not different  and  RL reads min <= 50), and a dump() on the ch_bam_for_metagenomics channel should report unmapped_other. Nr. of reads in dumped FASTQ should match unmapped reads as calculated from results/mapping/*.flagstat; and unmapped other fASTQ in bam_filtering directoryt.
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'unmapped' -dump-channels --bamfiltering_mappingquality 37

## Check what happens when we do paired-end merging and sending reads to metagenomics...
nextflow run ../main.nf -profile test,singularity --outdir ./results -ansi-log false --input data/samplesheet.tsv --fasta data/reference/Mammoth_MT_Krause.fasta --run_bamfiltering --bamfiltering_savefilteredbams --run_metagenomicscreening --metagenomicscreening_input 'unmapped' -dump-channels --bamfiltering_mappingquality 37 --preprocessing_skippairmerging
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

### CALCULATE DAMAGE

#### With DamageProfiler

```bash
## DamageProfiler with default parameters
## Expect:damageprofiler directory with txt, pdf, svg for each library (19 files total per library).
nextflow run main.nf -profile test,conda --outdir ./results -resume
```
