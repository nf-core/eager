# Adjustable parameters for nf-core/eager

This part of the readme contains a list of user-adjustable parameters in nf-core/eager. You can specify any of these parameters on the command line when calling the pipeline by simply prefixing the respective parameter with a double dash `--`.

Example:
```
nextflow run nf-core/eager -r 2.0 -profile standard,docker --singleEnd [...]
```
This would run the pipeline in single end mode, thus assuming that all entered `FastQ` files are sequenced following a single end sequencing protocol.


## General Pipeline Parameters
### `--name` false
### `--singleEnd` false
### `--genome <name>` "Custom"
### `--snpcapture` false
### `--bedfile` ''
### `--fasta` false
### `--bwa_index` false
### `--seq_dict` false
### `--fasta_index` false
### `--saveReference` false
### `--udg` false
### `--udg_type "Half"` (or it assumes fullUDG)
### `--multiqc_config` "$baseDir/conf/multiqc_config.yaml"
### `--email` false
### `--plaintext_email` false

## Step skipping parameters

### `--skip_preseq` false
### `--skip_damage_calculation` false
### `--skip_qualimap` false
### `--skip_deduplication` false
### `--complexity_filter` false

## Complexity Filtering Options

### `--complexity_filter_poly_g_min` 10

## Adapter Clipping and Merging Options

### `--clip_forward_adaptor` AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
### `--clip_reverse_adaptor` AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA
### `--clip_readlength` 30
### `--clip_min_read_quality` 20
### `--clip_min_adap_overlap` 1

## Read Mapping Parameters

### `--bwaalnn` 0.04
### `--bwaalnk` 2
### `--bwaalnl` 32

## Read Filtering and Conversion Parameters

### `--bam_keep_mapped_only` false
### `--bam_keep_all` true
### `--bam_filter_reads` false
### `--bam_mapping_quality_threshold` 0 

## DNA Damage Assessment Parameters

### `--damageprofiler_length` 100
### `--damageprofiler_threshold` 15

## Read DeDuplication Parameters

### `--dedupper` dedup

## Library Complexity Estimation Parameters

### `--preseq_step_size` 1000

## DNA Damage Filtering Methods

### `--run_pmdtools` false
### `--pmdtools_range` 10
### `--pmdtools_threshold ` 3
### `--pmdtools_reference_mask` ''
### `--pmdtools_max_reads` 1000
