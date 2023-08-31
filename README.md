# Configuration

Make sure you have the following dependencies installed:
   - Nextflow (tested on version 22.04.0.5697)
   - Conda (tested on version 23.3.1)

# Running

The script called `run`:

(1) If you have not run indexing, please input the fasta file list:

```sh
#input fasta files only
nextflow run com.nf \
	--fasta_info "fasta_info" \
	--reads "https://github.com/nf-core/test-datasets/raw/eager/testdata/Mammoth/fastq/JK2782_TGGCCGATCAACGA_L007_R1_001.fastq.gz.tengrand.fq.gz" \
	--threads 128 \
	--label "test" \
	-profile conda \
	-with-trace # -resume
```

(2) If you have run the indexing of concatenated fasta files, please input both fasta files list and bowtie2 index files

```sh
nextflow run com.nf \
        --fasta_info "fasta_info" \
        --reads "https://github.com/nf-core/test-datasets/raw/eager/testdata/Mammoth/fastq/JK2782_TGGCCGATCAACGA_L007_R1_001.fastq.gz.tengrand.fq.gz" \
        --threads 128 \
        --label "test" \
        -profile conda \
        -with-trace \
        --fasta_index "./results/index/test.*.bt2* #-resume
```

**Input**

`--fasta_info`: Specify the path to the fasta information file with the first column as the path to the fasta file and the second column as their name (reference ID).

`--reads`: Provide the reads to be mapped.

`--fasta_index`: Optional input. Input all the bowtie2 index files (.bt2l or .bt2) of the concatenated fasta file of all fasta files. With this input, the indexing process will be skipped.

`--threads`: Number of threads to use (for example the number of cores of your cluster). If your cluster has a partition of memory, please base the number of threads/cores on the memory needed during indexing, which is the most memory-intensive process, unless you skip indexing.

`--label`: Label used for output or naming bowtie2 indexes.

`-with-trace`: Reports time and memory usage of processes.

`-profile conda`: This makes nextflow take care of the environment. Otherwise, you can disable this and check out the `container/env.yml` file for which software you need to load in the module.

**Output**

All outputs are stored in a directory `results/`

`${label}_${reads_filename}.mapped_config_from` has 5 columns:
contig ID, reference ID, Contig length, the number of mapped read-segments, the number of unmapped read-segments.

`${label}_${reads_filename}_reads_sum.csv` output 3 columns: 
- `species`: reference ID given in fasta_info file.

- `n_reads_Sum`: the number of reads mapped to this fasta file

- `contigs_len_Sum`: the total genome size summed from all contigs of the fasta file

**Subdirectory of `results/`**

- `contigs`: contigs ID corresponding to the reference ID given in fasta_info file.

- `mapping`: bam and their index, also extracted bam file with only reads mapped to each fasta file, which has a `-ext.bam` as suffix.

- `index`: bowtie2 index files
