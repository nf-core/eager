## investigation notes for updating code to allow for PE inputs into metagenomics profiling (eg for kraken, malt)

see
https://github.com/nf-core/eager/issues/945

current issue is that the reads that go into mapping are not by default extracted as singletons and non-singletons, so we lose that information
Then downstream the inputs into the krakenuniq module (even if split correctly with meta vars) don't have the correct headers to parse the PE nature of the reads (since they have all been concatenated anyways, and just were ORIGINALLY PE)

So: needs to be fixed up higher (eg in bamfiltering.nf, likely with a new adjustment to the SAMTOOLS_FASTQ_UNMAPPED, SAMTOOLS_FASTQ_MAPPED, and SAMTOOLS_VIEW_BAM_FILTERING modules )

ISSUE FOUND: while the outputting of PE reads is OK in bamfiltering.nf (fastq_mapped & fastq_unmapped) when overlap merging is not done cat_fastq weirdly merges singletons to one PE file and other to the other PE file, so then everything gets fucked up
"""
cat input1/JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_unmapped_other.fastq.gz input3/JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_unmapped_1.fastq.gz > JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_1.merged.fastq.gz
cat input2/JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_unmapped_singleton.fastq.gz input4/JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_unmapped_2.fastq.gz > JK2782_JK2782_TGGCCGATCAACGA_Mammoth_MT_Krause_2.merged.fastq.gz
"""

Decision is needed on what behavior is wanted for unmapped singletons, other. and then likely remove the call to cat_fastq for PE reads
Possibly just split to also have the singletons parsed separately?
