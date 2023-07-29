#!/bin/bash

bam=$1 #the filename of the bam -- without .suffix
sorted_bam=${bam}_sorted.bam
bam_idxstats=${bam}_sorted.idxstats

samtools sort -o $sorted_bam ${bam}.bam
samtools index $sorted_bam 
samtools idxstats $sorted_bam > $bam_idxstats
