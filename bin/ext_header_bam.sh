#!/bin/bash

bam=$1
contigs=$2

samtools view -H $bam > header

first_line=$(cat header | LC_ALL=C fgrep -n "@SQ" | head -1 | cut -d ':' -f 1)
last_line=$(cat header | LC_ALL=C fgrep -n "@SQ" | tail -1 | cut -d ':' -f 1)

#before SQ
head -n "$((first_line - 1))" header > tmp
#extract SQ matching contig ID
sed -n "$((first_line - 1)), $last_line"'p' header | LC_ALL=C fgrep -w -f $contigs >> tmp
#after SQ
tail -n "+$((last_line + 1))" header >> tmp

cat tmp
