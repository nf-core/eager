#!/bin/bash

#input fasta
fa=$1
#fa="/crex/proj/snic2022-6-144/nobackup/CHENYU/check_dna_damage/ref/Artemisia_vulgaris_4220.fna"

taxon=$(basename "$fa" | rev | cut -d. -f2- | cut -f 2,3 -d'_' | rev)

#extract all contigs

if [[ $fa == *.gz ]]; then
  zcat "$fa" | grep ">" | sed 's/>//g' > "${taxon}.contigs"
else
  cat "$fa" | grep ">" | sed 's/>//g' > "${taxon}.contigs"
fi




#add a second column as the taxon, output ${taxon}.contigs
awk -v taxon="$taxon" '{ $2 = $2 "\t" taxon } 1' "${taxon}.contigs" > temp_file && mv temp_file "${taxon}.contigs"


