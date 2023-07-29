#!/usr/bin/env python3

import sys
idxstats_file = sys.argv[1]
contigs_file = sys.argv[2]
output_file = sys.argv[3]

#idxstats_file = '/crex/proj/snic2022-6-144/nobackup/CHENYU/check_dna_damage/map/Artemisia_all_sorted.idxstats'
#contigs_file = '/crex/proj/snic2022-6-144/nobackup/CHENYU/check_dna_damage/ref/Artemisia_vulgaris.contigs'
#output_file = 'joined_file.txt'

# Create a dictionary to store the mapping between the first column and the entire line of the idxstats file
mapping = {}
with open(idxstats_file, 'r') as idx_file:
    for line in idx_file:
        columns = line.strip().split('\t')
        if len(columns) >= 1:
            mapping[columns[0]] = line.strip()

# Open the contigs file and write the joined lines to the output file
with open(contigs_file, 'r') as contigs, open(output_file, 'w') as output:
    for line in contigs:
        columns = line.strip().split('\t')
        if len(columns) >= 1:
            contig_id = columns[0].replace(' ', '')
            if contig_id in mapping:
                output.write(f"{line.strip()}\t{mapping[contig_id]}\n")
            else:
                output.write(f"{line.strip()}\tNA\tNA\tNA\n")

print("Join operation completed successfully.")
