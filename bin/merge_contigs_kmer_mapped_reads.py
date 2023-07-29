#!/usr/bin/env python3

#allow an input named genus from command line
import sys
genus = sys.argv[1]

# Read the contents of the first file into a dictionary
file1_dict = {}
with open('map/'+genus+'_all_sorted.idxstats', 'r') as file1:
    for line in file1:
        columns = line.strip().split('\t')
        if len(columns) >= 1:
            file1_dict[columns[0]] = line.strip()

# Create a new file to store the joined data
output_file = open(genus+'_contig_info.txt', 'w')

# Read the contents of the second file and join the matching lines
with open('contigs/'+genus+'_all.contigs', 'r') as file2:
    for line in file2:
        columns = line.strip().split(' ')
        if len(columns) >= 2:
            key = columns[0]
            if key in file1_dict:
                joined_line = file1_dict[key] + '\t' + line.replace(' ','\t')
                output_file.write(joined_line)

# Close the output file
output_file.close()

