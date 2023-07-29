#!/usr/bin/env python3

#summing up the number of mapped reads, grouped by species

import pandas as pd

def sum_mapped_reads_by_species(data):

    # Use pandas.melt() to reshape the DataFrame and aggregate by sum
    melted_data = pd.melt(data, id_vars='species', value_vars='mapped_reads', value_name='n_reads_Sum')

    #print(melted_data.head(10))

    aggregated_data = melted_data.groupby('species')['n_reads_Sum'].sum().reset_index()

    return aggregated_data

def sum_contig_length_by_species(data):
    
    # Use pandas.melt() to reshape the DataFrame and aggregate by sum
    melted_data = pd.melt(data, id_vars='species', value_vars='contig_length', value_name='contigs_len_Sum')

    aggregated_data = melted_data.groupby('species')['contigs_len_Sum'].sum().reset_index()

    return aggregated_data



# Specify the file path of the data file
#file_path = 'Artemisia_contig_info.txt'
import sys
genus=sys.argv[1]
file_path = genus+'.mapped_config_from'
print(file_path)

# Read the data file into a pandas DataFrame
data = pd.read_csv(file_path, sep='\t', header=None, on_bad_lines = 'warn')

# Assign column names to the DataFrame
data.columns = ['contigID', 'species', 'contigID1', 'contig_length', 'mapped_reads', 'unmapped_reads']

# get the sum of the number of mapped reads
sum_mapped_reads = sum_mapped_reads_by_species(data)
print(sum_mapped_reads)

# get the sum of the contig lengths of a species
sum_contig_lens = sum_contig_length_by_species(data)
print(sum_contig_lens)

#combine the two dataframes
sum_all = pd.merge(sum_mapped_reads, sum_contig_lens, on = 'species')
print(sum_all)

sum_all.to_csv(genus+'_reads_sum.csv', sep='\t', index=False, header=True)
