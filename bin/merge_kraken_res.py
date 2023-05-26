#!/usr/bin/env python

# Written by Maxime Borry and released under the MIT license.
# Modifications for DSL2 compliance and integration into eager DSL2 release by @ilight1542
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import os
import pandas as pd
import numpy as np

def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='merge_kraken_res',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Merging csv count files in one table')
    parser.add_argument(
        '-or',
        dest="readout",
        default="kraken_read_count_table.csv",
        help="Read count output file. Default = kraken_read_count_table.csv")
    parser.add_argument(
        '-ok',
        dest="kmerout",
        default="kraken_kmer_unicity_table.csv",
        help="Kmer unicity output file. Default = kraken_kmer_unicity_table.csv")
    parser.add_argument(
        '-inr',
        dest="input_read_csvs",
        help="Paths to .read_ outputs from krakenparse")
    parser.add_argument(
        '-ink',
        dest="input_kmer_csvs",
        help="Paths to .kmer_ outputs from krakenparse"
    )

    args = parser.parse_args()

    read_paths_in = args.input_read_csvs
    kmer_paths_in = args.input_kmer_csvs
    readout       = args.readout
    kmerout       = args.kmerout

    return(readout, kmerout, read_paths_in, kmer_paths_in)


def get_csv(read_paths_in, kmer_paths_in):
    kmer = [i for i in kmer_paths_in]
    read = [i for i in read_paths_in]
    return(read, kmer)


def _get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def merge_csv(all_csv):
    df = pd.read_csv(all_csv[0], index_col=0)
    for i in range(1, len(all_csv)):
        df_tmp = pd.read_csv(all_csv[i], index_col=0)
        df = pd.merge(left=df, right=df_tmp, on='TAXID', how='outer')
    df.fillna(0, inplace=True)
    return(df)


def write_csv(pd_dataframe, outfile):
    pd_dataframe.to_csv(outfile)


if __name__ == "__main__":
    READOUT, KMEROUT, READ_PATHS_IN, KMER_PATHS_IN  = _get_args()
    reads, kmers = get_csv(READ_PATHS_IN,KMER_PATHS_IN)
    read_df = merge_csv(reads)
    kmer_df = merge_csv(kmers)
    write_csv(read_df, READOUT)
    write_csv(kmer_df, KMEROUT)
