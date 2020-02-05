#!/usr/bin/env python

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
        '-o',
        dest="output",
        default=None,
        help="Output file. Default = sources.csv")

    args = parser.parse_args()

    outfile = args.output

    return(outfile)


def get_csv():
    tmp = [i for i in os.listdir() if ".csv" in i]
    return(tmp)


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
    OUTFILE = _get_args()
    all_csv = get_csv()
    resdf = merge_csv(all_csv)
    write_csv(resdf, "kraken_otu_table.csv")
    print(resdf)
