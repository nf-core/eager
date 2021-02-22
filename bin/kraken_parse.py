#!/usr/bin/env python

# Written by Maxime Borry and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import csv

def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='kraken_parse',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='Parsing kraken')
    parser.add_argument('krakenReport', help="path to kraken report file")
    parser.add_argument(
        '-c',
        dest="count",
        default=50,
        help="Minimum number of hits on clade to report it. Default = 50")
    parser.add_argument(
        '-or',
        dest="readout",
        default=None,
        help="Read count output file. Default = <basename>.read_kraken_parsed.csv")
    parser.add_argument(
        '-ok',
        dest="kmerout",
        default=None,
        help="Kmer Output file. Default = <basename>.kmer_kraken_parsed.csv")

    args = parser.parse_args()

    infile = args.krakenReport
    countlim = int(args.count)
    readout = args.readout
    kmerout = args.kmerout

    return(infile, countlim, readout, kmerout)


def _get_basename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return(basename)


def parse_kraken(infile, countlim):
    '''
    INPUT:
        infile (str): path to kraken report file
        countlim (int): lowest count threshold to report hit
    OUTPUT:
        resdict (dict): key=taxid, value=readCount

    '''
    with open(infile, 'r') as f:
        read_dict = {}
        kmer_dict = {}
        csvreader = csv.reader(f, delimiter='\t')
        for line in csvreader:
            reads = int(line[1])
            if reads >= countlim:
                taxid = line[6]
                kmer = line[3]
                unique_kmer = line[4]
                try:
                    kmer_duplicity = float(kmer)/float(unique_kmer)
                except ZeroDivisionError:
                    kmer_duplicity = 0
                read_dict[taxid] = reads
                kmer_dict[taxid] = kmer_duplicity

        return(read_dict, kmer_dict)


def write_output(resdict, infile, outfile):
    with open(outfile, 'w') as f:
        basename = _get_basename(infile)
        f.write(f"TAXID,{basename}\n")
        for akey in resdict.keys():
            f.write(f"{akey},{resdict[akey]}\n")


if __name__ == '__main__':
    INFILE, COUNTLIM, readout, kmerout = _get_args()

    if not readout:
        read_outfile = _get_basename(INFILE)+".read_kraken_parsed.csv"
    else:
        read_outfile = readout
    if not kmerout:    
        kmer_outfile = _get_basename(INFILE)+".kmer_kraken_parsed.csv"
    else:
        kmer_outfile = kmerout

    read_dict, kmer_dict = parse_kraken(infile=INFILE, countlim=COUNTLIM)
    write_output(resdict=read_dict, infile=INFILE, outfile=read_outfile)
    write_output(resdict=kmer_dict, infile=INFILE, outfile=kmer_outfile)
