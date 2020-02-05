#!/usr/bin/env python


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
        '-o',
        dest="output",
        default=None,
        help="Output file. Default = <basename>.kraken_parsed.csv")

    args = parser.parse_args()

    infile = args.krakenReport
    countlim = int(args.count)
    outfile = args.output

    return(infile, countlim, outfile)


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
        countlim (int): lower count threshold to report hit
    OUTPUT:
        resdict (dict): key=taxid, value=readCount

    '''
    with open(infile, 'r') as f:
        resdict = {}
        csvreader = csv.reader(f, delimiter='\t')
        for line in csvreader:
            reads = int(line[1])
            if reads >= countlim:
                taxid = line[4]
                resdict[taxid] = reads
        return(resdict)


def write_output(resdict, infile, outfile):
    with open(outfile, 'w') as f:
        basename = _get_basename(infile)
        f.write(f"TAXID,{basename}\n")
        for akey in resdict.keys():
            f.write(f"{akey},{resdict[akey]}\n")


if __name__ == '__main__':
    INFILE, COUNTLIM, outfile = _get_args()

    if not outfile:
        outfile = _get_basename(INFILE)+".kraken_parsed.csv"

    tmp_dict = parse_kraken(infile=INFILE, countlim=COUNTLIM)
    write_output(resdict=tmp_dict, infile=INFILE, outfile=outfile)
