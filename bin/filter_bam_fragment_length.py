#!/usr/bin/env python3

# Written by Maxime Borry and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import pysam


def get_args():
    """This function parses and return arguments passed in"""
    parser = argparse.ArgumentParser(
        prog="bam_filter", description="Filter bam on fragment length"
    )
    parser.add_argument("bam", help="Bam aligment file")
    parser.add_argument(
        "-l",
        dest="fraglen",
        default=35,
        type=int,
        help="Minimum fragment length. Default = 35",
    )
    parser.add_argument(
        "-a",
        dest="all",
        default=False,
        action="store_true",
        help="Include all reads, even unmapped",
    )
    parser.add_argument(
        "-o",
        dest="output",
        default=None,
        help="Output bam basename. Default = {bam_basename}.filtered.bam",
    )

    args = parser.parse_args()

    bam = args.bam
    fraglen = args.fraglen
    allreads = args.all
    outfile = args.output

    return (bam, fraglen, allreads, outfile)


def getBasename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return basename


def filter_bam(infile, outfile, fraglen, allreads):
    """Write bam to file

    Args:
        infile (stream): pysam stream
        outfile (str): Path to output bam
        fraglen(int): Minimum fragment length to keep
        allreads(bool): Apply on all reads, not only mapped
    """
    bamfile = pysam.AlignmentFile(infile, "rb")
    bamwrite = pysam.AlignmentFile(outfile + ".filtered.bam", "wb", template=bamfile)

    for read in bamfile.fetch(until_eof=True):
        if allreads:
            if read.query_length >= fraglen:
                bamwrite.write(read)
        else:
            if read.is_unmapped == False and read.query_length >= fraglen:
                bamwrite.write(read)


if __name__ == "__main__":
    BAM, FRAGLEN, ALLREADS, OUTFILE = get_args()

    BAMFILE = pysam.AlignmentFile(BAM, "rb")

    if OUTFILE is None:
        OUTFILE = getBasename(BAM)

    filter_bam(BAM, OUTFILE, FRAGLEN, ALLREADS)

