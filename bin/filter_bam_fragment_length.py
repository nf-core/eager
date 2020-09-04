#!/usr/bin/env python3

import argparse
import multiprocessing
import pysam
from functools import partial


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
        "-p", dest="process", type=int, default=4, help="Number of parallel processes"
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
    process = args.process
    allreads = args.all
    outfile = args.output

    return (bam, fraglen, process, allreads, outfile)


def getBasename(file_name):
    if ("/") in file_name:
        basename = file_name.split("/")[-1].split(".")[0]
    else:
        basename = file_name.split(".")[0]
    return basename


def extract_mapped_chr(chr, fraglen, allreads):
    """
    Get mapped reads per chromosome
    Args:
        chr(str): chromosome
        fraglen(int): Minimum fragment length to keep
        allreads(bool): Apply on all reads, not only mapped
    Returns:
        res(list): list of filtered reads
    """
    res = []
    reads = BAMFILE.fetch(chr, multiple_iterators=True)
    for read in reads:
        if allreads:
            if read.query_length >= fraglen:
                res.append(read)
        else:
            if read.is_unmapped == False and read.query_length >= fraglen:
                res.append(read)
    return res


def extract_mapped(proc, fraglen, allreads):
    """Get mapped reads in parallel
    Returns:
        bamfile(stream): pysam stream
        fraglen(int): Minimum fragment length to keep
        allreads(bool): Apply on all reads, not only mapped
    """
    try:
        chrs = BAMFILE.references
    except ValueError as e:
        print(e)

    # Returns empty list if not reads mapped (because not ref match in bam)
    if len(chrs) == 0:
        return []

    # Checking that nb_process is not > nb_chromosomes
    proc = min(proc, len(chrs))

    extract_partial = partial(extract_mapped_chr, fraglen=fraglen, allreads=allreads)
    with multiprocessing.Pool(proc) as p:
        res = p.map(extract_partial, chrs)
    result = [i for ares in res for i in ares if len(i) > 0]
    return result


def write_bam(infile, outfile, reads):
    """Write bam to file

    Args:
        infile (stream): pysam stream
        outfile (str): Path to output bam
        reads (list): List of filtered reads
    """
    bamwrite = pysam.AlignmentFile(outfile + ".filtered.bam", "wb", template=infile)

    for read in reads:
        bamwrite.write(read)


if __name__ == "__main__":
    BAM, FRAGLEN, PROCESS, ALLREADS, OUTFILE = get_args()

    BAMFILE = pysam.AlignmentFile(BAM, "rb")

    if OUTFILE is None:
        OUTFILE = getBasename(BAM)

    reads = extract_mapped(PROCESS, FRAGLEN, ALLREADS)
    # won't work, see https://github.com/pysam-developers/pysam/issues/752
    write_bam(BAMFILE, OUTFILE, reads)

