#!/usr/bin/env python3

# Written by Maxime Borry and released under the MIT license.
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
from multiprocessing.pool import ThreadPool
import pysam
from xopen import xopen
from functools import partial
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from io import StringIO


def _get_args():
    """This function parses and return arguments passed in"""
    parser = argparse.ArgumentParser(
        prog="extract_mapped_reads",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Remove mapped in bam file from fastq files",
    )
    parser.add_argument("bam_file", help="path to bam file")
    parser.add_argument("fwd", help="path to forward fastq file")
    parser.add_argument(
        "-rev", dest="rev", default=None, help="path to reverse fastq file"
    )
    parser.add_argument(
        "-of", dest="out_fwd", default=None, help="path to forward output fastq file"
    )
    parser.add_argument(
        "-or", dest="out_rev", default=None, help="path to forward output fastq file"
    )
    parser.add_argument(
        "-m",
        dest="mode",
        default="remove",
        help="Read removal mode: remove reads (remove) or replace sequence by N (replace)",
    )
    parser.add_argument(
        "-t", dest="threads", default=4, help="Number of parallel threads"
    )

    args = parser.parse_args()

    bam = args.bam_file
    in_fwd = args.fwd
    in_rev = args.rev
    out_fwd = args.out_fwd
    out_rev = args.out_rev
    mode = args.mode
    threads = int(args.threads)

    return (bam, in_fwd, in_rev, out_fwd, out_rev, mode, threads)


def extract_mapped_chr(chr, BAMFILE):
    """
    Get mapped reads per chromosome
    Args:
        chr(str): chromosome
        BAMFILE(pysam alignement): PySam alignment file object
    Returns:
        res(list): list of mapped reads (str) name per chromosome
    """
    res = []
    reads = BAMFILE.fetch(chr, multiple_iterators=True)
    for read in reads:
        if read.is_unmapped == False:
            if read.query_name.startswith("M_"):
                read_name = read.query_name.replace("M_", "").split()[0].split("/")[0]
            else:
                read_name = read.query_name.split()[0].split("/")[0]
            res.append(read_name)
    return res


def extract_mapped(BAMFILE, threads):
    """Get mapped reads in parallel
    Args:
        BAMFILE(pysam alignement): PySam alignment file object
        threads(int): number of threads to use
    Returns:

    """
    try:
        chrs = BAMFILE.references
    except ValueError as e:
        print(e)

    # Returns empty list if not reads mapped (because not ref match in bam)
    if len(chrs) == 0:
        return []

    # Checking that nb_process is not > nb_chromosomes
    threads = min(threads, len(chrs))
    extract_mapped_chr_partial = partial(extract_mapped_chr, BAMFILE=BAMFILE)
    with ThreadPool(threads) as t:
        res = t.map(extract_mapped_chr_partial, chrs)
    result = [i for ares in res for i in ares if len(i) > 0]
    return result


def parse_write_fq(
    fq_in_path, fq_out_path, mapped_reads, write_mode, remove_mode, threads
):
    """Parse Fastq, compares to mapped reads, and writes cleaned fastq to disk

    Args:
        fq_in_path (str): path to input fastq file
        fq_out_path (str): path to output fastq file
        mapped_reads (list): list of the name of mapped reads
        write_mode (str): writing mode ('w' or 'wb')
        remove_mode (str): remove mode ('replace' or 'remove')
        threads (int): number of threads for parallel write compression by xopen
    """
    if write_mode == "wb":
        with pysam.FastqFile(fq_in_path) as fin, xopen(
            fq_out_path, mode="wb", threads=threads
        ) as fout:
            for entry in fin:
                if entry.name in mapped_reads:
                    if remove_mode == "remove":
                        continue
                    else:
                        entry.set_sequence(
                            sequence="N" * len(entry.sequence),
                            quality=entry.quality,
                        )
                        fout.write(f"{entry}\n".encode())
                else:
                    fout.write(f"{entry}\n".encode())
    else:
        with pysam.FastqFile(fq_in_path) as fin, open(fq_out_path, mode="w") as fout:
            for entry in fin:
                if entry.name in mapped_reads:
                    if remove_mode == "remove":
                        continue
                    else:
                        entry.set_sequence(
                            sequence="N" * len(entry.sequence),
                            quality=entry.quality,
                        )
                        fout.write(f"{entry}\n")
                else:
                    fout.write(f"{entry}\n")


def get_mapped_reads(fq_dict, mapped_reads):
    """Sort mapped reads from dictionary of fastq reads
    Args:
        fq_dict(dict) dictionary with read names as keys, seq and quality as values
        in a list
        mapped_reads(list) list of mapped reads
    Returns:
        fqd(dict) dictionary with read names as key, unmapped/mapped (u|m),
        seq and quality as values in a list
    """

    def intersection(list1, list2):
        return set(list1).intersection(list2)

    def difference(list1, list2):
        return set(list1).difference(list2)

    fqd = {}
    all_reads = list(fq_dict.keys())
    mapped = intersection(all_reads, mapped_reads)
    unmapped = difference(all_reads, mapped_reads)

    for rm in mapped:
        fqd[rm] = ["m"] + fq_dict[rm]
    for ru in unmapped:
        fqd[ru] = ["u"] + fq_dict[ru]

    return fqd


def check_remove_mode(mode):
    if mode.lower() not in ["replace", "remove"]:
        print(f"Mode must be {' or '.join(mode)}")
    return mode.lower()


if __name__ == "__main__":
    BAM, IN_FWD, IN_REV, OUT_FWD, OUT_REV, MODE, THREADS = _get_args()

    if OUT_FWD == None:
        out_fwd = f"{IN_FWD.split('/')[-1].split('.')[0]}.r1.fq.gz"
    else:
        out_fwd = OUT_FWD

    if out_fwd.endswith(".gz"):
        write_mode = "wb"
    else:
        write_mode = "w"

    remove_mode = check_remove_mode(MODE)
    BAMFILE = pysam.AlignmentFile(BAM, mode="r", threads=THREADS)

    # FORWARD OR SE FILE
    print(f"- Extracting mapped reads from {BAM}")
    mapped_reads = extract_mapped(BAMFILE, threads=THREADS)
    print(f"Comparing mapped reads with forward fastq files and writing to disk")
    parse_write_fq(
        IN_FWD,
        OUT_FWD,
        mapped_reads,
        write_mode=write_mode,
        remove_mode=remove_mode,
        threads=THREADS,
    )
    if IN_REV:
        print(f"Comparing mapped reads with forward fastq files and writing to disk")
        parse_write_fq(
            IN_REV,
            OUT_REV,
            mapped_reads,
            write_mode=write_mode,
            remove_mode=remove_mode,
            threads=THREADS,
        )
