#!/usr/bin/env python3

# Written by Maxime Borry and released under the MIT license.
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import pysam
from xopen import xopen
import logging
import os
from pathlib import Path


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
        "-merged",
        dest="merged",
        default=False,
        action="store_true",
        help="specify if bam file was created from merged fastq files",
    )
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
        help="Read removal mode: remove reads (remove) or replace sequence by N (replace). Default = remove",
    )
    parser.add_argument(
        "-t", dest="threads", default=4, help="Number of parallel threads"
    )

    args = parser.parse_args()

    bam = args.bam_file
    in_fwd = args.fwd
    merged = args.merged
    in_rev = args.rev
    out_fwd = args.out_fwd
    out_rev = args.out_rev
    mode = args.mode
    threads = int(args.threads)

    return (bam, in_fwd, merged, in_rev, out_fwd, out_rev, mode, threads)


def extract_mapped(bamfile, merged):
    """Get mapped reads in parallel
    Args:
        threads(int): number of threads to use
        bam(str): path to bamfile
    Returns:
        bamfile(str): path to bam alignment file
        result(set): list of mapped reads name (str)
    """
    if bamfile.endswith(".bam") or bamfile.endswith(".gz"):
        read_mode = "rb"
    else:
        read_mode = "r"
    mapped_reads = set()
    bamfile = pysam.AlignmentFile(bamfile, mode=read_mode)
    for read in bamfile.fetch():
        if read.flag != 4:
            if merged:
                if read.query_name.startswith("M_"):
                    mapped_reads.add(read.query_name[2:])
                elif read.query_name.startswith("MT_"):
                    mapped_reads.add(read.query_name[3:])
                else:
                    mapped_reads.add(read.query_name)
            else:
                mapped_reads.add(read.query_name)
    return mapped_reads


def read_write_fq(fq_in, fq_out, mapped_reads, mode, write_mode, proc):
    """
    Read and write fastq file with mapped reads removed
    Args:
        fq_in(str): path to input fastq file
        fq_out(str): path to output fastq file
        mapped_reads(set): set of mapped reads name (str)
        mode(str): read removal mode (remove or replace)
        write_mode(str): write mode (w or wb)
        proc(int): number of parallel processes
        merged(bool): True if bam file was created from merged fastq files
    """
    if write_mode == "w":
        cm = open(fq_out, write_mode)
    elif write_mode == "wb":
        cm = xopen(fq_out, mode=write_mode, threads=proc)
    with pysam.FastxFile(fq_in) as fh:
        with cm as fh_out:
            for read in fh:
                try:
                    if read.name in mapped_reads:
                        if mode == "replace":
                            read.sequence = "N" * len(read.sequence)
                            read = str(read) + "\n"
                            if write_mode == "w":
                                fh_out.write(read)
                            elif write_mode == "wb":
                                fh_out.write(read.encode())
                    else:
                        read = str(read) + "\n"
                        if write_mode == "w":
                            fh_out.write(read)
                        elif write_mode == "wb":
                            fh_out.write(read.encode())
                except Exception as e:
                    logging.error(f"Problem with {str(read)}")
                    logging.error(e)

def check_remove_mode(mode):
    if mode.lower() not in ["replace", "remove"]:
        logging.info(f"Mode must be {' or '.join(mode)}")
    return mode.lower()


if __name__ == "__main__":
    BAM, IN_FWD, MERGED, IN_REV, OUT_FWD, OUT_REV, MODE, PROC = _get_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")

    if OUT_FWD == None:
        out_fwd = os.path.join(os.getcwd(), Path(IN_FWD).stem + ".r1.fq.gz")
    else:
        out_fwd = OUT_FWD

    if out_fwd.endswith(".gz"):
        write_mode = "wb"
    else:
        write_mode = "w"

    remove_mode = check_remove_mode(MODE)

    # FORWARD OR SE FILE
    logging.info(f"- Extracting mapped reads from {BAM}")
    mapped_reads = extract_mapped(BAM, merged=MERGED)
    logging.info(f"- Checking forward fq file {IN_FWD}")
    read_write_fq(
        fq_in=IN_FWD,
        fq_out=out_fwd,
        mapped_reads=mapped_reads,
        mode=remove_mode,
        write_mode=write_mode,
        proc=PROC,
    )
    logging.info(f"- Cleaned forward FastQ file written to {out_fwd}")

    # REVERSE FILE
    if IN_REV:
        if OUT_REV == None:
            out_rev = os.path.join(os.getcwd(), Path(IN_REV).stem + ".r2.fq.gz")
        else:
            out_rev = OUT_REV
        logging.info(f"- Checking reverse fq file {IN_FWD}")
        read_write_fq(
            fq_in=IN_REV,
            fq_out=out_rev,
            mapped_reads=mapped_reads,
            mode=remove_mode,
            write_mode=write_mode,
            proc=PROC,
        )
        logging.info(f"- Cleaned reverse FastQ file written to {out_rev}")
