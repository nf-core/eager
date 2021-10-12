#!/usr/bin/env python3

# Written by Maxime Borry and released under the MIT license.
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
from multiprocessing import Pool
import pysam
from xopen import xopen
from functools import partial
import sys


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


def extract_mapped_chr(chr, bam):
    """
    Get mapped reads per chromosome
    Args:
        chr(str): chromosome
        bam(str): path to bamfile
    Returns:
        res(list): list of mapped reads (str) name per chromosome
    """
    res = []
    BAMFILE = pysam.AlignmentFile(bam, "rb")
    reads = BAMFILE.fetch(chr)
    for read in reads:
        if read.is_unmapped is False:
            if read.query_name.startswith("M_"):
                read_name = read.query_name.replace("M_", "").split()[0].split("/")[0]
            else:
                read_name = read.query_name.split()[0].split("/")[0]
            res.append(read_name)
    del BAMFILE
    return res


def extract_mapped(threads, bam):
    """Get mapped reads in parallel
    Args:
        threads(int): number of threads to use
        bam(str): path to bamfile
    Returns:

    """
    BAMFILE = pysam.AlignmentFile(bam, "rb", threads=threads)
    try:
        chrs = BAMFILE.references
    except ValueError as e:
        print(e)
    del BAMFILE

    # Returns empty list if not reads mapped (because not ref match in bam)
    if len(chrs) == 0:
        return []
    # Checking that nb_process is not > nb_chromosomes
    threads = min(threads, len(chrs))
    extract_mapped_chr_partial = partial(extract_mapped_chr, bam=bam)
    with Pool(threads) as p:
        res = p.map(extract_mapped_chr_partial, chrs)
    result = [i for ares in res for i in ares if len(i) > 0]
    return result


def compare_mapped_reads(all_read_dict, mapped_reads):
    """Sort mapped reads from dictionary of fastq reads
    Args:
        all_reads(dict) dict of all reads
        mapped_reads(list) list of mapped reads
    Returns:
        (dict) dictionary with read names as key, list of [unmapped/mapped (u|m),
        ,read entry] as value
    """

    fqd = dict()

    def intersection(list1, list2):
        return list(set(list1).intersection(set(list2)))

    def difference(list1, list2):
        return list(set(list1).difference(set(list2)))

    mapped = intersection(all_read_dict.keys(), mapped_reads)
    unmapped = difference(all_read_dict.keys(), mapped_reads)

    for rm in mapped:
        fqd[rm] = ["m", all_read_dict[rm]]
    for ru in unmapped:
        fqd[ru] = ["u", all_read_dict[ru]]

    return fqd


def parse_fq(fq_in_path):
    """Parse Fastq

    Args:
        fq_in_path (str): path to input fastq file
    Returns:
        (dict) dictionary with read names as key, read entry as value
    """
    all_read_dict = dict()
    with pysam.FastqFile(fq_in_path) as fin:
        for entry in fin:
            all_read_dict[entry.name] = entry
    return all_read_dict


def write_fq(fq_out_path, read_dict, write_mode, remove_mode, threads):
    """Writes cleaned fastq to disk

    Args:
        fq_out_path (str): path to output fastq file
        read_dict (dict): dictionary with read names as key, list of [unmapped/mapped (u|m),
        ,read entry] as value
        write_mode (str): writing mode ('w' or 'wb')
        remove_mode (str): remove mode ('replace' or 'remove')
        threads (int): number of threads for parallel write compression by xopen
    """
    if write_mode == "wb":
        with xopen(fq_out_path, mode="wb", threads=threads) as fout:
            for entry in read_dict:
                if read_dict[entry][0] == "u":
                    fout.write(f"{read_dict[entry][1]}\n".encode())
                else:
                    if remove_mode == "remove":
                        continue
                    else:
                        read_dict[entry][1].set_sequence(
                            sequence="N" * len(read_dict[entry][1].sequence),
                            quality=read_dict[entry][1].quality,
                        )
                        fout.write(f"{read_dict[entry][1]}\n".encode())
    else:
        with open(fq_out_path, mode="w") as fout:
            for entry in read_dict:
                if read_dict[entry][0] == "u":
                    fout.write(f"{read_dict[entry][1]}\n")
                else:
                    if remove_mode == "remove":
                        continue
                    else:
                        read_dict[entry][1].set_sequence(
                            sequence="N" * len(read_dict[entry][1].sequence),
                            quality=read_dict[entry][1].quality,
                        )
                        fout.write(f"{read_dict[entry][1]}\n")


def format_fastq_record(name, comment, sequence, quality):
    sequence = "N" * len(sequence)
    return f"{name} {comment}\n{sequence}\n+\n{quality}\n"


def check_remove_mode(mode):
    if mode.lower() not in ["replace", "remove"]:
        print(f"Mode must be {' or '.join(mode)}")
        sys.exit(1)
    return mode.lower()


if __name__ == "__main__":
    BAM, IN_FWD, IN_REV, OUT_FWD, OUT_REV, MODE, THREADS = _get_args()

    if OUT_FWD is None:
        out_fwd = f"{IN_FWD.split('/')[-1].split('.')[0]}.r1.fq.gz"
    else:
        out_fwd = OUT_FWD

    if out_fwd.endswith(".gz"):
        write_mode = "wb"
    else:
        write_mode = "w"

    remove_mode = check_remove_mode(MODE)

    if IN_REV:
        total_steps = 3
    else:
        total_steps = 2

    # FORWARD OR SE FILE
    print(f"* Step 1/{total_steps}: Extracting mapped reads from {BAM}")
    mapped_reads = extract_mapped(bam=BAM, threads=THREADS)
    print(
        f"* Step 2/{total_steps}: Comparing mapped reads with forward fastq files and writing {OUT_FWD} to disk "
    )
    print("\t - Step 2.1: Reading forward fastq")
    all_fwd_reads = parse_fq(IN_FWD)
    print("\t - Step 2.2: Comparing forward fastq")
    all_fwd_reads_checked = compare_mapped_reads(all_fwd_reads, mapped_reads)
    del all_fwd_reads
    print("\t - Step 2.3: Writing forward fastq")
    write_fq(
        OUT_FWD,
        all_fwd_reads_checked,
        write_mode=write_mode,
        remove_mode=remove_mode,
        threads=THREADS,
    )

    del all_fwd_reads_checked

    if IN_REV:
        print(
            f"* Step 3/{total_steps}: Comparing mapped reads with reverse fastq files and writing {OUT_REV} to disk"
        )
        print("\t - Step 3.1: Reading reverse fastq")
        all_rev_reads = parse_fq(IN_REV)
        print("\t - Step 3.2: Comparing reverse fastq")
        all_rev_reads_checked = compare_mapped_reads(all_rev_reads, mapped_reads)
        del all_rev_reads
        print("\t - Step 3.3: Writing reverse fastq")
        write_fq(
            OUT_REV,
            all_rev_reads_checked,
            write_mode=write_mode,
            remove_mode=remove_mode,
            threads=THREADS,
        )
        del all_rev_reads_checked
