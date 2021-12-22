#!/usr/bin/env python3

# Written by Maxime Borry and released under the MIT license.
# See git repository (https://github.com/nf-core/eager) for full license text.

import argparse
import multiprocessing
import pysam
from xopen import xopen
from functools import partial
import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from io import StringIO


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='extract_mapped_reads',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Remove mapped in bam file from fastq files")
    parser.add_argument('bam_file', help="path to bam file")
    parser.add_argument('fwd', help='path to forward fastq file')
    parser.add_argument(
        '-rev',
        dest="rev",
        default=None,
        help="path to reverse fastq file")
    parser.add_argument(
        '-of',
        dest="out_fwd",
        default=None,
        help="path to forward output fastq file")
    parser.add_argument(
        '-or',
        dest="out_rev",
        default=None,
        help="path to forward output fastq file")
    parser.add_argument(
        '-m',
        dest='mode',
        default='remove',
        help='Read removal mode: remove reads (remove) or replace sequence by N (replace)'
    )
    parser.add_argument(
        '-p',
        dest='process',
        default=4,
        help='Number of parallel processes'
    )

    args = parser.parse_args()

    bam = args.bam_file
    in_fwd = args.fwd
    in_rev = args.rev
    out_fwd = args.out_fwd
    out_rev = args.out_rev
    mode = args.mode
    proc = int(args.process)

    return(bam, in_fwd, in_rev, out_fwd, out_rev, mode, proc)


def extract_mapped_chr(chr):
    """
    Get mapped reads per chromosome
    Args:
        chr(str): chromosome
        bam(str): bamfile path
    Returns:
        res(list): list of mapped reads (str) name per chromosome
    """
    res = []
    reads = BAMFILE.fetch(chr, multiple_iterators=True)
    for read in reads:
        if read.is_unmapped == False:
            if read.query_name.startswith("M_"):
                read_name = read.query_name.replace(
                    "M_", "").split()[0].split("/")[0]
            else:
                read_name = read.query_name.split()[0].split("/")[0]
            res.append(read_name)
    return(res)


def extract_mapped(proc):
    """Get mapped reads in parallel
    Returns:
        bamfile(str): path to bam alignment file
        result(list): list of mapped reads name (str)
    """
    try:
        chrs = BAMFILE.references
    except ValueError as e:
        print(e)

    # Returns empty list if not reads mapped (because not ref match in bam)
    if len(chrs) == 0:
        return([])

    # Checking that nb_process is not > nb_chromosomes
    proc = min(proc, len(chrs))
    with multiprocessing.Pool(proc) as p:
        res = p.map(extract_mapped_chr, chrs)
    result = [i for ares in res for i in ares if len(i) > 0]
    return(result)


def parse_fq(fq):
    """Parse a FASTQ file
    Args:
        fq(str): path to fastq file
    Returns:
        fqd(dict): dictionary with read names as keys, seq and quality as values
        in a list
    """
    def get_fq_reads(allreads):
        read_dict = {}
        for title, seq, qual in FastqGeneralIterator(allreads):
            # NEED TO ONLY KEEP THE FIRST PART OF THE FASTQ READ NAME FOR CROSS
            # REFERENCING WITH BAM FILE: ONLY THIS INFORMATION IS KEPT WHEN
            # COLLAPSING READS WITH ADAPTERREMOVAL

            # Until fastq format 1.8
            # Split after slash
            # @HWUSI-EAS100R:6:73:941:1973#0/1
            suf_title = ""
            title_space = title.split()
            if len(title_space) > 1:
                title = title_space[0]
                suf_title = f" {title_space[1]}"

            # From fastq format 1.8
            # Split after space
            # @EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG
            title_slash = title.split("/")
            if len(title_slash) > 1:
                title = title_slash[0]
                suf_title = f"/{title_slash[1]}"

            read_dict[title] = [suf_title, seq, "+", qual]
        return(read_dict)

    if fq.endswith('.gz'):
        with xopen(fq, 'r') as allreads:
            fqd = get_fq_reads(allreads)
    else:
        with open(fq, 'r') as allreads:
            fqd = get_fq_reads(allreads)

    return(fqd)


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
        return(set(list1).intersection(list2))

    def difference(list1, list2):
        return(set(list1).difference(list2))

    fqd = {}
    all_reads = list(fq_dict.keys())
    mapped = intersection(all_reads, mapped_reads)
    unmapped = difference(all_reads, mapped_reads)

    for rm in mapped:
        fqd[rm] = ['m']+fq_dict[rm]
    for ru in unmapped:
        fqd[ru] = ['u']+fq_dict[ru]

    return(fqd)


def write_fq(fq_dict, fname, write_mode, remove_mode, proc):
    """Write to fastq file
    Args:
        fq_dict(dict): fq_dict with unmapped read names as keys,
            unmapped/mapped (u|m), seq, and quality as values in a list
        fname(string) Path to output fastq file
        write_mode (str): 'rb' or 'r'
        remove_mode (str): remove (remove read) or replace (replace read sequence) by Ns
        proc(int) number of processes
    """
    fq_dict_keys = list(fq_dict.keys())
    if write_mode == 'wb':
        with xopen(fname, mode='wb', threads=proc) as fw:
            for fq_dict_key in fq_dict_keys:
                wstring = ""
                if remove_mode == 'remove':
                    if fq_dict[fq_dict_key][0] == 'u':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        for i in fq_dict[fq_dict_key][2:]:
                            wstring += f"{i}\n"
                elif remove_mode == 'replace':
                    # if unmapped, write all the read lines
                    if fq_dict[fq_dict_key][0] == 'u':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        for i in fq_dict[fq_dict_key][2:]:
                            wstring += f"{i}\n"
                    # if mapped, write all the read lines, but replace sequence
                    # by N*(len(sequence))
                    elif fq_dict[fq_dict_key][0] == 'm':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        wstring += f"{'N'*len(fq_dict[fq_dict_key][2])}\n"
                        for i in fq_dict[fq_dict_key][3:]:
                            wstring += f"{i}\n"
                fw.write(wstring.encode())
    else:
        with open(fname, 'w') as fw:
            for fq_dict_key in fq_dict_keys:
                wstring = ""
                if remove_mode == 'remove':
                    if fq_dict[fq_dict_key][0] == 'u':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        for i in fq_dict[fq_dict_key][2:]:
                            wstring += f"{i}\n"
                elif remove_mode == 'replace':
                    # if unmapped, write all the read lines
                    if fq_dict[fq_dict_key][0] == 'u':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        for i in fq_dict[fq_dict_key][2:]:
                            wstring += f"{i}\n"
                    # if mapped, write all the read lines, but replace sequence
                    # by N*(len(sequence))
                    elif fq_dict[fq_dict_key][0] == 'm':
                        wstring += f"@{fq_dict_key+fq_dict[fq_dict_key][1]}\n"
                        wstring += f"{'N'*len(fq_dict[fq_dict_key][2])}\n"
                        for i in fq_dict[fq_dict_key][3:]:
                            wstring += f"{i}\n"
                fw.write(wstring)


def check_remove_mode(mode):
    if mode.lower() not in ['replace', 'remove']:
        print(f"Mode must be {' or '.join(mode)}")
    return(mode.lower())


if __name__ == "__main__":
    BAM, IN_FWD, IN_REV, OUT_FWD, OUT_REV, MODE, PROC = _get_args()

    if OUT_FWD == None:
        out_fwd = f"{IN_FWD.split('/')[-1].split('.')[0]}.r1.fq.gz"
    else:
        out_fwd = OUT_FWD

    if out_fwd.endswith(".gz"):
        write_mode = "wb"
    else:
        write_mode = "w"

    remove_mode = check_remove_mode(MODE)
    BAMFILE = pysam.AlignmentFile(BAM, 'r')

    # FORWARD OR SE FILE
    print(f"- Extracting mapped reads from {BAM}")
    mapped_reads = extract_mapped(proc=PROC)
    print(f"- Parsing forward fq file {IN_FWD}")
    fqd_fwd = parse_fq(IN_FWD)
    print("- Cross referencing mapped reads in forward fq")
    fq_dict_fwd = get_mapped_reads(fqd_fwd, mapped_reads)
    # print(fq_dict_fwd)
    print(f"- Writing forward fq to {out_fwd}")
    write_fq(fq_dict=fq_dict_fwd, fname=out_fwd,
            write_mode=write_mode, remove_mode=remove_mode, proc=PROC)

    # REVERSE FILE
    if IN_REV:
        if OUT_REV == None:
            out_rev = f"{IN_REV.split('/')[-1].split('.')[0]}.r2.fq.gz"
        else:
            out_rev = OUT_REV
        print(f"- Parsing reverse fq file {IN_REV}")
        fqd_rev = parse_fq(IN_REV)
        print("- Cross referencing mapped reads in reverse fq")
        fq_dict_rev = get_mapped_reads(fqd_rev, mapped_reads)
        print(f"- Writing reverse fq to {out_rev}")
        write_fq(fq_dict=fq_dict_rev, fname=out_rev,
                write_mode=write_mode, remove_mode=remove_mode, proc=PROC)
