#!/usr/bin/env python3

import argparse
import multiprocessing
import pysam
from functools import partial
import gzip
import sys


def _get_args():
    '''This function parses and return arguments passed in'''
    parser = argparse.ArgumentParser(
        prog='extract_mapped_reads',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=f'''
Remove mapped in bam file from fastq files
        ''')
    parser.add_argument('bam_file', help="path to bam file")
    parser.add_argument('fwd', help='path to forward fastq file')
    parser.add_argument(
        '-2',
        dest="rev",
        default=None,
        help="path to forward fastq file")
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
    proc = int(args.process)

    return(bam, in_fwd, in_rev, out_fwd, out_rev, proc)


def extract_mapped_chr(chr, bam):
    """
    Get mapped reads per chromosome
    INPUT:
    - chr(str): chromosome
    - bam(str): bamfile path
    OUTPUT:
    - res(list): list of mapped reads (str) name per chromosome
    """
    res = []
    bamfile = pysam.AlignmentFile(bam, "rb")
    reads = bamfile.fetch(chr, multiple_iterators=True)
    for read in reads:
        res.append(read.query_name)
    return(res)


def extract_mapped(bam, processes):
    """
    Get mapped reads in parallel
    INPUT:
    - bam(str): bamfile path
    OUTPUT:
    - result(list) list of mapped reads name (str)
    """
    try:
        bamfile = pysam.AlignmentFile(bam, "rb")
        chrs = bamfile.references
    except ValueError as e:
        print(e)
    extract_mapped_chr_partial = partial(extract_mapped_chr, bam=bam)
    p = multiprocessing.Pool(processes)
    res = p.map(extract_mapped_chr_partial, chrs)
    p.close()
    p.join()
    result = [i for ares in res for i in ares]
    return(result)


def parse_fq(fq):
    """
    Parse a FASTQ file
    INPUT:
    - fq(str): path to fastq file
    OUTPUT:
    - fqd(dict): dictionary with read names as keys, seq and quality as values
        in a list
    """

    def get_fq_reads(allreads):
        fqd = {}
        myflag = True
        for line in allreads:
            line = line.decode('utf-8').rstrip()
            if myflag == True:
                instrument = line.split()[0].split(":")[0]
                myflag = False
            if line.startswith(instrument):
                seqname = line[1:].split()[0]
                fqd[seqname] = []
                continue
            else:
                fqd[seqname].append(line)
        return(fqd)

    if fq.endswith('.gz'):
        with gzip.open(fq, 'rb') as allreads:
            fqd = get_fq_reads(allreads)
    else:
        with open(fq, 'r') as allreads:
            fqd = get_fq_reads(allreads)

    return(fqd)


def remove_mapped(fq_dict, mapped_reads):
    """
    Remove mapped reads from dictionary of fastq reads
    INPUT:
    - fq_dict(dict) dictionary with read names as keys, seq and quality as values
        in a list
    - mapped_reads(list) list of mapped reads
    OUTPUT:
    - ufqd(dict) dictionary with unmapped read names as keys, seq and quality as values
        in a list
    """
    ufqd = {}
    unmap = [i for i in list(fq_dict.keys()) if i not in mapped_reads]
    # print(unmap)
    for r in unmap:
        ufqd[r] = fq_dict[r]

    return(ufqd)


def write_fq(fq_dict, fname):
    """
    Write to fastq file
    INPUT:
    - fq_dict(dict) dictionary with unmapped read names as keys, seq and quality as values
        in a list
    - fname(string) Path to output fastq file
    """

    if fname.endswith('.gz'):
        with gzip.open(fname, 'wb') as f:
            for k in list(fq_dict.keys()):
                f.write(f"{k}\n".encode())
                for i in fq_dict[k]:
                    f.write(f"{i}\n".encode())

    else:
        with open(fname, 'w') as f:
            for k in list(fq_dict.keys()):
                f.write(f"{k}\n")
                for i in fq_dict[k]:
                    f.write(f"{i}\n")


if __name__ == "__main__":
    BAM, IN_FWD, IN_REV, OUT_FWD, OUT_REV, PROC = _get_args()

    if IN_REV and not OUT_REV:
        print('You specified an input reverse fastq, but no output reverse fastq')
        sys.exit(1)

    if OUT_FWD == None:
        out_fwd = f"{IN_FWD.split('/')[-1].split('.')[0]}.r1.fq.gz"
    else:
        out_fwd = OUT_FWD

    mapped_reads = extract_mapped(BAM, PROC)
    fwd_dict = parse_fq(IN_FWD)
    fwd_unmap = remove_mapped(fwd_dict, mapped_reads)
    write_fq(fwd_unmap, out_fwd)
    if IN_REV:
        if OUT_REV == None:
            out_rev = f"{IN_REV.split('/')[-1].split('.')[0]}.r2.fq.gz"
        else:
            out_rev = OUT_REV
        rev_dict = parse_fq(IN_REV)
        rev_unmap = remove_mapped(rev_dict, mapped_reads)
        write_fq(rev_unmap, out_rev)
