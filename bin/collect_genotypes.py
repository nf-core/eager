#!/usr/bin/env python

# MIT License (c) Thiseas C. Lamnidis (@TCLamnidis)

import argparse
import filecmp

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

## A function to return the number of genotypes per line in a .geno file.
def file_width(fname):
    with open(fname) as f:
        for i in f:
            return(len(i.strip()))
            break

## A function to check that there are no duplicate individual IDs across ind files.
def check_for_duplicate_ids(indf1, indf2):
    with open(indf1) as f:
        inds1 = [x.strip().split()[0] for x in f.readlines()]
    with open(indf2) as f:
        inds2 = [x.strip().split()[0] for x in f.readlines()]
    intersection = set(inds1).intersection(inds2)
    if len(intersection) > 0:
        raise IOError("Input .ind files contain duplicate individual IDs. Duplicates: {}".format(intersection))

## Function to check that the snp files are identical
def check_snp_files(snpf1, snpf2):
    if not filecmp.cmp(snpf1, snpf2):
        raise IOError("Input .snp files are not identical.")

## Function to check the consistency of an eigenstrat database
def validate_eigenstrat(genof, snpf, indf):
    dimsGeno = [file_len(genof), file_width(genof)]
    linesSnp = file_len(snpf)
    linesInd = file_len(indf)

    # print(dimsGeno,linesSnp,linesInd)
    ##Check geno and snp compatibility
    if dimsGeno[0] !=    linesSnp:
        raise IOError("Input .snp and .geno files do not match.")

    ##Check geno and ind compatibility
    if dimsGeno[1] !=    linesInd:
        raise IOError("Input .ind and .geno files do not match.")

VERSION = "1.0.0"

parser = argparse.ArgumentParser(usage="%(prog)s (-i <Input file prefix>) (-c <input ind file> | -R | -E) [-L <SAMPLE LIST> | -S Ind [-S Ind2]] [-o <OUTPUT FILE PREFIX>]" , description="A tool to check two different EingenStrat databses for shared individuals, and extract or remove individuals from an EigenStrat database.")
parser._optionals.title = "Available options"
parser.add_argument("-g1", "--genoFn1", type = str, metavar = "<GENO FILE 1 NAME>",    required = True, help = "The path to the input geno file of the first dataset.")
parser.add_argument("-s1", "--snpFn1",  type = str, metavar = "<SNP FILE 1 NAME>",     required = True, help = "The path to the input snp file of the first dataset.")
parser.add_argument("-i1", "--indFn1",  type = str, metavar = "<IND FILE 1 NAME>",     required = True, help = "The path to the input ind file of the first dataset.")
parser.add_argument("-g2", "--genoFn2", type = str, metavar = "<GENO FILE 2 NAME>",    required = True, help = "The path to the input geno file of the second dataset.")
parser.add_argument("-s2", "--snpFn2",  type = str, metavar = "<SNP FILE 2 NAME>",     required = True, help = "The path to the input snp file of the second dataset.")
parser.add_argument("-i2", "--indFn2",  type = str, metavar = "<IND FILE 2 NAME>",     required = True, help = "The path to the input ind file of the second dataset.")
parser.add_argument("-o", "--output",   type = str, metavar = "<OUTPUT FILES PREFIX>", required = True, help = "The desired output file prefix. Three output files are created, <OUTPUT FILES PREFIX>.geno , <OUTPUT FILES PREFIX>.snp and <OUTPUT FILES PREFIX>.ind .")
parser.add_argument("-v", "--version", action='version', version="{}".format(VERSION), help="Print the version and exit.")
args = parser.parse_args()

## Open input files
GenoFile1 = open(args.genoFn1, "r")
SnpFile1  = open(args.snpFn1,  "r")
IndFile1  = open(args.indFn1,  "r")

GenoFile2 = open(args.genoFn2, "r")
# SnpFile2  = open(args.snpFn2,  "r") ## Never actually read in line by line
IndFile2  = open(args.indFn2,  "r")

## open output files
GenoFileOut = open(args.output + ".geno", "w")
SnpFileOut  = open(args.output + ".snp",  "w")
IndFileOut  = open(args.output + ".ind",  "w")

## Perform basic validation on inputs
validate_eigenstrat(args.genoFn1, args.snpFn1, args.indFn1)
validate_eigenstrat(args.genoFn2, args.snpFn2, args.indFn2)
check_for_duplicate_ids(args.indFn1, args.indFn2)
check_snp_files(args.snpFn1, args.snpFn2)

## Now actually merge the data
## Geno
for line1, line2 in zip(GenoFile1, GenoFile2):
    geno_line="{}{}".format(line1.strip(),line2.strip())
    print(geno_line, file=GenoFileOut)

## Snp
##  Copying the file would be faster, but this way we do not rely on the os or external packages.
##  We already checked that the snp files are byte-identical, so we can just copy one of them.
for line in SnpFile1:
    print(line.strip(), file=SnpFileOut)

## Ind
##  The indfiles are simply concatenated in the same order as the geno file.
for line in IndFile1:
    print(line.strip(), file=IndFileOut)
for line in IndFile2:
    print(line.strip(), file=IndFileOut)
