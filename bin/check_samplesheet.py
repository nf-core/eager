#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/eager samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample_id	library_id	lane	colour_chemistry	pairment	strandedness	damage_treatment	r1	r2	bam
    Sample1	Sample1_Lib1	1	4	PE	double	full	Sample1_Lib1_L008_R1_001.fq.gz	Sample1_Lib1_L008_R2_001.fq.gz	NA
    Sample2	Sample2_Lib1	2	2	SE	double	full	Sample2_Lib1_L008_R1_001.fq.gz	NA	NA
    Sample3	Sample3_Lib1	9	4	SE	single	none	NA	NA	Sample3_Lib1.bam

    For an example see:
    https://github.com/nf-core/test-datasets/raw/eager/testdata/Mammoth/mammoth_design_fastq_bam_dsl2.tsv
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:

        ## Check header
        MIN_COLS = 10
        HEADER = ["sample_id", "library_id", "lane", "colour_chemistry", "pairment", "strandedness", "damage_treatment", "r1", "r2", "bam"]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format("\t".join(header), "\t".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample_id, library_id, lane, colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam = lspl[: len(HEADER)]

            sample_id = sample_id.replace(" ", "_")
            if not sample_id:
                print_error("[nf-core/eager] error: sample_id entry has not been specified!", "Line", line)

            library_id = library_id.replace(" ", "_")
            if not library_id:
                print_error("[nf-core/eager] error: library_id entry has not been specified!", "Line", line)

            ## Check input file extensions
            for reads in [r1, r2, bam]:
                if reads:
                    if reads.find(" ") != -1:
                        print_error("[nf-core/eager] error: FASTQ or BAM file(s) contains spaces!", "Line", line)
                    if not reads.endswith(".fastq.gz") and not reads.endswith(".fq.gz") and not reads.endswith(".bam") and not reads.endswith("NA"):
                        print_error(
                            "[nf-core/eager] error: FASTQ or BAM file(s) do not have extension '.fastq.gz', '.fq.gz', or '.bam'!",
                            "Line",
                            line,
                        )

            ## TODO UP TO HERE
            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fastq_1, fastq_2]

            if sample_id and fastq_1 and fastq_2:  ## Paired-end short reads
                sample_info = ["0", fastq_1, fastq_2]
            elif sample and fastq_1 and not fastq_2:  ## Single-end short reads
                sample_info = ["1", fastq_1, fastq_2]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ single_end, fastq_1, fastq_2 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "single_end", "fastq_1", "fastq_2"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):

                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error("Multiple runs of a sample must be of the same datatype!", "Sample: {}".format(sample))

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join(["{}_T{}".format(sample, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
