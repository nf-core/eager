#!/usr/bin/env python

# TODO nf-core: Update the script to check the samplesheet
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv

import os
import sys
import errno
import argparse

def isNAstr(var):
    x=False
    if isinstance(var, str) and var == "NA":
        x=True
    return x

def detect_multistrandedness(all_info_dict):
    for sample in all_info_dict.keys():
        lib_strands=[]
        for lib in all_info_dict[sample].keys():
            for lane in all_info_dict[sample][lib].keys():
                ## all_info_dict[sample][lib][lane] = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]
                lib_strands.append(all_info_dict[sample][lib][lane][2])
        if len(set(lib_strands)) > 1:
            print_error("Cannot have both single- and double-stranded libraries with the same sample_id.", "Sample", sample)


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
    error_str = "[check_samplesheet.py] error: Please check samplesheet. {}".format(error)
    if context != "" and context_str != "":
        error_str = "[check_samplesheet.py] error: Please check samplesheet. {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample_id	library_id	lane	colour_chemistry	pairment	strandedness	damage_treatment	r1	r2	bam
    Sample1	Sample1_Lib1	1	4	paired	double	full	Sample1_Lib1_L008_R1_001.fq.gz	Sample1_Lib1_L008_R2_001.fq.gz	NA
    Sample2	Sample2_Lib1	2	2	single	double	full	Sample2_Lib1_L008_R1_001.fq.gz	NA	NA
    Sample3	Sample3_Lib1	9	4	single	single	none	NA	NA	Sample3_Lib1.bam

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
            print("Please check samplesheet header: {} != {}".format("\t".join(header), "\t".join(HEADER)))
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
                print_error("sample_id entry has not been specified!", "Line:", line)

            library_id = library_id.replace(" ", "_")
            if not library_id:
                print_error("library_id entry has not been specified!", "Line:", line)

            if not lane.isnumeric():
                print_error("lane number is not numeric!", "Line:", line)

            if colour_chemistry not in ['2', '4']:
                print_error("colour_chemistry is not recognised (e.g. 2 for Illumina NextSeq/NovaSeq or 4 Illumina MiSeq/HiSeq/BGI)! Options: 2, 4.", "Line:", line)

            if pairment not in ['paired', 'single']:
                print_error("pairment is not recognised. Specify single for BAM input. Options: paired, single.", "Line", line)

            if strandedness not in ['double', 'single']:
                print_error("strandedness is not recognised. Options: double, single", "Line", line)

            if damage_treatment not in ['none', 'half', 'full']:
                print_error("damage_treatment is not recognised. Corresponds to UDG treatment. Options: none, half, full.", "Line", line)

            ## Check input file extensions
            for reads in [r1, r2, bam]:
                if reads.find(" ") != -1:
                    print_error("FASTQ or BAM file(s) contains spaces! Please rename.", "Line", line)
                if not reads.endswith(".fastq.gz") and not reads.endswith(".fq.gz") and not reads.endswith(".bam") and not isNAstr(reads):
                    print_error(
                        "FASTQ or BAM file(s) have unrecognised extension. Options: .fastq.gz, .fq.gz, or .bam!",
                        "Line",
                        line,
                    )

            if not isNAstr(bam) and not pairment == 'single':
                print_error("Pairment for BAM input can only be 'single'.", "Line", line)

            ## Prepare meta
            lane_info = []  ## [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]

            if sample_id and not isNAstr(r1) and isNAstr(r2) and isNAstr(bam): ## R1 only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]
            elif sample_id and not isNAstr(r1) and not isNAstr(r2) and isNAstr(bam): ## R1 and R2 only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]
            elif sample_id and isNAstr(r1) and isNAstr(r2) and not isNAstr(bam): ## BAM only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create a complex structure of dictionaries
            ## sample mapping dictionary = { sample1: [{ library1: [ { lane1: [ colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam ] }] }] }
            ## Each sample contains a dictionary that has library IDs as keys, and a dictionary of dictionaries as value. The library ID values are dictionaries with lanes as keys and the lane info as values
            sample_mapping_dict.setdefault(sample_id, {}) ## Add the sample id as key with an empty dictionary value if it doesnt exist, else do nothing. 
            sample_mapping_dict[sample_id].setdefault(library_id, {}) ## Add the library id as key with an empty dictionary value if it doesnt exist, else do nothing. 
            # sample_mapping_dict[sample_id][library_id].setdefault(lane, {}) ## Add the lane as key with an empty dictionary value if it doesnt exist, else do nothing. 
            
            ## Throw error if the sample_id/library_id/lane combination already exists.
            if lane in sample_mapping_dict[sample_id][library_id].keys():
                print_error("Samplesheet contains duplicate rows!", "Line", line)
            sample_mapping_dict[sample_id][library_id][lane] = lane_info ## Add lane_info to lane within library within sample.

    ## Ensure a single library strandedness per sample.
    detect_multistrandedness(sample_mapping_dict)

            # ## Create sample mapping dictionary = { sample: [ library_id, lane, colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam ] }
            # ## TODO Thiseas: check whether this can work with sample_id and library_id
            # if sample_id not in sample_mapping_dict:
            #     sample_mapping_dict[sample_id] = [sample_info]
            # else:
            #     if sample_info in sample_mapping_dict[sample_id]:
            #         print_error("Samplesheet contains duplicate rows!", "Line", line)
            #     else:
            #         sample_mapping_dict[sample_id].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write("\t".join(["sample_id", "library_id", "lane", "colour_chemistry", "pairment", "strandedness", "damage_treatment", "r1", "r2", "bam"]) + "\n")
            for sample_id in sorted(sample_mapping_dict.keys()):
                for library_id in sorted(sample_mapping_dict[sample_id].keys()):
                    for lane in sorted(sample_mapping_dict[sample_id][library_id].keys()):
                        fout.write("\t".join([sample_id, library_id, lane, *sample_mapping_dict[sample_id][library_id][lane]]) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
