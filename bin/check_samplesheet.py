#!/usr/bin/env python

# MIT License (c) The nf-core/eager community
# Modified for nf-core/eager by Thiseas C. Lamnidis (@TCLamnidis)

import os
import sys
import errno
import argparse


def isNAstr(var):
    x = False
    if isinstance(var, str) and var == "NA":
        x = True
    return x


def detect_multistrandedness(all_info_dict, error_counter):
    for sample in all_info_dict.keys():
        lib_strands = []
        for lib in all_info_dict[sample].keys():
            for lane in all_info_dict[sample][lib].keys():
                ## all_info_dict[sample][lib][lane] = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam]
                lib_strands.append(all_info_dict[sample][lib][lane][2])
        if len(set(lib_strands)) > 1:
            error_counter = print_error(
                "Cannot have both single- and double-stranded libraries with the same sample_id.",
                "Sample",
                sample,
                error_counter,
            )
    return error_counter


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


def print_error(error, context="Line", context_str="", error_counter=0):
    if isinstance(context_str, str):
        context_str = "'{}'".format(context_str.strip())
    error_str = "[check_samplesheet.py] Error in samplesheet: {}".format(error)
    if context != "" and context_str != "":
        error_str = "[check_samplesheet.py] Error in samplesheet @ {} {}: {}".format(
            context.strip(), context_str, error
        )
    print(error_str)
    error_counter += 1
    return error_counter


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample_id	library_id	lane	colour_chemistry	pairment	strandedness	damage_treatment	r1	r2	bam	bam_reference_id
    Sample1	Sample1_Lib1	1	4	paired	double	full	Sample1_Lib1_L008_R1_001.fq.gz	Sample1_Lib1_L008_R2_001.fq.gz	NA	NA
    Sample2	Sample2_Lib1	2	2	single	double	full	Sample2_Lib1_L008_R1_001.fq.gz	NA	NA	NA
    Sample3	Sample3_Lib1	9	4	single	single	none	NA	NA	Sample3_Lib1.bam	Reference

    For an example see:
    https://github.com/nf-core/test-datasets/raw/eager/testdata/Mammoth/mammoth_design_fastq_bam_dsl2.tsv
    """

    error_counter = 0
    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 11
        HEADER = [
            "sample_id",
            "library_id",
            "lane",
            "colour_chemistry",
            "pairment",
            "strandedness",
            "damage_treatment",
            "r1",
            "r2",
            "bam",
            "bam_reference_id",
        ]
        header = [x.strip('"') for x in fin.readline().strip().split("\t")]
        if header[: len(HEADER)] != HEADER:
            print("Please check samplesheet header: {} != {}".format("\t".join(header), "\t".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line_num, line in enumerate(fin):
            line_num += 2  ## From 0-based to 1-based. Add an extra 1 for the header line
            lspl = [x.strip().strip('"') for x in line.strip().split("\t")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                error_counter = print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)), "Line", line, error_counter
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                error_counter = print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS), "Line", line, error_counter
                )

            ## Check sample name entries
            (
                sample_id,
                library_id,
                lane,
                colour_chemistry,
                pairment,
                strandedness,
                damage_treatment,
                r1,
                r2,
                bam,
                bam_reference_id,
            ) = lspl[: len(HEADER)]

            sample_id = sample_id.replace(" ", "_")
            if not sample_id:
                error_counter = print_error("sample_id entry has not been specified!", "Line", line_num, error_counter)

            library_id = library_id.replace(" ", "_")
            if not library_id:
                error_counter = print_error("library_id entry has not been specified!", "Line", line_num, error_counter)

            if not lane.isnumeric():
                error_counter = print_error("lane number is not numeric!", "Line", line_num, error_counter)

            if colour_chemistry not in ["2", "4"]:
                error_counter = print_error(
                    "colour_chemistry is not recognised (e.g. 2 for Illumina NextSeq/NovaSeq or 4 Illumina MiSeq/HiSeq/BGI)! Options: 2, 4.",
                    "Line",
                    line_num,
                    error_counter,
                )

            if pairment not in ["paired", "single"]:
                error_counter = print_error(
                    "pairment is not recognised. Options: paired, single.", "Line", line_num, error_counter
                )

            if strandedness not in ["double", "single"]:
                error_counter = print_error(
                    "strandedness is not recognised. Options: double, single", "Line", line_num, error_counter
                )

            if damage_treatment not in ["none", "half", "full"]:
                error_counter = print_error(
                    "damage_treatment is not recognised. Corresponds to UDG treatment. Options: none, half, full.",
                    "Line",
                    line_num,
                    error_counter,
                )

            ## Check input file extensions
            for reads in [r1, r2, bam]:
                if reads.find(" ") != -1:
                    error_counter = print_error(
                        "FASTQ or BAM file(s) contains spaces! Please rename.", "Line", line_num, error_counter
                    )
                if (
                    not reads.endswith(".fastq.gz")
                    and not reads.endswith(".fq.gz")
                    and not reads.endswith(".bam")
                    and not isNAstr(reads)
                ):
                    error_counter = print_error(
                        "FASTQ or BAM file(s) have unrecognised extension. Options: .fastq.gz, .fq.gz, or .bam!",
                        "Line",
                        line,
                        error_counter,
                    )

            if not isNAstr(bam) and not pairment == "single":
                error_counter = print_error(
                    "Pairment for BAM input can only be 'single'.", "Line", line_num, error_counter
                )

            if (not isNAstr(bam) and isNAstr(bam_reference_id)) or (isNAstr(bam) and not isNAstr(bam_reference_id)):
                error_counter = print_error(
                    "A BAM and BAM reference id must always be provided together.", "Line", line_num, error_counter
                )

            if not isNAstr(r1) and not isNAstr(bam_reference_id):
                error_counter = print_error(
                    "FASTQ input cannot have a BAM reference id.", "Line", line_num, error_counter
                )

            ## Prepare meta
            lane_info = (
                []
            )  ## [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam, bam_reference_id]

            if (
                sample_id
                and pairment == "single"
                and not isNAstr(r1)
                and isNAstr(r2)
                and isNAstr(bam)
                and isNAstr(bam_reference_id)
            ):  ## SE: R1 only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam, bam_reference_id]
            elif (
                sample_id
                and pairment == "paired"
                and not isNAstr(r1)
                and not isNAstr(r2)
                and isNAstr(bam)
                and isNAstr(bam_reference_id)
            ):  ## PE: R1 and R2 only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam, bam_reference_id]
            elif (
                sample_id
                and pairment == "single"
                and isNAstr(r1)
                and isNAstr(r2)
                and not isNAstr(bam)
                and not isNAstr(bam_reference_id)
            ):  ## bam input(SE): BAM only
                lane_info = [colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam, bam_reference_id]
            ## Print errors only when pairment is valid but input files don't match pairment
            elif pairment in ["single", "paired"]:
                error_counter = print_error(
                    "Input files don't match pairment. 'single' pairment requires a valid r1 or bam (but not both). 'paired' pairment requires valid r1 and r2.",
                    "Line",
                    line_num,
                    error_counter,
                )

            ## Create a complex structure of dictionaries
            ## sample mapping dictionary = { sample1: [{ library1: [ { lane1: [ colour_chemistry, pairment, strandedness, damage_treatment, r1, r2, bam ] }] }] }
            ## Each sample contains a dictionary that has library IDs as keys, and a dictionary of dictionaries as value. The library ID values are dictionaries with lanes as keys and the lane info as values
            sample_mapping_dict.setdefault(
                sample_id, {}
            )  ## Add the sample id as key with an empty dictionary value if it doesnt exist, else do nothing.
            sample_mapping_dict[sample_id].setdefault(
                library_id, {}
            )  ## Add the library id as key with an empty dictionary value if it doesnt exist, else do nothing.
            # sample_mapping_dict[sample_id][library_id].setdefault(lane, {}) ## Add the lane as key with an empty dictionary value if it doesnt exist, else do nothing.

            ## Throw error if the sample_id/library_id/lane combination already exists.
            if lane in sample_mapping_dict[sample_id][library_id].keys():
                error_counter = print_error(
                    "Each combination of Sample_Id, Library_Id and Lane must be unique!",
                    "Line",
                    line_num,
                    error_counter,
                )
            sample_mapping_dict[sample_id][library_id][
                lane
            ] = lane_info  ## Add lane_info to lane within library within sample.

    ## If formatting errors have occurred print their number and fail.
    if error_counter > 0:
        print(
            "[Formatting check] {} formatting error(s) were detected in the input file. Please check samplesheet.".format(
                error_counter
            )
        )
        sys.exit(1)

    ## Ensure a single library strandedness per sample.
    error_counter = detect_multistrandedness(sample_mapping_dict, error_counter)

    ## If content validation errors have occurred print their number and fail. (e.g. same sample/lib/lane combination appears multiple times, or a line is duplicated.)
    if error_counter > 0:
        print(
            "[Strandedness validation] {} validation error(s) were detected in the input file. Please check samplesheet.".format(
                error_counter
            )
        )
        sys.exit(1)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                "\t".join(
                    [
                        "sample_id",
                        "library_id",
                        "lane",
                        "colour_chemistry",
                        "pairment",
                        "strandedness",
                        "damage_treatment",
                        "r1",
                        "r2",
                        "bam",
                        "bam_reference_id",
                    ]
                )
                + "\n"
            )
            for sample_id in sorted(sample_mapping_dict.keys()):
                for library_id in sorted(sample_mapping_dict[sample_id].keys()):
                    for lane in sorted(sample_mapping_dict[sample_id][library_id].keys()):
                        fout.write(
                            "\t".join([sample_id, library_id, lane, *sample_mapping_dict[sample_id][library_id][lane]])
                            + "\n"
                        )
    else:
        error_counter = print_error("No entries to process!", "", "", error_counter)


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
