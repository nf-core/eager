#!/usr/bin/env python3

# Written by Aida Andrades Valtueña and released under MIT license. 
# See git repository (https://github.com/aidaanva/endorS.py) for full license text.

"""Script to calculate the endogenous DNA in a sample from samtools flag stats.
It can accept up to two files: pre-quality and post-quality filtering. We recommend
to use both files but you can also use the pre-quality filtering.
"""
import re
import sys
import json
import argparse
import textwrap

parser = argparse.ArgumentParser(prog='endorS.py',
   usage='python %(prog)s [-h] [--version] <samplesfile>.stats [<samplesfile>.stats]',
   formatter_class=argparse.RawDescriptionHelpFormatter,
   description=textwrap.dedent('''\
   author:
     Aida Andrades Valtueña (aida.andrades[at]gmail.com)

   description:
     %(prog)s calculates endogenous DNA from samtools flagstat files and print to screen
     Use --output flag to write results to a file
   '''))
parser.add_argument('samtoolsfiles', metavar='<samplefile>.stats', type=str, nargs='+',
                    help='output of samtools flagstat in a txt file (at least one required). If two files are supplied, the mapped reads of the second file is divided by the total reads in the first, since it assumes that the <samplefile.stats> are related to the same sample. Useful after BAM filtering')
parser.add_argument('-v','--version', action='version', version='%(prog)s 0.4')
parser.add_argument('--output', '-o', nargs='?', help='specify a file format for an output file. Options: <json> for a MultiQC json output. Default: none')
parser.add_argument('--name', '-n', nargs='?', help='specify name for the output file. Default: extracted from the first samtools flagstat file provided')
args = parser.parse_args()

#Open the samtools flag stats pre-quality filtering:
try:
    with open(args.samtoolsfiles[0], 'r') as pre:
        contentsPre = pre.read()
    #Extract number of total reads
    totalReads = float((re.findall(r'^([0-9]+) \+ [0-9]+ in total',contentsPre))[0])
    #Extract number of mapped reads pre-quality filtering:
    mappedPre = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped ',contentsPre))[0])
    #Calculation of endogenous DNA pre-quality filtering:
    if totalReads == 0.0:
        endogenousPre = 0.000000
        print("WARNING: no reads in the fastq input, Endogenous DNA raw (%) set to 0.000000")
    elif mappedPre == 0.0:
        endogenousPre = 0.000000
        print("WARNING: no mapped reads, Endogenous DNA raw (%) set to 0.000000")
    else:
        endogenousPre = float("{0:.6f}".format(round((mappedPre / totalReads * 100), 6)))
except:
    print("Incorrect input, please provide at least a samtools flag stats as input\nRun:\npython endorS.py --help \nfor more information on how to run this script")
    sys.exit()
#Check if the samtools stats post-quality filtering have been provided:
try:
    #Open the samtools flag stats post-quality filtering:
    with open(args.samtoolsfiles[1], 'r') as post:
        contentsPost = post.read()
    #Extract number of mapped reads post-quality filtering:
    mappedPost = float((re.findall(r'([0-9]+) \+ [0-9]+ mapped',contentsPost))[0])
    #Calculation of endogenous DNA post-quality filtering:
    if totalReads == 0.0:
        endogenousPost = 0.000000
        print("WARNING: no reads in the fastq input, Endogenous DNA modified (%) set to 0.000000")
    elif mappedPost == 0.0:
        endogenousPost = 0.000000
        print("WARNING: no mapped reads, Endogenous DNA modified (%) set to 0.000000")
    else:
        endogenousPost = float("{0:.6f}".format(round((mappedPost / totalReads * 100),6)))
except:
    print("Only one samtools flagstat file provided")
    #Set the number of reads post-quality filtering to 0 if samtools
    #samtools flag stats not provided:
    mappedPost = "NA"

#Setting the name depending on the -name flag:
if args.name is not None:
    name = args.name
else:
    #Set up the name based on the first samtools flagstats:
    name= str(((args.samtoolsfiles[0].rsplit(".",1)[0]).rsplit("/"))[-1])
#print(name)


if mappedPost == "NA":
    #Creating the json file
    jsonOutput={
    "id": "endorSpy",
    "plot_type": "generalstats",
    "pconfig": {
        "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'}
    },
    "data": {
        name : { "endogenous_dna": endogenousPre}
    }
    }
else:
    #Creating the json file
    jsonOutput={
    "id": "endorSpy",
    "plot_type": "generalstats",
    "pconfig": {
        "endogenous_dna": { "max": 100, "min": 0, "title": "Endogenous DNA (%)", "format": '{:,.2f}'},
        "endogenous_dna_post": { "max": 100, "min": 0, "title": "Endogenous DNA Post (%)", "format": '{:,.2f}'}
    },
    "data": {
        name : { "endogenous_dna": endogenousPre, "endogenous_dna_post": endogenousPost}
    },
    }
#Checking for print to screen argument:
if args.output is not None:
   #Creating file with the named after the name variable:
   #Writing the json output:
   fileName = name + "_endogenous_dna_mqc.json"
   #print(fileName)
   with open(fileName, "w+") as outfile:
      json.dump(jsonOutput, outfile)
      print(fileName,"has been generated")
else:
   if mappedPost == "NA":
      print("Endogenous DNA (%):",endogenousPre)
   else:
      print("Endogenous DNA raw (%):",endogenousPre)
      print("Endogenous DNA modified (%):",endogenousPost)
