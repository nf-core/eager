#!/usr/bin/env python3

# Written by Thiseas C. Lamnidis and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import sys, json
from collections import OrderedDict

jsonOut = OrderedDict()
data = OrderedDict()


input = open(sys.argv[1], 'r')
for line in input:
  fields = line.strip().split()
  sample_id = fields[0]
  covered_snps = fields[1]
  total_snps = fields[2]
  if sample_id[0] == "#":
    continue
  
  data[sample_id] = {"Covered_Snps":covered_snps, "Total_Snps":total_snps}

jsonOut = {"plot_type": "generalstats", "id": "snp_coverage",
    "pconfig": {
        "Covered_Snps" : {"title" : "#SNPs Covered"},
        "Total_Snps" : {"title": "#SNPs Total"}
    }, 
    "data" : data
}

with open(sys.argv[1].rstrip('.txt')+'_mqc.json', 'w') as outfile:
    json.dump(jsonOut, outfile)
