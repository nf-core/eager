#!/usr/bin/env python3
import sys, re, json
from collections import OrderedDict

data=OrderedDict()

Input_files=sys.argv[1:]

output = open("nuclear_contamination.txt", 'w')
print ("Individual", "Method1_MOM_estimate", "Method1_MOM_SE", "Method1_ML_estimate", "Method1_ML_SE", "Method2_MOM_estimate", "Method2_MOM_SE", "Method2_ML_estimate", "Method2_ML_SE", sep="\t", file=output)
for fn in Input_files:
    ## For each file, reset the values to "N/A" so they don't carry over from last file.
    mom1, err_mom1= "N/A","N/A"
    ml1, err_ml1="N/A","N/A"
    mom2, err_mom2= "N/A","N/A"
    ml2, err_ml2="N/A","N/A"
    with open(fn, 'r') as f:
        Estimates={}
        Ind=re.sub('\.X.contamination.out$', '', fn)
        for line in f:
            fields=line.strip().split()
            if line.strip()[0:19] == "We have nSNP sites:":
                nSNPs=fields[4][:-1]
            elif line.strip()[0:7] == "Method1" and line.strip()[9:16] == 'new_llh':
                mom1=fields[3].split(":")[1]
                err_mom1=fields[4].split(":")[1]
                ml1=fields[5].split(":")[1]
                err_ml1=fields[6].split(":")[1]
                ## Sometimes angsd fails to run method 2, and the error is printed directly after the SE for ML. When that happens, exclude the first word in the error from the output. (Method 2 data will be shown as NA)
                if err_ml1.endswith("contamination"):
                    err_ml1 = err_ml1[:-13]
            elif line.strip()[0:7] == "Method2" and line.strip()[9:16] == 'new_llh':
                mom2=fields[3].split(":")[1]
                err_mom2=fields[4].split(":")[1]
                ml2=fields[5].split(":")[1]
                err_ml2=fields[6].split(":")[1]
        data[Ind]={ "Number_of_SNPs" : nSNPs, "Method1_MOM_estimate" : mom1, "Method1_MOM_SE" : err_mom1, "Method1_ML_estimate" : ml1, "Method1_ML_SE" : err_ml1, "Method2_MOM_estimate" : mom2, "Method2_MOM_SE" : err_mom2, "Method2_ML_estimate" : ml2, "Method2_ML_SE" : err_ml2 }
        print (Ind, nSNPs, mom1, err_mom1, ml1, err_ml1, mom2, err_mom2, ml2, err_ml2, sep="\t", file=output)

with open('nuclear_contamination.json', 'w') as outfile:
    json.dump(data, outfile)
