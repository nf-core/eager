#!/usr/bin/env python3

# Written by Thiseas C. Lamnidis and released under the MIT license. 
# See git repository (https://github.com/nf-core/eager) for full license text.

import sys, re, json
from collections import OrderedDict

jsonOut=OrderedDict()
data=OrderedDict()

## Function to convert a set of elements into floating point numbers, when possible, else leave them be.
def make_float(x):
    # print (x)
    output=[None for i in range(len(x))]
    ## If value for an estimate/error is -nan, replace with "NA". JSON does not accept NaN as a valid field.
    for i in range(len(x)):
        if x[i] == "-nan" or x[i] == "nan":
            output[i]="N/A"
            continue
        try:
            output[i]=float(x[i])
        except:
            output[i]=x[i]
    
    return(tuple(output))


Input_files=sys.argv[1:]

output = open("nuclear_contamination.txt", 'w')
print ("Individual", "Num_SNPs", "Method1_MOM_estimate", "Method1_MOM_SE", "Method1_ML_estimate", "Method1_ML_SE", "Method2_MOM_estimate", "Method2_MOM_SE", "Method2_ML_estimate", "Method2_ML_SE", sep="\t", file=output)
for fn in Input_files:
    ## For each file, reset the values to "N/A" so they don't carry over from last file.
    mom1, err_mom1= "N/A","N/A"
    ml1, err_ml1="N/A","N/A"
    mom2, err_mom2= "N/A","N/A"
    ml2, err_ml2="N/A","N/A"
    nSNPs="0"
    with open(fn, 'r') as f:
        Estimates={}
        Ind=re.sub('\.X.contamination.out$', '', fn).split("/")[-1]
        for line in f:
            fields=line.strip().split()
            if line.strip()[0:19] == "We have nSNP sites:":
                nSNPs=fields[4].rstrip(",")
            elif line.strip()[0:7] == "Method1" and line.strip()[9:16] == 'new_llh':
                mom1=fields[3].split(":")[1]
                err_mom1=fields[4].split(":")[1]
                ml1=fields[5].split(":")[1]
                err_ml1=fields[6].split(":")[1]
                ## Sometimes angsd fails to run method 2, and the error is printed directly after the SE for ML. When that happens, exclude the first word in the error from the output. (Method 2 jsonOut will be shown as NA)
                if err_ml1.endswith("contamination"):
                    err_ml1 = err_ml1[:-13]
            elif line.strip()[0:7] == "Method2" and line.strip()[9:16] == 'new_llh':
                mom2=fields[3].split(":")[1]
                err_mom2=fields[4].split(":")[1]
                ml2=fields[5].split(":")[1]
                err_ml2=fields[6].split(":")[1]
        ## Convert estimates and errors to floating point numbers
        (ml1, err_ml1, mom1, err_mom1, ml2, err_ml2, mom2, err_mom2) = make_float((ml1, err_ml1, mom1, err_mom1, ml2, err_ml2, mom2, err_mom2))
        data[Ind]={ "Num_SNPs" : int(nSNPs), "Method1_MOM_estimate" : mom1, "Method1_MOM_SE" : err_mom1, "Method1_ML_estimate" : ml1, "Method1_ML_SE" : err_ml1, "Method2_MOM_estimate" : mom2, "Method2_MOM_SE" : err_mom2, "Method2_ML_estimate" : ml2, "Method2_ML_SE" : err_ml2 }
        print (Ind, nSNPs, mom1, err_mom1, ml1, err_ml1, mom2, err_mom2, ml2, err_ml2, sep="\t", file=output)


jsonOut = {"plot_type": "generalstats", "id": "nuclear_contamination",
    "pconfig": {
        "Num_SNPs" : {"title" : "Number of SNPs"},
        "Method1_MOM_estimate" : {"title": "Contamination Estimate (Method1_MOM)"},
        "Method1_MOM_SE" : {"title": "Estimate Error (Method1_MOM)"},
        "Method1_ML_estimate" : {"title": "Contamination Estimate (Method1_ML)"},
        "Method1_ML_SE" : {"title": "Estimate Error (Method1_ML)"},
        "Method2_MOM_estimate" : {"title": "Contamination Estimate (Method2_MOM)"},
        "Method2_MOM_SE" : {"title": "Estimate Error (Method2_MOM)"},
        "Method2_ML_estimate" : {"title": "Contamination Estimate (Method2_ML)"},
        "Method2_ML_SE" : {"title": "Estimate Error (Method2_ML)"}
    }, 
    "data" : data
}
with open('nuclear_contamination_mqc.json', 'w') as outfile:
    json.dump(jsonOut, outfile)
