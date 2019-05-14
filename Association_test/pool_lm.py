#! /usr/bin/python3
# fit allele frequency of each SNP and phenotype to linear regression

import sys
from scipy import stats
import pandas as pd
import numpy as np
import math

if len(sys.argv)==1:
    sys.exit("python pool_lm.py pools.Neff.freq sample.info > pools.lm.out")

freq_input = sys.argv[1]
pheno = sys.argv[2]

# read phenotype data into a dict
ph = pd.read_csv(pheno, header=0, sep="\t")
ph = ph.set_index("ID")["Salinity"].to_dict()

# read frequency
# header
print("CHR","POS","R","P","logP",sep="\t")
with open(freq_input) as fr:
    phe_list = []
    header = fr.readline().strip().split("\t")
    for s in range(2,len(header)):
        if header[s] in ph.keys():
            phe_list.append(ph[header[s]])
    
    # go through each SNP
    for l in fr:
        l = l.strip().split("\t")
        if l.count("NA")-2 < 0.5*len(header):
            freq_list = l[2:]
            df = pd.DataFrame({"x": freq_list})
            df["y"] = phe_list
            df = df.replace("NA",np.NaN)
            df = df.dropna()
            df = df.astype(float)
        #print(df)
            slope, intercept, r_value, p_value, std_err = stats.linregress(df.x, df.y)
            if p_value < 0:
                pass
            else:
                print(l[0],l[1],round(r_value,3),round(p_value,3),round(-math.log10(p_value), 2),sep="\t")
    