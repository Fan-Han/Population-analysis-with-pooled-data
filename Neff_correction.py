#! /usr/bin/python3
# correct raw read depth in pooled sequencing data
# input: sample list including pool size; Allele depth of each pool
# output: allele frequency of each pool

import sys
import pandas as pd
import numpy as np

if len(sys.argv) == 1:
    print("python Neff_correction.py sample.list pool.AD.FORMAT > pool.Neff.freq")
    sys.exit(1)

pool_size = sys.argv[1]
AD = sys.argv[2]

# read pool size into a dictionary
size = pd.read_csv(pool_size, header=0, sep="\t")
size = size.set_index("ID")["Size"].to_dict()

# Neff corection
def Neff(AD, n):
    if AD == ".":
        AD = "0,0"
    AD = [int(x) for x in AD.split(",")]
    CT = sum(AD)
    Neff = ((n * CT) - 1)/(n + CT)
    ref = 0
    alt = 0
    if CT != 0:
        ref = int(Neff * (AD[0]/CT))
        alt = int(Neff * (AD[1]/CT))
    #print(AD, CT, n, Neff)
    return ref, alt


# Go through AD file
with open(AD) as f:
    header = f.readline().strip()
    print(header)
    header = header.split("\t")
    
    # handle the header
    size_dict = {}
    for i in range(len(header)):
        if header[i] in size.keys():
            size_dict[i] = size[header[i]]
    
    for line in f:
        line = line.strip()
        line = line.split("\t")
        new_line = []
        for s in range(len(line)):
            if s not in size_dict.keys():
                new_line.append(line[s])
            else:
                cor_ref, cor_alt = Neff(line[s], 2*size_dict[s])
                new_allele = ",".join([str(cor_ref), str(cor_alt)])
                #print(size_dict[s], new_allele)
                new_line.append(new_allele)
        new_line = "\t".join(new_line)
        print(new_line)
                
