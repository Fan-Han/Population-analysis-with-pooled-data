#! /usr/bin/python3
# Calculate the allele frequency for particular populations

import sys
import argparse

parser = argparse.ArgumentParser(description = "Calculate the minor allele frequency for particular group of samples")
parser.add_argument("--pop", type = str, help = "Population list that the MAF will calculated in. This flag can be specified multiple times", required = True, action = "append")
parser.add_argument("--AD", type = str, help = "File of allelic depth at each site of each pool", required = True)

args = vars(parser.parse_args())

AD = args["AD"]
pop = args["pop"]

# read all samples into a list
def read_sample(my_list):
    out_list = []
    for file in my_list:
        with open(file) as f:
            for l in f:
                l = l.strip().split("\t")
                out_list.append(l[0])
    return out_list

samples = read_sample(pop)

# header
print("CHR","POS","MAF", sep="\t")

# Go through AD file
with open(AD) as ad:
    header = ad.readline()
    header = header.strip().split("\t")
    sample_index = []
    for i in range(len(header)):
        if header[i] in samples:
            sample_index.append(i)
    
    # check if all samples are in the AD file
    #print(len(sample_index))
    if len(sample_index) != len(samples):
        sys.exit("Not all samples are in the AD file!")
    
    # read allelic depth
    for l in ad:
        l = l.strip().split("\t")
        CHR = l[0]
        POS = l[1]
        ref = 0
        alt = 0
        for x in range(len(l)):
            if x in sample_index and l[x] != ",":
                x = [int(n) for n in l[x].split(",")]
                ref += x[0]
                alt += x[1]
        
        # MAF
        if ref + alt != 0:
            MAF = round(min(ref, alt)/(ref + alt),3)
            if MAF > 0:
                print(CHR, POS, MAF, sep="\t")
