#! /usr/bin/python3
# calculate number of SNPs under each genomic category in each DeltaAF bin

import sys
import argparse
from collections import defaultdict
import math
import numpy as np

parser = argparse.ArgumentParser(description = "calculate number of SNPs under each genomic category in each DeltaAF bin")
parser.add_argument("--DAF", type = str, help = "File of delta allele frequency between two populations", required = True)
parser.add_argument("--Anno", type = str, help = "Annotation file of each SNP from snpEff)", required = True)
parser.add_argument("--Bin", type = float, help = "Size of bin for DAF. Default is 0.1", default = 0.1)

args = vars(parser.parse_args())

DAF = args["DAF"]
Anno = args["Anno"]
Bin = args["Bin"]

# read annotation file into a dictionary
dict_Anno = defaultdict(dict)
with open(Anno) as A:
    for l in A:
        if "WARNING" not in l:
            l = l.strip().split("\t")
            chrom = l[0]
            pos = l[1]
            anno = l[4].split("|")
            if anno[1] in ["5_prime_UTR_variant","3_prime_UTR_variant","missense_variant","synonymous_variant","downstream_gene_variant","upstream_gene_variant","intergenic_region","intron_variant"]:
                dict_Anno[chrom][pos] = anno[1]


# create a list of bins

# a function to find two closes value in a list
def two_closest(x):
    for i in np.arange(0, 1+Bin, Bin):
        if(x > i and x <= i+Bin):
            return i+Bin


# go through DAF file
dict_bin = defaultdict(lambda:defaultdict(int))
with open(DAF) as D:
    header = D.readline()
    for l in D:
        l = l.strip().split("\t")
        # find the bin of the dAF
        daf = abs(float(l[4]))
        bin_end = two_closest(daf)
        try:
            snp_anno = dict_Anno[l[0]][l[1]]
            dict_bin[bin_end][snp_anno] += 1
        except KeyError:
            pass


for a,b in dict_bin.items():
    for c,d in b.items():
        print(round(a,2),c,d, sep="\t")
