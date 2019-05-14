#! /usr/bin/python
# input: Neff correct allele count 60.Neff.AD
# input: popA.list
# input: popB.list

import sys
import argparse
from scipy.stats.distributions import chi2
import math

parser = argparse.ArgumentParser(description = "Perform genome-wide contrast between two pools by chi-square test on allelic counts")
parser.add_argument("--popA", type = str, help="Sample list of population A (no header)", required=True)
parser.add_argument("--popB", type = str, help="Sample list of population B (no header)", required=True)
parser.add_argument("--AD", type = str, help="Allelic depth file", required=True)

args = vars(parser.parse_args())

popA = args["popA"]
popB = args["popB"]
AD = args["AD"]

# print header
#h = ["CHR","POS", popA.replace(".list", ""), popB.replace(".list", ""), "chi", "log10P"]
h = ["CHR","POS","log10P"]
h = "\t".join(h)
print(h)

# read population list into a list
def read_pop(file):
    pop = []
    with open(file) as f:
        for l in f:
            l = l.strip()
            l = l.split("\t")
            l = l[0]
            pop.append(l)
    return pop

popA = read_pop(popA)
n_popA = len(popA)
popB = read_pop(popB)
n_popB = len(popB)

# chi-square test
def chi_square(ADA, ADB):
    ADA = [int(x) for x in ADA.split(",")]
    ADB = [int(x) for x in ADB.split(",")]
    refA, altA = ADA
    refB, altB = ADB
    
    total = refA + altA + refB + altB
    exp_refA = (refA + refB) * (refA + altA) / total
    exp_refB = (refA + refB) * (refB + altB) / total
    exp_altA = (altA + altB) * (refA + altA) / total
    exp_altB = (altA + altB) * (refB + altB) / total
    
    chi_value = (exp_refA - refA)**2/exp_refA + (exp_refB - refB)**2/exp_refB + (exp_altA - altA)**2/exp_altA + (exp_altB - altB)**2/exp_altB
    chi_value = round(chi_value, 2)
    
    p_value = chi2.sf(chi_value, 1)
    if p_value == 0:
        p_value += 0.1
    p_value = round(-math.log10(p_value), 2)
    return str(chi_value), str(p_value)

# Go through the AD file
with open(AD) as ad:
    header = ad.readline()
    header = header.strip().split("\t")
    popA_index = []
    popB_index = []
    for m in range(len(header)):
        if header[m] in popA:
            popA_index.append(m)
        elif header[m] in popB:
            popB_index.append(m)
            
    if len(popA_index) != n_popA or len(popB_index) != n_popB:
        sys.exit("Not all samples are in the AD file!")
    
    for l in ad:
        l = l.strip().split("\t")
        popA_ref = 0
        popA_alt = 0
        popB_ref = 0
        popB_alt = 0
        
        for i in range(len(l)):
            if i in popA_index:
                alle_A = [int(x) for x in l[i].split(",")]
                popA_ref += alle_A[0]
                popA_alt += alle_A[1]
            elif i in popB_index:
                alle_B = [int(x) for x in l[i].split(",")]
                popB_ref += alle_B[0]
                popB_alt += alle_B[1]
        
        if (popA_ref + popA_alt + popB_ref + popB_alt == 0) or (popA_ref + popB_ref == 0) or (popA_alt + popB_alt == 0) or (popA_ref + popA_alt == 0) or (popB_ref + popB_alt == 0) or (popA_ref + popA_alt < 20*n_popA) or (popB_ref + popB_alt < 20*n_popB) or (popA_ref + popA_alt > 200*n_popA) or (popB_ref + popB_alt > 200*n_popB):
            pass
        else:
            popA = ",".join([str(popA_ref), str(popA_alt)])
            popB = ",".join([str(popB_ref), str(popB_alt)])
            chi_value, p_value = chi_square(popA, popB)
            if float(chi_value) <= 0:
                pass
            else:
                #new_line = [l[0], l[1], popA, popB, chi_value, p_value]
                new_line = [l[0], l[1], p_value]
                new_line = "\t".join(new_line)
                print(new_line)
