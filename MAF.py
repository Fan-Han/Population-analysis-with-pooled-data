#! /usr/bin/python3
# Calculate the allele frequency for particular populations

import sys
import argparse

parser = argparse.ArgumentParser(description = "Calculate the minor allele frequency for particular group of samples")
parser.add_argument("--pop", type = str, help = "Population list that the MAF will calculated in. This flag can be specified multiple times", required = True, action = "append")
parser.add_argument("--AD", type = str, help = "File of allelic depth at each site of each pool", required = True)

args = parser.parse_arg()
print(args)