#! /usr/bin/python
# calculate delta allele frequency between two contrasting groups
# input: per pool allele frequency file
# input: popA list
# input popB list

import sys
import argparse

parser = argparse.ArgumentParser(description="Calculate average allele frequency across pools in each contrasting group and calculate delta allele frequency between two groups")
parser.add_argument("--popA", type = str, help="Sample list of population A (no header)", required=True)
parser.add_argument("--popB", type = str, help="Sample list of population B (no header)", required=True)
parser.add_argument("--freq", type = str, help="Per pool allele frequency file", required=True)

args = vars(parser.parse_args())

popA = args["popA"]
popB = args["popB"]
freq = args["freq"]

# print header
h = ["CHR","POS", popA.replace(".list", ""), popB.replace(".list", ""), "deltaAF"]
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

# Go through the freq file
with open(freq) as f:
	header = f.readline()
	header = header.strip().split("\t")
	popA_index = []
	popB_index = []
	for m in range(len(header)):
		if header[m] in popA:
			popA_index.append(m)
		elif header[m] in popB:
			popB_index.append(m)

	if len(popA_index) != n_popA or len(popB_index) != n_popB:
		sys.exit("Not all samples are in the freq file!")

	for l in f:
		new_line = []

		l = l.strip().split("\t")
		popA_list = []
		popB_list = []

		new_line.extend([l[0], l[1]])
		for i in range(2, len(l)):
			if i in popA_index and l[i] != "NA":
				popA_list.append(float(l[i]))
			elif i in popB_index and l[i] != "NA":
				popB_list.append(float(l[i]))

		popA_ave, popB_ave, delAF = [0, 0, 0]
		if len(popA_list) == n_popA and len(popB_list) == n_popB:
			popA_ave = round(sum(popA_list)/n_popA, 2)
			popB_ave = round(sum(popB_list)/n_popB, 2)
			delAF = round((popA_ave - popB_ave), 2)

			if delAF != 0:
				new_line.extend([ str(x) for x in [popA_ave, popB_ave, delAF]])
				new_line = "\t".join(new_line)
				print(new_line)
