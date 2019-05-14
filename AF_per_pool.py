#! /usr/bin/python
# calculate reference allele frequency for each pool
# input: file.AD.FORMAT

import sys

if len(sys.argv) <= 1:
	print("python AF_per_pool.py file.AD.FORMAT > file.freq")
	sys.exit(1)

AD = sys.argv[1]

with open(AD) as f:
	header = f.readline().strip()
	print(header)

	for l in f:
		l = l.strip().split("\t")
		new_l = []
		new_l.extend([l[0], l[1]])
		for i in range(2, len(l)):
			if  l[i] == ".":
				l[i] = "0,0"
			pool = [ int(x) for x in l[i].split(",")]
			if sum(pool) == 0 or sum(pool) < 15 or sum(pool) > 300:
				new_l.append("NA")
			else:
				new_l.append(str(round(pool[0]/sum(pool), 2)))
		new_l = "\t".join(new_l)
		print(new_l)
