# Population-analysis-with-pooled-data
This repository contains python scripts for analyzing population data from pooled sequencing.

## AF_per_pool.py
Calculate allele frequency for each sample given AD format out of vcftools.  
Allele frequency is based read count.  
python AF_per_pool.py file.AD.FORMAT > file.freq

## chi_square_contrast.py
Chi-square test of allele depth at each SNP.  
chi_square_contrast.py --popA POPA --popB POPB --AD AD  
The following arguments are required: --popA, --popB, --AD  
popA and popB are the lists of contrasting samples. One sample per row.   
AD is the allele depth after correction of population size.  

## DAF_Categorical_SNP.py
Count the number of SNPs in each genomic category of each binned delta allel frequency.  
DAF_Categorical_SNP.py --DAF DAF --Anno ANNO [--Bin BIN]  
The following arguments are required: --DAF, --Anno  
DAF is the file of delta allele frequency at each SNP.   
Anno is the annotate SNP file out of snpEff. The default bin value is 0.1

## deltaAF.py
Calculate delta allele frequency between the given populations.  
deltaAF.py --popA POPA --popB POPB --freq FREQ  
popA and popB are the lists of contrasting samples. One sample per row.  
freq is the corrected allele frequency file (For example, 60.Neff.freq)  

## Neff_correction.py
Correct allele depth based on sample size.  
Neff_correction.py sample.list pool.AD.FORMAT > pool.Neff.freq  
Example of sample list:  
ID      Size  
A   47  
B 50  
C  49  
D   35  
...  


