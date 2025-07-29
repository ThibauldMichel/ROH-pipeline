# Analysis modules
import numpy as np   
import pandas as pd
import gzip
import re
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns
#import sys
#import scipy

# Define the function used to save the data:
def write_values_to_file(file_path, value_a, value_b):
     with open(file_path, "a") as file:
         file.write(f"{value_a}\t{value_b}\n")




# Input genome length (in Kb):
gen = pd.read_csv('genofreq/size_genome.txt', header=None)

GEN_LENGTH =  gen.iloc[0, 0]
print(GEN_LENGTH)

# List populations:
pop = pd.read_csv('genofreq/list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()
print(ls_pop)

##################################################################################
# JELLYFISH UPDATE


# List populations:
pop = pd.read_csv('genofreq/list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()
print(ls_pop)

## Making the jellyfish heterozygosity table

# Load the data
table = pd.read_csv('jellyfish/table_het_jellyfish.csv')
popidx = pd.read_csv('config/table_reads.tsv', sep='\\s+', low_memory=False)

# Merge on sample names (case-sensitive)
merged = pd.merge(table, popidx, left_on='sample', right_on='Sample', how='inner')

# Drop redundant 'Sample' column if you don’t need it
jellhet = merged.drop(columns=['Sample'])
#######################################################################################


# Making the calculation of l, minimum number SNPs for each population:

for pop in ls_pop:

    # Upload genotype file:
##########################################################################
    df = pd.read_csv('genofreq/geno_freq_' + str(pop) + '.tsv', sep='\\s+', low_memory=False)
    dfi = df.drop("#CHROM", axis=1).drop("POS", axis=1)

#Convert all data to numeric, coercing errors to NaN
    df = dfi.apply(pd.to_numeric, errors='coerce')

# Sum homozygotes (0) and heterozygotes (1)
    sum_0s = (df == 0).sum()
    sum_1s = (df == 1).sum()

# Create new rows
    new_row_0 = pd.DataFrame([sum_0s], index=['homozygotes'])
    new_row_1 = pd.DataFrame([sum_1s], index=['heterozygotes'])

# Append new rows to the original DataFrame
    df = pd.concat([df, new_row_0, new_row_1])

# Transpose rows -> columns
    dfh = df.transpose()
#########################################################################
    # Calculate the mean heterozygosity in the population according to jellyfish
    subjell=jellhet[jellhet['Population'] == 1]
    meanhet=subjell['heterozygosity'].mean()
#########################################################################

    df=dfh

    # Calculation frequency heterozygotes
    df['q'] = df.loc[:,'heterozygotes']/(df.loc[:,'homozygotes'] + df.loc[:,'heterozygotes']) 


    # Calculate the mean
    # meanhet = df.loc[:,'q'].mean() JELLYFISH UPDATE IS NOT INTERESTED

    df.loc[:,'het_ratio'] = df.loc[:,'heterozygotes'] / GEN_LENGTH

    dfratio = df.loc[:,'het_ratio'].to_frame()

    # Average number of SNP per sample (will be used in the final equation)
    df.loc[:,'nber_SNPs'] = df.loc[:,'homozygotes'] + df.loc[:,'heterozygotes']

    meansnp =df.loc[:,'nber_SNPs'].mean() 

    # Now calculate l, the minimum number of SNPs that constitute a ROH

    # Number of individuals
    ni = len(df)
    # Number of SNPs (# mean number of SNPs)
    #ns = len(df.index)
    ns=meansnp
    # mean SNP heterozygosity
    het=meanhet

    # calculation l
    l = np.log(0.05/(ns*ni))/np.log(1-het)

    # print result
    print(l)

    # Save in a file:
    write_values_to_file("genofreq/list_min_number_het.txt", pop, l)

