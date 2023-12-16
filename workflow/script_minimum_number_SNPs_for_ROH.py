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



# Making the calculation of l, minimum number SNPs for each population:

for pop in ls_pop:

    # Upload genotype file:
    df = pd.read_csv('genofreq/geno_freq_' + pop + '.tsv', delim_whitespace=True, low_memory=False)

    # Removing non-values of the dataframe:
    dfedit=df.drop(['#CHROM', 'POS'], axis=1)

    # transpose rows -> columns:
    df=dfedit.transpose()

    dfdrop = df.dropna(axis=1)
    
    dfh = dfdrop

    # Make sum of the heterozygotes/homozygotes
    dfh['homozygotes'] = (dfdrop == 0).T.sum()
    dfh['heterozygotes'] = (dfdrop == 1).T.sum()

    df=dfh

    # Calculation frequency heterozygotes
    df['q'] = df['heterozygotes']/(df['homozygotes'] + df['heterozygotes']) 


    # Calculate the mean
    meanhet = df['q'].mean()

    df['het_ratio'] = df['heterozygotes'] / GEN_LENGTH

    dfratio = df['het_ratio'].to_frame()

    # Average number of SNP per sample (will be used in the final equation)
    df['nber_SNPs'] = df['homozygotes'] + df['heterozygotes']

    meansnp =df['nber_SNPs'].mean() 

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

