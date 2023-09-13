#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Analysis modules
import numpy as np   
import pandas as pd
#import matplotlib as mpl
#import matplotlib.pyplot as plt
#import seaborn as sns
#import sys
#import scipy

# Input genome length (in Kb):
GEN_LENGTH = 1182020

# List populations:
pop = pd.read_csv('list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()


for pop in ls_pop:

    # Upload genotype file:
    df = pd.read_csv('geno_freq_' + pop + '.tsv', delim_whitespace=True, low_memory=False)

    # Removing non-values of the dataframe:
    dfedit=df.drop(['#CHROM', 'POS', 'REF', 'ALT'], axis=1)

    # transpose rows -> columns:
    df=dfedit.transpose()

    # Make sum of the heterozygotes/homozygotes
    df['homozygotes'] = (df == 0).T.sum()
    df['heterozygotes'] = (df == 1).T.sum()

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

    with open("list_min_number_het.txt", "a") as f:
          print(pop, l, file=f)


# In[ ]:




