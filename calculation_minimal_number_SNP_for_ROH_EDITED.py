#!/usr/bin/env python
# coding: utf-8

# In[3]:


# Analysis modules
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy


# In[4]:


# Input path to to directory where are stored the header.tsv file and table_geno.tsv file:
#PATH = '/home/thibauld/Documents/Bioinformatics/Bioinformatics_lab_book/2023_08_14_long_reads_dataset/data'

# Input genome length (in Kb):
GEN_LENGTH = 1182020


# In[5]:


# Calculation of l, the minimum number of SNPs that constitute a ROH


# In[6]:


# Headers VCF file:
# CHROM %POS  %REF  %ALT [ %GT]\n'
#head = pd.read_csv(PATH + '/header.tsv', delim_whitespace=True, low_memory=False)

# Drop columns ID QUAL FILTER INFO FORMAT:
#headedit = head.drop(['ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], axis=1)

# Make a list of these columns
#listhead = headedit.columns.values.tolist()


# In[12]:


# After hand editing, upload the new file:
df = pd.read_csv('geno_freq_A12.tsv', delim_whitespace=True, low_memory=False)


# In[13]:


df


# In[14]:


# Put the header to the dataframe:
#df.columns = listhead


# In[15]:


# Removing non-values of the dataframe:
dfedit=df.drop(['#CHROM', 'POS', 'REF', 'ALT'], axis=1)


# In[33]:


# transpose rows -> columns:
df=dfedit.transpose()


# In[59]:


df


# In[79]:


# START FROM HERE
dfdrop = df.dropna(axis=1)


# In[80]:


# Make sum of the heterozygotes/homozygotes
dfh['homozygotes'] = (dfdrop == '0').T.sum()
dfh['heterozygotes'] = (dfdrop == '1').T.sum()


# In[45]:


# Make sum of the heterozygotes/homozygotes
#dftest['homozygotes'] = (df == 0).T.sum()
#dftest['heterozygotes'] = (df == 1).T.sum()


# In[89]:


df = dfh


# In[90]:


# Calculation frequency heterozygotes
df['q'] = df['heterozygotes']/(df['homozygotes'] + df['heterozygotes']) 

# Calculate the mean
meanhet = df['q'].mean()


# In[91]:


# We will calculate the ratio of heterozygoity on all the positions of the genome:

# Add FROH, the FROH isthe fraction of each genome in ROH > 0.5 Mb => ACTUALLY no ROH > 0.2 Mb!!!! We try this time without filtering.
# No filtering: all data included.
# In the stats_genome_file.txt
# Begonia_peltatifolia_scaffold	310462054  = 310,462,054 = 310462 Kb


# In[92]:


df['het_ratio'] = df['heterozygotes'] / GEN_LENGTH


# In[93]:


dfratio = df['het_ratio'].to_frame()


# In[94]:


dfratio.to_csv('matrix_het.csv')


# In[95]:


# Now, calculate the min number of SNP per sample


# In[96]:


# Average number of SNP per sample (will be used in the final equation)
df['nber_SNPs'] = df['homozygotes'] + df['heterozygotes']

meansnp =df['nber_SNPs'].mean() 


# In[97]:


df


# In[98]:


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
  print(l, file=f)


# In[ ]:




