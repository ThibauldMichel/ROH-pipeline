#!/usr/bin/env python

# Analysis modules
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy

# https://stackoverflow.com/questions/34693991/repel-annotations-in-matplotlib
from adjustText import adjust_text

# Input genome length (in Kb):
gen = pd.read_csv('genofreq/size_genome.txt', header=None)

GEN_LENGTH =  gen.iloc[0, 0]
print(GEN_LENGTH)

# List populations:
pop = pd.read_csv('genofreq/list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()
print(ls_pop)

# Load dataset:
#PATH = '/home/thibauld/Documents/Bioinformatics/Bioinformatics_lab_book/2022_07_18_ROH_coelocentrum'
df = pd.read_csv('plink/GVCF_SNPs.hom', delim_whitespace=True)

# Add a specimen column:
df['Specimen'] = df['FID'].astype(str) + '_' + df['IID'].astype(str)

# ## Filtering data with number SNPs

# Use the list number het for the filtering:
dffilt = pd.read_csv('list_min_number_het.txt', delim_whitespace=True, index_col=0)
dictionary=dffilt.to_dict()
dict=dictionary['Minimum_number_SNPs']

# Filtering of the data:
# Threshold value calculated in another script:

#open an empty datafame:
dffilt  = pd.DataFrame()

# Make a loop through the dictionnary of pop:threshold:
for f in dict:
# Select the population:
    dfsub = df[df['FID'] == f] ;
    dftemp = dfsub[dfsub['NSNP'] >= dict[f]] ;
    dffilt = dffilt.append(dftemp)
    #data = [dfint, dftemp]
    #dffinal= pd.concat(data, axis=1)



dffilt['Specimen'].unique()


# What % of ROH remaining?
a=len(df.index)
b=len(dffilt.index)
print(b/a)

# What is the number of samples remaining?
len(set(dffilt['Specimen'])) 


df = dffilt

#   # Calulate SROH, sum of ROH (Mb)

# Subset the df, and sum the KB column
dfsubset = df[["FID", 'KB', 'CHR', 'Specimen']]
dfsubset.groupby('FID')['KB'].transform('sum')

# Add the sum to the dataframe
dfsubset.loc[:,"SROH"]  = dfsubset.groupby(['Specimen', 'CHR'])['KB'].transform('sum')
dfsubset

# Now, removing the rows that have same ID, CHR, and SROH value. 
dfsu = dfsubset.drop_duplicates(subset=['Specimen', 'CHR'], keep='last').drop('KB', axis=1)
dfsup = dfsu.drop('FID', axis=1)

# Pivot the dataset to have one 
dffinal = dfsup.pivot(index='CHR', columns='Specimen')

dfna = dffinal.fillna(0)

sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)

dfg = dfr.droplevel(0, axis=0)

dfh = dfg.rename({0: 'SROH'}, axis='columns')

dfsumROH = pd.DataFrame(dfh['SROH'])


# # Calculate the NROH (total number of ROH)

# Add a count column to the dataframe

df['IID'] = df['Specimen']

cols=['IID', 'CHR']

df['count']  = df.groupby(cols)['IID'].transform('count')
df.sort_values(by='count')

dfsubset = df[["IID", "FID", "CHR", 'KB', 'count']]
dfsubset.sort_values(by='count')


# Remove the duplicates from the dataframe
dfd = dfsubset.drop_duplicates(['IID', 'CHR']).sort_values(by='count')
dfkl = dfd.drop('KB', axis=1)

# transpose the dataframe
dffinal = dfkl.pivot(index='CHR', columns='IID')

# Remove NA in table
dfna = dffinal.fillna(0)

# check level of indexing
nblevels = dfna.index.nlevels


# In[245]:


# sum of the dataframe for each species
sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)


# In[246]:


# Remove one level of the index
dfg = dfr.droplevel(0, axis=0)


# In[247]:





# In[248]:


# Rename the ROH column
dfNROH = dfg.rename({0: 'NROH'}, axis='columns')


# In[249]:





# In[250]:


dfNROH.append(dfsumROH)


# In[251]:


dfNROH['SROH'] = dfsumROH



# In[252]:


# Add FROH, the FROH isthe fraction of each genome in ROH > 0.5 Mb => ACTUALLY no ROH > 0.2 Mb!!!! We try this time without filtering.
# No filtering: all data included.
# In the stats_genome_file.txt
# Begonia_peltatifolia_scaffold	310462054  = 310,462,054 = 310462 Kb



dfNROH['FROH'] = dfNROH['SROH'] / GEN_LENGTH


# In[253]:


dfNROH['temp'] = dfNROH.index
dfNROH['Population'] = dfNROH["temp"].str[:-2]


# In[254]:


dfNROH.to_csv('matrix_ROH.csv')


