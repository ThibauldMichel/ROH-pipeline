#!/usr/bin/env python
# coding: utf-8

# In[62]:


import pandas as pd


# In[63]:


##########################################################################################################################################################


# # Load data

# In[64]:


# Input genome length (in Kb):
gen = pd.read_csv('genofreq/size_genome.txt', header=None)

GEN_LENGTH =  gen.iloc[0, 0]
f"Reference genome is {GEN_LENGTH} bases long."


# In[65]:


# List populations:
pop = pd.read_csv('genofreq/list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()
f"Populations studied are {ls_pop}."


# In[66]:

# Make a population dictionary to input Population column in the .hom file:
dfpop = pd.read_csv(
    'genofreq/list_specimen_population.tsv',
    header=None,
    sep=r'\s+',     # replaces delim_whitespace=True
    index_col=0
)

# Reset and set new index
df1 = dfpop.reset_index()
df2 = df1.set_index(df1.columns[1])

# Get simple dictionary (index -> value)
dictionary = df2[0].to_dict()
print(dictionary)

# In[67]:


# Load dataset:
#PATH = '/home/thibauld/Documents/Bioinformatics/Bioinformatics_lab_book/2022_07_18_ROH_coelocentrum'
df = pd.read_csv('plink/all.hom', delim_whitespace=True)

# Add a specimen column:
df['Specimen'] = df['FID'].astype(str) + '_' + df['IID'].astype(str)


# In[70]:


# Longest/shortest ROH in this dataset (in Kb)?
df['KB'].describe()


# In[ ]:





# In[71]:


##########################################################################################################################################################


# # Input population tag

# In[72]:


## The Population tag is missing in the .hom file:
## Recover the populations label for each ROH row in the .hom file:
# Empty list:
list_roh_pop=[]

# Put the population label in a list:
ID=df['IID']
for f in ID:
    temp_pop=dict_pop_inv[f] ;
    list_roh_pop.append(temp_pop)


# In[73]:


# Eventually, add the Population column to the .hom file:
df['Population']= list_roh_pop


# In[74]:


# Add a specimen column:
df['Specimen'] = df['IID'].astype(str)


# In[75]:


df


# In[76]:


##########################################################################################################################################################


# # Filter the data with min number of SNPs per ROH

# In[77]:


# Use the list number het for the filtering:
dffilt = pd.read_csv('genofreq/list_min_number_het.txt', header=None, delim_whitespace=True, index_col=0)
dictionary=dffilt.to_dict()
dict=dictionary[1]


# In[78]:


# Filtering of the data:
# Threshold value calculated in another script:

#open an empty datafame to concatenate the filtered dataframes:
dffilt  = pd.DataFrame()

# Make a loop through the dictionnary of 'pop:threshold' for filtering:
for f in dict:
# Select the population, filter, and concat the output to get final filtered dataframe:
    dfsub = df[df['Population'] == f] ;
    dftemp = dfsub[dfsub['NSNP'] >= dict[f]] ;
    dffilt = pd.concat([dffilt, dftemp])


# In[79]:


# What % of ROH remaining?
a=len(df.index)
b=len(dffilt.index)
ratio_ROH_remaining=b/a*100
f"There are {ratio_ROH_remaining:.1f}% ROH remaining after min SNPs filtration."


# In[80]:


# What is the number of samples remaining?
samples_remaining=len(set(dffilt['Specimen'])) 
f"There are {samples_remaining} samples remaining after min SNPs filtration."

# save these stats in a .txt document
# Calculate values
a = len(df.index)
b = len(dffilt.index)
ratio_ROH_remaining = b / a * 100
samples_remaining = len(set(dffilt['Specimen']))

# Prepare threshold text
line_threshold = f"The filtration SNP thresholds are {dict}\n"

# Prepare other lines
line1 = f"There are {ratio_ROH_remaining:.1f}% ROH remaining after min SNPs filtration.\n"
line2 = f"There are {samples_remaining} samples remaining after min SNPs filtration.\n"

# Write all to file
with open("plink/filtration_summary.txt", "w") as f:
    f.write(line_threshold)
    f.write(line1)
    f.write(line2)






# Save the dataframe in the plink directory
dffilt
dffilt.to_csv('plink/filtered_ROH.csv')

# In[81]:


##########################################################################################################################################################


# # Calculate SROH

# In[82]:


# Subset the df, and sum the KB column
dfsubset = dffilt[['KB', 'CHR', 'Specimen', 'Population']]


# In[83]:


# To manipulate the dataframe in order to get sum, we used the command 'copy' explained here:
# https://stackoverflow.com/questions/49728421/pandas-dataframe-settingwithcopywarning-a-value-is-trying-to-be-set-on-a-copy
dfsubset = dfsubset[dfsubset[['KB', 'CHR', 'Specimen', 'Population']].notnull()].copy()


# In[84]:


# Add the sum to the dataframe:
dfsubset.loc[:,"SROH"]  = dfsubset.groupby(['Specimen', 'CHR'])['KB'].transform('sum')


# In[85]:


# Now, removing the rows that have same ID, CHR, and SROH value. 
dfsu = dfsubset.drop_duplicates(subset=['Specimen', 'CHR'], keep='last').drop('KB', axis=1)
dfsup = dfsu.drop('Population', axis=1)


# In[86]:


# Pivot the dataset:
dffinal = dfsup.pivot(index='CHR', columns='Specimen')


# In[87]:


# Fill NaN:
dfna = dffinal.fillna(0)


# In[88]:


# Total length of ROH sum per samples:
sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)
dfg = dfr.droplevel(0, axis=0)
dfh = dfg.rename({0: 'SROH'}, axis='columns')
dfsumROH = pd.DataFrame(dfh['SROH'])
dfsumROH


# In[89]:


# To see if the values make sense, calculate the ratio compared to the reference genome:


# In[90]:


##########################################################################################################################################################


# # Calculate NROH

# In[91]:


# Now, calculate the NROH (total number of ROH)

# Add a count column to the dataframe
df['IID'] = df['Specimen']
# Group by specimen and chromosome, and add the count of ROH on same chr per specimen in the column "count"
cols=['IID', 'CHR']
df['count']  = df.groupby(cols)['IID'].transform('count')
df.sort_values(by='count')
# subset with values of interest.
dfsubset = df[["IID", "FID", "CHR", 'KB', 'count']]


# In[92]:


# Remove the duplicates counts of ROH, we want only one line per chromosome giving the number of ROH:
dfd = dfsubset.drop_duplicates(['IID', 'CHR']).sort_values(by='count')
dfk = dfd.drop('FID', axis=1)
dfkl = dfk.drop('KB', axis=1)


# In[93]:


# transpose the dataframe
dffinal = dfkl.pivot(index='CHR', columns='IID')
# Remove NA in table
dfna = dffinal.fillna(0)
# check level of indexing
nblevels = dfna.index.nlevels
# sum of the dataframe for each species
sumroh = dfna.sum(axis=0)


# In[94]:


##########################################################################################################################################################


# # Join NROH and SROH df

# In[95]:


dfr = pd.DataFrame(sumroh)
# Remove one level of the index
dfg = dfr.droplevel(0, axis=0)
# Rename the ROH column
dfNROH = dfg.rename({0: 'NROH'}, axis='columns')
dfNROH.append(dfsumROH)
dfNROH['SROH'] = dfsumROH


# # Add FROH

# In[96]:


# Add FROH, the FROH isthe fraction of each genome in ROH > 0.5 Mb => ACTUALLY no ROH > 0.2 Mb!!!! We try this time without filtering.
# No filtering: all data included.
# In the stats_genome_file.txt
# Begonia_peltatifolia_scaffold	310462054  = 310,462,054 = 310462 Kb


# In[97]:


# Add FROH value (calculation SROH * 1000, as it is in Kb and GEN_LENGTH in bases)
dfNROH['FROH'] = dfNROH['SROH']*1000 / GEN_LENGTH


# In[98]:


# Add population value

# Create column with ID:
dfNROH['Population'] =  dfNROH.index

# Use dictionary to put pops:
dfNROH['Population'] = dfNROH['Population'].map(dict_pop_inv)


# # Add heterozygosity

# In[99]:


dfhet = pd.read_csv('heterozygosity/matrix_heterozygosity.tsv',  delim_whitespace=True, index_col=0, header=None)
dfhet['sum'] = dfhet.iloc[:, 0] + dfhet.iloc[:, 1] + dfhet.iloc[:, 2]
dfhet['heterozygosity'] = dfhet.iloc[:, 2] / dfhet['sum']
dfhetero=pd.DataFrame(dfhet['heterozygosity'])


# In[100]:


dffinal=pd.concat([dfNROH, dfhetero], axis=1)


# In[101]:


dffinal.to_csv('plink/matrix_ROH.csv')


# In[ ]:




