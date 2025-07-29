#!/usr/bin/env python
# coding: utf-8

# In[142]:


import pandas as pd


# In[143]:


# In[3]:


# Input genome length (in Kb):# Make a population dictionary to input Population column in the .hom file:

gen = pd.read_csv('genofreq/size_genome.txt', header=None)

GEN_LENGTH =  gen.iloc[0, 0]
f"Reference genome is {GEN_LENGTH} bases long."


# In[144]:


# List populations:
pop = pd.read_csv('genofreq/list_populations.tsv', header=None)
ls_pop = pop.iloc[:, 0].to_list()
f"Populations studied are {ls_pop}."


# In[145]:


# List specimens:
dfspe = pd.read_csv('genofreq/list_specimen_population.tsv', header=None, sep='\\s+', index_col=0)
list_specimens=dfspe[1].to_list()


# In[146]:


# Make a population dictionary to input Population column in the .hom file:
dfpop = pd.read_csv('genofreq/list_specimen_population.tsv', header=None, sep='\\s+', index_col=0)
# Convert the DataFrame to a dictionary of lists
dict_pop_inv = dfpop.groupby(1).apply(lambda x: list(x.index)).to_dict()

# Convert lists with single elements to individual integers
dict_pop_inv = {k: v[0] for k, v in dict_pop_inv.items()}

print(dict_pop_inv)


# In[147]:


dict_pop = dfpop.groupby(1).apply(lambda x: list(x.index)).to_dict()

# Convert lists with single elements to individual integers
dict_pop_inv = {k: v[0] for k, v in dict_pop.items()}

print(dict_pop_inv)


# In[148]:


# Load dataset:
#PATH = '/home/thibauld/Documents/Bioinformatics/Bioinformatics_lab_book/2022_07_18_ROH_coelocentrum'
df = pd.read_csv('plink/all.hom', sep='\\s+')

# Add a specimen column:
df['Specimen'] = df['FID'].astype(str) + '_' + df['IID'].astype(str)


# In[149]:


# Longest/shortest ROH in this dataset (in Kb)?
df['KB'].describe()


# In[150]:


## The Population tag is missing in the .hom file:
## Recover the populations label for each ROH row in the .hom file:
# Empty list:
list_roh_pop=[]

# Put the population label in a list:
ID=df['IID']
for f in ID:
    temp_pop=dict_pop_inv[f] ;
    list_roh_pop.append(temp_pop)


# In[151]:


# Eventually, add the Population column to the .hom file:
df['Population']= list_roh_pop


# In[152]:


# Add a specimen column:
df['Specimen'] = df['IID'].astype(str)


# In[153]:


df


# In[154]:


# Use the list number het for the filtering:
dffilt = pd.read_csv('genofreq/list_min_number_het.txt', header=None, sep='\\s+', index_col=0)
dictionary=dffilt.to_dict()
dict=dictionary[1]


# In[155]:


dffilt


# In[156]:


# In[15]:


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


# In[157]:


# In[16]:


# What % of ROH remaining?
a=len(df.index)
b=len(dffilt.index)
ratio_ROH_remaining=b/a*100
f"There are {ratio_ROH_remaining:.1f}% ROH remaining after min SNPs filtration."


# In[158]:


# In[17]:


# What is the number of samples remaining?
samples_remaining=len(set(dffilt['Specimen'])) 
f"There are {samples_remaining} samples remaining after min SNPs filtration."


# In[159]:


# In[19]:


# Subset the df, and sum the KB column
dfsubset = dffilt[['KB', 'CHR', 'Specimen', 'Population']]


# In[160]:


# In[20]:


# To manipulate the dataframe in order to get sum, we used the command 'copy' explained here:
# https://stackoverflow.com/questions/49728421/pandas-dataframe-settingwithcopywarning-a-value-is-trying-to-be-set-on-a-copy
dfsubset = dfsubset[dfsubset[['KB', 'CHR', 'Specimen', 'Population']].notnull()].copy()


# In[161]:


# In[21]:


# Add the sum to the dataframe:
dfsubset.loc[:,"SROH"]  = dfsubset.groupby(['Specimen', 'CHR'])['KB'].transform('sum')


# In[162]:


# In[22]:


# Now, removing the rows that have same ID, CHR, and SROH value. 
dfsu = dfsubset.drop_duplicates(subset=['Specimen', 'CHR'], keep='last').drop('KB', axis=1)
dfsup = dfsu.drop('Population', axis=1)



# In[163]:


# In[23]:


# Pivot the dataset:
dffinal = dfsup.pivot(index='CHR', columns='Specimen')



# In[164]:


# In[24]:


# Fill NaN:
dfna = dffinal.fillna(0)


# In[165]:


# In[25]:


# Total length of ROH sum per samples:
sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)
dfg = dfr.droplevel(0, axis=0)
dfh = dfg.rename({0: 'SROH'}, axis='columns')
dfsumROH = pd.DataFrame(dfh['SROH'])
dfsumROH


# In[166]:


# In[28]:


# Now, calculate the NROH (total number of ROH)

# Add a count column to the dataframe
df['IID'] = df['Specimen']
# Group by specimen and chromosome, and add the count of ROH on same chr per specimen in the column "count"
cols=['IID', 'CHR']
df['count']  = df.groupby(cols)['IID'].transform('count')
df.sort_values(by='count')
# subset with values of interest.
dfsubset = df[["IID", "FID", "CHR", 'KB', 'count']]


# In[167]:


# In[29]:


# Remove the duplicates counts of ROH, we want only one line per chromosome giving the number of ROH:
dfd = dfsubset.drop_duplicates(['IID', 'CHR']).sort_values(by='count')
dfk = dfd.drop('FID', axis=1)
dfkl = dfk.drop('KB', axis=1)


# In[168]:


# In[30]:


# transpose the dataframe
dffinal = dfkl.pivot(index='CHR', columns='IID')
# Remove NA in table
dfna = dffinal.fillna(0)
# check level of indexing
nblevels = dfna.index.nlevels
# sum of the dataframe for each species
sumroh = dfna.sum(axis=0)


# In[169]:


# In[32]:


dfr = pd.DataFrame(sumroh)
# Remove one level of the index
dfg = dfr.droplevel(0, axis=0)
# Rename the ROH column
dfNROH = dfg.rename({0: 'NROH'}, axis='columns')
dfNROH = pd.concat([dfNROH, dfsumROH], axis=1)


# In[170]:


# # Add FROH

# In[33]:


# Add FROH, the FROH isthe fraction of each genome in ROH > 0.5 Mb => ACTUALLY no ROH > 0.2 Mb!!!! We try this time without filtering.
# No filtering: all data included.
# In the stats_genome_file.txt
# Begonia_peltatifolia_scaffold	310462054  = 310,462,054 = 310462 Kb


# In[171]:


# In[34]:


# Add FROH value (calculation SROH * 1000, as it is in Kb and GEN_LENGTH in bases)
dfNROH['FROH'] = dfNROH['SROH']*1000 / GEN_LENGTH


# In[172]:


# In[35]:


# Add population value

# Create column with ID:
dfNROH['Population'] =  dfNROH.index

# Use dictionary to put pops:
dfNROH['Population'] = dfNROH['Population'].map(dict_pop_inv)


# In[173]:


# In[61]:


dffinal = dfNROH


# In[174]:


# # Add heterozygosity calculated by hand

# In[65]:


df = pd.read_csv('genofreq/geno_freq_all.tsv', low_memory=False, header=0, index_col=0, sep='\t') ;
dfd=df.drop('POS', axis=1) ;

df = dfd.apply(pd.to_numeric, errors='coerce')

idxe = df[df.iloc[:, 0] == 1]


# In[175]:


idxe


# In[176]:


# In[55]:


list_specimens


# In[177]:


# In[56]:


het_list = []
for spe in list_specimens:
   temp = idxe[idxe[spe] == 1].index ;
   heterozygotes = len(temp) ;
   homozygotes = GEN_LENGTH - heterozygotes ;
   heteroygosity = heterozygotes/(heterozygotes+homozygotes) ;
   het_list.append(heteroygosity)


# In[178]:


# In[57]:


het_list


# In[179]:


# In[63]:


#dffinal['het_hand'] = het_list


# In[180]:


# In[64]:


dffinal


# In[181]:
###########################################################################################################
## JELLYFISH UPDATE
## Making the jellyfish heterozygosity table

# Load the data
table = pd.read_csv('jellyfish/table_het_jellyfish.csv')
popidx = pd.read_csv('config/table_reads.tsv', sep='\\s+', low_memory=False)

# Merge on sample names (case-sensitive)
merged = pd.merge(table, popidx, left_on='sample', right_on='Sample', how='inner')

# Drop redundant 'Sample' column if you don’t need it
jellhet = merged.drop(columns=['Sample'])

# Rename columns
jellhet1 = jellhet.rename(columns={
    'sample': 'Sample',
    'heterozygosity': 'Heterozygosity',
    'Population': 'Population_het'
})

# Set 'Sample' as index
jellhetfinal = jellhet1.set_index('Sample')

# angsdhet = pd.read_csv('heterozygosity/heterozygosity_results_angsd.csv', header=0, index_col=0) # JELLYFISH UPDATE ANGSD OUTDATED
dffinal = dffinal.join(jellhetfinal)
###############################################################################################################
# In[66]:


dffinal.to_csv('output/matrix_ROH.csv')


# In[ ]:





# In[ ]:





# In[ ]:





















