#!/usr/bin/env python
# coding: utf-8

# In[48]:


# Analysis modules
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy

# https://stackoverflow.com/questions/34693991/repel-annotations-in-matplotlib
from matplotlib import pyplot as plt
from adjustText import adjust_text


# In[49]:


# Load dataset:
#PATH = '/home/thibauld/Documents/Bioinformatics/Bioinformatics_lab_book/2022_07_18_ROH_coelocentrum'
df = pd.read_csv('both.hom', delim_whitespace=True)
df


# ## Filtering data with number SNPs

# In[50]:


# Use the list number het for the filtering:
dffilt = pd.read_csv('list_min_number_het.txt', delim_whitespace=True, index_col=0)
dffilt


# In[51]:


dictionary=dffilt.to_dict()
dict=dictionary['Minimum_number_SNPs']
dict


# In[52]:


# Filtering of the data:
# Threshold value calculated in another script:

#open an empty 
dffilt  = pd.DataFrame()

# Make a loop through the dictionnary of pop:threshold:
for f in dict:
# Select the population:
    dfsub = df[df['IID'] == f] ;
    dftemp = dfsub[dfsub['NSNP'] >= dict[f]] ;
    dffilt = dffilt.append(dftemp)
    #data = [dfint, dftemp]
    #dffinal= pd.concat(data, axis=1)


# In[53]:


dffilt['IID'].unique()


# In[54]:


# What % of ROH remaining?
a=len(df.index)
b=len(dffilt.index)
print(b/a)


# In[55]:


# What is the number of samples remaining?
len(set(dffilt['IID'])) 


# In[56]:


df = dffilt


# 
# 
#   # Calulate SROH, sum of ROH (Mb)

# In[58]:


# Subset the df, and sum the KB column
dfsubset = df[['FID', "IID", 'KB', 'CHR']]
dfsubset.groupby('FID')['KB'].transform('sum')


# In[59]:


# Add the sum to the dataframe
dfsubset.loc[:,"SROH"]  = dfsubset.groupby(['FID', 'CHR'])['KB'].transform('sum')
dfsubset


# In[60]:


# Now, removing the rows that have same ID, CHR, and SROH value. 
dfsu = dfsubset.drop_duplicates(subset=['FID', 'CHR'], keep='last').drop('KB', axis=1)
dfsu


# In[61]:


# Pivot the dataset to have one 
dffinal = dfsu.pivot(index='CHR', columns='FID')
dffinal


# In[62]:


dfna = dffinal.fillna(0)


# In[63]:


sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)


# In[64]:


dfg = dfr.droplevel(0, axis=0)


# In[65]:


dfh = dfg.rename({0: 'SROH'}, axis='columns')


# In[66]:


#baits_size = 1988179


# In[67]:


dfsumROH = pd.DataFrame(dfh['SROH'])
dfsumROH


# # Calculate the NROH (total number of ROH)

# In[68]:


# Add a count column to the dataframe

cols=['FID', 'CHR']

df['count']  = df.groupby(cols)['FID'].transform('count')
df.sort_values(by='count')


# In[70]:


dfsubset = df[["FID", "IID", "CHR", 'KB', 'count']]
dfsubset.sort_values(by='count')


# In[71]:


# Remove the duplicates from the dataframe
dfd = dfsubset.drop_duplicates(['FID', 'CHR']).sort_values(by='count')
dfkl = dfd.drop('KB', axis=1)


# In[72]:


# transpose the dataframe
dffinal = dfkl.pivot(index='CHR', columns='FID')


# In[73]:


dffinal


# In[74]:


# Remove NA in table
dfna = dffinal.fillna(0)


# In[75]:


# check level of indexing
nblevels = dfna.index.nlevels


# In[76]:


# sum of the dataframe for each species
sumroh = dfna.sum(axis=0)
dfr = pd.DataFrame(sumroh)


# In[77]:


# Remove one level of the index
dfg = dfr.droplevel(0, axis=0)


# In[78]:


# Rename the ROH column
dfNROH = dfg.rename({0: 'NROH'}, axis='columns')


# In[79]:


dfNROH


# In[80]:


dfNROH.append(dfsumROH)


# In[81]:


dfNROH['SROH'] = dfsumROH
dfNROH


# In[82]:


# Add FROH, the FROH isthe fraction of each genome in ROH > 0.5 Mb => ACTUALLY no ROH > 0.2 Mb!!!! We try this time without filtering.
# No filtering: all data included.
# In the stats_genome_file.txt
# Begonia_peltatifolia_scaffold	310462054  = 310,462,054 = 310462 Kb

# Input genome length (in Kb):
GEN_LENGTH = 1182020

dfNROH['FROH'] = dfNROH['SROH'] / GEN_LENGTH


# In[83]:


dfNROH.to_csv('matrix_ROH.csv')


# In[87]:


dfNROH['Population'] = dfNROH.index.str[:-2]
dfNROH['Species'] = dfNROH.index.str[:-3]


# In[88]:


dfNROH


# In[89]:


#dfNROH.to_csv('matrix_plot.csv')


# In[90]:


dfNROH = pd.read_csv('matrix_plot.csv')
dfNROH


# ### Full samples, hue by sections:

# In[91]:


# Analysis modules
import numpy as np
import pandas as pd
import matplotlib as plt
import seaborn as sns
import sys
import scipy
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

# Analysis modules
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import scipy

# https://stackoverflow.com/questions/34693991/repel-annotations-in-matplotlib
from matplotlib import pyplot as plt
from adjustText import adjust_text


# In[92]:


df= dfNROH

# Draw the scatterplot:
plt.figure(figsize=(14,8))
sns.set_theme(style="ticks")
ax = sns.scatterplot(x='SROH', y='NROH', data=df, alpha=0.8, hue='Population', palette='colorblind', style='Species', edgecolor= 'black', s=70)
ax.set_xlabel("SROH (kb)")
#ax.set_title("Estimator of Run Of Homozygosity for all Coelocentrum specimens")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


dflabel = df[['NROH', 'SROH']]
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
#for k, v in dflabel.iterrows():
#    ax.annotate(k, v, xytext=(0, 0), textcoords='offset points', family='sans-serif', fontsize=12, color='darkslategrey')

x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()

#texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))


fig = ax.get_figure()
fig.savefig('./plots/plot_SROH_NROH.png', dpi=400, bbox_inches='tight')


# In[93]:


df= dfNROH

# Draw the scatterplot:
plt.figure(figsize=(15,8))
sns.set_theme(style="ticks")
ax = sns.scatterplot(x='SROH', y='NROH', data=df, alpha=1)
ax.set_xlabel("SROH (kb)")
#ax.set_title("Estimator of Run Of Homozygosity for all Coelocentrum specimens")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


dflabel = df[['NROH', 'SROH']]
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
#for k, v in dflabel.iterrows():
#    ax.annotate(k, v, xytext=(0, 0), textcoords='offset points', family='sans-serif', fontsize=12, color='darkslategrey')

x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()

texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))


fig = ax.get_figure()
fig.savefig('./plots/plot_SROH_NROH_annots.png', dpi=400, bbox_inches='tight')


# In[94]:


df= dfNROH

# Draw the scatterplot:
plt.figure(figsize=(15,8))
sns.set_theme(style="ticks")
ax = sns.scatterplot(x='SROH', y='NROH', data=df, alpha=1)
ax.set_xlabel("SROH (kb)")
#ax.set_title("Estimator of Run Of Homozygosity for all Coelocentrum specimens")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


dflabel = df[['NROH', 'SROH']]
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
#for k, v in dflabel.iterrows():
#    ax.annotate(k, v, xytext=(0, 0), textcoords='offset points', family='sans-serif', fontsize=12, color='darkslategrey')

x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()

texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))


fig = ax.get_figure()
fig.savefig('./plots/plot_SROH_NROH_annots.png', dpi=400, bbox_inches='tight')


# In[98]:


plt.figure(figsize=(15,8))
sns.set_theme(style='ticks')

dfplot = df

dfplot['IID'] = df.index

dfpl = dfplot.sort_values(by='FROH')




ax = sns.barplot(x=dfpl['Population'], y=dfpl['FROH'], palette='colorblind')
plt.xticks(rotation=90, ha='right');

fig = ax.get_figure()
fig.savefig('./plots/plot_FROH.png', dpi=400, bbox_inches='tight')


# In[ ]:


dfpl


# In[ ]:


# Use the key to put sample ID in dfroh:
dfROH


# In[ ]:


# Modified dataframe with FST:

df = pd.read_csv('matrix_ROH_FIS.csv', index_col=0)
df


# In[ ]:


# Draw the scatterplot:
plt.figure(figsize=(15,8))
sns.set_theme(style="ticks")
ax = sns.scatterplot(x='SROH', y='NROH', data=df, alpha=1, hue='FIS', s=500)
ax.set_xlabel("SROH (kb)")
#ax.set_title("Estimator of Run Of Homozygosity for all Coelocentrum specimens")
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


dflabel = df[['NROH', 'SROH']]
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
#for k, v in dflabel.iterrows():
#    ax.annotate(k, v, xytext=(0, 0), textcoords='offset points', family='sans-serif', fontsize=12, color='darkslategrey')

x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()

texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))


fig = ax.get_figure()
fig.savefig('./plots/plot_SROH_NROH_annots_fIS.png', dpi=400, bbox_inches='tight')


# In[ ]:




