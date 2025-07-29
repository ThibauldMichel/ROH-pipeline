#!/usr/bin/env python
# coding: utf-8

# In[79]:


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


# # Download dataset:

# In[81]:


# Load dataset:
df = pd.read_csv('plink/all.hom', delim_whitespace=True)



# In[82]:


# Load the config file to replace samples names:
nm = pd.read_csv('config/table_reads.tsv', delim_whitespace=True)
dict_name = nm.set_index('read_ID')['Sample'].to_dict()


# In[83]:


# replace samples name in plink dataframe:
df['IID'] = df['IID'].replace(dict_name)



# # Density of ROH size:

# In[84]:


# List of the specimens:
list_iid=df['IID'].drop_duplicates().to_list()


# In[85]:


# Subset the plink output:
dft=df[['IID', 'KB']]
dfe=dft.set_index('IID')


# In[86]:


# Separate values per samples:
pivot_table = df.pivot(columns='IID', values='KB')


# In[87]:


# make a x limit for the plot:
x_lim = pivot_table.mean().max()


# In[88]:


# Make an histogram of the distribution of ROH related to their size in kb
# Draw the scatterplot:
plt.figure(figsize=(8,4))
sns.set_theme(style="ticks")

ax=pivot_table.plot.hist(bins=8000, histtype='step')

ax.set_xlabel("Length ROH (kb)")
ax.set_ylabel("Number of ROH")

ax.set_xlim(0, x_lim)

fig=ax.get_figure()
fig.savefig('plots/plot_ROH_length_hist.png', dpi=400, bbox_inches='tight')


# In[89]:


# Density plot:
# Draw the scatterplot:
plt.figure(figsize=(8,4))
sns.set_theme(style="ticks")


ax = pivot_table.plot.kde()

ax.set_xlim(0, x_lim)
ax.set_xlabel("Length ROH (kb)")

fig=ax.get_figure()
fig.savefig('plots/plot_ROH_length_density.png') 


# # Density number of SNPs:

# In[90]:


# Subset the plink output:
dft=df[['IID', 'NSNP']]
dfe=dft.set_index('IID')

# Separate values per samples:
pivot_table = df.pivot(columns='IID', values='NSNP')

# make a x limit for the plot:
x_lim = pivot_table.mean().max()


# In[91]:


# make a x limit for the plot:
x_lim = pivot_table.mean().max()


# In[92]:


# Make an histogram of the distribution of ROH related to their size in kb
# Draw the scatterplot:
plt.figure(figsize=(8,4))
sns.set_theme(style="ticks")

ax=pivot_table.plot.hist(bins=500, histtype='step')

ax.set_xlabel("Number of SNPs per ROH")
ax.set_ylabel("Number of ROH")

ax.set_xlim(0, x_lim)

fig=ax.get_figure()
fig.savefig('plots/plot_ROH_SNPs_hist.png', dpi=400, bbox_inches='tight')


# In[93]:


# Density plot:
# Draw the scatterplot:
plt.figure(figsize=(8,4))
sns.set_theme(style="ticks")


ax = pivot_table.plot.kde()

ax.set_xlim(0, x_lim)
ax.set_xlabel("Number SNPs per ROH")

fig=ax.get_figure()
fig.savefig('plots/plot_ROH_SNPs_density.png') 

