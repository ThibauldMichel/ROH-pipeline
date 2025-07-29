#!/usr/bin/env python
# coding: utf-8

# In[69]:


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
from matplotlib import pyplot
# https://stackoverflow.com/questions/34693991/repel-annotations-in-matplotlib
from adjustText import adjust_text




# In[70]:


df = pd.read_csv('output/matrix_ROH.csv', index_col=0)


# In[71]:


# Load the config file to replace samples names:
nm = pd.read_csv('config/table_reads.tsv', delim_whitespace=True)
dict_name = nm.set_index('read_ID')['Sample'].to_dict()

# replace samples name in plink dataframe:
df['IID']=df.index
df['IID'] = df['IID'].replace(dict_name)
dfi=df.set_index('IID')
df = dfi


# In[72]:


# creating X-Y Plots With a Regression Line
sns.set_theme(style="ticks")
sns.color_palette('colorblind')

# slope, intersept, and correlation coefficient calculation 
slope, intercept, r, p, stderr = scipy.stats.linregress(df['FROH'], df['Heterozygosity'])
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'

# plotting
fig, ax = pyplot.subplots()
ax.plot(df['FROH'], df['Heterozygosity'], linewidth=0, marker='o', alpha=0.7)
ax.plot(df['FROH'], intercept + slope * df['FROH'], label=line, color='orange')
ax.set_xlabel('FROH')
ax.set_ylabel('Heterozygosity_ratio')
ax.legend(facecolor='white')
pyplot.show()


# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
#dflabel = df
#x = dflabel['SROH'].to_list()
#y = dflabel['NROH'].to_list()
#s = dflabel.index.to_list()
#texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
#adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))


fig = ax.get_figure()
fig.savefig('plots/plot_heterozygosity_FROH_regline.png', dpi=400, bbox_inches='tight')



# In[73]:


# Draw the scatterplot NROH, SROH, FROH:

# Seaborn parameters:
sns.set_theme(style="ticks")
sns.color_palette("bright")

# The plot:
ax = sns.scatterplot(x='SROH',
                     y='NROH', 
                     hue='FROH', 
                     data=df, 
                     alpha=0.8,
                     palette='colorblind',
                     edgecolor= 'black')
ax.set_xlabel("SROH (kb)")
ax.legend()
ax.legend(title='FROH' ,loc='center left', bbox_to_anchor=(1, 0.5),
         labels=df['FROH'].round(3))

# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
dflabel = df
x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()
texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))

# Save figure:
ax.get_figure()
fig = ax.get_figure()
fig.savefig('plots/plot_SROH_NROH_FROH.png', dpi=400, bbox_inches='tight')


# In[74]:


# Draw the scatterplot NROH, SROH, Populations:

# Seaborn parameters:
sns.set_theme(style="ticks")
sns.color_palette("bright")

# The plot:
ax = sns.scatterplot(x='SROH',
                     y='NROH', 
                     hue='Population', 
                     data=df, 
                     alpha=0.8,
                     palette='colorblind',
                     edgecolor= 'black')
ax.set_xlabel("SROH (kb)")
ax.legend()
ax.legend(title='Populations' ,loc='center left', bbox_to_anchor=(1, 0.5))

# Input the annotation:
# https://stackoverflow.com/questions/15910019/annotate-data-points-while-plotting-from-pandas-dataframe
dflabel = df
x = dflabel['SROH'].to_list()
y = dflabel['NROH'].to_list()
s = dflabel.index.to_list()
texts = [plt.text(x[i], y[i], s[i]) for i in range(len(x))]
adjust_text(texts, arrowprops=dict(arrowstyle='-', color='grey'))

# Save figure:
ax.get_figure()
fig = ax.get_figure()
fig.savefig('plots/plot_SROH_NROH_Population.png', dpi=400, bbox_inches='tight')


# In[ ]:




