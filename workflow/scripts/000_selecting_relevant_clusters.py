#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd


# In[ ]:


clust_size = snakemake.params.cluster_size
input_file = snakemake.input[0]
output_file = snakemake.output[0]


# clust_size = 10

# In[18]:


print('Filtering clusters created by MMseq2')
print('-' * 50)
print(f'Removing clusters smaller than {clust_size} members')
print('-' * 50)


# In[15]:


df = pd.read_csv(input_file, sep='\t', header=None)


# In[17]:


print(f'Amount of clusters before filtering: {df[0].nunique()}')


# In[4]:


df0 = df.groupby(0).size()


# In[23]:


df0 = df0[df0 >= clust_size].index


# In[26]:


df = df[df[0].isin(df0)]


# In[28]:


print(f'Amount of clusters after filtering: {df[0].nunique()}')


# In[27]:


df.to_csv(output_file, header=None, index=None, sep='\t')


# In[ ]:




