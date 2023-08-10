#!/usr/bin/env python
# coding: utf-8

# In[ ]:





# In[99]:


import pandas as pd


# ## funtion for pipe 

# In[100]:



def get_only_uniprot_accession(df, colname):
    
    colname_out = colname + '_uniprot_accession'
    
    df[colname_out] = df[colname].str.split('-', expand=True)[1]
    
    return df


def calc_coverage(df, which_one ,start, end, lenght):
    
    col_name = 'COV_' + which_one
    
    df[col_name] = (df.loc[: , end] - df.loc[: , start]) / df.loc[: , lenght]
    
    df[col_name] = df[col_name].round(2)
    
    return df


# ## concatenating all tsv data frames 

# In[101]:


list_of_df = []
for file in snakemake.input:
    
    df = pd.read_csv(file, sep='\t', names=['query','target','alnlen','fident', 'evalue', 'qstart','qend', 'qlen','tstart','tend','tlen','aligmentinfo'])
    
    #info from file name
    info = file.split('_')[-4:-1]
    df['proteome'] = info[0]
    df['spp'] = info[-1]
    
    #appending_list
    list_of_df.append(df)
    #print(file)


# In[102]:


df = pd.concat(list_of_df)


# In[106]:


df = (df
 .pipe(get_only_uniprot_accession, 'query')
 .pipe(get_only_uniprot_accession, 'target')
 .pipe(calc_coverage, 'query', 'qstart', 'qend','qlen')
 .pipe(calc_coverage, 'target', 'tstart', 'tend','tlen'))


# In[ ]:


df.to_csv(snakemake.output[0], sep='\t')

