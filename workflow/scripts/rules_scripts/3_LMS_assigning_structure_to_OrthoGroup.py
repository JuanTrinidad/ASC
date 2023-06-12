#!/usr/bin/env python
# coding: utf-8

# In[28]:


import pandas as pd




pLDDTmean_file = snakemake.input.file1
orthoG_to_geneID = snakemake.input.file2
output =snakemake.output[0]



df = pd.read_csv(pLDDTmean_file, sep='\t', names=['geneID','pLDDT_mean','file_name'])

df_og = pd.read_csv(orthoG_to_geneID, sep='\t', names=['Ortholog_Group','geneID'])

df_merge = df_og.merge(df, how='left')

df_merge_best_structure = df_merge.sort_values(['pLDDT_mean', 'file_name'], ascending=[False, True]).drop_duplicates(subset=['Ortholog_Group'], keep='first')

df_merge_best_structure.to_csv(output, sep='\t', index=False)






