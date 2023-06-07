#!/usr/bin/env python
# coding: utf-8

# In[83]:


import pandas as pd


# In[ ]:





# 
# file1 = '../report_files/protein_structure_pLDDT_mean.tsv'
# file2 = '../mandatory_files/fasta_header_to_uniprot.tsv'
# file3 = '../mandatory_files/Ortholog_group_to_geneID.tsv'
# 
# #file_fasta = '../genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa'
# #file_fasta_output = '../genome_data_sets/query_proteomes/fasta_files/orthoMCL_groups_to_be_modelated.fa'
# 
# output_info = '../report_files/ortholog_groups_structure_stats.tsv'
# output = '../report_files/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'

# In[85]:


df_pLDDT = pd.read_csv(snakemake.input.file1, 
                       sep='\t', 
                       header=None, 
                       names = ['uniprot', 'pLDDT_mean'])

#print(df_pLDDT.shape)


# In[86]:


df_geneId_uniprot = pd.read_csv(snakemake.input.file2, 
                       sep='\t', 
                       names = ['Gene ID', 'uniprot'])

#print(df_geneId_uniprot.shape)


# In[87]:


df_orthologG = pd.read_csv(snakemake.input.file3, 
                       sep='\t', 
                       names = ['Ortholog_Group', 'Gene ID'])

#print(df_orthologG.shape)
#print(df_orthologG.Ortholog_Group.nunique())


# In[ ]:





# In[ ]:





# In[88]:


# uno los DF para agregar a los cluster el codigo uniprot y con eso agregar el valor promedio de pLDDT
df_merged = (df_geneId_uniprot
             .merge(df_pLDDT, left_on='uniprot', right_on='uniprot', how='left')
)
#print(df_merged.shape)


# In[ ]:





# In[89]:


df_orthologG_plus_structure = df_orthologG.merge(df_merged, left_on='Gene ID', right_on='Gene ID', how='left')
df_orthologG_plus_structure = df_orthologG_plus_structure.drop_duplicates(keep='first')


# In[ ]:





# In[ ]:





# In[90]:


print('Total number of otholog groups:', df_orthologG_plus_structure.Ortholog_Group.nunique())
print('Amount of otholog groups with protein structure assigned:', df_orthologG_plus_structure.dropna(subset=['uniprot']).Ortholog_Group.nunique())
print('Fraction:',  round(df_orthologG_plus_structure.dropna(subset=['uniprot']).Ortholog_Group.nunique() / df_orthologG.Ortholog_Group.nunique(), 3) )


# In[ ]:





# In[ ]:





# In[91]:


#csv output info
(df_orthologG_plus_structure
 .groupby('Ortholog_Group')['pLDDT_mean']
 .agg(['mean', 'median', 'min', 'max', 'std', 'skew'])
 .to_csv(snakemake.output.outstats, sep='\t')
)


# In[92]:



#csv output


df_orthologG_plus_structure = (

    df_orthologG_plus_structure
    .sort_values('pLDDT_mean', ascending=False)
    .drop_duplicates('Ortholog_Group', keep='first')
)


df_orthologG_plus_structure.to_csv(snakemake.output.out, sep='\t',index=False)


# In[ ]:





# In[ ]:




