#!/usr/bin/env python
# coding: utf-8

# In[158]:


import pandas as pd
import os


# In[ ]:





# 
# file1 = '../report/protein_structure_pLDDT_mean.tsv'
# file2 = '../../config/mandatory_files/fasta_header_to_uniprot.tsv'
# file3 = '../../config/mandatory_files/Ortholog_group_to_geneID.tsv'
# 
# 
# #file_fasta = '../genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa'
# #file_fasta_output = '../genome_data_sets/query_proteomes/fasta_files/orthoMCL_groups_to_be_modelated.fa'
# 
# output_info = '../report/ortholog_groups_structure_stats.tsv'
# output = '../report/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'

# In[160]:


df_pLDDT = pd.read_csv(snakemake.input.file1 , #snakemake.input.file1 
                       sep='\t', 
                       header=None, 
                       names = ['uniprot', 'pLDDT_mean'])

#print(df_pLDDT.shape)


# In[161]:


df_geneId_uniprot = pd.read_csv(snakemake.input.file2 , #snakemake.input.file2 
                       sep='\t', 
                       names = ['Gene ID', 'uniprot'])

#print(df_geneId_uniprot.shape)


# In[162]:


df_orthologG = pd.read_csv(snakemake.input.file3, #snakemake.input.file3 
                       sep='\t', 
                       names = ['Ortholog_Group', 'Gene ID'])

#print(df_orthologG.shape)
#print(df_orthologG.Ortholog_Group.nunique())


# In[163]:


# uno los DF para agregar a los cluster el codigo uniprot y con eso agregar el valor promedio de pLDDT
df_merged = df_geneId_uniprot.merge(df_pLDDT, left_on='uniprot', right_on='uniprot', how='left')
#print(df_merged.shape)


# In[167]:


#corrected structure assigment

file4 = snakemake.input.file4
if os.path.exists(file4):
    print('PDB files have been checked to ensure same sequence')
    print('-' * 50)
    df_PDBcheking = pd.read_csv(file4, sep='\t')
    df_merged = df_merged.merge(df_PDBcheking, left_on=['Gene ID', 'uniprot'], right_on=['genID', 'UNIPROT'], how='left')
    df_merged['fasta_div_pdb'] = df_merged['len_fasta'] / df_merged['lenPDB']
    df_merged['score_div_fasta'] = df_merged['score'] / df_merged['len_fasta']
    df_merged['score_div_pdb'] = df_merged['score'] / df_merged['lenPDB']
    
    
    df_merged['min_of_both'] = df_merged[['score_div_fasta', 'score_div_pdb']].min(axis=1) 
    
    df_merged = df_merged[((df_merged.fasta_div_pdb >= .9) & (df_merged.fasta_div_pdb <= 1.1)) & ((df_merged.min_of_both >= .9) & (df_merged.min_of_both <= 1))]
    df_merged = df_merged.iloc[:,:3]
else:
    print('PDB files are not been checked to ensure same sequence')


# In[ ]:





# In[ ]:





# In[165]:


df_orthologG_plus_structure = df_orthologG.merge(df_merged, left_on='Gene ID', right_on='Gene ID', how='left')
df_orthologG_plus_structure = df_orthologG_plus_structure.drop_duplicates(keep='first')


# In[ ]:





# In[ ]:





# In[166]:


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




