#!/usr/bin/env python
# coding: utf-8

# In[48]:


import pandas as pd
import argparse


# In[ ]:



parser = argparse.ArgumentParser()

parser.add_argument('--file1', type=argparse.FileType('r'), help='Path to input file')
parser.add_argument('--file2', type=argparse.FileType('r'), help='Path to input file')
parser.add_argument('--file3', type=argparse.FileType('r'), help='Path to input file')
parser.add_argument('--initialfasta', type=argparse.FileType('r'), help='Path to input fasta file')
parser.add_argument('--outputfasta', type=argparse.FileType('w'), help='Path to input fasta file')
parser.add_argument('--output', type=argparse.FileType('w'), help='Path to output file')
parser.add_argument('--outputstats', type=argparse.FileType('w'), help='Path to output file')

args = parser.parse_args()



# In[ ]:



file2 = args.file1.name
file3 = args.file2.name
file4 = args.file3.name

file_fasta = args.initialfasta.name
file_fasta_output = args.outputfasta.name


output = args.output.name
output_info = args.outputstats.name


# In[ ]:





# 
# file2 = '../report_files/protein_structure_pLDDT_mean.tsv'
# file3 = '../mandatory_files/fasta_header_to_uniprot.tsv'
# file4 = '../mandatory_files/Ortholog_group_to_geneID.tsv'
# 
# file_fasta = '../genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa'
# file_fasta_output = '../genome_data_sets/query_proteomes/fasta_files/orthoMCL_groups_to_be_modelated.fa'
# 
# output_info = '../report_files/ortholog_groups_structure_stats.tsv'
# output = '../report_files/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'

# In[67]:


df_pLDDT = pd.read_csv(file2, 
                       sep='\t', 
                       header=None, 
                       names = ['uniprot', 'pLDDT_mean'])

#print(df_pLDDT.shape)


# In[68]:


df_geneId_uniprot = pd.read_csv(file3, 
                       sep='\t', 
                       names = ['Gene ID', 'uniprot'])

#print(df_geneId_uniprot.shape)


# In[69]:


df_orthologG = pd.read_csv(file4, 
                       sep='\t', 
                       names = ['Ortholog_Group', 'Gene ID'])

#print(df_orthologG.shape)
#print(df_orthologG.Ortholog_Group.nunique())


# In[ ]:





# In[ ]:





# In[70]:


# uno los DF para agregar a los cluster el codigo uniprot y con eso agregar el valor promedio de pLDDT
df_merged = (df_geneId_uniprot
             .merge(df_pLDDT, left_on='uniprot', right_on='uniprot', how='left')
)
#print(df_merged.shape)


# In[ ]:





# In[76]:


df_orthologG_plus_structure = df_orthologG.merge(df_merged, left_on='Gene ID', right_on='Gene ID', how='left')
df_orthologG_plus_structure = df_orthologG_plus_structure.drop_duplicates(keep='first')


# In[ ]:





# In[ ]:





# In[83]:


print('Total number of otholog groups:', df_orthologG_plus_structure.Ortholog_Group.nunique())
print('The number of otholog groups with protein structure assigned are:', df_orthologG_plus_structure.dropna(subset=['uniprot']).Ortholog_Group.nunique())
print('Fraction:',  round(df_orthologG_plus_structure.dropna(subset=['uniprot']).Ortholog_Group.nunique() / df_orthologG.Ortholog_Group.nunique(), 3) )


# In[ ]:





# In[ ]:





# In[84]:


#csv output info
(df_orthologG_plus_structure
 .groupby('Ortholog_Group')['pLDDT_mean']
 .agg(['mean', 'median', 'min', 'max', 'std', 'skew'])
 .to_csv(output_info, sep='\t')
)


# In[85]:



#csv output


df_orthologG_plus_structure = (

    df_orthologG_plus_structure
    .sort_values('pLDDT_mean', ascending=False)
    .drop_duplicates('Ortholog_Group', keep='first')
)


df_orthologG_plus_structure.to_csv(output, sep='\t',index=False)


# # Creating fasta file for no modelated proteins 

# In[86]:


from Bio import SeqIO


# In[87]:


to_be_modelated = df_orthologG_plus_structure[df_orthologG_plus_structure['uniprot'].isna()]['Gene ID'].unique()


# In[89]:


with open(file_fasta, "r") as input_file:
    
    with open(file_fasta_output, "w") as output_file:
        
        for record in SeqIO.parse(input_file, "fasta"):
            
            if record.id in to_be_modelated:
                
                SeqIO.write(record, output_file, "fasta")


# In[ ]:




