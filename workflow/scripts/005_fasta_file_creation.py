#!/usr/bin/env python
# coding: utf-8

# In[33]:


import pandas as pd
from Bio import SeqIO
import os


# In[ ]:





# In[24]:


#funtion to search into fasta file and create cluster fasta files

def looking_in_fasta_file(name_geid_list, dest_dir, name, original_fasta_file):

    with open(original_fasta_file, "r") as input_file:
        
        #dest_dir = '../../tmp/report_files/sequences_to_be_modelated/'
        
        
        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)
            
        with open(name, "w") as output_file:

            for record in SeqIO.parse(input_file, "fasta"):

                if record.id in name_geid_list:

                    SeqIO.write(record, output_file, "fasta")


# ## file import 

# fasta_file_input = '../genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa'
# report_ortho_g = '../report_files/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
# ortho_info = '../mandatory_files/Ortholog_group_to_geneID.tsv'
# output_fasta_file = 'algo.fasta'
# ortho_group_size = 1

# In[ ]:



fasta_file_input = snakemake.input.original_fasta
report_ortho_g = snakemake.input.report_ortho_g
ortho_info = snakemake.input.ortho_info


output_fasta_file = snakemake.output.output_fasta_file
ortho_group_size = snakemake.params.OG_size


# ## creating df 

# In[26]:


df_OGSCU = pd.read_csv(report_ortho_g, 
            sep='\t')


# In[27]:


df_Oinfo = pd.read_csv(ortho_info, 
            sep='\t', 
            header=None, 
            names= ['OG','geneID'])


# In[28]:


#clusters without uniprot accession in the report of downloaded structures from AFDB
OG_nan = df_OGSCU[df_OGSCU['uniprot'].isna()].Ortholog_Group.unique()


# In[29]:


#from those count amount of times as proxi of cluster members
df_count = df_Oinfo[df_Oinfo.OG.isin(OG_nan)].value_counts('OG')


# In[30]:


#selecting the cluster member size to filter before creting fasta file
df_count_gt = df_count[df_count.gt(ortho_group_size)]


# In[31]:


#creating the list of gene ID to extract from fasta
df_Oinfo_cluster_of_interest = df_Oinfo[df_Oinfo.OG.isin(df_count_gt.index)]['geneID'].unique()


# In[34]:



#creating fasta file
looking_in_fasta_file(df_Oinfo_cluster_of_interest, 'report_files/fasta_files/', output_fasta_file, fasta_file_input)


# In[ ]:




