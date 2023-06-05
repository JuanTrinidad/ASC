#!/usr/bin/env python
# coding: utf-8

# In[78]:


import pandas as pd 
from Bio.PDB import *
import argparse
import glob


# In[79]:


def Average(lst):
    return round(sum(lst) / len(lst))


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--output', type=argparse.FileType('w'), help='Path to first input file')
args = parser.parse_args()

output = args.output.name


# In[80]:


print('Calculating pLDDT mean of protein structures')


# In[81]:


df_UNIPROT = pd.read_csv('mandatory_files/fasta_header_to_uniprot.tsv', sep='\t', header=None, names=['GeneID', 'UNIPROT'])

lista_archivos = [f'genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df_UNIPROT.UNIPROT.unique()]


#output = '../report_files/protein_structure_pLDDT_mean.tsv'


# In[ ]:





# In[73]:


PDBparser = PDBParser(PERMISSIVE = True, QUIET = True)

#dict to storage data
uniprot_pLDDT = {}

for file in lista_archivos:
    
    uniprot = file.split('/')[-1].split('-')[1]
    
    data = PDBparser.get_structure(uniprot, file)
    
    bfactor = []
    
    for atom in data.get_atoms():
    
        bfactor.append(atom.get_bfactor())
    
    uniprot_pLDDT[uniprot] = Average(bfactor)


# In[75]:


df_output = pd.DataFrame.from_dict(uniprot_pLDDT, orient='index')


df_output.to_csv(output, sep='\t', header=None)


# In[ ]:




