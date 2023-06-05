#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd 
from Bio.PDB import *
import argparse
import glob
import multiprocessing


# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--output', type=argparse.FileType('w'), help='Path to output file')
parser.add_argument('--threads', type=int, help='Number of threads INT')

args = parser.parse_args()

output = args.output.name
num_threads = args.threads


# In[2]:


print('Calculating pLDDT mean of protein structures')


# In[ ]:


df_UNIPROT = pd.read_csv('mandatory_files/fasta_header_to_uniprot.tsv', sep='\t', header=None, names=['GeneID', 'UNIPROT'])

lista_archivos = [f'genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df_UNIPROT.UNIPROT.unique()]


#output = '../report_files/protein_structure_pLDDT_mean_multi_prueba.tsv'
#num_threads = 20


# In[4]:


def Average(lst):
    return round(sum(lst) / len(lst))


# In[5]:




#dict to storage data

def extract_values_from_PDB_files(PDBfile):
    
    PDBparser = PDBParser(PERMISSIVE = True, QUIET = True)
    
    uniprot = PDBfile.split('/')[-1].split('-')[1]
    
    data = PDBparser.get_structure(uniprot, PDBfile)
    
    bfactor = []
    
    for atom in data.get_atoms():
    
        bfactor.append(atom.get_bfactor())
    
    return uniprot, Average(bfactor)


# In[6]:


def main(num_threads):
    # Get a list of all the file names.
    file_names = lista_archivos

    # Create a pool of processes.
    pool = multiprocessing.Pool(num_threads)

    # Map the extract_values function to each file name in the list.
    results = pool.map(extract_values_from_PDB_files, file_names)

    # Close the pool.
    pool.close()
    pool.join()

    # Create a dictionary with the key and values from each file.
    final_dictionary = {}
    for key, value in results:
        final_dictionary[key] = value

    # return the final dictionary.
    df_output = pd.DataFrame.from_dict(final_dictionary, orient='index')
    df_output.to_csv(output, sep='\t', header=None)
    
    print('All calculations done')


# In[ ]:





# In[ ]:


if __name__ == "__main__":
  main(num_threads)


# In[ ]:




