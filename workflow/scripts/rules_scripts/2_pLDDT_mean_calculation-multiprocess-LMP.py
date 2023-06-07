#!/usr/bin/env python
# coding: utf-8

# In[26]:


import pandas as pd 
from Bio.PDB import *
import glob
import multiprocessing



output = snakemake.output[0]
num_threads = snakemake.threads


path = snakemake.params.path
# output = '../report/probandopLDDT.tsv'
# num_threads = 30

# In[46]:


print('Calculating pLDDT mean of protein structures')


# In[47]:



lista_archivos = glob.glob(f'{path}/*.pdb')


# In[48]:


def Average(lst):
    return round(sum(lst) / len(lst))


# In[49]:




#dict to storage data

def extract_values_from_PDB_files(PDBfile):
    
    PDBparser = PDBParser(PERMISSIVE = True, QUIET = True)
    
    uniprot = PDBfile.split('/')[-1]
    
    uniprot = uniprot[:uniprot.find('_unrelaxed_rank_')]
    
    data = PDBparser.get_structure(uniprot, PDBfile)
    
    bfactor = []
    
    for atom in data.get_atoms():
    
        bfactor.append(atom.get_bfactor())
    
    return uniprot, Average(bfactor)


# In[50]:


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





# In[51]:


if __name__ == "__main__":
  main(num_threads)


# In[ ]:




