#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os
import argparse


# In[ ]:





# In[ ]:


parser = argparse.ArgumentParser()

parser.add_argument('--file1', type=argparse.FileType('r'), help='Path to first input file')
parser.add_argument('--threads', type=int, help='Number of threads INT')

args = parser.parse_args()



file1 = args.file1.name
threads = args.threads


# file1 = '../report_files/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
# threads = 30

# In[ ]:





# In[13]:


df = pd.read_csv(file1, sep='\t')
df = df.dropna()


# In[ ]:





# In[ ]:





# In[11]:


import os
import shutil
import concurrent.futures

# Get the list of files to copy
files = [f'genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df.uniprot.unique()]

# Create the destination directory
destination_dir = os.path.join(os.getcwd(), "genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers")
if not os.path.exists(destination_dir):
    os.mkdir(destination_dir)

# Copy the files in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
    future_to_file = {executor.submit(shutil.copy, file, destination_dir): file for file in files}

# Wait for all the copies to finish
for future in concurrent.futures.as_completed(future_to_file):
    future.result()

# Print a message to indicate that the copies are finished
print("All files have been copied.")


# In[ ]:




