#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import os


# In[13]:


df = pd.read_csv(snakemake.input[0], sep='\t')
df = df.dropna()

path = snakemake.params.path
threads = snakemake.threads
# In[11]:


import os
import shutil
import concurrent.futures

# Get the list of files to copy
files = [f'genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df.uniprot.unique()]

# Create the destination directory
destination_dir = os.path.join(os.getcwd(), path)
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




