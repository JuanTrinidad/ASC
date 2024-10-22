#!/usr/bin/env python
# coding: utf-8

# In[2]:


import pandas as pd


# In[20]:


#'../../config/mandatory_files/fasta_header_to_uniprot.tsv'
df_clusters = pd.read_csv(snakemake.input[0], 
                          sep='\t', 
                          names=['gene_name', 'UNIPROT_accession'])



df_clusters = df_clusters.dropna()
df_clusters = df_clusters[~df_clusters['UNIPROT_accession'].str.contains(snakemake.params.prefix_of_added_pdbs)]


# In[21]:


df_mmseq2_relevant_clusters = pd.read_csv(snakemake.input[1],
                                         sep='\t',
                                         names=['cluster_representer', 'cluster_members'])


# In[26]:


#avoid downloading structures from non important clusters
df_clusters = df_clusters[df_clusters.gene_name.isin(df_mmseq2_relevant_clusters.cluster_members)]


# # Download files from AlphaFold DB 

# In[10]:


import os
import requests
from multiprocessing.pool import ThreadPool
from tqdm import tqdm

# Definir una función para descargar un archivo desde una URL y guardarla en un directorio
def download_file(url, dest_dir):
    # Obtener el nombre del archivo a partir de la URL
    file_name = os.path.basename(url)
    # Construir la ruta completa del archivo de destino
    dest_path = os.path.join(dest_dir, file_name)
    
    # Descargar el archivo y escribirlo en disco
    if os.path.exists(dest_path):
        #print(f"Skipping {url} - file already exists in {output_dir}")
        return
    with open(dest_path, "wb") as file:
        response = requests.get(url)
        file.write(response.content)
        
# Lista de enlaces que se van a descargar
links = [f'https://alphafold.ebi.ac.uk/files/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df_clusters.UNIPROT_accession.unique()]

# Directorio de destino para las descargas
dest_dir = "genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB"

# Crear el directorio de destino si no existe
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)
    
# Crear un objeto ThreadPool con un número de hilos adecuado
pool = ThreadPool(snakemake.threads) # aquí se utiliza 4 hilos, pero esto se puede ajustar según el tamaño de su máquina

#
print('Downloading structures from AlphaFold Data Base.')

# Mapa los enlaces a la función de descarga en paralelo
for _ in tqdm(pool.imap_unordered(lambda url: download_file(url, dest_dir), links), total=len(links)):
    pass

# Cerrar el objeto ThreadPool para liberar recursos
pool.close()
pool.join()

print('Download Finished')
        


# In[ ]:




