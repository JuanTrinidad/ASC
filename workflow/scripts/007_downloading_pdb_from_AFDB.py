# %%
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#pd.set_option('display.max_columns', 500)

import shutil
import multiprocessing

import os
import requests
from multiprocessing.pool import ThreadPool
from tqdm import tqdm


# %%
file1 = snakemake.input.file1 

outputfilepath = snakemake.output.tsv_to_match_pdbs_names

#here we will download pdbs from AFDB
first_target_directory = 'tmp/FATCAT_target_pdb_files/'

# destination directory
destination_dir = snakemake.params.destination_dir

# cluster_structure_representers path
cluster_str_rep_path = 'genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers/'

#cores
num_cores = snakemake.threads

#top hits
top_hits = snakemake.params.top_hits

#prefix added pdbs
prefix = snakemake.params.prefix_of_added_pdbs

# %%
#removing and creating destination dir

if os.path.exists(destination_dir):
    shutil.rmtree(destination_dir)
os.makedirs(destination_dir)

# %%
#RBH results
df_rbh = pd.read_csv(file1, sep='\t', index_col=0)

# %%

df_rbh = df_rbh.drop_duplicates(subset=['query_uniprot_accession', 'target_uniprot_accession'])


# %%
df_rbh_top_hits = (
    df_rbh
    .sort_values('evalue')
    .groupby('query')
    .head(top_hits)
    
)

# %%
#simplifiying names of pdb files
df_rbh_top_hits['new_simple_name'] =  (

    df_rbh_top_hits['query_uniprot_accession']
    .str.replace(prefix,'')
    .str.split('_', expand=True)[0]
    .str.replace('.','')
    + '.pdb'
)



#simplifiying names of pdb files
df_rbh_top_hits['new_simple_name_target'] =  (

    df_rbh_top_hits['target_uniprot_accession']
    .str.replace(prefix,'')
    .str.split('_', expand=True)[0]
    .str.replace('.','')
    + '.pdb'
)

# %% [markdown]
# ### creating path to copy query structures and simplify file name

# %%
#dropping duplicates
df_rbh_km_tophits_to_move_files = df_rbh_top_hits[['query', 'query_uniprot_accession', 'new_simple_name']].drop_duplicates()
#adding path to file
df_rbh_km_tophits_to_move_files['path_to_file'] = cluster_str_rep_path + df_rbh_km_tophits_to_move_files['query']
#creating the list of tuples
files_list = list(zip(df_rbh_km_tophits_to_move_files.path_to_file, df_rbh_km_tophits_to_move_files.new_simple_name))

# %%
def copy_file(file, destination_dir):
    """Copy a file to a destination directory."""
    file_path, new_name = file
    shutil.copy(file_path, os.path.join(destination_dir, new_name))


def copy_files(files, destination_dir, num_cores=4):
    """Copy a list of files to a destination directory in parallel."""
    if os.path.exists(destination_dir):
        shutil.rmtree(destination_dir)
    os.makedirs(destination_dir)
    
    with multiprocessing.Pool(num_cores) as pool:
        pool.starmap(copy_file, [(file, destination_dir) for file in files])
        

# %%
# copy the files to the destination directory in parallel
copy_files(files_list, destination_dir, num_cores= num_cores)

# %% [markdown]
# ## Downloading target structures from AFDB

# %%
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
links = [f'https://alphafold.ebi.ac.uk/files/AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df_rbh_top_hits.target_uniprot_accession.unique()]

# Directorio de destino para las descargas
dest_dir = first_target_directory

# Crear el directorio de destino si no existe
if not os.path.exists(dest_dir):
    os.makedirs(dest_dir)
    
# Crear un objeto ThreadPool con un número de hilos adecuado
pool = ThreadPool(num_cores) # aquí se utiliza 4 hilos, pero esto se puede ajustar según el tamaño de su máquina

#
print('Downloading structures from AlphaFold Data Base.')

# Mapa los enlaces a la función de descarga en paralelo
for _ in tqdm(pool.imap_unordered(lambda url: download_file(url, dest_dir), links), total=len(links)):
    pass

# Cerrar el objeto ThreadPool para liberar recursos
pool.close()
pool.join()

print('Download Finished')

# %% [markdown]
# ### move files and delete directory

# %%
def copy_files_wo_rmDir(files, destination_dir, num_cores=4):
    """Copy a list of files to a destination directory in parallel."""
    
    with multiprocessing.Pool(num_cores) as pool:
        pool.starmap(copy_file, [(file, destination_dir) for file in files])

# %%


# List all files in the directory
files = os.listdir(first_target_directory)


target_file_list = [(first_target_directory + file,  file.split('-')[1] + '.pdb') for file in files]


# %%
# copy the files to the destination directory in parallel
copy_files_wo_rmDir(target_file_list, destination_dir, num_cores= num_cores)



# %%
if os.path.exists(first_target_directory):
    shutil.rmtree(first_target_directory)

# %% [markdown]
# ## creating comparison file

# %%
def creating_file_OLD_NEW_names(df, outputfilepath):
    
    df = df.sort_values('query_uniprot_accession')

    df[['query','target', 'query_uniprot_accession', 'target_uniprot_accession']].to_csv(outputfilepath, sep='\t', index=False)
    
    
creating_file_OLD_NEW_names(df_rbh_top_hits, outputfilepath)

# %%



