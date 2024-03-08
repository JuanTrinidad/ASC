# %%
import pandas as pd
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#pd.set_option('display.max_columns', 500)

import shutil
import multiprocessing

import concurrent.futures
import requests

# %%
file1 = snakemake.input.file1 #'../../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv'

outputfilepath = snakemake.output.tsv_to_match_pdbs_names

# destination directory
destination_dir = snakemake.params.destination_dir
# cluster_structure_representers path
cluster_str_rep_path = 'genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers/'
#path to model unziped
path_to_model_unzip = 'genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/'
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

#avoid downloading cif files because they are not compatible with fatcat

structures_to_download = [ 'https://alphafold.ebi.ac.uk/files/' + structure.replace('.gz','').replace('.cif', '.pdb') for structure in df_rbh_top_hits.target.unique()]

# %%
#this will load all result in the ram and after that write it in a file. For large dataset can be a problem. 
print('Downloading from AFDB structures:')

def download_files(links, num_cores):
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
        # Submit download tasks to the executor
        futures = [executor.submit(requests.get, link) for link in links]
        
        # Wait for the tasks to complete and retrieve the responses
        responses = [future.result() for future in concurrent.futures.as_completed(futures)]
        
        # Save the downloaded files
        for i, response in enumerate(responses):
            if response.status_code == 200:
                # Extract the filename from the link
                filename = destination_dir + links[i].split('/')[-1].split('-')[1] + '.pdb'
                with open(filename, 'wb') as file:
                    file.write(response.content)
                    print(f'Downloaded {filename}')
            else:
                print(f'Failed to download file_{i}.txt')



download_files(structures_to_download, num_cores)


# %%
def creating_file_OLD_NEW_names(df, outputfilepath):
    df['new_simple_name'] = df['new_simple_name'].str.split('.', expand=True)[0]
    df.sort_values('new_simple_name', inplace=True)

    df[['query','target', 'new_simple_name', 'target_uniprot_accession']].to_csv(outputfilepath, sep='\t', index=False)
    
    
creating_file_OLD_NEW_names(df_rbh_top_hits, outputfilepath)

# %%



