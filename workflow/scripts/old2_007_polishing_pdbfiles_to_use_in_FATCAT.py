# %%
import pandas as pd
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
#pd.set_option('display.max_columns', 500)

# %%
file1 = snakemake.input.file1 #'../report/TriTrypDB-65_All_species_clean_Ortholog_Group_Full_of_hypotetical_boolean.tsv'
file2 = snakemake.input.file2 #'../../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv'
file3 = snakemake.input.file3 #'../../results/Gene_annotation_info_from_uniprot_model_spp.tsv'
file4 = snakemake.input.file4 #'../report/TriTrypDB-65_All_species_clean_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
outputfilepath2 = snakemake.output.tsv_to_match_pdbs_names

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

# %% [markdown]
# file1 = '../report/TriTrypDB-65_All_species_clean_Ortholog_Group_Full_of_hypotetical_boolean.tsv'
# file2 = '../../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv'
# file3 = '../../results/Gene_annotation_info_from_uniprot_model_spp.tsv'
# file4 = '../report/TriTrypDB-65_All_species_clean_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
# 
# outputfilepath2 = '../tmp/quert_taget_accesion_to_fatcat_list.tsv'
# 
# # destination directory
# destination_dir = '../tmp/FATCAT_pdb_files/'
# # cluster_structure_representers path
# cluster_str_rep_path = '../genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers/'
# #path to model unziped
# path_to_model_unzip = '../genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/'
# #cores
# num_cores = 40
# #top hits
# top_hits = 5
# #prefix added pdbs
# prefix = 'wheelerlab_'

# %%
#hypothetical OG data
df_hypothetical_OG = pd.read_csv(file1, sep='\t')

#RBH results
df_rbh = pd.read_csv(file2, sep='\t', index_col=0)

#model organims annotation data
df_MO_annotation = pd.read_csv(file3, sep='\t', low_memory=False, index_col=0)


#estructure to genID to OG

df_OG_genID_uniprot  = pd.read_csv(file4, sep='\t')




#this are the db that i want to use as reference of annotation plus some usefull info. But removing not relevant columns for this analysis.
columns_for_df =  ['Entry','Entry Name','Protein names', 'Organism', 'Organism (ID)', 'Proteomes','Taxonomic lineage' ,'EC number', 'Protein families', 'OrthoDB','PANTHER', 'InterPro', 'Pfam', 'eggNOG', 'Gene Ontology (cellular component)'] #'KEGG'


#removing low annotation proteins from model org
df_MO_annotation = df_MO_annotation[columns_for_df].dropna(thresh=10)

#keeping only relevant columns
#df_kineto_annotation = df_kineto_annotation[columns_for_df]

# %%

#aca uso inner porque quiero evitar las proteinas sin buena anotacion en el resultado
df_rbh_km = df_rbh.merge(df_MO_annotation, left_on='target_uniprot_accession', right_on='Entry', how='left', suffixes=['_kineto', '_model'])
#print(df_rbh_km.shape)


# adding cluster representer
df_rbh_km = df_rbh_km.merge(df_OG_genID_uniprot, left_on='query_uniprot_accession', right_on='uniprot', how='left')
#print(df_rbh_km.shape)

df_rbh_km = df_rbh_km.merge(df_hypothetical_OG, left_on='Ortholog_Group', right_on='cluster_representer_or_OG', how='left')
#print(df_rbh_km.shape)


# %%


# %%


# %%
df_rbh_km = (
    df_rbh_km
    .sort_values('target', ascending=False)
    .drop_duplicates(subset=['query_uniprot_accession', 'target_uniprot_accession']
                    )
)

# %%


# %%
df_rbh_km_tophits =(
    df_rbh_km
    .sort_values(['evalue'])
    .groupby('query')
    .head(top_hits)
)

# %%


# %%


# %%
#simplifiying names of pdb files
df_rbh_km_tophits['new_simple_name'] =  (

    df_rbh_km_tophits['query_uniprot_accession']
    .str.replace(prefix,'')
    .str.split('_', expand=True)[0]
    .str.replace('.','')
    + '.pdb'
)



#dropping duplicates
df_rbh_km_tophits_to_move_files = df_rbh_km_tophits[['query', 'query_uniprot_accession', 'new_simple_name']].drop_duplicates()
#adding path to file
df_rbh_km_tophits_to_move_files['path_to_file'] = cluster_str_rep_path + df_rbh_km_tophits_to_move_files['query']
#creating the list of tuples
files_list = list(zip(df_rbh_km_tophits_to_move_files.path_to_file, df_rbh_km_tophits_to_move_files.new_simple_name))

# %% [markdown]
# #### Funtion to copie pdb files from cluster structure representer to use FATCAT. Also name simplified.

# %%
import os
import shutil
import multiprocessing

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

def copy_files_wo_rm(files, destination_dir, num_cores=4):
    """Copy a list of files to a destination directory in parallel."""
    with multiprocessing.Pool(num_cores) as pool:
        pool.starmap(copy_file, [(file, destination_dir) for file in files])



# copy the files to the destination directory in parallel
copy_files(files_list, destination_dir, num_cores= num_cores)


# %% [markdown]
# #### same approach to model org files

# %%
def filter_df(df, proteome, path):
    
    df = df[df.proteome == proteome]
    df = df.drop_duplicates(subset=['target'])
    df['path_to_file'] = path + df['target']
    df['new_simple_name'] = df.target.str.split('-',expand=True)[1]  + '.' +  df.target.str.split(".", 1, expand=True)[1]

    return list(zip(df.path_to_file, df.new_simple_name))

# %%
for proteome in os.listdir(path_to_model_unzip):
    UPproteome = proteome.split('_')[0]
    
    if UPproteome in df_rbh_km_tophits['proteome'].unique():
        print(f'Coping files from {proteome} to FACTCAT_folder')

        path = path_to_model_unzip + proteome + '/'
        
        files_list = filter_df(df_rbh_km_tophits, UPproteome, path)
        
        copy_files_wo_rm(files_list, destination_dir, num_cores= num_cores)

            



# %% [markdown]
# #### Final output list of files to compare with FATCAT

# %%
def creating_file_for_FATCATQue(df, outputfilepath2):
    df['new_simple_name'] = df['new_simple_name'].str.split('.', expand=True)[0]
    df.sort_values('new_simple_name', inplace=True)

    
    df[['query','target', 'new_simple_name', 'target_uniprot_accession']].to_csv(outputfilepath2, sep='\t', index=False)


creating_file_for_FATCATQue(df_rbh_km_tophits, outputfilepath2)

# %% [markdown]
# 


