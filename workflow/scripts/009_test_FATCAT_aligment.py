# %%
import pandas as pd
import subprocess
import multiprocessing
import glob
import os
import shutil

# %%

folder_path = "tmp/FATCAT_aligments/"

# Check if the folder exists
if os.path.exists(folder_path):
    # Remove the folder and its contents
    shutil.rmtree(folder_path)
    
# Create the folder
os.makedirs(folder_path)


# %%
num_cores = snakemake.threads # Specify the number of cores to use

df = pd.read_csv(snakemake.input[0], sep='\t') #'../tmp/TriTrypDB-65_All_species_clean_query_taget_accesion_to_fatcat_list.tsv'

original_shape = df.shape

# %%
#to ensure that we will continue with pdb files that are in the correct folder.
uniprot_names = []
for file in glob.glob('tmp/FATCAT_pdb_files/*.pdb'):

    uniprot = file.split('/')[-1][:-4]

    uniprot_names.append(uniprot)

# %%
#checking if the uniprot names are in the dataframe
df = (
df[
    df['new_simple_name'].isin(uniprot_names) &
    df['target_uniprot_accession'].isin(uniprot_names)
]
)

# %%
#checking if the number of rows is the same as the original dataframe
if original_shape == df.shape:
    print('All of the .pdb files are present in the FATCAT folder.')
else:
    pdb_not_present = original_shape[0] - df.shape[0]
    print(f'{pdb_not_present} .pdb files are not present in the FATCAT folder. Probably their are not in AFDB with the same UNIPROT accession than in the model organisms proteome files. Please check. The script will continue with the files that are present.')

# %%
command_args = []

for index, row in df.iterrows():
    pdb1 = row['target_uniprot_accession']
    pdb2 = row['new_simple_name'] 

    command_args.append(["git_repo_cloned/FATCAT/FATCATMain/FATCAT", 
                         "-p1", f"tmp/FATCAT_pdb_files/{pdb1}.pdb",
                         "-p2", f"tmp/FATCAT_pdb_files/{pdb2}.pdb",
                         "-o", f"tmp/FATCAT_aligments/{pdb1}_{pdb2}", 
                         "-m", "-t"])


# %%

def run_command(args):

    try:
        result = subprocess.run(args, capture_output=True, text=True, check=True)
        
    except:
        pass

            
    #except subprocess.CalledProcessError as e:

    #    if len(e.stdout) != 0:
    #        print(f"No output for {args[2]} vs {args[4]} comparison by FATCAT.\n")
    #        print("Output:")
    #        print(result.stdout)

        

# Run the commands in parallel using multiple cores
with multiprocessing.Pool(num_cores) as pool:
    
    pool.map(run_command, command_args)


# %%


# %% [markdown]
# def run_command(args):
#     result = None  # Initialize result
#     try:
#         result = subprocess.run(args, capture_output=True, text=True, check=True)
#     except subprocess.CalledProcessError as e:
#         print(f"Command failed: {e}")
#         if result:  # Check if result was assigned
#             print(result.stdout)
#     return result
# 
# # Run the commands in parallel using multiple cores
# with multiprocessing.Pool(num_cores) as pool:
#     pool.map(run_command, command_args)


