# %%
import pandas as pd
import glob
import multiprocessing
import subprocess

path = 'tmp/FATCAT_aligments_twisted/'
files = glob.glob(path + '*.pdb')


# Number of cores to use for parallel execution
num_cores = snakemake.threads #10 #multiprocessing.cpu_count()


# Creating a DataFrame with a column named 'File_Name' containing the file paths
df = pd.DataFrame({'File_Name': files})

# Splitting the file paths by '/' and extracting the last part of the path
df['File_Name_Split'] = df['File_Name'].str.split('/').str[-1]

# Splitting the 'File_Name_Split' column by '_' and extracting the second and third parts
df[['pdb_1', 'pdb_2']] = df['File_Name_Split'].str.split('_', expand=True)[[1,2]]


# %%
# Sorting the dataframe 'df' based on the columns 'pdb_1' and 'pdb_2'
df = df.sort_values(['pdb_1', 'pdb_2'])

# Dropping duplicate rows in the dataframe 'df' based on the columns 'pdb_1' and 'pdb_2',
# and keeping the first occurrence of each unique combination
df_first = df.drop_duplicates(subset=['pdb_1', 'pdb_2'], keep='first')

# Dropping duplicate rows in the dataframe 'df' based on the columns 'pdb_1' and 'pdb_2',
# and keeping the last occurrence of each unique combination
df_last =  df.drop_duplicates(subset=['pdb_1', 'pdb_2'], keep='last')


# %%
# Merging the 'df_first' and 'df_last' dataframes based on the columns 'pdb_1' and 'pdb_2'
merged_df = df_first.merge(df_last, on=['pdb_1', 'pdb_2'])

# Selecting only the 'File_Name_x' and 'File_Name_y' columns from the merged dataframe
merged_df = merged_df[['File_Name_x', 'File_Name_y']]

# %%
list_of_lists = merged_df.values.tolist()

# %%
command_args = []

for pdb1, pdb2 in list_of_lists:

    pdb1_name, pdb2_name = pdb1.split('/')[-1].split('.')[0].split('_')[1:3]

    command_args.append(["TMalign", 
                         f"{pdb1}",
                         f"{pdb2}",
                         f"{path}{pdb1_name}_{pdb2_name}_TMalgin_results.txt"])


# %%
def run_command(command):
    # Extract the file name from the command_args list
    file_name = command[3]
    command = command[:3]

    # Execute the command in the terminal using subprocess
    subprocess.run(command, stdout=open(file_name, 'w'))


with multiprocessing.Pool(num_cores) as pool:
    pool.map(run_command, command_args)



