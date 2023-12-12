# %%
import glob
import pandas as pd
import re

# %%
path = 'tmp/FATCAT_aligments_twisted/'
files = glob.glob(path + '*.txt')

# %%
#extracting relevant information from the files

lines_with_tmscore = []

for file in files:
    file_and_data = [file]

    with open(file, 'r') as f:
        for line in f:
            if line.startswith('TM-score=') | line.startswith('Aligned length=') | line.startswith('Length of Chain'):
                
                file_and_data.append( line.strip() )
        
    lines_with_tmscore.append(file_and_data)



# %%
#creating dataframe with the information
df = pd.DataFrame(lines_with_tmscore, columns=['File', 'Length_Chain1', 'Length_Chain2', 'Aligned length', 'TM-score_Chain1', 'TM-score_Chain2'])

# %%
#extracting columns and converting to numeric
df[['pdb1', 'pdb2']] = (df
 .File
 .str.split('/', expand=True).iloc[:, -1]
 .str.split('_', expand=True).iloc[:,:2]
 )

# %%
df[['Aligned_length', 'RMSD', 'Seq_ID']] = df['Aligned length'].str.split(',', expand=True)

# %%
#regluar expression for each column to avoid problems
#extracting length of chains
df['Length_Chain1'] = df['Length_Chain1'].str.extract(r'(\s+\d+\.?\d*)')
df['Length_Chain2'] = df['Length_Chain2'].str.extract(r'(\s+\d+\.?\d*)')

#extracting aligned length, RMSD, Seq_ID
df['Aligned_length'] = df['Aligned_length'].str.extract(r'(\s+\d+\.?\d*)')
df['RMSD'] = df['RMSD'].str.extract(r'(\s+\d+\.?\d*)')
df['Seq_ID'] = df['Seq_ID'].str.extract(r'(\s+\d+\.?\d*)')

#extracting TMscores
df['TM-score_Chain1'] = df['TM-score_Chain1'].str.extract(r'(\s+\d+\.?\d*\s+)')
df['TM-score_Chain2'] = df['TM-score_Chain2'].str.extract(r'(\s+\d+\.?\d*\s+)')


# %%
df = df.drop('Aligned length', axis=1)

# %%
df = df.apply(pd.to_numeric, errors='ignore')

# %%
df.to_csv(snakemake.output[0], index=False)

# %%



