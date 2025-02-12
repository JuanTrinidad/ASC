# %%
import pandas as pd

# %%
#load SRBH results
df_RBH = pd.read_csv(snakemake.input.SRBH, sep='\t', index_col=0)

#load FATCAT & TMscore
df_TMscore_FATCAT = pd.read_csv(snakemake.input.FATCATTMscore, sep=',')

#load Reference Organism Annotation
df_reference_org_annotation = pd.read_csv(snakemake.input.annotation, sep='\t', low_memory=False, index_col=0)

#load cluster info
df_cluster_info = pd.read_csv(snakemake.input.pLDDT, sep='\t', index_col=0)
df_cluster_info.columns = ['query_GeneID', 'pdb2', 'query_struct_pLDDTmean']

cluster_file_input = snakemake.input.cluster_file

output_FinalTable = snakemake.output[0]
output_ExtensionFile = snakemake.output[1]

# %%
# Merge three datasets

df_merge = df_RBH.merge(df_TMscore_FATCAT, left_on='target_uniprot_accession', right_on='pdb1', how='left')

df_merge = df_merge.merge(df_reference_org_annotation, left_on='target_uniprot_accession', right_on='Entry', how='left')

# %%
print(df_merge.shape)

df_merge = df_cluster_info.merge(df_merge, right_on='query_uniprot_accession', left_on='pdb2', how='right').drop(columns=['pdb2_x'])



# %%
df_merge.to_csv(output_FinalTable, sep='\t')

# %% [markdown]
# # Cluster extention file

# %%
#OGroup assignation
OG_and_members = pd.read_csv(cluster_file_input, 
                             sep='\t', 
                             header=None, 
                             names=['cluster_representer', 'members'])

# %%
OG_and_members_tmp = OG_and_members[OG_and_members['members'].isin(df_merge['query_GeneID'].unique())]

OG_and_members_tmp.columns = ['members', 'cluster_representer']

# %%
OG_and_members_tmp = OG_and_members_tmp.merge(OG_and_members, left_on='members', right_on='cluster_representer')[['cluster_representer_x', 'members_y']]
OG_and_members_tmp.columns = ['cluster_representer', 'members']
OG_and_members_tmp.to_csv(output_ExtensionFile)

# %%



