
import pandas as pd


#load SRBH results
df_RBH = pd.read_csv(snakemake.input.SRBH, sep='\t', index_col=0)

#load FATCAT & TMscore
df_TMscore_FATCAT = pd.read_csv(snakemake.input.FATCATTMscore, sep=',')

#load Reference Organism Annotation
df_reference_org_annotation = pd.read_csv(snakemake.input.annotation, sep='\t', low_memory=False, index_col=0)


# Merge three datasets

df_merge = df_RBH.merge(df_TMscore_FATCAT, left_on='target_uniprot_accession', right_on='pdb1', how='left')

df_merge = df_merge.merge(df_reference_org_annotation, left_on='target_uniprot_accession', right_on='Entry', how='left')

df_merge.to_csv(snakemake.output[0], sep='\t')




