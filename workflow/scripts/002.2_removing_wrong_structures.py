#!/usr/bin/env python
# coding: utf-8

# In[34]:


import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.PDB import *
import multiprocessing
import glob


# In[62]:


#snakemake.input.fh_to_uniprot

df_fasta_header_to_uniprot = pd.read_csv(snakemake.input.fh_to_uniprot, 
                                         sep='\t', 
                                         names=['gene_name', 'UNIPROT_accession'])


# In[63]:


list_for_df = []
for file in glob.glob('genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF*.pdb'):
    
    uniprot_accession = file.split('/')[-1].split('-')[1]
    
    list_for_df.append([file, uniprot_accession])
    
df_available_structures = pd.DataFrame(list_for_df, columns=['file_path', 'UNIPROT_accession'])




df_fasta_header_to_uniprot = df_fasta_header_to_uniprot.merge(df_available_structures, how='inner')


# In[ ]:





# In[65]:


genID_aaSequence = []
#'../genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa'
#snakemake.input.fasta_query_db
with open( snakemake.input.fasta_query_db , 'r') as fasta_file:


    for record in SeqIO.parse(fasta_file, "fasta"):
        genID_aaSequence.append([record.id, record.seq])


# In[66]:



df1 = pd.DataFrame(genID_aaSequence, columns=['gene_name', 'sequence'])

df_allinfo = df_fasta_header_to_uniprot.merge(df1, left_on='gene_name', right_on='gene_name', how='left').dropna()


# In[67]:


df_allinfo_aslist = df_allinfo.to_numpy().tolist()



def align_fasta_sequence(parameters_list):
    """Aligns a FASTA sequence against a PDB or CIFF file.

    Args:
    fasta_sequence: The FASTA sequence to align.
    pdb_file: The PDB or CIFF file to align against.

    Returns:
    The alignment of the FASTA sequence against the PDB or CIFF file.
    """

    
    pdb_file = parameters_list[2]
    fasta_sequence = parameters_list[3]

    
    p = PDBParser()
    structure = p.get_structure("XXX", pdb_file)
    ppb=PPBuilder()
    PDBsequence = ppb.build_peptides(structure)[0].get_sequence()
    
    len_fasta_seq = len(fasta_sequence)
    len_PDBsequence = len(PDBsequence)
    
    if len_fasta_seq == len_PDBsequence:
        
        score = len_fasta_seq
    
    else:
        alignment = pairwise2.align.globalxx(fasta_sequence, PDBsequence)
        score = alignment[0].score
    
    return [parameters_list[0], parameters_list[1], score, len_fasta_seq, len_PDBsequence]


def align_fasta_sequence_parallel(num_threads):
    
    # Create a pool of processes.
    pool = multiprocessing.Pool(num_threads) 


    # Map the extract_values function to each file name in the list.
    results = pool.map(align_fasta_sequence, df_allinfo_aslist)
    
    # Close the pool.
    pool.close()
    pool.join()    
    
    return results




salida = align_fasta_sequence_parallel(snakemake.threads)




# In[9]:


df = pd.DataFrame(salida, columns=['genID', 'UNIPROT', 'score', 'len_fasta', 'lenPDB'])

df.to_csv(snakemake.output[0], sep='\t', index=None)

#'../tmp/sequence_comparison_between_fasta_and_downloadedPDBfile.tsv'





