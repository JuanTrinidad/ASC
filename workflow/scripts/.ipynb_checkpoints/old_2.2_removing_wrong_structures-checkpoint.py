#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.PDB import *
import multiprocessing



# In[ ]:





# In[ ]:





# In[2]:


genID_UNIPROT = []
#'../../config/mandatory_files/fasta_header_to_uniprot.tsv'
with open(snakemake.input.fh_to_uniprot, 'r') as file:
    
    
    for line in file.readlines():
        
        line = line.split()
        
        
        
        genID = line[0]
        UNIPROTaccession = line[1]
        
        pars = [genID, f'genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/AF-{UNIPROTaccession}-F1-model_v4.pdb', UNIPROTaccession]
        
        #print(pars)
        genID_UNIPROT.append(pars)
        
    


# In[3]:


genID_aaSequence = []

with open(snakemake.input.fasta_TriTrypDB , 'r') as fasta_file:


    for record in SeqIO.parse(fasta_file, "fasta"):
        genID_aaSequence.append([record.id, record.seq])


# In[4]:


df0 = pd.DataFrame(genID_UNIPROT)
df1 = pd.DataFrame(genID_aaSequence)

df_allinfo = df0.merge(df1, right_on=0, left_on=0, how='left').dropna()


# In[5]:


df_allinfo_aslist = df_allinfo.to_numpy().tolist()


# In[ ]:





# In[6]:


#df_allinfo_aslist = df_allinfo_aslist[:1000]


# In[ ]:





# In[ ]:





# In[ ]:





# In[7]:


def align_fasta_sequence(parameters_list):
    """Aligns a FASTA sequence against a PDB or CIFF file.

    Args:
    fasta_sequence: The FASTA sequence to align.
    pdb_file: The PDB or CIFF file to align against.

    Returns:
    The alignment of the FASTA sequence against the PDB or CIFF file.
    """

    
    pdb_file = parameters_list[1]
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
    
    return [parameters_list[0], parameters_list[2], score, len_fasta_seq, len_PDBsequence]


def align_fasta_sequence_parallel(num_threads):
    
    # Create a pool of processes.
    pool = multiprocessing.Pool(num_threads) 


    # Map the extract_values function to each file name in the list.
    results = pool.map(align_fasta_sequence, df_allinfo_aslist)
    
    # Close the pool.
    pool.close()
    pool.join()    
    
    return results


# In[ ]:





# In[ ]:





# In[ ]:





# In[8]:


salida = align_fasta_sequence_parallel(30)


# In[ ]:





# In[9]:


df = pd.DataFrame(salida, columns=['genID', 'UNIPROT', 'score', 'len_fasta', 'lenPDB'])

df.to_csv(snakemake.output[0], sep='\t', index=None)

#'../tmp/sequence_comparison_between_fasta_and_downloadedPDBfile.tsv'


# In[ ]:




