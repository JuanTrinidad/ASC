__author__ = "Juan Trinidad"
__copyright__ = "Copyright 2024, Juan Trinidad"
__email__ = "jtrindad@fcien.edu.uy"
__license__ = "MIT"


###################
# LIBRARIES
###################
import pandas as pd
import glob
import pickle



##############
# CONFIG
##############

configfile: "../config/config.yaml"


################
# STARTING FILES
################

#-------------------------------------------------------
# starting fasta file
#-------------------------------------------------------
initial_fasta_file_name = config['input_files']['all_sequence_fasta']
print('Fasta file:\n' , initial_fasta_file_name, '\n')

initial_fasta_file_name_clean = initial_fasta_file_name.split('/')[-1][:-3]

print(initial_fasta_file_name_clean)
#-------------------------------------------------------




#-------------------------------------------------------
# model organisms  
#-------------------------------------------------------
model_organisms_files = glob.glob('../config/mandatory_files/link_to_download_from_AFDB/*.tar')
model_organisms_files_final = [file.split('/')[-1] for file in model_organisms_files]

#individual org db files
model_organisms_files_final_db = [file.split('/')[-1][:-4] for file in model_organisms_files]




print('This are the model organisms that will be used from AFDB v4:')
for file in model_organisms_files_final:
    print(file)

#-------------------------------------------------------
    


optional_file = f'report/fasta_files/{initial_fasta_file_name_clean}_protein_sequences_from_cluster_wo_structure_in_AFDB.fasta' if config['fasta_file_of_non_modelated_clusters']['create_fasta_file'] else '../config/mandatory_files/fasta_header_to_uniprot.tsv'



include:'rules/testing.smk'
include:'rules/00_MMseq2_sequence_clustering.smk'
include:'rules/01_Downloading_Selecting_and_Filtering_PDBs.smk'
include:'rules/02_Foldseek_rules-SingleOrgApproach.smk'
include:'rules/03_FATCAT.smk'
include:'rules/04_Downloading_annotation_from_uniprot.smk'

rule all: 
  input: f'tmp/{initial_fasta_file_name_clean}_extract_twisted_structures.out', f'report/{initial_fasta_file_name_clean}_TMscores_from_TMalign_twisted_structure.tsv'
  
  
  

