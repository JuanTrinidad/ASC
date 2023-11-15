###################
# LIBRARIES
###################
import pandas as pd
import glob

##############
# CONFIG
##############

configfile: "../config/config.yaml"


#############
# FILES
############

#-------------------------------------------------------
#starting fasta file
#-------------------------------------------------------
initial_fasta_file_name = config['input_files']['all_sequence_fasta'].split('/')[-1]
print('Fasta file:\n' ,initial_fasta_file_name, '\n')
#-------------------------------------------------------




#-------------------------------------------------------
# model organisms  
#-------------------------------------------------------
model_organisms_files = glob.glob('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/*')
model_organisms_files_final = [file.split('/')[-1] for file in model_organisms_files]

#individual org db files
model_organisms_files_final_db = [file.split('/')[-1][:-4] for file in model_organisms_files]

print('This are the model organisms files provided:')
for file in model_organisms_files:
    print(file.split('/')[-1] )
    
#-------------------------------------------------------



rule all:
  input: '../results/foldseek_search_TSV/cluster_representer_locally_modelated_proteins_vs_SwissProt.tsv', expand('../results/reciprocal_best_hit_TSV/cluster_representer_locally_modelated_proteins_vs_{organisms}.tsv', organisms = [file[:-4] for file in model_organisms_files_final])


#############################################################################################
# CALCULATION OF pLDDT MEAN to use as filter to select representitative structure in cluster
#############################################################################################

rule LMP_pLDDT_mean_calculation:
  input: 
    'tmp/download_finished.out'
  output: 
    'report/protein_structure_pLDDT_mean_locally_modelated_prot.tsv'
  threads: 
    config['pLDDT_mean_calculation']['threads']
  conda: 
    '../envs/pLDDT_mean_calc.yaml'
  params:
    path = 'genome_data_sets/query_proteomes/pdb_files/locally_modelated_proteins'
  script: 
    '../scripts/rules_scripts/2_pLDDT_mean_calculation-multiprocess-LMP.py'



#############################
# CREATING REPORTS pLDDT DATA
#############################

#add info, stats?
rule LMP_ortho_MCL_modelated_report:
  input: 
    file1 = 'report/protein_structure_pLDDT_mean_locally_modelated_prot.tsv',
    file2 = config['input_files']['ortholog_group_to_geneID']
  output: 
    'report/ortholog_groups_x_sequence_clustering_x_UNIPROT_locally_modelated_prot.tsv'
  conda:
    '../envs/pLDDT_mean_calc.yaml'
  script:
    '../scripts/rules_scripts/3_LMS_assigning_structure_to_OrthoGroup.py'



####################################
# SELECTED STRUCTURES TO CREATE DB
####################################

rule cluster_representer_protein_structures:
  input: 
    'report/protein_structure_pLDDT_mean_locally_modelated_prot.tsv'
  output: 
    touch('tmp/files_copied_LMP.done')
  conda:
    '../envs/env_nameconverter_and_AFDBdownload.yaml'
  params: 
    path ='genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers_locally_modelated'
  threads: workflow.cores
  script:
    '../scripts/rules_scripts/4_cluster_representative_protein_structures_selection_to_db.py'
    
    
########################################
#
#
#     FOLDSEEK
#
#
########################################


###################################################################################
# CREATING FOLDSEEK DB FROM DOWNLOADED PROTEIN STRUCTURES - CLUSTER REPRESENTERS
####################################################################################


rule foldseek_db_query_proteins:
  input: 
    'tmp/files_copied_LMP.done' #protein_files
  output:
    'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_locally_modelated_proteins'
  conda:
    '../envs/env_foldseek.yaml'
  shell:
    'foldseek createdb genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers_locally_modelated/ {output}'
    
    
    
##############################
# RECIPROCAL BEST HIT
##############################


rule foldseek_reciprocal_best_hit:
  input: 
    query = 'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_locally_modelated_proteins', 
    subject = 'genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}.source'
  output:
    '../results/reciprocal_best_hit/cluster_representer_locally_modelated_proteins_vs_{organisms}'  
  conda:
    '../envs/env_foldseek.yaml'
  threads: 40
  params:
    sensitivity = config['reciprocal_best_hit_parameters']['sensitivity'],
    cov_mode = config['reciprocal_best_hit_parameters']['cov_mode'],
    coverage = config['reciprocal_best_hit_parameters']['coverage']
  shell:
    'foldseek rbh {input.query} genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{wildcards.organisms} {output} ./tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -a'    


##############################
# EXTRACT TSV FROM RBH
##############################


rule foldseek_reciprocal_best_hit_extract_result_tsv:
  input: 
    query = 'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_locally_modelated_proteins',
    #subject = 'genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}.source',
    table = expand('../results/reciprocal_best_hit/cluster_representer_locally_modelated_proteins_vs_{organisms}', organisms = [file[:-4] for file in model_organisms_files_final])
  output:
    '../results/reciprocal_best_hit_TSV/cluster_representer_locally_modelated_proteins_vs_{organisms}.tsv'  
  conda:
    '../envs/env_foldseek.yaml'
  threads: 40
  shell:
    'foldseek createtsv {input.query} genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{wildcards.organisms} ../results/reciprocal_best_hit/cluster_representer_locally_modelated_proteins_vs_{wildcards.organisms} {output}'
    
    
    
##############################
# FOLDSEEK SEARCH
##############################


rule foldseek_search_ SwissProt:
  input: 
    query = 'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_locally_modelated_proteins',
    subject = 'protein_data_bases/Swiss-Prot'
  output:
    '../results/foldseek_search/cluster_representer_locally_modelated_proteins_vs_SwissProt.index'  
  conda:
    '../envs/env_foldseek.yaml'
  threads: 40
  params:
    sensitivity = config['foldseek_search_against_SwissProt_parameters']['sensitivity'],
    cov_mode = config['foldseek_search_against_SwissProt_parameters']['cov_mode'],
    coverage = config['foldseek_search_against_SwissProt_parameters']['coverage'],
    evalue = config['foldseek_search_against_SwissProt_parameters']['evalue']
  shell:
    'foldseek search \
    {input.query} \
    {input.subject} \
    ../results/foldseek_search/cluster_representer_locally_modelated_proteins_vs_SwissProt \
    ./tmp \
    --threads {threads} \
    -s {params.sensitivity} \
    -c {params.coverage} \
    -a \
    -e {params.evalue}'
    
    
##############################
# EXTRACT TSV FROM FOLDSEEK SEARCH
##############################
# esto estaria bueno generalizarlo para que se pueda poner otra DB

rule foldseek_search_ SwissProt_extract_result_tsv:
  input: 
    query = 'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_locally_modelated_proteins',
    subject = 'protein_data_bases/Swiss-Prot',
    table = '../results/foldseek_search/cluster_representer_locally_modelated_proteins_vs_SwissProt.index'  
  output:
     '../results/foldseek_search_TSV/cluster_representer_locally_modelated_proteins_vs_SwissProt.tsv'
  conda:
    '../envs/env_foldseek.yaml'
  threads: 40
  shell:
    'foldseek createtsv {input.query} {input.subject} ../results/foldseek_search/cluster_representer_locally_modelated_proteins_vs_SwissProt {output} ./tmp'