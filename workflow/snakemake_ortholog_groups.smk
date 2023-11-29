__author__ = "Juan Trinidad"
__copyright__ = "Copyright 2023, Juan Trinidad"
__email__ = "jtrindad@fcien.edu.uy"
__license__ = "MIT"


###################
# LIBRARIES
###################
import pandas as pd
import glob



##############
# CONFIG
##############

configfile: "../config/config.yaml"


################
# STARTING FILES
################

#-------------------------------------------------------
#starting fasta file
#-------------------------------------------------------
initial_fasta_file_name = config['input_files']['all_sequence_fasta'].split('/')[-1]
print('Fasta file:\n' , initial_fasta_file_name, '\n')
#-------------------------------------------------------




#-------------------------------------------------------
# model organisms  
#-------------------------------------------------------
model_organisms_files = glob.glob('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/*.tar')
model_organisms_files_final = [file.split('/')[-1] for file in model_organisms_files]

#individual org db files
model_organisms_files_final_db = [file.split('/')[-1][:-4] for file in model_organisms_files]


print('This are the model organisms files provided:')
for file in model_organisms_files_final:
    print(file)


#-------------------------------------------------------
    
#### CIF ###
#cif_files = [cif for cif in glob.glob('tmp/FATCAT_pdb_files/*.cif')]
  

#-------------------------------------------------------
#creating path to control data download
#-------------------------------------------------------
#df_UNIPROT = pd.read_csv('../config/mandatory_files/fasta_header_to_uniprot.tsv', sep='\t', header=None, names=['GeneID', 'UNIPROT'])
#pdb_files = [f'AF-{UNIPROTaccession}-F1-model_v4.pdb' for UNIPROTaccession in df_UNIPROT.UNIPROT.unique()]

#full path
#pdb_files_full_path = expand('genome_data_sets/query_proteomes/pdb_files/prot_structure_download_from_AlphaFoldDB/{pdb_files}', pdb_files = pdb_files)
#-------------------------------------------------------



# name="some/file.txt" if config["condition"] else "other/file.txt"

#print('protein_sequences_from_cluster_wo_structure_in_AFDB.fasta' if config['fasta_file_of_non_modelated_clusters']['create_fasta_file'] else 'config/mandatory_files/fasta_header_to_uniprot.tsv')
#optional_file = 'protein_sequences_from_cluster_wo_structure_in_AFDB.fasta' if config['fasta_file_of_non_modelated_clusters']['create_fasta_file'] else '../config/mandatory_files/fasta_header_to_uniprot.tsv'



include:'rules/testing.smk'
include:'rules/00_MMseq2_sequence_clustering.smk'
include:'rules/01_Downloading_Selecting_and_Filtering_PDBs.smk'
include:'rules/02_Foldseek_rules-SingleOrgApproach.smk'
include:'rules/03_FATCAT.smk'
include:'rules/04_Downloading_annotation_from_uniprot.smk'

rule all: 
  input: 'tmp/TriTrypDB-65_All_species_clean_extract_twisted_structures.out'
  
  #expand('genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}', organisms = model_organisms_files_final_db)
  #expand('tmp/FATCAT_pdb_files/{file}.cif', file = cif_files)
  #expand('genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}', organisms = model_organisms_files_final_db)
  #'tmp/to_compare.list' #'../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv'
  #'genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB'
  #'../results/foldseek_search_TSV/cluster_representer_vs_Uniprot50.tsv', '../results/reciprocal_best_hit_TSV/rbh_all_in_one_file.tsv' #, '../results/Gene_annotation_info_from_uniprot_model_spp.tsv'        
    #'../results/foldseek_search_TSV/cluster_representer_vs_Uniprot50.tsv'
    #expand('../results/reciprocal_best_hit_TSV/cluster_representer_vs_{organisms}.tsv', organisms = model_organisms_files_final_db),
    #'../results/foldseek_search_TSV/cluster_representer_vs_SwissProt.tsv' #
    #'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer' #
    #expand('../results/reciprocal_best_hit/cluster_representer_vs_{organisms}', organisms = model_organisms_files_final_db)
    #'../results/foldseek_search_TSV/cluster_representer_vs_SwissProt.tsv'
    #expand('../results/reciprocal_best_hit/cluster_representer_vs_{organisms}', organisms = [file[:-4] for file in model_organisms_files_final]),
    #'../results/foldseek_search/cluster_representer_vs_SwissProt.index',
    #'../results/foldseek_search_TSV/cluster_representer_vs_SwissProt.tsv',
    #'../results/Gene_annotation_info_from_uniprot_model_spp.tsv'
    #expand('genome_data_sets/subject_proteomes/annotation_info/{organisms}.tsv', organisms = [file[:-4] for file in model_organisms_files_final])
    #optional_file,
    #'genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB'
    
    
    
   
    

















###########################################
#
#
#                TESTING
#
#
###########################################





 
###############
#
#
#   OLD RULES
#
#
##############


'''

#########################################
# CREATE FASTA FOR NON MODELATED CLUSTERS
#########################################

rule fasta_of_orthoGroups_without_structure_in_AFDB:
  input: 
    original_fasta = 'genome_data_sets/query_proteomes/fasta_files/TriTrypDB-63_All_species_clean.fa',
    report_ortho_g = 'report/ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
    ortho_info = '../config/mandatory_files/Ortholog_group_to_geneID.tsv'    
  output:
    output_fasta_file = 'report_files/fasta_files/protein_sequences_from_cluster_wo_structure_in_AFDB.fasta'
  conda:
    'envs/env_pLDDT_mean_calc.yaml'
  params:
    OG_size = config['fasta_file_of_non_modelated_clusters']['ortholog_group_above_this_num_of_members']
  script:
    'scripts/5_fasta_file_creation.py'



##############################################
# CLUSTERING QUERY DB BY STRUCTURE
##############################################

# SEARCH
rule clustering_foldseek_db_query_proteins_STEP1:
  input:
    multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  output:
    multiext('tmp/clustering_intermediate_files/aln', '.dbtype', '.index')
  conda:
    'envs/env_foldseek.yaml'
  params: 
    path = lambda wildcards, output: output[0][:-7],
    sensitivity = config['clustering_foldseek_db_query_proteins']['sensitivity'],
    cov_mode = config['clustering_foldseek_db_query_proteins']['cov_mode'],
    coverage = config['clustering_foldseek_db_query_proteins']['coverage'],
    evalue = config['clustering_foldseek_db_query_proteins']['evalue'], 
    minseqid = config['clustering_foldseek_db_query_proteins']['minseqid']
  threads: workflow.cores
  shell:
    'foldseek search {input[0]} {input[0]} {params.path} /tmp -c {params.coverage} -e {params.evalue} -s {params.sensitivity} --min-seq-id {params.minseqid} --threads {threads}'


# CLUSTER 
rule clustering_foldseek_db_query_proteins_STEP2:
  input:
        DB = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
        aln = multiext('tmp/clustering_intermediate_files/aln', '.dbtype', '.index')
  output:
    multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_clustered', '','.dbtype', '.index')
  conda:
    'envs/env_foldseek.yaml'  
  params: 
    path1 = lambda wildcards, input: input.aln[0][:-7],
    path2 = lambda wildcards, output: output[0]
  threads: workflow.cores
  shell:
    'foldseek clust {input.DB[0]} {params.path1} {params.path2} --threads {threads}'

# TSV
rule clustering_foldseek_db_query_proteins_STEP3:
  input:
    DB = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    clustered = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_clustered', '','.dbtype', '.index')
  output:
    'genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer_clustered.tsv'
  conda:
    'envs/env_foldseek.yaml'      
  threads: workflow.cores
  shell:
    'foldseek createtsv {input.DB[0]} {input.DB[0]} {input.clustered[0]} {output}'
    
    
rule foldseek_search_ SwissProt:
  input: 
    query = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    subject = multiext('protein_data_bases/Swiss-Prot', '' ,'.dbtype','.index','.lookup','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index') #_mapping, _taxonomy, not mandatory
  output:
    multiext('../results/foldseek_search/cluster_representer_vs_SwissProt', '.index', '.dbtype') 
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  params:
    sensitivity = config['foldseek_search_against_SwissProt_parameters']['sensitivity'],
    cov_mode = config['foldseek_search_against_SwissProt_parameters']['cov_mode'],
    coverage = config['foldseek_search_against_SwissProt_parameters']['coverage'],
    evalue = config['foldseek_search_against_SwissProt_parameters']['evalue'],
    path =  lambda wildcards, output: output[0][:-6] 
    
  shell:
    'foldseek search {input.query[0]} {input.subject[0]} {params.path} ./tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -a -e {params.evalue}'.
    
    
    
    
rule clustering_Uniprot50DB:
  input: 
    multiext('protein_data_bases/UniProt50', '' ,'.dbtype','.index','.lookup','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  output: 
    multiext('protein_data_bases/UniProt50_clustered', '.index','.dbtype')
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  params:
    path = lambda wildcards, output: output[0][:-6],
    sensitivity = config['foldseek_clustering_Uniprot50_parameters']['sensitivity'],
    cov_mode = config['foldseek_clustering_Uniprot50_parameters']['cov_mode'],
    coverage = config['foldseek_clustering_Uniprot50_parameters']['coverage'],
    evalue = config['foldseek_clustering_Uniprot50_parameters']['evalue']
  shell:
    'foldseek search {input[0]} {input[0]} {params.path} tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -e {params.evalue}'
    
'''

#-------------------------------------------------------------------------------------------- 
'''    
rule clustering_single_DB_for_model_organisms_STEP1:
  input: 
    multiext('genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  output: 
    multiext('tmp/all_model_organisms_DB_searched', '.index','.dbtype')
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  params:
    path = lambda wildcards, output: output[0][:-6],
    sensitivity = config['foldseek_clustering_model_org_DB_parameters']['sensitivity'],
    cov_mode = config['foldseek_clustering_model_org_DB_parameters']['cov_mode'],
    coverage = config['foldseek_clustering_model_org_DB_parameters']['coverage'],
    evalue = config['foldseek_clustering_model_org_DB_parameters']['evalue'],
    fident = config['foldseek_clustering_model_org_DB_parameters']['fident']
  shell:
    'foldseek search {input[0]} {input[0]} {params.path} tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -e {params.evalue} --min-seq-id {params.fident}' 


rule clustering_single_DB_for_model_organisms_SETP2:
  input: 
    db = multiext('genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    aln = multiext('tmp/all_model_organisms_DB_searched', '.index','.dbtype')
  output:
    multiext('genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB_clustered', '.index' , '.dbtype')
  conda:
    'envs/env_foldseek.yaml'  
  params: 
    path1 = lambda wildcards, input: input.aln[0][:-6],
    path2 = lambda wildcards, output: output[0][:-6]
  threads: workflow.cores
  shell:
    'foldseek clust {input.db[0]} {params.path1} {params.path2} --threads {threads}'

'''
#--------------------------------------------------------------------------------------------