####################################################
#
#
#                    FOLDSEEK
#
#
####################################################





####################################################################################
# CREATING FOLDSEEK DB FROM DOWNLOADED PROTEIN STRUCTURES - CLUSTER REPRESENTERS
####################################################################################


rule foldseek_db_query_proteins:
  input: 
    'report/{all_sequence_fasta}_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
    'tmp/{all_sequence_fasta}_files_copied.done' 
  output:
    temp(multiext('genome_data_sets/query_proteomes/foldseek_data_base/{all_sequence_fasta}_DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'))
  conda:
    '../envs/env_foldseek.yaml'
  params: path =  'genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers/'
  threads: workflow.cores
  shell:
    'foldseek createdb {params.path} {output[0]} --threads {threads}'




#############################################################################
# DOWNLOADING MODEL SPECIES AND CREATING FOLDSEEK DB FOR EACH MODEL ORGANISMS
#############################################################################


rule download_model_organisms:
  input:
    '../config/mandatory_files/link_to_download_from_AFDB/{organism}.tar'
  output:
    temp('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/{organism}.tar')
  shell:
    'wget -q -P genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/ https://ftp.ebi.ac.uk/pub/databases/alphafold/v4/{wildcards.organism}.tar'



rule create_separates_DB_for_model_organisms:
  input: 
    'genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/{organisms}.tar'
  output:
    multiext('genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  conda:
    '../envs/env_foldseek.yaml'
  threads: workflow.cores
  shell:
    'foldseek createdb {input} {output[0]} --threads {threads}'



##############################
# RECIPROCAL BEST HIT
##############################


rule foldseek_reciprocal_best_hit_single_organisms:
  input: 
    query = multiext('genome_data_sets/query_proteomes/foldseek_data_base/{all_sequence_fasta}_DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'), 
    subject = multiext('genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  output: 
    multiext('../results/reciprocal_best_hit_SingleOrgApproach/{all_sequence_fasta}_cluster_representer_vs_{organisms}','','.dbtype','.index')
  conda:
    '../envs/env_foldseek.yaml'
  threads: workflow.cores
  params:
    sensitivity = config['reciprocal_best_hit_parameters']['sensitivity'],
    cov_mode = config['reciprocal_best_hit_parameters']['cov_mode'],
    coverage = config['reciprocal_best_hit_parameters']['coverage']
  shell:
    'foldseek rbh {input.query[0]} {input.subject[0]} {output[0]} /tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -a'



##############################
# EXTRACT TSV FROM RBH
##############################


rule foldseek_reciprocal_best_hit_extract_result_tsv:
  input: 
    query = multiext('genome_data_sets/query_proteomes/foldseek_data_base/{all_sequence_fasta}_DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    subject = multiext('genome_data_sets/subject_proteomes/foldseek_data_base/individual_org_DB/{organisms}', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    table = multiext('../results/reciprocal_best_hit_SingleOrgApproach/{all_sequence_fasta}_cluster_representer_vs_{organisms}', '', '.dbtype','.index')
  output:
    '../results/reciprocal_best_hit_SingleOrgApproach_TSV/{all_sequence_fasta}_cluster_representer_vs_{organisms}.tsv'  
  conda:
    '../envs/env_foldseek.yaml'
  threads: workflow.cores
  shell:
    'foldseek createtsv {input.query[0]} {input.subject[0]} {input.table[0]} {output} --threads {threads}'






##############################
# CONCAT TSV FILES IN ONE FILE
##############################

rule foldseek_reciprocal_best_hit_concat_TSV:
  input:
    expand('../results/reciprocal_best_hit_SingleOrgApproach_TSV/{all_sequence_fasta}_cluster_representer_vs_{organisms}.tsv', all_sequence_fasta = initial_fasta_file_name_clean , organisms = model_organisms_files_final_db)
  output:
    '../results/reciprocal_best_hit_SingleOrgApproach_TSV/{all_sequence_fasta}_rbh_all_in_one_file.tsv'
  conda:
    '../envs/env_pLDDT_mean_calc.yaml'
  script:
    '../scripts/006_concat_files_and_add_info.py'

































###################
#
# Rules cementery
#
###################









##############################
# FOLDSEEK SEARCH
##############################
'''
rule download_Uniprot50DB:
  output: 
    multiext('protein_data_bases/UniProt50', '' ,'.dbtype','.index','.lookup','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  shell:
    'foldseek databases Alphafold/UniProt50 protein_data_bases/UniProt50 /tmp --threads {threads}'
'''

'''
rule foldseek_search_Uniprot50DB:
  input: 
    query = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    subject = multiext('protein_data_bases/UniProt50', '' ,'.dbtype','.index','.lookup','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index') #_mapping, _taxonomy, not mandatory
  output:
    multiext('../results/foldseek_search/cluster_representer_vs_UniProt50', '.index', '.dbtype') 
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  params:
    path = lambda wildcards, output: output[0][:-6],
    sensitivity = config['foldseek_search_against_Uniprot50_parameters']['sensitivity'],
    cov_mode = config['foldseek_search_against_Uniprot50_parameters']['cov_mode'],
    coverage = config['foldseek_search_against_Uniprot50_parameters']['coverage'],
    evalue = config['foldseek_search_against_Uniprot50_parameters']['evalue']
  shell:
    'foldseek search {input.query[0]} {input.subject[0]} {params.path} /tmp --threads {threads} -s {params.sensitivity} -c {params.coverage} -e {params.evalue}'    
'''  
    
    


##############################
# EXTRACT TSV FROM FOLDSEEK SEARCH
##############################
# esto estaria bueno generalizarlo para que se pueda poner otra DB
'''
rule foldseek_search_ Uniprot50DB_extract_result_tsv:
  input: 
    query = multiext('genome_data_sets/query_proteomes/foldseek_data_base/DB_cluster_representer', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'),
    subject = multiext('protein_data_bases/UniProt50', '' ,'.dbtype','.index','.lookup', '_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index'), #_mapping, _taxonomy, not mandatory
    table =  multiext('../results/foldseek_search/cluster_representer_vs_UniProt50','.index', '.dbtype')
  output:
     '../results/foldseek_search_TSV/cluster_representer_vs_Uniprot50.tsv'
  conda:
    'envs/env_foldseek.yaml'
  threads: workflow.cores
  params: path =  lambda wildcards, input: input.table[0][:-6]
  shell:
    'foldseek createtsv {input.query[0]} {input.subject[0]} {params.path} {output} --threads {threads}'
'''

'''
rule clustering_single_DB_for_model_organisms:
    input:
        multiext('genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB/all_model_organisms_DB', '', '.dbtype','.index','.lookup','.source','_ca', '_ca.dbtype', '_ca.index','_h','_h.dbtype','_h.index', '_ss', '_ss.dbtype','_ss.index')
    output:
        multiext('genome_data_sets/subject_proteomes/foldseek_data_base/all_org_DB_clustered/all_model_organisms_DB_clustered', '.dbtype', '.index')
    conda:
        '../envs/env_foldseek.yaml'  
    params: 
        path = lambda wildcards, output: output[1][:-6],
        sensitivity = config['foldseek_clustering_model_org_DB_parameters']['sensitivity'],
        cov_mode = config['foldseek_clustering_model_org_DB_parameters']['cov_mode'],
        coverage = config['foldseek_clustering_model_org_DB_parameters']['coverage'],
        evalue = config['foldseek_clustering_model_org_DB_parameters']['evalue']
    threads: workflow.cores
    shell:
        'foldseek cluster {input[0]} {params.path} /tmp --threads {threads} -s {params.sensitivity} --cov-mode {params.cov_mode} -c {params.coverage} -e {params.evalue}'
        
'''
