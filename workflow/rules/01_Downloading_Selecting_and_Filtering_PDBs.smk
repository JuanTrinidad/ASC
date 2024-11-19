####################################################
#
#
#    Obtaining all available query structure data
#
#
####################################################

################################
# DOWNLOADING FROM ALPHA FOLD DB
################################

rule download_structure_files_from_alphafold_DB:
  input:
    config['input_files']['header_to_uniprot'],
    'report/{all_sequence_fasta}_Ortholog_group_to_geneID.tsv'
  output: 
    #pdb_files_full_path, #esto es necesario despues descomentar
    touch('tmp/{all_sequence_fasta}_download_finished.out')
  retries: 
    config['downloading_from_AFDB']['retries']
  threads:
    workflow.cores
  conda:
    '../envs/env_nameconverter_and_AFDBdownload.yaml'
  params:
    prefix_of_added_pdbs = config['structure_from_outside_AFDB']['prefix_in_PDB_name']
  script:
    '../scripts/001_download_prot_struct_from_AFDB.py'
    

#############################################################################################
# Calculation of pLDDT mean to use as filter to select representitative structure in cluster
#############################################################################################

rule pLDDT_mean_calculation:
  input: 
    tmp_download = 'tmp/{all_sequence_fasta}_download_finished.out'
  output: 
    'report/{all_sequence_fasta}_protein_structure_pLDDT_mean.tsv'
  threads: 
    workflow.cores
  conda:
    '../envs/env_pLDDT_mean_calc.yaml'
  script:
    '../scripts/002_pLDDT_mean_calculation-multiprocess.py'


###################################
#  Removing wrong structures files
###################################


rule checking_if_SameSimilar_fastaseq_vs_pdbseq:
  input:
    fh_to_uniprot = config['input_files']['header_to_uniprot'],
    fasta_query_db = config['input_files']['all_sequence_fasta'],
    tmp_download = 'tmp/{all_sequence_fasta}_download_finished.out'
  output:
    'tmp/{all_sequence_fasta}_sequence_comparison_between_fasta_and_downloadedPDBfile.tsv'
  threads:
    workflow.cores
  conda:
    '../envs/env_pLDDT_mean_calc.yaml'
  script:
    '../scripts/002.2_removing_wrong_structures.py'


#############################
# CREATING REPORTS pLDDT DATA
#############################

rule ortho_MCL_modelated_report:
  input: 
    file1 = 'report/{all_sequence_fasta}_protein_structure_pLDDT_mean.tsv',
    file2 = config['input_files']['header_to_uniprot'],
    file3 = 'report/{all_sequence_fasta}_Ortholog_group_to_geneID.tsv',
    file4 = 'tmp/{all_sequence_fasta}_sequence_comparison_between_fasta_and_downloadedPDBfile.tsv' #puede ser opcional, es largo el script que lo genera
  output: 
    out = 'report/{all_sequence_fasta}_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
    outstats= 'report/{all_sequence_fasta}_ortholog_groups_structure_stats.tsv'
  conda: 
    '../envs/env_pLDDT_mean_calc.yaml'
  script:
    '../scripts/003_assigning_structure_to_ortholog_group.py'
    
    

####################################
# SELECTED STRUCTURES TO CREATE DB
####################################

rule cluster_representer_protein_structures:
  input: 
    'report/{all_sequence_fasta}_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
  output: 
    temp(touch('tmp/{all_sequence_fasta}_files_copied.done'))
  conda:
    '../envs/env_nameconverter_and_AFDBdownload.yaml'
  params: 
    path = 'genome_data_sets/query_proteomes/pdb_files/cluster_structure_representers'
  threads: workflow.cores
  script:
    '../scripts/004_cluster_representative_protein_structures_selection_to_db.py'



#########################################
# CREATE FASTA FOR NON MODELATED CLUSTERS
#########################################

rule fasta_of_orthoGroups_without_structure_in_AFDB:
  input: 
    original_fasta = config['input_files']['all_sequence_fasta'],
    report_ortho_g = 'report/{all_sequence_fasta}_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
    ortho_info = 'report/{all_sequence_fasta}_Ortholog_group_to_geneID.tsv'    
  output:
    output_fasta_file = 'report/fasta_files/{all_sequence_fasta}_protein_sequences_from_cluster_wo_structure_in_AFDB.fasta'
  conda:
    '../envs/env_pLDDT_mean_calc.yaml'
  params:
    OG_size = config['fasta_file_of_non_modelated_clusters']['ortholog_group_above_this_num_of_members']
  script:
    '../scripts/005_fasta_file_creation.py'
