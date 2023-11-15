print(workflow.basedir)


rule moving_files_from_clusters_to_FACTCATfolder:
    input: 
        file1 = 'report/TriTrypDB-65_All_species_clean_Ortholog_Group_Full_of_hypotetical_boolean.tsv',
        file2 = '../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv',
        file3 = '../results/Gene_annotation_info_from_uniprot_model_spp.tsv',
        file4 = 'report/TriTrypDB-65_All_species_clean_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv'
    output: touch('tmp/cluster_structure_files_copied.out')
    conda: '../envs/env_pLDDT_mean_calc.yaml'
    params: 
        destination_dir = 'tmp/FATCAT_pdb_files/',
        top_hits = 5,
        prefix_of_added_pdbs = config['structure_from_outside_AFDB']['prefix_in_PDB_name']
    threads: workflow.cores 
    script: '../scripts/007_polishing_pdbfiles_to_use_in_FATCAT.py'
        



rule unzip_uniprot_proteomes:
    input: 'genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/{proteome}.tar'
    output: 'genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/'
    shell: 'tar -xvf {input} -C {output}'









################
# CLONING GITHUB
################


rule clone_FATCAT_repository:
    output: 'git_repo_cloned/FATCAT/Install'
    shell: 'git clone https://github.com/GodzikLab/FATCAT-dist.git git_repo_cloned/FATCAT'
    
    
###################
# INSTALLING FATCAT
###################

rule installing_FATCAT_repository:
    input: 'git_repo_cloned/FATCAT/Install'
    output: touch('tmp/FATCAT_installed_sucesfully.out')
    shell: '''
    cd git_repo_cloned/FATCAT/
    ./Install
    
    '''
    
    
rule run_FATCAT:
    input: 'tmp/FATCAT_installed_sucesfully.out'
    output: '/tmp/probando.out'
    params: workflow.basedir
    shell: 'echo {params}'
  
    
    
    
    
    
    