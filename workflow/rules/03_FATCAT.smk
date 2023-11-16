





rule unzip_uniprot_proteomes:
    input: 
        'genome_data_sets/subject_proteomes/pdb_files/model_organisms_files/{proteome}.tar'
    output: 
        temp(directory('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/'))
    shell: 
        '''
        echo unziping {wildcards.proteome}
        mkdir -p {output}
        tar -xf {input} -C {output}
        cd {output}
        find . -mindepth 2 -type f -exec mv -f -- {{}} . \;
        '''


rule moving_files_from_clusters_to_FACTCATfolder:
    input: 
        file1 = 'report/TriTrypDB-65_All_species_clean_Ortholog_Group_Full_of_hypotetical_boolean.tsv',
        file2 = '../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv',
        file3 = '../results/Gene_annotation_info_from_uniprot_model_spp.tsv',
        file4 = 'report/TriTrypDB-65_All_species_clean_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
        dirs = expand('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/', proteome = model_organisms_files_final_db),
        files_directories = expand('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/', proteome= model_organisms_files_final_db)
    output: 
        list_file_to_fatcat = 'tmp/to_compare.list',
        tsv_to_match_pdbs_names = 'tmp/query_taget_accesion_to_fatcat_list.tsv'
    conda: 
        '../envs/env_pLDDT_mean_calc.yaml'
    params: 
        destination_dir = 'tmp/FATCAT_pdb_files/',
        top_hits = 5,
        prefix_of_added_pdbs = config['structure_from_outside_AFDB']['prefix_in_PDB_name']
    threads: workflow.cores 
    script: '../scripts/007_polishing_pdbfiles_to_use_in_FATCAT.py'
        

rule unzip_pdb_files:
    input:
        'tmp/query_taget_accesion_to_fatcat_list.tsv',
        'tmp/FATCAT_pdb_files/{file}.gz'
    output:
        'tmp/FATCAT_pdb_files/{file}'
    shell:
        'gunzip {input[1]}'



#primero correr el script 007
#rule unzip_pdb_files:
#despues convertir cif a pdb





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
  
    
    
    
    
    
    