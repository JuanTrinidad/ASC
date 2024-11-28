






###############################
#
#
#
################################


rule downloading_and_coping_files_to_FACTCATfolder:
    input: 
        file1 = '../results/reciprocal_best_hit_SingleOrgApproach_TSV/{all_sequence_fasta}_rbh_all_in_one_file.tsv',
    output: 
        tsv_to_match_pdbs_names = 'tmp/{all_sequence_fasta}_query_taget_accesion_to_fatcat_list.tsv'#, directory('tmp/FATCAT_pdb_files/')
    conda: 
        '../envs/env_pLDDT_mean_calc.yaml'
    params: 
        destination_dir = 'tmp/FATCAT_pdb_files/',
        top_hits = config['FATCAT_top_N_hits_to_align']['top_N_hits'],
        prefix_of_added_pdbs = config['structure_from_outside_AFDB']['prefix_in_PDB_name']
    threads: workflow.cores 
    script: '../scripts/007_downloading_pdb_from_AFDB.py'


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



#################
# USING FATCAT
#################

rule FATCAT_aligment:
    input:
        'tmp/{all_sequence_fasta}_query_taget_accesion_to_fatcat_list.tsv',
        'tmp/FATCAT_installed_sucesfully.out'
    output:
        touch('tmp/{all_sequence_fasta}_FATCAT_aligment.out')
    conda: '../envs/env_pLDDT_mean_calc.yaml'
    threads: workflow.cores
    script:
        '../scripts/009_test_FATCAT_aligment.py'



rule extract_twisted_structures:
    input: 'tmp/{all_sequence_fasta}_FATCAT_aligment.out'
    output: touch('tmp/{all_sequence_fasta}_extract_twisted_structures.out')
    conda: '../envs/env_pLDDT_mean_calc.yaml'
    script: '../scripts/011_prot_structure_sup_to_twisted2pdbs.py'



rule TMalign_twisted_structures:
    input: 'tmp/{all_sequence_fasta}_extract_twisted_structures.out'
    output: touch('tmp/{all_sequence_fasta}_TMalign_twisted_structures.out')
    conda: '../envs/env_TMalign.yaml'
    threads: workflow.cores
    script: '../scripts/012_FATCAT_twisted_pdbs_TMalign.py'



rule extract_TMscore_from_TMalign_results:
    input: 'tmp/{all_sequence_fasta}_TMalign_twisted_structures.out'
    output: 'report/{all_sequence_fasta}_TMscores_from_TMalign_twisted_structure.tsv'
    conda: '../envs/env_pLDDT_mean_calc.yaml'
    script: '../scripts/013_extract_TMscore_TMalign.py'

















rule run_FATCAT:
    input: 'tmp/FATCAT_installed_sucesfully.out'
    output: '/tmp/probando.out'
    params: workflow.basedir
    shell: 'echo {params}'
  
    
    
    
    
    
### solucion a cuando hay que levantar archivos sin tenerlos facil prar imput

def get_files(wildcards):
    return glob.glob("tmp/FATCAT_pdb_files/*.pdb")


rule use_files:
    input:
        get_files
    output: touch('tmp/using_files.out')
    shell:
        """
        echo {input}
        """


