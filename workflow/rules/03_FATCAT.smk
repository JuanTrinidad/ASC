


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
        file1 = 'report/{all_sequence_fasta}_Ortholog_Group_Full_of_hypotetical_boolean.tsv',
        file2 = '../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv',
        file3 = '../results/Gene_annotation_info_from_uniprot_model_spp.tsv',
        file4 = 'report/{all_sequence_fasta}_ortholog_groups_x_sequence_clustering_x_UNIPROT.tsv',
        dirs = expand('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/', proteome = model_organisms_files_final_db),
        files_directories = expand('genome_data_sets/subject_proteomes/pdb_files/model_organisms_files_unzip/{proteome}/', proteome= model_organisms_files_final_db)  
    output: 
        tsv_to_match_pdbs_names = 'tmp/{all_sequence_fasta}_query_taget_accesion_to_fatcat_list.tsv'#, directory('tmp/FATCAT_pdb_files/')
    conda: 
        '../envs/env_pLDDT_mean_calc.yaml'
    params: 
        destination_dir = 'tmp/FATCAT_pdb_files/',
        top_hits = 5,
        prefix_of_added_pdbs = config['structure_from_outside_AFDB']['prefix_in_PDB_name']
    threads: workflow.cores 
    script: '../scripts/007_polishing_pdbfiles_to_use_in_FATCAT.py'
        



#ver de usar la funcion para levantar archivos y que esto sea mas rapdio
def get_gz_files(wildcards):
    return glob.glob("tmp/FATCAT_pdb_files/*.gz")

#da algunos problemas cuando se corre por segunda vez. Esto se debe a no usar correctamente snakemake. Despues solucionar, ej. usando temp en la regla anterior.
rule unzip_pdb_files:
    input:
        pdb_gz = get_gz_files, 
        tsv_to_match_pdbs_names = 'tmp/{all_sequence_fasta}_query_taget_accesion_to_fatcat_list.tsv'
    output:
        touch('tmp/{all_sequence_fasta}_unziped_pdbcif_files.out')
    #params: files = [file for file in glob.glob('tmp/FATCAT_pdb_files/*.gz')]
    run:
        for file in input.pdb_gz:
            shell("gunzip {file}")




rule cif_to_pdb:
    input:
        'tmp/{all_sequence_fasta}_unziped_pdbcif_files.out'
    output:
        touch('tmp/{all_sequence_fasta}_cif_to_pdb.out')
    conda:
        '../envs/env_cif2pdb.yaml'
    threads: workflow.cores
    script:
        '../scripts/008_cif2pdb.py'



rule FATAT_aligment:
    input:
        'tmp/{all_sequence_fasta}_query_taget_accesion_to_fatcat_list.tsv'
    output:
        touch('tmp/{all_sequence_fasta}_FATCAT_aligment.out')
    conda: '../envs/env_pLDDT_mean_calc.yaml'
    threads: workflow.cores
    script:
        '../scripts/009_FATCAT_aligment.py'



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