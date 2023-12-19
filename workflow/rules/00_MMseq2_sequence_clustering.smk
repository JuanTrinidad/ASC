##################################################
#
#
#                     MMseq2
# 
# 
##################################################


####################################################################
# creating query database to cluster from initial protein fasta file
####################################################################

rule mmseqs2_DB_query_protein_sequences:
    input: 
        'genome_data_sets/query_proteomes/fasta_files/{all_sequence_fasta}.fa' #config['input_files']['all_sequence_fasta']
    output: 
        'genome_data_sets/query_proteomes/mmseqs2_data_base/{all_sequence_fasta}.db'
    conda: 
        '../envs/env_mmseq2.yaml'
    shell: 
        'mmseqs createdb {input} {output}' 


###############################
# clustering the query database
###############################
        
rule mmseqs2_clustering_query_protein_sequences:
    input:
        'genome_data_sets/query_proteomes/mmseqs2_data_base/{all_sequence_fasta}.db' #'genome_data_sets/query_proteomes/mmseqs2_data_base/' + initial_fasta_file_name + '.db'
    output:
        multiext('genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}', '.index', '.dbtype')
    conda:
        '../envs/env_mmseq2.yaml'
    params:
        mmseqs2_cluster_module = config['mmseq2_query_sequences_clustering']['mmseqs2_cluster_module'],
        cluster_mode = config['mmseq2_query_sequences_clustering']['cluster_mode'],
        similarity_type = config['mmseq2_query_sequences_clustering']['similarity_type'],
        minseqID = config['mmseq2_query_sequences_clustering']['minimum_sequence_identity'],
        mincoverage = config['mmseq2_query_sequences_clustering']['minimum_coverage'],
        coveragemode = config['mmseq2_query_sequences_clustering']['coverage_mode'],
        output_path = 'genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}',
        evalue = config['mmseq2_query_sequences_clustering']['evalue']
        
    shell:
        'mmseqs {params.mmseqs2_cluster_module} {input} {params.output_path} /tmp --cluster-mode {params.cluster_mode} --similarity-type {params.similarity_type} --min-seq-id {params.minseqID} -c {params.mincoverage} --cov-mode {params.coveragemode} -e {params.evalue}'
        

##################################
# extract TSV file from clustering
##################################

rule extract_TSV_from_mmseqs2_clustering_query_protein_sequences:
  input: 
    DB = 'genome_data_sets/query_proteomes/mmseqs2_data_base/{all_sequence_fasta}.db', 
    DB_clu = multiext('genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}', '.index', '.dbtype')
  output: 
    'genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}.tsv'
  conda:
    '../envs/env_mmseq2.yaml'
  params:
    output_path = 'genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}'
  shell: 
    'mmseqs createtsv {input.DB} {input.DB} {params.output_path} {output}'


################################
#  Extracting important clusters
################################
# the idea is to remove not relevant clusters, such us 1 member clusters

rule mmseqs2_cluster_selection_to_annotate:
    input: 
      'genome_data_sets/query_proteomes/mmseqs2_seq_clusters/{all_sequence_fasta}.tsv'
    output:
      'report/{all_sequence_fasta}_Ortholog_group_to_geneID.tsv'
    conda:
      '../envs/env_pLDDT_mean_calc.yaml'
    params:
      cluster_size = config['cluster_filtering']['cluster_size']
    script:
      '../scripts/000_selecting_relevant_clusters.py'


