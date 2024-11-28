




########################################
#
# PROTEOMES ANNOTATION DATA FROM UNIPROT
#
########################################

rule download_proteomes_annotation_info:
  input: '../config/mandatory_files/link_to_download_from_AFDB/{UP}_{num}_{spp}_{v}.tar'
  output: 'genome_data_sets/subject_proteomes/annotation_info/{UP}_{num}_{spp}_{v}.tsv'
  shell: 
    'curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Corganism_id%2Cxref_proteomes%2Clineage%2Cec%2Ccc_function%2Ccc_pathway%2Creviewed%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id%2Ccc_domain%2Cft_domain%2Cft_motif%2Cprotein_families%2Cxref_kegg%2Cxref_orthodb%2Cxref_eggnog%2Cxref_panther%2Cxref_interpro%2Cxref_pfam&format=tsv&query=%28{wildcards.UP}%29" | gzip -c -d > {output}'


rule concat_TSV_in_one_file:
  input: expand('genome_data_sets/subject_proteomes/annotation_info/{organisms}.tsv', organisms = [file[:-4] for file in model_organisms_files_final])
  output: 'protein_data_bases/annotation_info/Gene_annotation_info_from_uniprot_model_spp.tsv'
  run:
    pd.concat([pd.read_csv(file, sep="\t",low_memory=False) for file in input]).to_csv(output[0], sep='\t')
    

#########################################
#
# DATA BASES ANNOTATION DATA FROM UNIPROT
#
#########################################

rule download_kinetoplastea_annotation_info:
  input: 'tmp/TriTrypDB-65_All_species_clean_download_finished.out' #modificar
  output: tsv = 'protein_data_bases/annotation_info/kinetoplastea_taxid5653_annotation_info_from_uniprot.tsv'
  shell: 
    'curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Corganism_id%2Cxref_proteomes%2Clineage%2Cec%2Ccc_function%2Ccc_pathway%2Creviewed%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id%2Ccc_domain%2Cft_domain%2Cft_motif%2Cprotein_families%2Cxref_kegg%2Cxref_orthodb%2Cxref_eggnog%2Cxref_panther%2Cxref_interpro%2Cxref_pfam&format=tsv&query=%28%28taxonomy_id%3A5653%29%29" | gzip -c -d > {output.tsv}' 




#########################################
#
# FINAL TABLE
#
#########################################


rule final_table:
  input: 
    SRBH = '../results/reciprocal_best_hit_SingleOrgApproach_TSV/rbh_all_in_one_file.tsv',
    FATCATTMscore = 'report/toy_file_TMscores_from_TMalign_twisted_structure.tsv',
    annotation = 'protein_data_bases/annotation_info/Gene_annotation_info_from_uniprot_model_spp.tsv',
  output:
    '../results/FINAL_TABLE.tsv'
  script:
    '../scripts/014_creating_final_table.py'





################################
#
# ALL ANNOTATION DATA FROM UNIPROT
#
#################################

#esta regla esta funcional pero tarda mucho la corro una vez y despues la comento para que no me rompa las pelotas
#rule download_annotation_data_from_UNIPROT:
#  output: 'protein_data_bases/annotation_info/ALL_UNIPROT_ANNOTATION.gz'
#  conda:
#    'envs/env_nameconverter_and_AFDBdownload.yaml'  
#  script:
#    'scripts/7_download_from_uniprot.py'
  

#rule unzip_annotation_data_from_UNIPROT:
#  input:  'protein_data_bases/annotation_info/ALL_UNIPROT_ANNOTATION.gz'
#  output: 'protein_data_bases/annotation_info/ALL_UNIPROT_ANNOTATION.tsv'
#  shell:
#    'gzip -c -d {input} > {output}'




#rule download_SwissProt_annotation_info:
#  input: 'tmp/download_finished.out' #modificar
#  output: tsv = 'protein_data_bases/annotation_info/Swiss-Prot_annotation_data.tsv'
#  shell: 
#    'curl -H "Accept: text/plain; format=tsv" "https://rest.uniprot.org/uniprotkb/stream?compressed=true&fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Corganism_id%2Cxref_proteomes%2Clineage%2Cec%2Ccc_function%2Ccc_pathway%2Creviewed%2Cgo_p%2Cgo_c%2Cgo%2Cgo_f%2Cgo_id%2Ccc_domain%2Cft_domain%2Cft_motif%2Cprotein_families%2Cxref_kegg%2Cxref_orthodb%2Cxref_eggnog%2Cxref_panther%2Cxref_interpro%2Cxref_pfam&format=tsv&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29" | gzip -c -d > {output.tsv}' 
 



