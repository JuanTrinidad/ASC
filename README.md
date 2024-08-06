# ASC Snakemake workflow

Welcome to the ASC (Annotation by Structural Comparisons) Pipeline!

This repository contains the implementation of the ASC pipeline, which leverages structural comparisons for genome annotation. The workflow is built using the Snakemake workflow management system, providing an automated and reproducible approach to bioinformatics analyses.

### Repository Contents:
* Workflow directory: The main execution file "snakemake_ASC.smk" and all necessary scripts to run the pipeline.
* Config Directory: Contains the config.yaml file where you can adjust parameters according to your requirements.


### Prerequisites:
Input Files:
* A FASTA file containing the protein sequences to be annotated.
* A TSV file linking the FASTA headers with UniProt accessions for AFDB protein structure download.



### Pipeline Overview:
* Clustering Sequences: 
    * Protein sequences are clustered using MMseqs2, with user-defined parameters for the minimum number of proteins per cluster. 
    * Clusters are filtered based on the number of sequences.

Structure Download:
The pipeline downloads all available structures for the sequences within the clusters from the AlphaFold Database (AFDB) via FTP.
Inconsistencies between sequences from different databases are prevented by comparing each protein structure with the sequences in the FASTA file.
The best structure prediction for each cluster is chosen based on the average pLDDT.
Database Creation:

A database for the Foldseek program is created using the representative structures.
Clusters without available protein structures are excluded from the analysis.
Model Organism Databases:

Snakemake downloads all available model organisms from AFDB and creates independent databases for each organism.
Structural Comparisons:

Using Foldseek reciprocal-best-hit (SRBH), comparisons are performed between the query database and each target database.
The pipeline selects the top SRBH results as specified by the user and performs flexible structural alignments using the FATCAT algorithm.
TM-scores are calculated for the aligned structures using the TM-align algorithm.