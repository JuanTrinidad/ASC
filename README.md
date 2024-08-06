# Snakemake workflow

Welcome to the ASC (Annotation by Structural Comparisons) Pipeline!

This repository contains the implementation of the ASC pipeline, which leverages structural comparisons for genome annotation. The workflow is built using the Snakemake workflow management system, providing an automated and reproducible approach to bioinformatics analyses.

### Repository Contents:
* Workflow directory: The main execution file "snakemake_ASC.smk" and all necessary scripts to run the pipeline.
* Config Directory: Contains the config.yaml file where you can adjust parameters according to your requirements.


### Prerequisites:
Input Files:
* A FASTA file containing the protein sequences to be annotated.
* A TSV file linking the FASTA headers with UniProt accessions for AFDB protein structure download.