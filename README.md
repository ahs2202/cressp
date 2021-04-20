# Cross-Reactive-Epitope-Search-using-Structural-Properties-of-proteins (CRESSP)
 a program to find cross-reactive epitopes with structural information from known protein structures.



## Installation 

CRESSP requires BLASTP from BLAST+ and HMMER3, which can be installed through conda with the following commands:

```sh
    conda install -c bioconda blast
    conda install -c bioconda hmmer
```

CRESSP can be installed from PyPI (https://pypi.org/project/cressp/)

```sh
    pip install cressp
```





## Usage

```sh
usage: cressp [-h] [-t DIR_FILE_PROTEIN_TARGET] [-q DIR_FILE_PROTEIN_QUERY] [-o DIR_FOLDER_OUTPUT] [-c CPU] [-w WINDOW_SIZE] [-s FLOAT_THRES_AVG_SCORE_BLOSUM_WEIGHTED] [-e FLOAT_THRES_E_VALUE] [-H]
              [-d DIR_FILE_QUERY_HMMDB] [-Q]


arguments:
  -h, --help            show this help message and exit
  -t DIR_FILE_PROTEIN_TARGET, --dir_file_protein_target DIR_FILE_PROTEIN_TARGET
                        (Required) an input FASTA file containing target protein sequences.
  -q DIR_FILE_PROTEIN_QUERY, --dir_file_protein_query DIR_FILE_PROTEIN_QUERY
                        (Default: UniProt human proteins) an input FASTA file containing query protein sequences.
  -o DIR_FOLDER_OUTPUT, --dir_folder_output DIR_FOLDER_OUTPUT
                        (Default: a subdirectory of the current directory) an output directory
  -c CPU, --cpu CPU     (Default: 1) Number of logical CPUs (threads) to use in the current compute node.
  -w WINDOW_SIZE, --window_size WINDOW_SIZE
                        (Default: 30) list of window sizes separated by comma. Example: 15,30,45
  -s FLOAT_THRES_AVG_SCORE_BLOSUM_WEIGHTED, --float_thres_avg_score_blosum_weighted FLOAT_THRES_AVG_SCORE_BLOSUM_WEIGHTED
                        (Default: 0.15) threshold for average weighted BLOSOM62 alignment score for filtering aligned sequences
  -e FLOAT_THRES_E_VALUE, --float_thres_e_value FLOAT_THRES_E_VALUE
                        (Default: 1e-20) threshold for the global alignment e-value in a scientific notation Example: 1e-3
  -H, --flag_use_HMM_search
                        (Default: False) Set this flag to perform HMM search in addition to BLASTP search. HMM profile search is performed with HMMER3. The search usually takes several hours for
                        metagenome-assembled genomes
  -d DIR_FILE_QUERY_HMMDB, --dir_file_query_hmmdb DIR_FILE_QUERY_HMMDB
                        (Default: a HMM profile database of 1012 human proteins searched against UniProt Pan Proteomes. These proteins consist of experimentally validated human autoantigens) a file
                        containing HMM DB of query proteins aligned against pan-proteomes
  -Q, --flag_skip_struc_prop_for_protein_target
                        (Default: False) Set this flag to skip the estimation of structural properties of target proteins. Only structural properties of query proteins will be used to calculate
                        accessibility-weighted-similarity scores
```





## Tutorial

Download the proteome of SARS-CoV-2 from UniProt ([UP000464024](https://www.uniprot.org/proteomes/UP000464024)) as a fasta sequence



run the following command to 

1. Search SARS-CoV-2 protein sequences with human protein sequences with BLASTP and HMMER3, and
2. Calculate similarity scores based on relative-surface-availability of residues of human proteins

```sh
    cressp -t UP000464024_2697049.fasta.gz --flag_use_HMM_search --FLOAT_THRES_E_VALUE 5e-2 --cpu 2
```



