#!/usr/bin/env python


from biobookshelf.main import *
from biobookshelf import *


from cressp.structural_property_estimation import Estimate_structural_property
from cressp.cross_reactivity_prediction import Predict_T_cell_cross_reactivity, Predict_B_cell_cross_reactivity
from cressp.web_application import Prepare_data_for_web_application

import argparse
import traceback
import os, sys, getopt
from io import StringIO
import time
import math

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining
import warnings
warnings.filterwarnings( action = 'ignore' )
# from scarab.quality_control import Estimate_structural_property

# retrieve a function for logging
Info = multiprocessing.get_logger( ).info

#     # read dict_blosum62 from the tsv file
#     df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
#     dict_blosum62 = dict( )
#     for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
#         dict_blosum62[ aa_0, aa_1 ] = score

def cressp( dir_file_protein_target = None, dir_file_protein_query = 'human', dir_folder_output = 'default', n_threads = 1, l_window_size = [ 30 ], float_thres_e_value = 30, flag_use_HMM_search = False, dir_file_query_hmmdb = 'human', flag_use_rcsb_pdb_only = False, int_number_of_proteins_in_a_batch_during_dnn_prediction = 1000, flag_only_use_structural_properties_of_query_proteins = False, float_thres_avg_score_blosum_weighted__b_cell = 0.15, float_thres_avg_score_blosum__b_cell = 0.0, float_thres_rsa_correlation = 0.0, float_thres_avg_blosum62_score_for_mhc = 2, float_thres_min_mhc_allele_frequency = 0.5, float_thres_binding_affinities_in_nM = 500, flag_replace_unconventional_acid_code = False, flag_use_all_gpu_devices = False, flag_deduplicate_based_on_aligned_subsequences_for_visualization = False ) :
    """
    The main function of Cross-Reactive-Epitope-Search-using-Structural-Properties-of-proteins (CRESSP)
    
    
    
    dir_file_protein_target = None : (Required) an input FASTA file containing target protein sequences.
            
    dir_file_protein_query = 'human' : (Default: UniProt human proteins) an input FASTA file containing query protein sequences.
            
    dir_folder_output = 'default' : (Default: a subdirectory 'cressp_out/' of the current directory) an output directory.
            
    n_threads = 1 : (Default: 1) Number of logical CPUs (threads) to use in the current compute node.
            
    l_window_size = [ 30 ] : list of window sizes for searching cross-reactive epitopes from the alignments between query and target proteins.
            
    float_thres_e_value = 30 : (Default: 30, meaning no filtering) threshold for the global alignment e-value.
            
    flag_use_HMM_search = False : (Default: False) Set this flag to perform HMM search in addition to BLASTP search. HMM profile search is performed with HMMER3. The search usually takes several hours for metagenome-assembled genomes.
            
    dir_file_query_hmmdb = 'human' : (Default: a HMM profile database of 1012 human proteins searched against UniProt Pan Proteomes. These proteins consist of experimentally validated human autoantigens) a file containing HMM DB of query proteins aligned against pan-proteomes.
            
    flag_use_rcsb_pdb_only = False : (Default: False) When calculating consensus structural properties of input proteins, do not use homology-based modeled structures from SWISS-MODEL repositories and only use experimental protein structures from RCSB PDB.
    
    int_number_of_proteins_in_a_batch_during_dnn_prediction = 1000 : (Default: 1000) number of proteins in a batch when predicting structural properties of protein residues not covered by known/predict structures.
            
    flag_only_use_structural_properties_of_query_proteins = False : (Default: False) Only use estimated structural properties of the query proteins when predicting cross-reactivity between target and query proteins (skip the estimation of structural properties of target proteins). When 'dir_file_protein_query' == 'human' (default value), it will significantly reduce computation time by skipping estimation and prediction steps of structural properties of target and query proteins.
            
    float_thres_avg_score_blosum_weighted__b_cell = 0.15 : (Default: 0.15) threshold for average weighted BLOSOM62 alignment score for filtering predicted cross-reactive epitopes.
            
    float_thres_avg_score_blosum__b_cell = 0.0 : (Default: 0) threshold for average BLOSOM62 alignment score for filtering predicted cross-reactive epitopes.
            
    float_thres_rsa_correlation = 0 : (Default: 0) threshold for correlation coefficient of Relative Surface Area (RSA) values for filtering predicted cross-reactive epitopes.
    
    float_thres_avg_blosum62_score_for_mhc = 2.0 : (Default: 2.0) threshold for average BLOSOM62 alignment score for filtering predicted cross-reactive T-cell epitopes (cross-reactive MHC epitopes)
    
    float_thres_min_mhc_allele_frequency = 0.5 : (Default: 0.5) a threshold for filtering out MHC alleles with low allele frequencies. MHC alleles with allele frequency above the threshold for at least one population will be used for cross-reactive T-cell epitope prediction
            
    float_thres_binding_affinities_in_nM = 500 : (Default: 500) a threshold predicted IC50 values for filtering predicted T-cell cross-reactive epitopes. A pair of peptides whose geometric average of predicted binding affinities (IC50) values above this threshold will be removed.
    
    flag_deduplicate_based_on_aligned_subsequences_for_visualization = False : (Default: False) perform the deduplication step for removing redundant predicted cross-reactive epitopes based on aligned subsequences
            
    """
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package


    """
    Parse arguments
    """
    # define logging levels
    logger = multiprocessing.log_to_stderr( )
    logger.setLevel( logging.INFO )
    Info( 'CRESSP: comparative analysis of two sets of proteomes for searching potentially cross-reactive B-cell and T-cell epitopes' )
    Info( '[Pipeline Start] Pipeline Started at ' + TIME_GET_timestamp( True ) )
    
    # fixed arguments
    float_search_thres_e_value = 30 # set e-value threshold for search (fixed for maximum sensitivity)

    
    ''' check whether the program is called from the command-line interface or from an interactive Python programming environment '''
    str_name_program = sys.argv[ 0 ]
    if '/' in str_name_program :
        str_name_program = str_name_program.rsplit( '/', 1 )[ 1 ]
    flag_usage_from_command_line_interface = str_name_program[ : len( 'cressp' ) ] == 'cressp'
    
    ''' parse arguments when the function was called from the command-line interface '''
    if flag_usage_from_command_line_interface :
        # command line arguments
        parser = argparse.ArgumentParser( description = "Cross-Reactive-Epitope-Search-using-Structural-Properties-of-proteins (cressp): a program to find cross-reactive epitopes with structural information from known protein structures." )
        parser.add_argument( "-t", "--dir_file_protein_target", help = "(Required) an input FASTA file containing target protein sequences." )
        parser.add_argument( "-q", "--dir_file_protein_query", help = "(Default: UniProt human proteins) an input FASTA file containing query protein sequences.", default = 'human' )
        parser.add_argument( "-o", "--dir_folder_output", help = "(Default: a subdirectory 'cressp_out/' of the current directory) an output directory", default = "default" )
        parser.add_argument( "-c", "--cpu", help = "(Default: 1) Number of logical CPUs (threads) to use in the current compute node.", default = '1' )
        parser.add_argument( "-w", "--window_size", help = "(Default: 30) list of window sizes separated by comma. Example: 15,30,45", default = "30" )
        parser.add_argument( "-e", "--float_thres_e_value", help = "(Default: 30, meaning no filtering) threshold for the global alignment e-value in a scientific notation Example: 1e-3", default = "30" )
        parser.add_argument( "-H", "--flag_use_HMM_search", help = "(Default: False) Set this flag to perform HMM search in addition to BLASTP search. HMM profile search is performed with HMMER3. The search usually takes several hours for metagenome-assembled genomes", action = 'store_true' )
        parser.add_argument( "-d", "--dir_file_query_hmmdb", help = "(Default: a HMM profile database of 1012 human proteins searched against UniProt Pan Proteomes. These proteins consist of experimentally validated human autoantigens) a file containing HMM DB of query proteins aligned against pan-proteomes", default = "human" )
        parser.add_argument( "-R", "--flag_use_rcsb_pdb_only", help = "(Default: False) When calculating consensus structural properties of input proteins, do not use homology-based modeled structures from SWISS-MODEL repositories and only use experimental protein structures from RCSB PDB", action = 'store_true' )
        parser.add_argument( "-B", "--int_number_of_proteins_in_a_batch_during_dnn_prediction", help = "(Default: 1000) number of proteins in a batch when predicting structural properties of protein residues not covered by known/predict structures.", default = '1000' )
        parser.add_argument( "-Q", "--flag_only_use_structural_properties_of_query_proteins", help = "(Default: False) Only use estimated structural properties of the query proteins when predicting cross-reactivity between target and query proteins (skip the estimation of structural properties of target proteins). When 'dir_file_protein_query' == 'human' (default value), it will significantly reduce computation time by skipping estimation and prediction steps of structural properties of target and query proteins.", action = 'store_true' )
        # for filtering predicted cross-reactive epitopes
        parser.add_argument( "-s", "--float_thres_avg_score_blosum_weighted__b_cell", help = "(Default: 0.15) threshold for average rsa-weighted BLOSOM62 alignment score for filtering predicted cross-reactive b-cell epitopes", default = '0.15' )
        parser.add_argument( "-S", "--float_thres_avg_score_blosum__b_cell", help = "(Default: 0) threshold for average BLOSOM62 alignment score for filtering predicted cross-reactive b-cell epitopes", default = '0.0' )
        parser.add_argument( "-C", "--float_thres_rsa_correlation", help = "(Default: 0) threshold for correlation coefficient of Relative Surface Area (RSA) values for filtering predicted cross-reactive b-cell epitopes", default = '0.0' )
        parser.add_argument( "-b", "--float_thres_avg_blosum62_score_for_mhc", help = "(Default: 2.0) threshold for average BLOSOM62 alignment score for filtering predicted cross-reactive T-cell epitopes (cross-reactive MHC epitopes)", default = '2.0' )
        parser.add_argument( "-m", "--float_thres_min_mhc_allele_frequency", help = "(Default: 0.5) a threshold for filtering out MHC alleles with low allele frequencies. MHC alleles with allele frequency above the threshold for at least one population will be used for cross-reactive T-cell epitope prediction", default = '0.5' )
        parser.add_argument( "-a", "--float_thres_binding_affinities_in_nM", help = "(Default: 500) a threshold predicted IC50 values for filtering predicted T-cell cross-reactive epitopes. A pair of peptides whose geometric average of predicted binding affinities (IC50) values above this threshold will be removed.", default = '500' )
        parser.add_argument( "-U", "--flag_replace_unconventional_acid_code", help = "(Default: False) If this flag is set, unconventional amino acids in the input protein sequences will be replaced with chemically similar amino acid. Specifically, Selenocysteine (U) to Cysteine (C), Pyrrolysine (O) to Tyrosine (Y)", action = 'store_true' )
        parser.add_argument( "-G", "--flag_use_all_gpu_devices", help = "(Default: False) Use all available GPU devices during prediction of RSA values. When this flag is set to True, the RSA prediction might be completed faster, but all GPU memories will be occupied by Tensorflow", action = 'store_true' )
        parser.add_argument( "-D", "--flag_deduplicate_based_on_aligned_subsequences_for_visualization", help = "(Default: False) perform the deduplication step for removing redundant predicted cross-reactive epitopes based on aligned subsequences", action = 'store_true' )
        args = parser.parse_args( )
        if args.dir_file_protein_target is None :
            Info( '[Error] Required argument(s) is missing. to view help message, use -h or --help flag' )
            sys.exit( )

        # parse arguments # no further processing required
        flag_use_HMM_search = args.flag_use_HMM_search
        flag_use_rcsb_pdb_only = args.flag_use_rcsb_pdb_only
        flag_only_use_structural_properties_of_query_proteins = args.flag_only_use_structural_properties_of_query_proteins
        flag_use_all_gpu_devices = args.flag_use_all_gpu_devices
        flag_replace_unconventional_acid_code = args.flag_replace_unconventional_acid_code
        int_number_of_proteins_in_a_batch_during_dnn_prediction = int( args.int_number_of_proteins_in_a_batch_during_dnn_prediction )
        n_threads = min( int( args.cpu ), int( multiprocessing.cpu_count( ) ) ) # max number of n_threads is the CPU number of the current local machine
        float_thres_avg_blosum62_score_for_mhc = float( args.float_thres_avg_blosum62_score_for_mhc )
        float_thres_min_mhc_allele_frequency = float( args.float_thres_min_mhc_allele_frequency )
        l_window_size = list( int( e ) for e in args.window_size.split( ',' ) ) # set default window size
        float_thres_e_value = float( args.float_thres_e_value )
        float_thres_avg_score_blosum_weighted__b_cell = float( args.float_thres_avg_score_blosum_weighted__b_cell ) 
        float_thres_avg_score_blosum__b_cell = float( args.float_thres_avg_score_blosum__b_cell ) 
        float_thres_rsa_correlation = float( args.float_thres_rsa_correlation ) 
        float_thres_binding_affinities_in_nM = float( args.float_thres_binding_affinities_in_nM )

        # parse directory arguments
#         dir_file_protein_target_representative
        dir_file_protein_target = args.dir_file_protein_target
        dir_file_protein_query = args.dir_file_protein_query
        dir_folder_output = args.dir_folder_output
        dir_file_query_hmmdb = args.dir_file_query_hmmdb
    else :
        ''' parse arguments when the function was called from an interactive Python interpreter '''
        if dir_file_protein_target is None :
            Info( "[Error] required input 'dir_file_protein_target' was not given" )
            if flag_usage_from_command_line_interface : sys.exit( )
            else : return - 1

    # handle default settings for input datafiles
    flag_default_protein_query_was_used = False # a flag indicating whether the default protein_query was used
    if dir_file_protein_query == 'human' :
        flag_default_protein_query_was_used = True
        PKG.Download_Data( "data/human/uniprot.tsv.gz", dir_remote, name_package ) # download data
        dir_file_protein_query = f'{dir_folder_cressp}data/human/human_uniprot.fa' # set default 'dir_file_protein_query' # file_name will be used to refer to the given query protein, and thus 'human_uniprot' is used as a file_name.
        if not os.path.exists( dir_file_protein_query ) : # if default human protein fasta file is not available, extract the sequence from the dataframe
            df_protein_query = pd.read_csv( f'{dir_folder_cressp}data/human/uniprot.tsv.gz', sep = '\t' )
            dict_fasta_protein_query = df_protein_query.set_index( 'fasta_header' ).seq.to_dict( )
            FASTA_Write( dir_file_protein_query, dict_fasta = dict_fasta_protein_query )

        if flag_use_HMM_search : 
            if dir_file_query_hmmdb == 'human' :
                PKG.Download_Data( "data/human/hmmdb_autoantigen.hmm.gz", dir_remote, name_package ) # download data
                PKG.Gunzip_Data( "data/human/hmmdb_autoantigen.hmm.gz", name_package ) # unzip data
                dir_file_query_hmmdb = pkg_resources.resource_filename( name_package, 'data/human/hmmdb_autoantigen.hmm' ) # set default 'dir_file_query_hmmdb'
    else :
        if flag_use_HMM_search and dir_file_query_hmmdb == 'human' :
            Info( "[Error] exiting since query proteins other than default human proteins were given, the default HMM profile database cannot be used." )
            if flag_usage_from_command_line_interface : sys.exit( )
            else : return - 1
    # get absolute paths for file arguments        
    dir_file_protein_target = os.path.abspath( dir_file_protein_target )
    dir_file_protein_query = os.path.abspath( dir_file_protein_query )
    dir_file_query_hmmdb = os.path.abspath( dir_file_query_hmmdb )

    # handle default setting for the output folder
    if dir_folder_output == 'default' :
        dir_folder_output = f"{os.getcwd( )}/cressp_out/" # set default output folder
    # get absolute paths for folder arguments        
    dir_folder_output = os.path.abspath( dir_folder_output )

    # handle output folder
    if dir_folder_output[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_output += '/'
    # define folder directories for each task
    dir_folder_pipeline = f"{dir_folder_output}pipeline/"
    dir_folder_pipeline_temp = f'{dir_folder_pipeline}temp/' 
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # create a working directory of estimating structural properties
    dir_folder_pipeline_web = f'{dir_folder_pipeline}web_application/' # a working directory for exporting data for web applications
    dir_folder_web = f'{dir_folder_output}web_application/'
    # create folders
    for dir_folder in [ dir_folder_output, dir_folder_pipeline, dir_folder_pipeline_temp, dir_folder_pipeline_struc, dir_folder_pipeline_web, dir_folder_web ] :
        os.makedirs( dir_folder, exist_ok = True )


    """
    Read and move input protein Fasta files
    """
    def __Get_File_Name_Without_Extension__( dir_file ) :
        ''' get file_name of a given directory to a fasta file (gzipped or uncompressed) and return it '''
        name_file_without_extension = dir_file.rsplit( '/', 1 )[ 1 ]
        # when a file does not have extension
        if '.' not in name_file_without_extension :
            return name_file_without_extension
        name_file_without_extension, str_file_extension = name_file_without_extension.rsplit( '.', 1 )
        # when a file is gzipped, perform extension spliting one more time
        if str_file_extension.lower( ) == 'gz' :
            name_file_without_extension, str_file_extension = name_file_without_extension.rsplit( '.', 1 )
        return name_file_without_extension
    
    # retrieve file names
    name_file_protein_query = __Get_File_Name_Without_Extension__( dir_file_protein_query ) # retrieve file name containing the query proteins
    name_file_protein_target = __Get_File_Name_Without_Extension__( dir_file_protein_target ) # retrieve file name containing the query proteins
    
    """ check flag """
    dir_file_flag = f"{dir_folder_pipeline}copying_input_files_completed.flag"
    if not os.path.exists( dir_file_flag ) :
        
        try :
            dict_fasta_protein_query = FASTA_Read( dir_file_protein_query )
            FASTA_Write( f"{dir_folder_pipeline}protein_query.fasta", dict_fasta = dict_fasta_protein_query )
            dir_file_protein_query = f"{dir_folder_pipeline}protein_query.fasta" # set directory of fasta file to the new file location
        except :
            print( f"exiting due to an error while reading and moving 'dir_file_protein_query' {dir_file_protein_query}" )
            if flag_usage_from_command_line_interface : sys.exit( )
            else : return - 1

        try :
            dict_fasta_protein_target = FASTA_Read( dir_file_protein_target )
            FASTA_Write( f"{dir_folder_pipeline}protein_target.fasta", dict_fasta = dict_fasta_protein_target )
            dir_file_protein_target = f"{dir_folder_pipeline}protein_target.fasta" # set directory of fasta file to the new file location
        except :
            print( f"exiting due to an error while reading and moving 'dir_file_protein_target' {dir_file_protein_target}" )
            if flag_usage_from_command_line_interface : sys.exit( )
            else : return - 1
            
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    else :
        # define directories of input protein sequences and load query protein sequences
        dir_file_protein_query = f"{dir_folder_pipeline}protein_query.fasta" # set directory of fasta file to the new file location
        dir_file_protein_target = f"{dir_folder_pipeline}protein_target.fasta" # set directory of fasta file to the new file location
        dict_fasta_protein_query = FASTA_Read( dir_file_protein_query ) # load protein sequences of query (required during BLASTP output processing)
    Info( "[Task Completion] copying input protein sequence files was completed at " + TIME_GET_timestamp( True ) )
        
    """
    Report external and internal program settings
    """
    # compile setting (metadata)
    dict_cressp_setting = {  
            # program setting
        'n_threads' : n_threads,
        'flag_use_HMM_search' : flag_use_HMM_search,
        'flag_use_rcsb_pdb_only' : flag_use_rcsb_pdb_only,
        'flag_only_use_structural_properties_of_query_proteins' : flag_only_use_structural_properties_of_query_proteins,
        'flag_deduplicate_based_on_aligned_subsequences_for_visualization' : flag_deduplicate_based_on_aligned_subsequences_for_visualization,
        'float_thres_avg_blosum62_score_for_mhc' : float_thres_avg_blosum62_score_for_mhc,
        'float_thres_min_mhc_allele_frequency' : float_thres_min_mhc_allele_frequency,
        'l_window_size' : l_window_size,
        'float_thres_e_value' : float_thres_e_value,
        'float_thres_avg_score_blosum_weighted__b_cell' : float_thres_avg_score_blosum_weighted__b_cell,
        'float_thres_avg_score_blosum__b_cell' : float_thres_avg_score_blosum__b_cell,
        'float_thres_rsa_correlation' : float_thres_rsa_correlation,
        'float_thres_binding_affinities_in_nM' : float_thres_binding_affinities_in_nM,
        'dir_file_protein_target' : dir_file_protein_target,
        'dir_file_protein_query' : dir_file_protein_query,
        'dir_folder_output' : dir_folder_output,
        'dir_file_query_hmmdb' : dir_file_query_hmmdb,
            # internal setting
        'dir_folder_pipeline' : dir_folder_pipeline,
        'dir_folder_pipeline_temp' : dir_folder_pipeline_temp,
        'dir_folder_pipeline_struc' : dir_folder_pipeline_struc,
        'dir_folder_pipeline_web' : dir_folder_pipeline_web,
        'dir_folder_web' : dir_folder_web,
        'name_file_protein_query' : name_file_protein_query,
        'name_file_protein_target' : name_file_protein_target }
    
    Info( f"[Setting] cressp will be run with the following setting: {str( dict_cressp_setting )}" )
    
    ''' export CRESSP setting '''
    dir_file_json_setting_cressp = f"{dir_folder_pipeline}cressp_setting.json"
    if os.path.exists( dir_file_json_setting_cressp ) :
        with open( dir_file_json_setting_cressp, 'r' ) as file :
            j = json.load( file )
        if j != dict_cressp_setting :
            Info( f"[Warning] the current CRESSP setting is different from the previous CRESSP setting recorded in the pipeline folder. The previous setting will be used." )
            with open( dir_file_json_setting_cressp, 'r' ) as file :
                dict_cressp_setting = json.load( file ) # override current CRESSP setting with previous CRESSP setting
    with open( dir_file_json_setting_cressp, 'w' ) as newfile :
        json.dump( dict_cressp_setting, newfile )    
    
    """
    Perform BLASTP alignment
    """
    """ check flag """
    dir_file_flag = f"{dir_folder_pipeline}blastp_completed.flag"
    if not os.path.exists( dir_file_flag ) :

        # create blastp_db using query_protein sequences
        dir_prefix_blastdb_protein_query = f"{dir_folder_pipeline}makeblastdb_out/protein_query"
        os.makedirs( f"{dir_folder_pipeline}makeblastdb_out/", exist_ok = True )  
        OS_Run( [ "makeblastdb", "-in", f"{dir_folder_pipeline}protein_query.fasta", '-dbtype', 'prot', '-parse_seqids', '-max_file_sz', '1GB', '-out', dir_prefix_blastdb_protein_query ], dir_file_stdout = f"{dir_prefix_blastdb_protein_query}.makeblastdb.stdout.txt", dir_file_stderr = f"{dir_prefix_blastdb_protein_query}.makeblastdb.stderr.txt", return_output = False ) # make blast db for protein_query

        # run blastp
        dir_file_blastp_output = f'{dir_folder_pipeline}blastp.tsv'
        OS_Run( [ 'blastp', '-query', dir_file_protein_target, '-db', dir_prefix_blastdb_protein_query, '-out', dir_file_blastp_output, '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop', '-num_threads', f'{n_threads}', '-evalue', f'{float_search_thres_e_value}' ], dir_file_stdout = f"{dir_folder_pipeline}blastp.stdout.txt", dir_file_stderr = f"{dir_folder_pipeline}blastp.stderr.txt", return_output = False ) # run blastp
        OS_Run( [ 'gzip', dir_file_blastp_output ], dir_file_stdout = f"{dir_file_blastp_output}.gzip.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.gzip.stderr.txt", return_output = False ) # compress blastp output
        dir_file_blastp_output += '.gz'
        
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    else :
        # define blastp output file
        dir_file_blastp_output = f'{dir_folder_pipeline}blastp.tsv.gz'
    Info( "[Task Completion] BLASTP search was completed at " + TIME_GET_timestamp( True ) )
            
    """
    Perform HMMER alignment using a given HMM profile DB (using a couple of HMM profile DBs would be also helpful in the near future)
    """
    # run hmmsearch according to 'flag_use_HMM_search' flag
    if flag_use_HMM_search : 
        """ check flag """
        dir_file_flag = f"{dir_folder_pipeline}hmmsearch_completed.flag"
        dir_file_hmmsearch_output = f'{dir_folder_pipeline}hmmsearch.out'
        if not os.path.exists( dir_file_flag ) :
            OS_Run( [ 'hmmsearch', '-o', dir_file_hmmsearch_output, '--acc', '--notextw', '--cpu', f'{n_threads}', '-E', f'{float_search_thres_e_value}', dir_file_query_hmmdb, dir_file_protein_target ], dir_file_stdout = f"{dir_folder_pipeline}hmmsearch.stdout.txt", dir_file_stderr = f"{dir_folder_pipeline}hmmsearch.stderr.txt", return_output = False ) # run hmmsearch
            """ set flag """
            with open( dir_file_flag, 'w' ) as newfile :
                newfile.write( 'completed\n' )
        Info( "[Task Completion] HMMER search was completed at " + TIME_GET_timestamp( True ) )
    
    """
    Combine BLASTP and HMMSEARCH outputs 
    """
    dir_file_matched = f'{dir_folder_pipeline}matched.tsv.gz'
    
    """ check flag """
    dir_file_flag = f"{dir_file_matched}.write_completed.flag"
    if not os.path.exists( dir_file_flag ) :
        # load blastp result
        dict_qacc_to_seq = FASTA_Read( dir_file_protein_target ) # read query protein sequences
        dict_qacc_to_seq = dict( ( header.split( ' ', 1 )[ 0 ], dict_qacc_to_seq[ header ] ) for header in list( dict_qacc_to_seq ) )
        set_acc_target = set( dict_qacc_to_seq ) # retrieve a set of acc_target
        df_blastp = BLAST_Read( dir_file_blastp_output, dict_qaccver_to_seq = dict_qacc_to_seq ) 
        df_blastp = df_blastp[ [ 'saccver', 'qaccver', 'sstart', 'send', 'qstart', 'qend', 'subject_seq_aligned', 'query_seq_aligned', 'evalue', 'pident' ] ]
        df_blastp.pident = df_blastp.pident / 100
        df_blastp.columns = [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'identity' ]
        df_blastp[ 'source' ] = 'blastp'

        # load hmmer result according to 'flag_use_HMM_search' flag
        if flag_use_HMM_search : 
            dict_qacc_to_seq = dict( ( header.split( ' ', 1 )[ 0 ], dict_fasta_protein_query[ header ] ) for header in list( dict_fasta_protein_query ) )
            # ignore if records from profiles that do not match input query protein sequences 
            df = PD_Select( HMMER_HMMSEARCH_Read_output( dir_file_hmmsearch_output ), query_accession = set( dict_qacc_to_seq ) ) 
            l_query_alignment = list( ) # replace query consensus sequence with query sequence
            for query_accession, query_alignment, query_start, query_end in df[ [ 'query_accession', 'query_alignment', 'query_start', 'query_end' ] ].values :
                query_seq = dict_qacc_to_seq[ query_accession ][ query_start - 1 : query_end ] # retrive a subsequence of query sequence
                l_subsequence = list( )
                int_start = 0
                for subsequence in query_alignment.split( '.' ) :
                    int_subsequence_length = len( subsequence )
                    l_subsequence.append( query_seq[ int_start : int_start + int_subsequence_length ] )
                    int_start += int_subsequence_length
                l_query_alignment.append( '-'.join( l_subsequence ) )
            df[ 'query_alignment' ] = l_query_alignment
            df[ 'target_alignment' ] = df.target_alignment.str.upper( ) # target alignment string contains amino acids in lower cases when alignment confidence is low, and it is convenient to convert them to upper characters
            df_hmmer = df
            df_hmmer = df_hmmer[ [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'conditional_Evalue', 'accuracy' ] ] # subset common columns
            df_hmmer.columns = [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'identity' ] # rename columns
            df_hmmer[ 'source' ] = 'hmmer' 
            set_acc_query = set( header.split( ' ', 1 )[ 0 ] for header in dict_fasta_protein_query ) # retrieve a set of acc_query
            df_hmmer = PD_Select( df_hmmer, query_accession = set_acc_query ) # exclude records with accessions that does not exist in the input query protein sequences

        df_matched = pd.concat( [ df_hmmer, df_blastp ], ignore_index = True ) if flag_use_HMM_search else df_blastp
        df_matched.to_csv( dir_file_matched, sep = '\t', index = False )

        """ set flag """
        with open( dir_file_flag, 'w' ) as file :
            file.write( f"search results were written at {TIME_GET_timestamp( )}" )
    Info( "[Task Completion] Combining HMMER and BLASTP search results was completed at " + TIME_GET_timestamp( True ) )
    
    """
    Estimate structural properties of proteins 
    """
    
    """ check flag """
    dir_file_flag = f"{dir_folder_pipeline}rsa_estimation_completed.flag"
    if not os.path.exists( dir_file_flag ) :
        # use previously calculated structural properties when the default query proteins were used
        if flag_default_protein_query_was_used :
            shutil.copyfile( f'{dir_folder_cressp}data/human/uniprot.tsv.gz', f'{dir_folder_pipeline}protein_query.tsv.gz' ) 
            # use previously calculated alignment between human proteins and RCSB_PDB structures
            PKG.Download_Data( "data/human/uniprot.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz", dir_remote, name_package ) # download data (alignment between human proteins and RCSB_PDB structures)
            shutil.copyfile( f'{dir_folder_cressp}data/human/uniprot.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz', f'{dir_folder_pipeline}struc/protein_query.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz' ) 

        # estimate structural properties
        for name_file in [ 'protein_target' ] if flag_default_protein_query_was_used else [ 'protein_target', 'protein_query' ] : # skip prediction of query proteins if default query proteins are used
            """ check flag_2 """
            dir_file_flag_2 = f"{dir_folder_pipeline}{name_file}.tsv.gz.completed.flag"
            if not os.path.exists( dir_file_flag_2 ) :
                Estimate_structural_property( f'{dir_folder_pipeline}{name_file}.fasta', n_threads, dir_folder_pipeline, dir_folder_pipeline_temp, flag_use_rcsb_pdb_only, int_number_of_proteins_in_a_batch_during_dnn_prediction, flag_use_all_gpu_devices )

                """ set flag_2 """
                with open( dir_file_flag_2, 'w' ) as newfile :
                    newfile.write( 'completed\n' )
                    
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    Info( "[Task Completion] Estimation of structural properties of input protein sequences completed at " + TIME_GET_timestamp( True ) )

    """
    Calculate B-cell epitope similarity scores based on structural properties of proteins 
    """
    
    """ check flag """
    dir_file_flag = f"{dir_file_matched}.predicting_B_cell_cross_reactivity_completed.flag"
    if not os.path.exists( dir_file_flag ) :
        
        Predict_B_cell_cross_reactivity( dir_folder_pipeline, dir_folder_pipeline_temp, n_threads, l_window_size, float_thres_e_value, flag_only_use_structural_properties_of_query_proteins, float_thres_avg_score_blosum_weighted__b_cell, float_thres_avg_score_blosum__b_cell, float_thres_rsa_correlation )
        
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    Info( "[Task Completion] Prediction of B-cell cross-reactivity was completed at " + TIME_GET_timestamp( True ) )
    
    """ 
    Calculate T-cell epitope similarity scores based on BLOSUM62 scores and predicted binding affinity scores.
    """
    
    """ check flag """
    dir_file_flag = f"{dir_file_matched}.predicting_T_cell_cross_reactivity_completed.flag"
    if not os.path.exists( dir_file_flag ) :
        
        Predict_T_cell_cross_reactivity( dir_folder_pipeline, float_thres_avg_blosum62_score_for_mhc, float_thres_min_mhc_allele_frequency, float_thres_binding_affinities_in_nM, flag_replace_unconventional_acid_code )
        
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    Info( "[Task Completion] Prediction of T-cell cross-reactivity was completed at " + TIME_GET_timestamp( True ) )
    
    """
    Further process and export data for visualization using a web application
    """
    
    """ check flag """
    dir_file_flag = f"{dir_file_matched}.prepare_data_for_web_application_completed.flag"
    if not os.path.exists( dir_file_flag ) :

        # combine results of all 'window_size' values
        OS_FILE_Combine_Files_in_order( glob.glob( f"{dir_folder_pipeline}b_cell.subsequence__window_size_*.tsv.gz" ), f"{dir_folder_pipeline}b_cell.subsequence.tsv.gz", remove_n_lines = 1, flag_use_header_from_first_file = True )
        # prepare data for web application using the combined subsequence
        # copy data for web application and encode using base64 encoding, and write metadata
        Prepare_data_for_web_application( f"{dir_folder_pipeline}b_cell.subsequence.tsv.gz", f"{dir_folder_pipeline}t_cell.mhc_binding.tsv.gz", dict_cressp_setting )
        
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )
    Info( "[Task Completion] Exporting data for CRESSP web application was completed at " + TIME_GET_timestamp( True ) )


def Parse_Structural_Properties( dir_file_df_sp, int_datatype = None, name_col_for_identifying_protein = 'id_protein' ) :
    """
    # 2021-05-31 20:29:02 
    parse structural properties with typical parameters for parsing ascii-encoded structural properties and typical column names  
    
    'dir_file_df_sp' : a file directory to TSV file containing structural properties estimated from CRESSP. Alternatively, a dataframe of the TSV file can be given as an input.
    'int_datatype' : used for initializing 'datatype_acc' if 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' column does not exist
    'name_col_for_identifying_protein' : columne name of the given tabular data for identifying a unique protein. The unique value in the column will be used as a 'key' in the returned dicitonaries
    
    """
    df_sp = pd.read_csv( dir_file_df_sp, sep = '\t' ) if isinstance( dir_file_df_sp, ( str ) ) else dir_file_df_sp # read database of structural properties # if an object that is not a string datatype is given, assumes 'dir_file_df_sp' is a dataframe containing structural properties
    df_sp.set_index( name_col_for_identifying_protein, inplace = True )
    
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    
    dict_sp = dict( )
    dict_sp[ 'acc' ] = ASCII_Decode( df_sp.rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1.to_dict( ), ** dict_kw_rsa ) # decode mkdssp outputs of RCSB PDB data
    dict_sp[ 'phi' ] = ASCII_Decode( df_sp[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
    dict_sp[ 'psi' ] = ASCII_Decode( df_sp[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
    dict_sp[ 'ss8' ] = ASCII_Decode( df_sp[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].to_dict( ), ** dict_kw_ss8 )
    if 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' in df_sp.columns.values :
        dict_sp[ 'datatype_acc' ] = ASCII_Decode( df_sp[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ].to_dict( ), ** dict_kw_datatype )
    else :
        if int_datatype is not None :
            return -1
            
        # initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
        dict_datatype_acc = dict( )
        for h in dict_sp[ 'acc' ] :
            arr = np.zeros_like( dict_sp[ 'acc' ][ h ] ) # default = 0 (predicted structural property)
            arr[ ~ np.isnan( dict_sp[ 'acc' ][ h ] ) ] = int_datatype
            dict_datatype_acc[ h ] = arr
        dict_sp[ 'datatype_acc' ] = dict_datatype_acc
    if 'structure_id___redundancy_reduced' in df_sp.columns.values :
        dict_fasta = df_sp[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
        dict_sp[ 'structure_id' ] = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )
        dict_sp[ 'seq' ] = df_sp[ 'seq' ].to_dict( )
    
    return dict_sp

def Encode_Structural_Properties( dict_sp, name_col_for_identifying_protein = 'id_protein' ) :
    """
    # 2021-05-31 20:28:58 
    Compose a dataframe containing encoded structural properties using given 'dict_sp'

    'dict_sp' : dictionary containing structural properties and protein sequences
    """
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    ''' initialize dataframe with protein sequences '''
    df_sp = pd.Series( dict_sp[ 'seq' ], name = 'seq' ).reset_index( ).rename( columns = { 'index' : name_col_for_identifying_protein } ).set_index( name_col_for_identifying_protein )

    # encode combined structural properties into ASCII strings using ASCII encoding
    df_sp[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ] = pd.Series( ASCII_Encode( dict_sp[ 'acc' ], ** dict_kw_rsa ) )
    df_sp[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'phi' ], ** dict_kw_torsion_angle ) )
    df_sp[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'psi' ], ** dict_kw_torsion_angle ) )
    df_sp[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ] = pd.Series( ASCII_Encode( dict_sp[ 'ss8' ], ** dict_kw_ss8 ) )
    df_sp[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ] = pd.Series( ASCII_Encode( dict_sp[ 'datatype_acc' ], ** dict_kw_datatype ) )
    df_sp[ 'structure_id___redundancy_reduced' ] = pd.Series( dict( ( h, Encode_List_of_Strings( dict_sp[ 'structure_id' ][ h ] ) ) for h in dict_sp[ 'structure_id' ] ) )
    df_sp.reset_index( drop = False, inplace = True )
    return df_sp
    
    
if __name__ == "__main__" :
    cressp( ) # run CRESSP at the top-level environment
    
