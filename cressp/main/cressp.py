#!/usr/bin/env python


from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

from cressp.structural_property_estimation import Estimate_structural_property

from cressp.cross_reactivity_prediction import Calculate_Similarity_Scores_in_Aligned_Sequences, Combine_result_files_for_each_window_size

import argparse
import traceback
import os, sys, getopt
from io import StringIO
import time
import math


#     # read dict_blosum62 from the tsv file
#     df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
#     dict_blosum62 = dict( )
#     for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
#         dict_blosum62[ aa_0, aa_1 ] = score

def main( ) :
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp2/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package


    """
    Parse arguments
    """
    # fixed arguments
    float_search_thres_e_value = 30 # set e-value threshold for search (fixed for maximum sensitivity)

    # command line arguments
    parser = argparse.ArgumentParser( description = "Cross-Reactive-Epitope-Search-using-Structural-Properties-of-proteins (cressp): a program to find cross-reactive epitopes with structural information from known protein structures." )
    parser.add_argument( "-t", "--dir_file_protein_target", help = "(Required) an input FASTA file containing target protein sequences." )
    parser.add_argument( "-q", "--dir_file_protein_query", help = "(Default: UniProt human proteins) an input FASTA file containing query protein sequences.", default = 'human' )
    parser.add_argument( "-o", "--dir_folder_output", help = "(Default: a subdirectory of the current directory) an output directory", default = "default" )
    parser.add_argument( "-c", "--cpu", help = "(Default: 1) Number of logical CPUs (threads) to use in the current compute node.", default = '1' )
    parser.add_argument( "-w", "--window_size", help = "(Default: 30) list of window sizes separated by comma. Example: 15,30,45", default = "30" )
    parser.add_argument( "-s", "--float_thres_avg_score_blosum_weighted", help = "(Default: 0.15) threshold for average weighted BLOSOM62 alignment score for filtering aligned sequences", default = '0.15' )
    parser.add_argument( "-e", "--float_thres_e_value", help = "(Default: 1e-20) threshold for the global alignment e-value in a scientific notation Example: 1e-3", default = "1e-20" )
    parser.add_argument( "-H", "--flag_use_HMM_search", help = "(Default: False) Set this flag to perform HMM search in addition to BLASTP search. HMM profile search is performed with HMMER3. The search usually takes several hours for metagenome-assembled genomes", action = 'store_true' )
    parser.add_argument( "-d", "--dir_file_query_hmmdb", help = "(Default: a HMM profile database of 1012 human proteins searched against UniProt Pan Proteomes. These proteins consist of experimentally validated human autoantigens) a file containing HMM DB of query proteins aligned against pan-proteomes", default = "human" )
    parser.add_argument( "-S", "--flag_use_rcsb_pdb_only", help = "When calculating consensus structural properties of input proteins, do not use homology-based modeled structures from SWISS-MODEL repositories and only use experimental protein structures from RCSB PDB", action = 'store_true' )
    parser.add_argument( "-Q", "--flag_only_use_structural_properties_of_query_proteins", help = "Only use estimated structural properties of the query proteins when predicting cross-reactivity between target and query proteins (skip the estimation of structural properties of target proteins). When 'dir_file_protein_query' == 'human' (default value), it will significantly reduce computation time by skipping estimation and prediction steps of structural properties of target & query proteins.", action = 'store_true' )

    args = parser.parse_args( )
    if args.dir_file_protein_target is None :
        print( "Required argument(s) is missing. to view help message, use -h or --help flag" )
        sys.exit( )

    # parse arguments # no further processing required
    flag_use_HMM_search = args.flag_use_HMM_search
    flag_use_rcsb_pdb_only = args.flag_use_rcsb_pdb_only
    flag_only_use_structural_properties_of_query_proteins = args.flag_only_use_structural_properties_of_query_proteins
    n_threads = min( int( args.cpu ), int( multiprocessing.cpu_count( ) ) ) # max number of n_threads is the CPU number of the current local machine
    l_window_size = list( int( e ) for e in args.window_size.split( ',' ) ) # set default window size
    float_thres_e_value = float( args.float_thres_e_value )
    float_thres_avg_score_blosum_weighted = float( args.float_thres_avg_score_blosum_weighted ) 

    # parse directory arguments
    dir_file_protein_target = args.dir_file_protein_target
    dir_file_protein_query = args.dir_file_protein_query
    dir_folder_output = args.dir_folder_output
    dir_file_query_hmmdb = args.dir_file_query_hmmdb

    # handle default settings for input datafiles
    flag_default_protein_query_was_used = False # a flag indicating whether the default protein_query was used
    if dir_file_protein_query == 'human' :
        flag_default_protein_query_was_used = True
        PKG.Download_Data( "data/human/uniprot.tsv.gz", dir_remote, name_package ) # download data
        dir_file_protein_query = f'{dir_folder_cressp}data/human/uniprot.fa' # set default 'dir_file_protein_query'
        if not os.path.exists( dir_file_protein_query ) : # if default human protein fasta file is not available, extract the sequence from the dataframe
            df_protein_query = pd.read_csv( f'{dir_folder_cressp}data/human/uniprot.tsv.gz', sep = '\t' )
            dict_fasta_protein_query = df_protein_query.set_index( 'id_protein' ).seq.to_dict( )
            FASTA_Write( dir_file_protein_query, dict_fasta = dict_fasta_protein_query )

        if flag_use_HMM_search : 
            if dir_file_query_hmmdb != 'human' :
                print( "since query proteins is the default human proteins, the default HMM profile database will be used instead." )
            PKG.Download_Data( "data/human/hmmdb_autoantigen.hmm.gz", dir_remote, name_package ) # download data
            PKG.Gunzip_Data( "data/human/hmmdb_autoantigen.hmm.gz", name_package ) # unzip data
            dir_file_query_hmmdb = pkg_resources.resource_filename( name_package, 'data/human/hmmdb_autoantigen.hmm' ) # set default 'dir_file_query_hmmdb'
    else :


        if flag_use_HMM_search and dir_file_query_hmmdb == 'human' :
            print( "exiting since query proteins other than default human proteins were given, the default HMM profile database cannot be used" )
            sys.exit( )
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
    if os.path.exists( dir_folder_output ) :
        print( "exiting since the given output directory already exists" )
        sys.exit( )
    else :
        os.makedirs( dir_folder_output, exist_ok = True )
    # create a temporary output folder inside the output folder
    dir_folder_pipeline = f"{dir_folder_output}pipeline/"
    os.makedirs( dir_folder_pipeline, exist_ok = True )
    dir_folder_pipeline_temp = f'{dir_folder_pipeline}temp/' 
    os.makedirs( dir_folder_pipeline_temp, exist_ok = True )


    """
    Read and move input protein Fasta files
    """
    try :
        dict_fasta_protein_query = FASTA_Read( dir_file_protein_query )
        FASTA_Write( f"{dir_folder_pipeline}protein_query.fasta", dict_fasta = dict_fasta_protein_query )
        dir_file_protein_query = f"{dir_folder_pipeline}protein_query.fasta" # set directory of fasta file to the new file location
    except :
        print( f"exiting due to an error while reading and moving 'dir_file_protein_query' {dir_file_protein_query}" )
        sys.exit( )

    try :
        dict_fasta_protein_target = FASTA_Read( dir_file_protein_target )
        FASTA_Write( f"{dir_folder_pipeline}protein_target.fasta", dict_fasta = dict_fasta_protein_target )
        dir_file_protein_target = f"{dir_folder_pipeline}protein_target.fasta" # set directory of fasta file to the new file location
    except :
        print( f"exiting due to an error while reading and moving 'dir_file_protein_target' {dir_file_protein_target}" )
        sys.exit( )


    """
    Perform BLASTP alignment
    """
    # create blastp_db using query_protein sequences
    dir_prefix_blastdb_protein_query = f"{dir_folder_pipeline}makeblastdb_out/protein_query"
    os.makedirs( f"{dir_folder_pipeline}makeblastdb_out/", exist_ok = True )  
    OS_Run( [ "makeblastdb", "-in", f"{dir_folder_pipeline}protein_query.fasta", '-dbtype', 'prot', '-parse_seqids', '-max_file_sz', '1GB', '-out', dir_prefix_blastdb_protein_query ], dir_file_stdout = f"{dir_prefix_blastdb_protein_query}.makeblastdb.stdout.txt", dir_file_stderr = f"{dir_prefix_blastdb_protein_query}.makeblastdb.stderr.txt", return_output = False ) # make blast db for protein_query

    # run blastp
    dir_file_blastp_output = f'{dir_folder_pipeline}blastp.tsv'
    OS_Run( [ 'blastp', '-query', dir_file_protein_target, '-db', dir_prefix_blastdb_protein_query, '-out', dir_file_blastp_output, '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop', '-num_threads', f'{n_threads}', '-evalue', f'{float_search_thres_e_value}' ], dir_file_stdout = f"{dir_folder_pipeline}blastp.stdout.txt", dir_file_stderr = f"{dir_folder_pipeline}blastp.stderr.txt", return_output = False ) # run blastp
    OS_Run( [ 'gzip', dir_file_blastp_output ], dir_file_stdout = f"{dir_file_blastp_output}.gzip.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.gzip.stderr.txt", return_output = False ) # compress blastp output
    dir_file_blastp_output += '.gz'


    """
    Perform HMMER alignment using a given HMM profile DB (using a couple of HMM profile DBs would be also helpful in the near future)
    """
    # run hmmsearch according to 'flag_use_HMM_search' flag
    if flag_use_HMM_search : 
        dir_file_hmmsearch_output = f'{dir_folder_pipeline}hmmsearch.out'
        OS_Run( [ 'hmmsearch', '-o', dir_file_hmmsearch_output, '--acc', '--notextw', '--cpu', f'{n_threads}', '-E', f'{float_search_thres_e_value}', dir_file_query_hmmdb, dir_file_protein_target ], dir_file_stdout = f"{dir_folder_pipeline}hmmsearch.stdout.txt", dir_file_stderr = f"{dir_folder_pipeline}hmmsearch.stderr.txt", return_output = False ) # run hmmsearch

    """
    Combine BLASTP and HMMSEARCH outputs 
    """
    dir_file_matched = f'{dir_folder_pipeline}matched.tsv.gz'
    dir_file_matched_write_complete_flag = f"{dir_file_matched}.write_completed.flag" # a flag indicating the write operation of the file has been completed.

    # load blastp result
    dict_qacc_to_seq = FASTA_Read( dir_file_protein_target ) # read query protein sequences
    dict_qacc_to_seq = dict( ( header.split( ' ', 1 )[ 0 ], dict_qacc_to_seq[ header ] ) for header in list( dict_qacc_to_seq ) )
    df_blastp = BLAST_Read( dir_file_blastp_output, dict_qaccver_to_seq = dict_qacc_to_seq ) 
    df_blastp = df_blastp[ [ 'saccver', 'qaccver', 'sstart', 'send', 'qstart', 'qend', 'subject_seq_aligned', 'query_seq_aligned', 'evalue', 'pident' ] ]
    df_blastp.pident = df_blastp.pident / 100
    df_blastp.columns = [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'identity' ]
    df_blastp[ 'source' ] = 'blastp'

    # load hmmer result according to 'flag_use_HMM_search' flag
    if flag_use_HMM_search : 
        df = HMMER_HMMSEARCH_Read_output( dir_file_hmmsearch_output )
        dict_qacc_to_seq = dict_fasta_protein_query
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

    df_matched = pd.concat( [ df_hmmer, df_blastp ], ignore_index = True ) if flag_use_HMM_search else df_blastp
    df_matched.to_csv( dir_file_matched, sep = '\t', index = False )

    # write a flag (file) indicating the writing operation was completed.
    with open( dir_file_matched_write_complete_flag, 'w' ) as file :
        file.write( f"search results were written at {TIME_GET_timestamp( )}" )

    """
    Estimate structural properties of proteins 
    """
    # use previously calculated structural properties when the default query proteins were used
    if flag_default_protein_query_was_used :
        shutil.copyfile( f'{dir_folder_cressp}data/human/uniprot.tsv.gz', f'{dir_folder_pipeline}protein_query.tsv.gz' ) 
    # estimate structural properties
    for name_file in [ 'protein_target', 'protein_query' ] :
        if not os.path.exists( f'{dir_folder_pipeline}{name_file}.tsv.gz' ) :
            Estimate_structural_property( f'{dir_folder_pipeline}{name_file}.fasta', n_threads, dir_folder_output, dir_folder_pipeline, dir_folder_pipeline_temp, flag_use_rcsb_pdb_only )


    """
    Calculate similarity scores based on structural properties of proteins 
    """
    df = pd.read_csv( dir_file_matched, sep = '\t' ) # read alignments between query and target protein sequences
    df.index.name = 'id_alignment' # retrieve id_alignment (index of df_matched) during retrieving subsequences
    df.reset_index( drop = False, inplace = True ) # add id_alignment column to to the dataframe
    print( f"number of records: {len( df )}" )
    df = df[ df.e_value <= float_thres_e_value ] # drop entries with too low global similarity
    print( f"number of records after filtering: {len( df )}" )


    # predict cross-reactivity
    l_uuid_process = Multiprocessing( df, Calculate_Similarity_Scores_in_Aligned_Sequences, n_threads, dir_temp = dir_folder_pipeline_temp, global_arguments = [ float_thres_avg_score_blosum_weighted, l_window_size, dir_folder_cressp, dir_folder_pipeline, dir_folder_pipeline_temp, flag_only_use_structural_properties_of_query_proteins ] ) # process similarity search result with multiple processes, and collect uuid of the processes

    # combine output files for each window size
    Multiprocessing( l_window_size, Combine_result_files_for_each_window_size, n_threads = min( len( l_window_size ), n_threads ), dir_temp = dir_folder_pipeline_temp, global_arguments = [ dir_folder_pipeline, dir_folder_pipeline_temp ] ) # combine result files for each window_size



    #     # Bin similarity scores by acc_query and a given binning size for analysis
    #     size_window_binning = 100
    #     size_overlap_binning = 50

    #     Multiprocessing( l_window_size, Bin_Similarity_Scores, n_threads = min( len( l_window_size ), int( OS_Memory( )[ 'MemAvailable' ] / 1e7 ) ), dir_temp = dir_folder_pipeline_temp ) # combine result files for each window_size



    
if __name__ == "__main__" :
    main( )
    