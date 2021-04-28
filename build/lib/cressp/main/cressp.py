#!/usr/bin/env python
# 20210419

from biobookshelf.main import *
from biobookshelf import *

from cressp.structural_property_estimation import Estimate_structural_property

import argparse
import traceback
import os, sys, getopt
from io import StringIO
import time
import math

# search for local similarities with all human genes # 2020-10-23 17:35:34 
def Calculate_similarity_score_using_structural_data_of_query_protein_only( dir_file_input, float_thres_avg_score_blosum_weighted, l_window_size, dir_folder_cressp, dir_folder_pipeline, dir_folder_pipeline_temp ) :
    """
    # Calculate_similarity_score_using_structural_data_of_query_protein_only
    """
    uuid_process = UUID( ) # uuid of current process
    np.warnings.filterwarnings( 'ignore' )

    # weighted similarity scoring of aligned sequences of a given window size
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

     
    df_protein_query = pd.read_csv( f'{dir_folder_pipeline}protein_query.tsv.gz', sep = '\t' ) # # load structural property data for query sequences
    dict_acc_to_arr_acc_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ].dropna( ).to_dict( ), ** dict_kw_rsa )
    dict_acc_to_arr_datatype_acc_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ].dropna( ).to_dict( ), ** dict_kw_datatype )
    dict_acc_to_arr_phi_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
    dict_acc_to_arr_psi_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
    dict_acc_to_arr_ss8_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].dropna( ).to_dict( ), ** dict_kw_ss8 )
    dict_acc_to_arr_datatype_ss8_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'ss8_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ].dropna( ).to_dict( ), ** dict_kw_datatype )
    dict_fasta = df_protein_query.set_index( 'id_protein' )[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
    dict_acc_to_arr_structure_id_query = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )

    # prepare Blosum62 score matrix including all non-conventional annotations of amino acids # 2020-10-29 02:15:32 
    # read dict_blosum62 from the tsv file
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score

    # list of ss8 classifications
    l_ss8 = [ 'G', 'H', 'I', 'E', 'B', 'T', 'S', 'C' ]
    
    df_matched = pd.read_csv( dir_file_input, sep = '\t' )
    l_l_output_file = list( )
    for int_window_size in l_window_size :
        dir_file_output = f"{dir_folder_pipeline_temp}{uuid_process}_window_size_{int_window_size}.tsv.gz"
        file_output = gzip.open( dir_file_output, 'wb' )
        l_col = [ 'id_alignment', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'score_blosum', 'score_blosum_weighted', 'sum_of_weights', 'n_residues_acc', 'prop_pdb_evidence_query', 'structure_id_query', 'count_structure_id_query', 'most_frequent_ss8_query', 'count_most_frequent_ss8_query' ]
        file_output.write( ( '\t'.join( l_col ) + '\n' ).encode( ) ) # write header to the output file
        n_errors = 0
        df = df_matched
        df = df[ df.query_alignment.apply( len ) >= int_window_size ] # drop entries with alignment shorter than the 'int_window_size'
        l_l_value = list( )
        for id_alignment, query_accession, target_accession, query_start, query_end, target_start, target_end, query_alignment, target_alignment, e_value, identity, source in df[ [ 'id_alignment', 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'identity', 'source' ] ].values :
            try :
                if query_accession not in dict_acc_to_arr_acc_query : continue # if accession does not exist, skip the entry
                arr_acc_query = dict_acc_to_arr_acc_query[ query_accession ]
                arr_acc_datatype_query = dict_acc_to_arr_datatype_acc_query[ query_accession ]
                arr_ss8_query = dict_acc_to_arr_ss8_query[ query_accession ]
                arr_ss8_datatype_query = dict_acc_to_arr_datatype_ss8_query[ query_accession ]

                len_seq_query, len_alignment = len( arr_acc_query ), len( query_alignment ) # retrieve lengths of sequences and alignment
                arr_phi_query = dict_acc_to_arr_phi_query[ query_accession ] if query_accession in dict_acc_to_arr_phi_query else np.full( len_seq_query, np.nan ) # phi and psi data might not be available for all query and target sequences. if the data is not available, use an empty array (array filled with np.nan) as phi and psi data
                arr_psi_query = dict_acc_to_arr_psi_query[ query_accession ] if query_accession in dict_acc_to_arr_psi_query else np.full( len_seq_query, np.nan )
                arr_structure_id_query = dict_acc_to_arr_structure_id_query[ query_accession ] if query_accession in dict_acc_to_arr_structure_id_query else np.full( len_seq_query, np.nan )

                arr_alignment_score_blosum, arr_alignment_acc_query, arr_alignment_acc_datatype_query, arr_alignment_ss8_query, arr_alignment_ss8_datatype_query, arr_alignment_phi_query, arr_alignment_psi_query = np.full( ( 6 + 1, len_alignment ), np.nan, dtype = float ) # initialize arrays for containing scores of alignment and structural properties of aligned target and query sequences # np.nan represents invalid values.
                arr_alignment_structure_id_query = np.full( len_alignment, np.nan, dtype = object )
                arr_alignment_gap = np.zeros( len_alignment, dtype = bool ) # True when the position contain a gap
                index_pos_in_query, index_pos_in_target = query_start - 1, target_start - 1 # 0-based coordinates        
                for index_pos_in_alignment, residue_query, residue_target in zip( np.arange( len_alignment ), query_alignment.upper( ), target_alignment.upper( ) ) : 
                    arr_alignment_score_blosum[ index_pos_in_alignment ] = dict_blosum62[ residue_query, residue_target ]
                    if residue_query == '-' : # ignore gaps in the alignment for retrieving aligned structural properties, except for RSA values and datatypes for RSA values (RSA values in the gaps are used as weights of the gap penalty score)
                        arr_alignment_acc_query[ index_pos_in_alignment ] = 1 # gaps were filled with 1 to use blosom gap score without weight 
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = np.nan # put invalid datatype value to a position with gap
                        arr_alignment_gap[ index_pos_in_alignment ] = True
                        index_pos_in_target += 1 
                    elif residue_target == '-' : 
                        arr_alignment_acc_query[ index_pos_in_alignment ] = arr_acc_query[ index_pos_in_query ]
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = arr_acc_datatype_query[ index_pos_in_query ]
                        arr_alignment_gap[ index_pos_in_alignment ] = True
                        index_pos_in_query += 1
                    else :
                        arr_alignment_acc_query[ index_pos_in_alignment ] = arr_acc_query[ index_pos_in_query ]
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = arr_acc_datatype_query[ index_pos_in_query ]
                        arr_alignment_ss8_query[ index_pos_in_alignment ] = arr_ss8_query[ index_pos_in_query ]
                        arr_alignment_ss8_datatype_query[ index_pos_in_alignment ] = arr_ss8_datatype_query[ index_pos_in_query ]
                        arr_alignment_phi_query[ index_pos_in_alignment ] = arr_phi_query[ index_pos_in_query ]
                        arr_alignment_psi_query[ index_pos_in_alignment ] = arr_psi_query[ index_pos_in_query ]
                        arr_alignment_structure_id_query[ index_pos_in_alignment ] = arr_structure_id_query[ index_pos_in_query ]
                        index_pos_in_query += 1
                        index_pos_in_target += 1
                arr_alignment_acc_query_scaled = arr_alignment_acc_query * 2 # scale relative accessible surface area values (multiply with 2 and replace values > 1 with 1)
                arr_alignment_acc_query_scaled[ arr_alignment_acc_query_scaled > 1 ] = 1
                arr_alignment_acc_weight = arr_alignment_acc_query_scaled
                arr_alignment_score_blosum_weighted = arr_alignment_score_blosum * arr_alignment_acc_weight
                mask_gap = np.isnan( arr_alignment_score_blosum_weighted )
                if mask_gap.sum( ) > 1 : # detect residues where surface availability data is unavailable
#                     print( f'surface availability data unavailable for {mask_gap.sum( )} residues' )
#                     print( 'sample_id\t', sample_id, 'id_alignment\t', id_alignment, 'query_accession\t', query_accession, 'target_accession\t', target_accession, 'query_start\t', query_start, 'query_end\t', query_end, 'target_start\t', target_start, 'target_end\t', target_end, 'query_alignment\t', query_alignment, 'target_alignment\t', target_alignment, 'e_value\t', e_value, 'identity\t', identity, 'source\t', source )
                    continue
                arr_alignment_acc_query_scaled[ arr_alignment_gap ] = np.nan # reassign RSA value of the positions containing gaps to np.nan (representing invalid value) for calculation of similarity score based on RSA values

                n_windows = len_alignment + 1 - int_window_size
                for index_pos_start in range( n_windows ) : 
                    slice_window = slice( index_pos_start, index_pos_start + int_window_size )
                    float_score_blosum_for_window = arr_alignment_score_blosum[ slice_window ].sum( )
                    float_score_blosum_weighted_for_window = arr_alignment_score_blosum_weighted[ slice_window ].sum( )
                    arr_alignment_acc_weight_for_window = arr_alignment_acc_weight[ slice_window ]
                    float_sum_of_acc_weights = arr_alignment_acc_weight_for_window[ ~ np.isnan( arr_alignment_acc_weight_for_window ) ].sum( ) # calculate total sum of weights

                    arr_alignment_acc_datatype_query_for_window = arr_alignment_acc_datatype_query[ slice_window ]
                    arr_alignment_acc_datatype_query_for_window_valid = arr_alignment_acc_datatype_query_for_window[ ~ np.isnan( arr_alignment_acc_datatype_query_for_window ) ] # drop invalid datatype values
                    float_prop_pdb_evidence_query = ( arr_alignment_acc_datatype_query_for_window_valid == 2 ).sum( ) / len( arr_alignment_acc_datatype_query_for_window_valid ) if len( arr_alignment_acc_datatype_query_for_window_valid ) > 0 else np.nan # calculate proportion of residues supported by PDB structures

                    n_aligned_residues_for_window = int_window_size - np.sum( arr_alignment_gap[ slice_window ] ) # retrieve number of aligne residue in the alignment

                    structure_id_query, count_structure_id_query = DICTIONARY_Find_Max( COUNTER( arr_alignment_structure_id_query[ slice_window ] ) )
                    # find major secondary structure classification for the current window
                    arr_alignment_ss8_query_for_window = arr_alignment_ss8_query[ slice_window ]
                    arr_alignment_ss8_query_for_window_valid = arr_alignment_ss8_query_for_window[ ~ np.isnan( arr_alignment_ss8_query_for_window ) ].astype( int )
                    id_ss8_most_frequent_query, count_id_ss8_most_frequent_query = DICTIONARY_Find_Max( COUNTER( arr_alignment_ss8_query_for_window_valid, ignore_float = True ) )

                    if structure_id_query is None and count_structure_id_query is None : structure_id_query, count_structure_id_query = np.nan, np.nan # retrieve structure_id of most frequence PDB structure (either from RCSB PDB or structural model)
                    if id_ss8_most_frequent_query is None and count_id_ss8_most_frequent_query is None :
                        str_ss8_query, id_ss8_most_frequent_query, count_id_ss8_most_frequent_query = np.nan, np.nan, np.nan
                    else :
                        id_ss8_most_frequent_query = int( id_ss8_most_frequent_query )
                        str_ss8_query = l_ss8[ id_ss8_most_frequent_query ]

                    query_alignment_subsequence = query_alignment[ slice_window ]
                    target_alignment_subsequence = target_alignment[ slice_window ]
                    int_gap_count_subsequence_query, int_gap_count_subsequence_target = query_alignment_subsequence.count( '-' ), target_alignment_subsequence.count( '-' )
                    int_gap_count_before_subsequence_query, int_gap_count_before_subsequence_target = query_alignment[ : index_pos_start ].count( '-' ), target_alignment[ : index_pos_start ].count( '-' )

                    if float_score_blosum_for_window < 0 and float_score_blosum_weighted_for_window < 0 : continue # do not write subsequence whose blosum score or weighted blosum scores is below 0
                    if float_score_blosum_weighted_for_window < float_thres_avg_score_blosum_weighted * int_window_size : continue # filter records with low avg alignment score weighted with accessiblities
                    l = [ id_alignment, index_pos_start + 1, index_pos_start + int_window_size, query_start + index_pos_start - int_gap_count_before_subsequence_query, query_start + index_pos_start + int_window_size - 1 - int_gap_count_subsequence_query - int_gap_count_before_subsequence_query, target_start + index_pos_start - int_gap_count_before_subsequence_target, target_start + index_pos_start + int_window_size - 1 - int_gap_count_subsequence_target - int_gap_count_before_subsequence_target, np.round( float_score_blosum_for_window, 3 ), np.round( float_score_blosum_weighted_for_window, 3 ), np.round( float_sum_of_acc_weights, 3 ), n_aligned_residues_for_window, float_prop_pdb_evidence_query, structure_id_query, count_structure_id_query, str_ss8_query, count_id_ss8_most_frequent_query ]
                    file_output.write( ( '\t'.join( list( map( str, l ) ) ) + '\n' ).encode( ) )
            except Exception as e : # count errors
                n_errors += 1
                print( f"error number {n_errors}", traceback.format_exc( ) )
        file_output.close( )
    return uuid_process
    
def Combine_result_files_for_each_window_size( dir_file_input, dir_folder_pipeline, dir_folder_pipeline_temp ) :
    """
    combine result files for each window_size
    """
    for arr_input in pd.read_csv( dir_file_input, sep = '\t', header = None ).values :
        int_window_size = arr_input[ 0 ]
        l_dir_file = glob.glob( f"{dir_folder_pipeline_temp}*_window_size_{int_window_size}.tsv.gz" ) # retrieve list of output files to combine into a single output file
        dir_file_output_combining = f'{dir_folder_pipeline}subsequence__window_size_{int_window_size}__struc_data_query_only.combining.tsv.gz'
        dir_file_output_combining_completed = f'{dir_folder_pipeline}subsequence__window_size_{int_window_size}__struc_data_query_only.tsv.gz'
        l_col = [ 'id_alignment', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'score_blosum', 'score_blosum_weighted', 'sum_of_weights', 'n_residues_acc', 'prop_pdb_evidence_query', 'structure_id_query', 'count_structure_id_query', 'most_frequent_ss8_query', 'count_most_frequent_ss8_query' ]
        OS_FILE_Combine_Files_in_order( l_dir_file, dir_file_output_combining, header = '\t'.join( l_col ) + '\n', remove_n_lines = 1 )
        print( f"combining output files for window size {int_window_size} is completed" )
        os.rename( dir_file_output_combining, dir_file_output_combining_completed ) # rename the file once completed
    
    
# def Bin_Similarity_Scores( dir_file_input ) :
#     """
#     # 2021-03-31 13:01:07 
#     bin similarity scores 
#     """
#     for arr_input in pd.read_csv( dir_file_input, header = None, sep = '\t' ).values :
#         int_window_size = int( arr_input[ 0 ] )
#         df_matched = pd.read_csv( dir_file_matched, sep = '\t', low_memory = False )
#         df_subsequence = pd.read_csv( f'{dir_folder_pipeline}subsequence__window_size_{int_window_size}__struc_data_human_only.tsv.gz', sep = '\t', low_memory = False )
#         df_subsequence.index.name = 'id_subsequence'
#         df_subsequence.reset_index( drop = False, inplace = True )
#         df = df_subsequence[ [ 'id_subsequence', "query_start", "query_end", 'score_blosum_weighted' ] ]
#         m = MAP.Map( df_matched.query_accession.to_dict( ) )
#         df[ 'acc_query' ] = df_subsequence.id_alignment.apply( m.a2b )
        
#         # Bin similarity scores by acc_query for analysis
#         s = df[ [ 'score_blosum_weighted', 'acc_query' ] ].groupby( 'acc_query' ).max( ).score_blosum_weighted
#         dir_file_output = f"{dir_folder_pipeline}similarity_score__window_size_{int_window_size}.max.bin_accession.tsv.gz"
#         s.to_csv( dir_file_output, sep = '\t' )
        
#         # Binning by 100 amino acid long region of acc_query with 50 amino acid overlaps
#         # build interval trees for score binning
#         df_uniprot_humans = pd.read_csv( dir_file_human_protein_db, sep = '\t' ) # load human proteins used for sequence similarity search
#         dict_acc_to_it = dict( )
#         for acc, length in df_uniprot_humans[ [ 'id_protein', 'seq_length' ] ].values :
#             dict_acc_to_it[ acc ] = intervaltree.IntervalTree( )
#             start = 0
#             while True :
#                 dict_acc_to_it[ acc ].addi( start, start + size_window_binning, [ f"{acc}_from_{start + 1}_to_{min( start + size_window_binning, length )}" ] )
#                 if start + size_window_binning >= length :  
#                     break
#                 else :
#                     start += size_window_binning - size_overlap_binning
#         # retrieve max_score and id_subsequence with the max_score for each id_feature (interval)
#         dict_id_feature_to_max_score = dict( )
#         dict_id_feature_to_id_subsequence_with_max_score = dict( )
#         for id_subsequence, query_start, query_end, score_blosum_weighted, acc_query in df[ [ 'id_subsequence', "query_start", "query_end", 'score_blosum_weighted', 'acc_query' ] ].values :
#             if acc_query not in dict_acc_to_it : continue 
#             it = dict_acc_to_it[ acc_query ]
#             id_feature = list( it.at( query_start - 1 ).intersection( it.at( query_end - 1 ) ) )[ 0 ][ 2 ][ 0 ]
#             if id_feature not in dict_id_feature_to_max_score :
#                 dict_id_feature_to_max_score[ id_feature ] = score_blosum_weighted
#                 dict_id_feature_to_id_subsequence_with_max_score[ id_feature ] = id_subsequence
#             else :
#                 if dict_id_feature_to_max_score[ id_feature ] < score_blosum_weighted : # update id_subsequence with max score
#                     dict_id_feature_to_max_score[ id_feature ] = score_blosum_weighted
#                     dict_id_feature_to_id_subsequence_with_max_score[ id_feature ] = id_subsequence
#         s_max_score = pd.Series( dict_id_feature_to_max_score )
#         s_max_score.index.name = 'id_feature'
#         s_max_score.name = "score_blosum_weighted"
#         s_id_subsequence_with_max_score = pd.Series( dict_id_feature_to_id_subsequence_with_max_score )
#         s_id_subsequence_with_max_score.index.name = 'id_feature'
#         s_id_subsequence_with_max_score.name = "id_subsequence"
        
#         dir_file_output_binning_by_acc_query_region__score = f"{dir_folder_pipeline}similarity_score__window_size_{int_window_size}.max.bin_{size_window_binning}_overlap_{size_overlap_binning}.score.tsv.gz"
#         dir_file_output_binning_by_acc_query_region__id = f"{dir_folder_pipeline}similarity_score__window_size_{int_window_size}.max.bin_{size_window_binning}_overlap_{size_overlap_binning}.id.tsv.gz"
#         s_max_score.to_csv( dir_file_output_binning_by_acc_query_region__score, sep = '\t' )
#         s_id_subsequence_with_max_score.to_csv( dir_file_output_binning_by_acc_query_region__id, sep = '\t' )

def main( ) :
    
#     # read dict_blosum62 from the tsv file
#     df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
#     dict_blosum62 = dict( )
#     for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
#         dict_blosum62[ aa_0, aa_1 ] = score


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
    ''' 
    currently not used
    '''
    parser.add_argument( "-Q", "--flag_skip_struc_prop_for_protein_target", help = "(Default: False) Set this flag to skip the estimation of structural properties of target proteins. Only structural properties of query proteins will be used to calculate accessibility-weighted-similarity scores", action = 'store_true' )
    
    args = parser.parse_args( )
    if args.dir_file_protein_target is None :
        print( "Required argument(s) is missing. to view help message, use -h or --help flag" )
        sys.exit( )
        
    # parse arguments # no further processing required
    flag_use_HMM_search = args.flag_use_HMM_search
    n_threads = int( args.cpu )
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
        
        '''
        <TEMPORARY
        '''
        print( "currently only the default human proteins are available, exiting" )
        sys.exit( )
        '''
        TEMPORARY>
        '''
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
    Estimate structural properties of proteins (to be implemented)
    """
    # use previously calculated structural properties when the default query proteins were used
    if flag_default_protein_query_was_used :
        shutil.copyfile( f'{dir_folder_cressp}data/human/uniprot.tsv.gz', f'{dir_folder_pipeline}protein_query.tsv.gz' ) 
        
    Estimate_structural_property( f'{dir_folder_pipeline}protein_target.fasta', n_threads, dir_folder_output, dir_folder_pipeline, dir_folder_pipeline_temp ) # test 
        
    """
    Calculate similarity scores based on structural properties of proteins 
    """
    df = pd.read_csv( dir_file_matched, sep = '\t' ) # read alignments between query and target protein sequences
    df.index.name = 'id_alignment' # retrieve id_alignment (index of df_matched) during retrieving subsequences
    df.reset_index( drop = False, inplace = True ) # add id_alignment column to to the dataframe
    print( f"number of records: {len( df )}" )
    df = df[ df.e_value <= float_thres_e_value ] # drop entries with too low global similarity
    print( f"number of records after filtering: {len( df )}" )
    
    l_uuid_process = Multiprocessing( df, Calculate_similarity_score_using_structural_data_of_query_protein_only, n_threads, dir_temp = dir_folder_pipeline_temp, global_arguments = [ float_thres_avg_score_blosum_weighted, l_window_size, dir_folder_cressp, dir_folder_pipeline, dir_folder_pipeline_temp ] ) # process similarity search result with multiple processes, and collect uuid of the processes
    
    
    # combine output files for each window size
    Multiprocessing( l_window_size, Combine_result_files_for_each_window_size, n_threads = min( len( l_window_size ), n_threads ), dir_temp = dir_folder_pipeline_temp, global_arguments = [ dir_folder_pipeline, dir_folder_pipeline_temp ] ) # combine result files for each window_size

    
    
#     # Bin similarity scores by acc_query and a given binning size for analysis
#     size_window_binning = 100
#     size_overlap_binning = 50
    
#     Multiprocessing( l_window_size, Bin_Similarity_Scores, n_threads = min( len( l_window_size ), int( OS_Memory( )[ 'MemAvailable' ] / 1e7 ) ), dir_temp = dir_folder_pipeline_temp ) # combine result files for each window_size




    
if __name__ == "__main__" :
    main( )
    