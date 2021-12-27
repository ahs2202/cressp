from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

# define global read-only variables
dict_acc_to_arr_acc_query, dict_acc_to_arr_datatype_acc_query, dict_acc_to_arr_phi_query, dict_acc_to_arr_psi_query, dict_acc_to_arr_ss8_query, dict_acc_to_arr_structure_id_query, dict_acc_to_arr_acc_target, dict_acc_to_arr_datatype_acc_target, dict_acc_to_arr_phi_target, dict_acc_to_arr_psi_target, dict_acc_to_arr_ss8_target, dict_acc_to_arr_structure_id_target = [ dict( ) ] * 12

def BCell_Calculate_Similarity_Scores_in_Aligned_Sequences( dir_file_input, float_thres_avg_score_blosum_weighted, float_thres_avg_score_blosum, float_thres_rsa_correlation, l_window_size, dir_folder_cressp, dir_folder_pipeline, dir_folder_pipeline_temp, flag_only_use_structural_properties_of_query_proteins ) :
    """
    Calculate_Similarity_Scores_in_Aligned_Sequences
    using structural properties of both query and target proteins
    """
    uuid_process = UUID( ) # uuid of current process
    np.warnings.filterwarnings( 'ignore' )
    
    float_scaling_factor_for_diff_torsion_angle = 20 # scaling of absolute difference of torsion angle (np.abs): 'float_scaling_factor_for_diff_torsion_angle' / np.abs( angle diff. ) - 'float_scaling_factor_for_diff_torsion_angle' / ( max np.abs( angle diff. ) or 180 )
    float_scaling_factor_for_diff_acc = 0.25 

    # read input alignment between query and target proteins
    df_matched = pd.read_csv( dir_file_input, sep = '\t' )

    # prepare Blosum62 score matrix including all non-conventional annotations of amino acids # 2020-10-29 02:15:32 
    # read dict_blosum62 from the tsv file
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score
        
    # list of ss8 classifications
    l_ss8 = [ 'G', 'H', 'I', 'E', 'B', 'T', 'S', 'C' ]
    
    # for each window_size
    for int_window_size in l_window_size :
        ''' initialize the output file '''
        dir_file_output = f"{dir_folder_pipeline_temp}{uuid_process}_window_size_{int_window_size}.tsv.gz"
        file_output = gzip.open( dir_file_output, 'wb' )
        l_col = [ 'window_size', 'id_alignment', 'source', 'query_accession', 'target_accession', 'e_value', 'identity', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'score_blosum', 'score_blosum_weighted', 'sum_of_weights', 'score_similarity_acc', 'score_similarity_phi', 'score_similarity_psi', 'score_similarity_ss8', 'n_residues_acc', 'n_residues_phi', 'n_residues_psi', 'n_residues_ss8', 'correl_coeffi_acc', 'correl_p_value_acc', 'correl_coeffi_phi', 'correl_p_value_phi', 'correl_coeffi_psi', 'correl_p_value_psi', 'prop_pdb_evidence_query', 'prop_pdb_evidence_target', 'structure_id_query', 'count_structure_id_query', 'structure_id_target', 'count_structure_id_target', 'most_frequent_ss8_query', 'count_most_frequent_ss8_query', 'most_frequent_ss8_target', 'count_most_frequent_ss8_target' ]
        file_output.write( ( '\t'.join( l_col ) + '\n' ).encode( ) ) # write header to the output file
        n_errors = 0
        
        df = df_matched # copy reference to the alignment file
        df = df[ df.query_alignment.apply( len ) >= int_window_size ] # drop entries smaller than the 'int_window_size'
        l_l_value = list( )
        ''' for each alignment between query and target protein '''
        for id_alignment, query_accession, target_accession, query_start, query_end, target_start, target_end, query_alignment, target_alignment, e_value, identity, source in df[ [ 'id_alignment', 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'identity', 'source' ] ].values :
            try :
                ''' retrieve structural properties of query and target proteins '''
                arr_acc_query = dict_acc_to_arr_acc_query[ query_accession ]
                arr_acc_datatype_query = dict_acc_to_arr_datatype_acc_query[ query_accession ]
                arr_ss8_query = dict_acc_to_arr_ss8_query[ query_accession ]
                arr_acc_target = dict_acc_to_arr_acc_target[ target_accession ]
                arr_acc_datatype_target = dict_acc_to_arr_datatype_acc_target[ target_accession ]
                arr_ss8_target = dict_acc_to_arr_ss8_target[ target_accession ]

                len_seq_query, len_seq_target, len_alignment = len( arr_acc_query ), len( arr_acc_target ), len( query_alignment ) # retrieve lengths of sequences and alignment
                arr_phi_query = dict_acc_to_arr_phi_query[ query_accession ] if query_accession in dict_acc_to_arr_phi_query else np.full( len_seq_query, np.nan ) # phi and psi data might not be available for all query and target sequences. if the data is not available, use an empty array (array filled with np.nan) as phi and psi data
                arr_psi_query = dict_acc_to_arr_psi_query[ query_accession ] if query_accession in dict_acc_to_arr_psi_query else np.full( len_seq_query, np.nan )
                arr_structure_id_query = dict_acc_to_arr_structure_id_query[ query_accession ] if query_accession in dict_acc_to_arr_structure_id_query else np.full( len_seq_query, np.nan )
                arr_phi_target = dict_acc_to_arr_phi_target[ target_accession ] if target_accession in dict_acc_to_arr_phi_target else np.full( len_seq_target, np.nan )
                arr_psi_target = dict_acc_to_arr_psi_target[ target_accession ] if target_accession in dict_acc_to_arr_psi_target else np.full( len_seq_target, np.nan )
                arr_structure_id_target = dict_acc_to_arr_structure_id_target[ target_accession ] if target_accession in dict_acc_to_arr_structure_id_target else np.full( len_seq_target, np.nan )

                ''' initialize '''
                arr_alignment_score_blosum, arr_alignment_acc_query, arr_alignment_acc_datatype_query, arr_alignment_ss8_query, arr_alignment_phi_query, arr_alignment_psi_query, arr_alignment_acc_target, arr_alignment_acc_datatype_target, arr_alignment_ss8_target, arr_alignment_phi_target, arr_alignment_psi_target = np.full( ( 10 + 1, len_alignment ), np.nan, dtype = float ) # initialize arrays for containing scores of alignment and structural properties of aligned target and query sequences # np.nan represents invalid values.
                arr_alignment_structure_id_query, arr_alignment_structure_id_target = np.full( ( 2, len_alignment ), np.nan, dtype = object )
                arr_alignment_gap = np.zeros( len_alignment, dtype = bool ) # True when the position contain a gap
                index_pos_in_query, index_pos_in_target = query_start - 1, target_start - 1 # 1-based to 0-based coordinates        
                ''' map structural properties of query and target proteins to the alignment '''
                # loop through each position in the alignment
                for index_pos_in_alignment, residue_query, residue_target in zip( np.arange( len_alignment ), query_alignment.upper( ), target_alignment.upper( ) ) : 
                    arr_alignment_score_blosum[ index_pos_in_alignment ] = dict_blosum62[ residue_query, residue_target ]
                    if residue_query == '-' : # ignore gaps in the alignment for retrieving aligned structural properties, except for RSA values and datatypes for RSA values (RSA values in the gaps are used as weights of the gap penalty score)
                        arr_alignment_acc_target[ index_pos_in_alignment ] = arr_acc_target[ index_pos_in_target ]
                        arr_alignment_acc_datatype_target[ index_pos_in_alignment ] = arr_acc_datatype_target[ index_pos_in_target ]
                        arr_alignment_acc_query[ index_pos_in_alignment ] = 1
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = np.nan # put invalid datatype value to a position with gap
                        arr_alignment_gap[ index_pos_in_alignment ] = True
                        index_pos_in_target += 1 
                    elif residue_target == '-' : 
                        arr_alignment_acc_target[ index_pos_in_alignment ] = 1
                        arr_alignment_acc_datatype_target[ index_pos_in_alignment ] = np.nan # put invalid datatype value to a position with gap
                        arr_alignment_acc_query[ index_pos_in_alignment ] = arr_acc_query[ index_pos_in_query ]
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = arr_acc_datatype_query[ index_pos_in_query ]
                        arr_alignment_gap[ index_pos_in_alignment ] = True
                        index_pos_in_query += 1
                    else :
                        arr_alignment_acc_query[ index_pos_in_alignment ] = arr_acc_query[ index_pos_in_query ]
                        arr_alignment_acc_datatype_query[ index_pos_in_alignment ] = arr_acc_datatype_query[ index_pos_in_query ]
                        arr_alignment_ss8_query[ index_pos_in_alignment ] = arr_ss8_query[ index_pos_in_query ]
                        arr_alignment_phi_query[ index_pos_in_alignment ] = arr_phi_query[ index_pos_in_query ]
                        arr_alignment_psi_query[ index_pos_in_alignment ] = arr_psi_query[ index_pos_in_query ]
                        arr_alignment_structure_id_query[ index_pos_in_alignment ] = arr_structure_id_query[ index_pos_in_query ]
                        arr_alignment_acc_target[ index_pos_in_alignment ] = arr_acc_target[ index_pos_in_target ]
                        arr_alignment_acc_datatype_target[ index_pos_in_alignment ] = arr_acc_datatype_target[ index_pos_in_target ]
                        arr_alignment_ss8_target[ index_pos_in_alignment ] = arr_ss8_target[ index_pos_in_target ]
                        arr_alignment_phi_target[ index_pos_in_alignment ] = arr_phi_target[ index_pos_in_target ]
                        arr_alignment_psi_target[ index_pos_in_alignment ] = arr_psi_target[ index_pos_in_target ]
                        arr_alignment_structure_id_target[ index_pos_in_alignment ] = arr_structure_id_target[ index_pos_in_target ]
                        index_pos_in_query += 1
                        index_pos_in_target += 1
                """ predict cross-reactivity """
                arr_alignment_acc_query_scaled = arr_alignment_acc_query * 2 # scale relative accessible surface area values (multiply with 2 and replace values > 1 with 1)
                arr_alignment_acc_target_scaled = arr_alignment_acc_target * 2
                arr_alignment_acc_query_scaled[ arr_alignment_acc_query_scaled > 1 ] = 1
                arr_alignment_acc_target_scaled[ arr_alignment_acc_target_scaled > 1 ] = 1
                arr_alignment_acc_weight = arr_alignment_acc_query_scaled * arr_alignment_acc_target_scaled
                arr_alignment_score_blosum_weighted = arr_alignment_score_blosum * arr_alignment_acc_weight
                mask_gap = np.isnan( arr_alignment_score_blosum_weighted )
                if mask_gap.sum( ) > 1 : 
                    pass
                arr_alignment_acc_query_scaled[ arr_alignment_gap ] = np.nan # reassign RSA value of the positions containing gaps to np.nan (representing invalid value) for calculation of similarity score based on RSA values
                arr_alignment_acc_target_scaled[ arr_alignment_gap ] = np.nan # reassign RSA value of the positions containing gaps to np.nan (representing invalid value) for calculation of similarity score based on RSA values
                arr_alignment_score_similarity_acc = ( 1 + float_scaling_factor_for_diff_acc ) - np.abs( arr_alignment_acc_query_scaled - arr_alignment_acc_target_scaled ) * ( 1 + float_scaling_factor_for_diff_acc ) # calculate rsa similarity scores using scaled rsa values after alignment
                arr_alignment_score_similarity_acc[ arr_alignment_score_similarity_acc > 1 ] = 1 

                arr_alignment_phi_diff = np.abs( arr_alignment_phi_query - arr_alignment_phi_target ) # calculate torsion angle similarity scores using scaled torsion angle values after alignment
                mask = arr_alignment_phi_diff > 180 # max difference of torsion angles should be 180
                arr_alignment_phi_diff[ mask ] = 360 - arr_alignment_phi_diff[ mask ]
                arr_alignment_psi_diff = np.abs( arr_alignment_psi_query - arr_alignment_psi_target )
                mask = arr_alignment_psi_diff > 180
                arr_alignment_psi_diff[ mask ] = 360 - arr_alignment_psi_diff[ mask ]
                float_min_diff_torsion_angle = float_scaling_factor_for_diff_torsion_angle / 2
                arr_alignment_phi_diff[ arr_alignment_phi_diff < float_min_diff_torsion_angle ] = float_min_diff_torsion_angle # replace 0 with 'float_min_diff_torsion_angle' to avoid zero division
                arr_alignment_psi_diff[ arr_alignment_psi_diff < float_min_diff_torsion_angle ] = float_min_diff_torsion_angle # replace 0 with 'float_min_diff_torsion_angle' to avoid zero division
                arr_alignment_score_similarity_phi = float_scaling_factor_for_diff_torsion_angle / arr_alignment_phi_diff - float_scaling_factor_for_diff_torsion_angle / 180 # max difference of torsion angles is 180
                arr_alignment_score_similarity_psi = float_scaling_factor_for_diff_torsion_angle / arr_alignment_psi_diff - float_scaling_factor_for_diff_torsion_angle / 180 # max difference of torsion angles is 180
                arr_alignment_score_similarity_phi[ arr_alignment_score_similarity_phi > 1 ] = 1
                arr_alignment_score_similarity_psi[ arr_alignment_score_similarity_psi > 1 ] = 1 

                arr_alignment_score_similarity_ss8 = ( arr_alignment_ss8_query == arr_alignment_ss8_target ).astype( float ) # calculate ss8 classification similarity scores using ss8 classifications after alignment
                arr_alignment_score_similarity_ss8[ np.isnan( arr_alignment_ss8_query ) | np.isnan( arr_alignment_ss8_target ) ] = np.nan

                n_windows = len_alignment + 1 - int_window_size
                for index_pos_start in range( n_windows ) : 
                    slice_window = slice( index_pos_start, index_pos_start + int_window_size )
                    float_score_blosum_for_window = arr_alignment_score_blosum[ slice_window ].sum( )
                    float_score_blosum_weighted_for_window = arr_alignment_score_blosum_weighted[ slice_window ].sum( )
                    arr_alignment_acc_weight_for_window = arr_alignment_acc_weight[ slice_window ]
                    float_sum_of_acc_weights = arr_alignment_acc_weight_for_window[ ~ np.isnan( arr_alignment_acc_weight_for_window ) ].sum( ) # calculate total sum of weights

                    float_score_similarity_acc_for_window, n_valid_residues_score_similarity_acc_for_window = NUMPY_Calculate_Average( arr_alignment_score_similarity_acc[ slice_window ] )
                    float_score_similarity_phi_for_window, n_valid_residues_score_similarity_phi_for_window = NUMPY_Calculate_Average( arr_alignment_score_similarity_phi[ slice_window ] )
                    float_score_similarity_psi_for_window, n_valid_residues_score_similarity_psi_for_window = NUMPY_Calculate_Average( arr_alignment_score_similarity_psi[ slice_window ] )
                    float_score_similarity_ss8_for_window, n_valid_residues_score_similarity_ss8_for_window = NUMPY_Calculate_Average( arr_alignment_score_similarity_ss8[ slice_window ] )
                    float_correl_coeffi_acc, float_correl_p_value_acc, _ = PEARSONR( arr_alignment_acc_query_scaled[ slice_window ], arr_alignment_acc_target_scaled[ slice_window ] )
                    float_correl_coeffi_phi, float_correl_p_value_phi, _ = PEARSONR( arr_alignment_phi_query[ slice_window ], arr_alignment_phi_target[ slice_window ] )
                    float_correl_coeffi_psi, float_correl_p_value_psi, _ = PEARSONR( arr_alignment_psi_query[ slice_window ], arr_alignment_psi_target[ slice_window ] )

                    arr_alignment_acc_datatype_query_for_window = arr_alignment_acc_datatype_query[ slice_window ]
                    arr_alignment_acc_datatype_target_for_window = arr_alignment_acc_datatype_target[ slice_window ]
                    arr_alignment_acc_datatype_query_for_window_valid = arr_alignment_acc_datatype_query_for_window[ ~ np.isnan( arr_alignment_acc_datatype_query_for_window ) ] # drop invalid datatype values
                    arr_alignment_acc_datatype_target_for_window_valid = arr_alignment_acc_datatype_target_for_window[ ~ np.isnan( arr_alignment_acc_datatype_target_for_window ) ] # drop invalid datatype values
                    float_prop_pdb_evidence_query = ( arr_alignment_acc_datatype_query_for_window_valid == 2 ).sum( ) / len( arr_alignment_acc_datatype_query_for_window_valid ) if len( arr_alignment_acc_datatype_query_for_window_valid ) > 0 else np.nan # calculate proportion of residues supported by PDB structures
                    float_prop_pdb_evidence_target = ( arr_alignment_acc_datatype_target_for_window_valid == 2 ).sum( ) / len( arr_alignment_acc_datatype_target_for_window_valid ) if len( arr_alignment_acc_datatype_target_for_window_valid ) > 0 else np.nan # calculate proportion of residues supported by PDB structures

                    ''' retrieve most frequent structure_id at the aligned region (for query and target proteins) '''
                    structure_id_query, count_structure_id_query = DICTIONARY_Find_Max( COUNTER( arr_alignment_structure_id_query[ slice_window ] ) )
                    structure_id_target, count_structure_id_target = DICTIONARY_Find_Max( COUNTER( arr_alignment_structure_id_target[ slice_window ] ) )
                    if structure_id_query is None and count_structure_id_query is None : structure_id_query, count_structure_id_query = np.nan, np.nan # retrieve structure_id of most frequence PDB structure (either from RCSB PDB or structural model)
                    if structure_id_target is None and count_structure_id_target is None : structure_id_target, count_structure_id_target = np.nan, np.nan # retrieve structure_id of most frequence PDB structure (either from RCSB PDB or structural model)
                    
                    ''' find most frequent secondary structure classification (ss8) at the aligned region (for query and target proteins) '''
                    # find major secondary structure classification for the current window
                    # query protein
                    arr_alignment_ss8_query_for_window = arr_alignment_ss8_query[ slice_window ]
                    arr_alignment_ss8_query_for_window_valid = arr_alignment_ss8_query_for_window[ ~ np.isnan( arr_alignment_ss8_query_for_window ) ].astype( int )
                    id_ss8_most_frequent_query, count_id_ss8_most_frequent_query = DICTIONARY_Find_Max( COUNTER( arr_alignment_ss8_query_for_window_valid, ignore_float = True ) )
                    if id_ss8_most_frequent_query is None and count_id_ss8_most_frequent_query is None :
                        str_ss8_query, id_ss8_most_frequent_query, count_id_ss8_most_frequent_query = np.nan, np.nan, np.nan
                    else :
                        id_ss8_most_frequent_query = int( id_ss8_most_frequent_query )
                        str_ss8_query = l_ss8[ id_ss8_most_frequent_query ]
                    # target protein
                    arr_alignment_ss8_target_for_window = arr_alignment_ss8_target[ slice_window ]
                    arr_alignment_ss8_target_for_window_valid = arr_alignment_ss8_target_for_window[ ~ np.isnan( arr_alignment_ss8_target_for_window ) ].astype( int )
                    id_ss8_most_frequent_target, count_id_ss8_most_frequent_target = DICTIONARY_Find_Max( COUNTER( arr_alignment_ss8_target_for_window_valid, ignore_float = True ) )
                    if id_ss8_most_frequent_target is None and count_id_ss8_most_frequent_target is None :
                        str_ss8_target, id_ss8_most_frequent_target, count_id_ss8_most_frequent_target = np.nan, np.nan, np.nan
                    else :
                        id_ss8_most_frequent_target = int( id_ss8_most_frequent_target )
                        str_ss8_target = l_ss8[ id_ss8_most_frequent_target ]
                    
                    ''' retrieve subsequence of query and target proteins and count gaps '''
                    query_alignment_subsequence = query_alignment[ slice_window ]
                    target_alignment_subsequence = target_alignment[ slice_window ]
                    int_gap_count_subsequence_query, int_gap_count_subsequence_target = query_alignment_subsequence.count( '-' ), target_alignment_subsequence.count( '-' )
                    int_gap_count_before_subsequence_query, int_gap_count_before_subsequence_target = query_alignment[ : index_pos_start ].count( '-' ), target_alignment[ : index_pos_start ].count( '-' )
                    
                    ''' perform filtering using given threshold values '''
                    if float_score_blosum_for_window < float_thres_avg_score_blosum * int_window_size : continue # do not write subsequence whose average blosum score is below a given threshold
                    if float_correl_coeffi_acc < float_thres_rsa_correlation : continue # do not write subsequence whose correlation coefficient of rsa values is below a given threshold
                    if float_score_blosum_weighted_for_window < float_thres_avg_score_blosum_weighted * int_window_size : continue # filter records with low avg alignment score weighted with accessiblities
                        
                    ''' write a record to the output file '''
                    l = [ int_window_size, id_alignment, source, query_accession, target_accession, e_value, identity, index_pos_start + 1, index_pos_start + int_window_size, query_start + index_pos_start - int_gap_count_before_subsequence_query, query_start + index_pos_start + int_window_size - 1 - int_gap_count_subsequence_query - int_gap_count_before_subsequence_query, target_start + index_pos_start - int_gap_count_before_subsequence_target, target_start + index_pos_start + int_window_size - 1 - int_gap_count_subsequence_target - int_gap_count_before_subsequence_target, query_alignment_subsequence, target_alignment_subsequence, float_score_blosum_for_window, float_score_blosum_weighted_for_window, float_sum_of_acc_weights, float_score_similarity_acc_for_window, float_score_similarity_phi_for_window, float_score_similarity_psi_for_window, float_score_similarity_ss8_for_window, n_valid_residues_score_similarity_acc_for_window, n_valid_residues_score_similarity_phi_for_window, n_valid_residues_score_similarity_psi_for_window, n_valid_residues_score_similarity_ss8_for_window, float_correl_coeffi_acc, float_correl_p_value_acc, float_correl_coeffi_phi, float_correl_p_value_phi, float_correl_coeffi_psi, float_correl_p_value_psi, float_prop_pdb_evidence_query, float_prop_pdb_evidence_target, structure_id_query, count_structure_id_query, structure_id_target, count_structure_id_target, str_ss8_query, count_id_ss8_most_frequent_query, str_ss8_target, count_id_ss8_most_frequent_target ]
                    file_output.write( ( '\t'.join( list( map( str, l ) ) ) + '\n' ).encode( ) ) # write a record to the output file
            except Exception as e : # count errors
                n_errors += 1
                print( f"error number {n_errors}", traceback.format_exc( ) )
        file_output.close( ) # close output file
    return uuid_process # return the UUID of the current process 
def BCell_Combine_result_files_for_each_window_size( dir_file_input, dir_folder_pipeline, dir_folder_pipeline_temp ) :
    """
    combine result files for each window_size
    """
    for arr_input in pd.read_csv( dir_file_input, sep = '\t', header = None ).values :
        int_window_size = arr_input[ 0 ]
        l_dir_file = glob.glob( f"{dir_folder_pipeline_temp}*_window_size_{int_window_size}.tsv.gz" ) # retrieve list of output files to combine into a single output file
        dir_file_output_combining = f'{dir_folder_pipeline}b_cell.subsequence__window_size_{int_window_size}.combining.tsv.gz'
        dir_file_output_combining_completed = f'{dir_folder_pipeline}b_cell.subsequence__window_size_{int_window_size}.tsv.gz'
        OS_FILE_Combine_Files_in_order( l_dir_file, dir_file_output_combining, flag_use_header_from_first_file = True, remove_n_lines = 1, delete_input_files = True )
        print( f"combining output files for window size {int_window_size} is completed" )
        os.rename( dir_file_output_combining, dir_file_output_combining_completed ) # rename the file once completed
def Predict_B_cell_cross_reactivity( dir_folder_pipeline = None, dir_folder_pipeline_temp = None, n_threads = 1, l_window_size = [ 30 ], float_thres_e_value = 30, flag_only_use_structural_properties_of_query_proteins = False, float_thres_avg_score_blosum_weighted__b_cell = 0.15, float_thres_avg_score_blosum__b_cell = 0, float_thres_rsa_correlation = 0 ) :
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package

    
    """
    Predict_B_cell_cross_reactivity
    """
    dir_file_matched = f'{dir_folder_pipeline}matched.tsv.gz'
    df_matched = pd.read_csv( dir_file_matched, sep = '\t' ) # read alignments between query and target protein sequences
    df_matched.index.name = 'id_alignment' # retrieve id_alignment (index of df_matched) during retrieving subsequences
    df_matched.reset_index( drop = False, inplace = True ) # add id_alignment column to to the dataframe
    print( f"number of records: {len( df_matched )}" )
    df_matched = df_matched[ df_matched.e_value <= float_thres_e_value ] # drop entries with too low global similarity
#     print( f"number of records after filtering: {len( df_matched )}" )
    
    # update global variables
    global dict_acc_to_arr_acc_query, dict_acc_to_arr_datatype_acc_query, dict_acc_to_arr_phi_query, dict_acc_to_arr_psi_query, dict_acc_to_arr_ss8_query, dict_acc_to_arr_structure_id_query, dict_acc_to_arr_acc_target, dict_acc_to_arr_datatype_acc_target, dict_acc_to_arr_phi_target, dict_acc_to_arr_psi_target, dict_acc_to_arr_ss8_target, dict_acc_to_arr_structure_id_target

    # setting for decoding structural properties
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    # load structural property data for query protein sequences
    df_protein_query = PD_Select( pd.read_csv( f'{dir_folder_pipeline}protein_query.tsv.gz', sep = '\t' ), id_protein = df_matched.query_accession.unique( ) ) # subset structural property database with id_proteins in the alignments
    dict_acc_to_arr_acc_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ].dropna( ).to_dict( ), ** dict_kw_rsa )
    dict_acc_to_arr_datatype_acc_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ].dropna( ).to_dict( ), ** dict_kw_datatype )
    dict_acc_to_arr_phi_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
    dict_acc_to_arr_psi_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
    dict_acc_to_arr_ss8_query = ASCII_Decode( df_protein_query.set_index( 'id_protein' )[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].dropna( ).to_dict( ), ** dict_kw_ss8 )
    dict_fasta = df_protein_query.set_index( 'id_protein' )[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
    dict_acc_to_arr_structure_id_query = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )
    del dict_fasta, df_protein_query
    
    # load structural property data for target protein sequences
    if not flag_only_use_structural_properties_of_query_proteins :
        df_protein_target = PD_Select( pd.read_csv( f'{dir_folder_pipeline}protein_target.tsv.gz', sep = '\t' ), id_protein = df_matched.target_accession.unique( ) ) # subset structural property database with id_proteins in the alignments
        dict_acc_to_arr_acc_target = ASCII_Decode( df_protein_target.set_index( 'id_protein' )[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ].dropna( ).to_dict( ), ** dict_kw_rsa )
        dict_acc_to_arr_datatype_acc_target = ASCII_Decode( df_protein_target.set_index( 'id_protein' )[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ].dropna( ).to_dict( ), ** dict_kw_datatype )
        dict_acc_to_arr_phi_target = ASCII_Decode( df_protein_target.set_index( 'id_protein' )[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
        dict_acc_to_arr_psi_target = ASCII_Decode( df_protein_target.set_index( 'id_protein' )[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle )
        dict_acc_to_arr_ss8_target = ASCII_Decode( df_protein_target.set_index( 'id_protein' )[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].dropna( ).to_dict( ), ** dict_kw_ss8 )
        dict_fasta = df_protein_target.set_index( 'id_protein' )[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
        dict_acc_to_arr_structure_id_target = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )
        del dict_fasta, df_protein_target
    else :
        ''' retrieve target protein sequence exist in the alignments '''
        set_id_protein_target = set( df_matched.target_accession.unique( ) )
        dict_fasta_protein_target = FASTA_Read( f"{dir_folder_pipeline}protein_target.fasta" )
        dict_fasta_protein_target = dict( ( h.split( ' ', 1 )[ 0 ], dict_fasta_protein_target[ h ] ) for h in dict_fasta_protein_target if h.split( ' ', 1 )[ 0 ] in set_id_protein_target ) # retrieve sequence_id by spliting the header at the first space (to makes sequence_id consistent with that used with blastp) # subset protein sequences with id_proteins in the alignments
        dict_arr_float_empty = dict( ( h, np.full( len( dict_fasta_protein_target[ h ] ), np.nan ) ) for h in dict_fasta_protein_target )
        
        ''' initialize structural properties of target proteins with empty arrays '''
        dict_acc_to_arr_acc_target = deepcopy( dict_arr_float_empty )
        dict_acc_to_arr_datatype_acc_target = deepcopy( dict_arr_float_empty )
        dict_acc_to_arr_phi_target = deepcopy( dict_arr_float_empty )
        dict_acc_to_arr_psi_target = deepcopy( dict_arr_float_empty )
        dict_acc_to_arr_ss8_target = deepcopy( dict_arr_float_empty )
        dict_acc_to_arr_structure_id_target = dict( ( h, np.full( len( dict_fasta_protein_target[ h ] ), np.nan, dtype = object ) ) for h in dict_fasta_protein_target )
        del set_id_protein_target, dict_fasta_protein_target, dict_arr_float_empty
    
    """ predict b-cell cross-reactivity """
    l_uuid_process = Multiprocessing( df_matched, BCell_Calculate_Similarity_Scores_in_Aligned_Sequences, n_threads, dir_temp = dir_folder_pipeline_temp, global_arguments = [ float_thres_avg_score_blosum_weighted__b_cell, float_thres_avg_score_blosum__b_cell, float_thres_rsa_correlation, l_window_size, dir_folder_cressp, dir_folder_pipeline, dir_folder_pipeline_temp, flag_only_use_structural_properties_of_query_proteins ] ) # process similarity search result with multiple processes, and collect uuid of the processes
    # combine output files for each window size
    Multiprocessing( l_window_size, BCell_Combine_result_files_for_each_window_size, n_threads = min( len( l_window_size ), n_threads ), dir_temp = dir_folder_pipeline_temp, global_arguments = [ dir_folder_pipeline, dir_folder_pipeline_temp ] ) # combine result files for each window_size
    

    #     # Bin similarity scores by acc_query and a given binning size for analysis
    #     size_window_binning = 100
    #     size_overlap_binning = 50

    #     Multiprocessing( l_window_size, Bin_Similarity_Scores, n_threads = min( len( l_window_size ), int( OS_Memory( )[ 'MemAvailable' ] / 1e7 ) ), dir_temp = dir_folder_pipeline_temp ) # combine result files for each window_size
    
    
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


    