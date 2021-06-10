from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

# prepare B-cell CrossReactivity data for Web-based visualization application (step1: retrieve aligned positions to structures) # 2020-08-13 22:42:15 
# retrieve id_structure_alignment and structure_alignment_position for record with structure evidence
def Retrieve_Alignment( arr_data, query_start, query_end ) :
    """ 'arr_data' contains folloing columns [ _qaccver, _saccver, _pident, _length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, _bitscore, _btop, _query_seq_aligned, _subject_seq_aligned, _id_alignment ]
    pick 'arr' in 'arr_data' that overlap most significantly with the given pair of 'query_start', 'query_end' and return 1-based positions of corrected 'query_start', 'query_end' (so that 'query_start', 'query_end' is covered by the alignment) and 1-based position of corrected 'query_start', 'query_end' in the alignment along with _id_alignment and max_overlap.
    return: query_start, query_end, max_overlap, _id_alignment, alignment_start, alignment_end """
    l_value = list( )
    for arr in arr_data : # retrieve alignment to structure that overlaps with the current alignment 
        _qaccver, _saccver, _pident, _length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, _bitscore, _btop, _query_seq_aligned, _subject_seq_aligned, _id_alignment = arr
        l_value.append( INTERVAL_Overlap( ( _qstart, _qend ), ( query_start, query_end ), flag_0_based_coordinate_system = False, flag_sort_to_retrieve_start_and_end = False ) )
    index_max_overlap = np.argmax( l_value ) # index of record with maximum overlap with the query/target subsequence in 'arr_data'
    max_overlap = l_value[ index_max_overlap ] # maximum overlap with the query/target subsequence
    _qaccver, _saccver, _pident, _length, _mismatch, _gapopen, _qstart, _qend, _sstart, _send, _evalue, _bitscore, _btop, _query_seq_aligned, _subject_seq_aligned, _id_alignment = arr_data[ index_max_overlap ] # pick the alignment record with the largest overlap with the query subsequence
    if max_overlap < query_end - query_start + 1 : query_start, query_end = max( query_start, _qstart ), min( query_end, _qend ) # if max_overlap is smaller than the length of subsequence of the query sequence, reassign 'query_start', 'query_end' based on the query subsequence in the retrieved alignment record.
    len_alignment = len( _query_seq_aligned ) # length of the allignment 
    if _gapopen == 0 : 
        alignment_start, alignment_end = max( 1, query_start - _qstart + 1 ), min( len_alignment, query_end - _qstart + 1 ) # Retrieve 1-based coordinates of the query/target subsequence on the alignment # when there is no gap in the alignment.
        subject_start, subject_end = _sstart + alignment_start - 1, _sstart + alignment_end - 1 # retrieve 1-based coordinates of the subject subsequence on the alignment, when there is no gap in the alignment
    else : # Retrieve 1-based coordinates of the query/target subsequence on the alignment # when there are gaps in the alignment.
        dict_query_seq_1_based_to_alignment_1_based, dict_query_seq_1_based_to_subject_seq_1_based = dict( ), dict( ) # build dictionary for converting a 1-based coordinate of query_seq to a 1-based coordinate of the alignment and the subject sequence
        int_pos_query, int_pos_subject = _qstart - 1, _sstart - 1 # initialize the positions of query and subject sequences
        for index_alignment, residue_query, residue_subject in zip( np.arange( len( _query_seq_aligned ) ), _query_seq_aligned, _subject_seq_aligned ) : # iterate through the alignment
            if residue_query == '-' :
                int_pos_subject += 1
            elif residue_subject == '-' :
                int_pos_query += 1
                dict_query_seq_1_based_to_alignment_1_based[ int_pos_query ] = index_alignment + 1 # 1-based coordinate
                dict_query_seq_1_based_to_subject_seq_1_based[ int_pos_query ] = int_pos_subject
            else :
                int_pos_subject += 1
                int_pos_query += 1
                dict_query_seq_1_based_to_alignment_1_based[ int_pos_query ] = index_alignment + 1 # 1-based coordinate
                dict_query_seq_1_based_to_subject_seq_1_based[ int_pos_query ] = int_pos_subject
        alignment_start, alignment_end, subject_start, subject_end = dict_query_seq_1_based_to_alignment_1_based[ query_start ], dict_query_seq_1_based_to_alignment_1_based[ query_end ], dict_query_seq_1_based_to_subject_seq_1_based[ query_start ], dict_query_seq_1_based_to_subject_seq_1_based[ query_end ]
    return max_overlap, _id_alignment, query_start, query_end, subject_start, subject_end, alignment_start, alignment_end

def Retrieve_Overlapping_Structures( dir_file_input, name_file, dir_folder_pipeline ) :
    """ retrieve id_structure of most widely overlapping structures """
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # a working directory of estimating structural properties
    
    # read input file
    df = pd.read_csv( dir_file_input, sep = '\t' )
    
    # load pdb aligned to target proteins
    df_blastp_pdb_target = pd.read_csv( f'{dir_folder_pipeline_struc}protein_target.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz', sep = '\t' )
    df_blastp_pdb_target = PD_Select( df_blastp_pdb_target, qaccver = set( df.target_accession.values ) ) # subset for the target accession in the current input
    df_blastp_pdb_target[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_target ) ) # retrieve integer index (id) of each alignment
    # load pdb aligned to query proteins
    df_blastp_pdb_query = pd.read_csv( f'{dir_folder_pipeline_struc}protein_query.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz', sep = '\t' )
    df_blastp_pdb_query = PD_Select( df_blastp_pdb_query, qaccver = set( df.query_accession.values ) ) # subset for the query accession in the current input
    df_blastp_pdb_query[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_query ) ) # retrieve integer index (id) of each alignment
   
    # retrieve dictionary indices of the blastp alignment dataframes for more efficient accessing of rows
    dict_index_blastp_pdb_query = DF_Build_Index_Using_Dictionary( df_blastp_pdb_query, [ 'qaccver', 'saccver' ] )
    arr_data_blastp_pdb_query = df_blastp_pdb_query.values
    dict_index_blastp_pdb_target = DF_Build_Index_Using_Dictionary( df_blastp_pdb_target, [ 'qaccver', 'saccver' ] )
    arr_data_blastp_pdb_target = df_blastp_pdb_target.values

    
    ''' search maximum overlap with RCSB_PDB structures for each record in the input data '''
    l_col = [ 'window_size', 'id_alignment', 'source', 'query_accession', 'target_accession', 'e_value', 'identity', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'score_blosum', 'score_blosum_weighted', 'sum_of_weights', 'score_similarity_acc', 'score_similarity_phi', 'score_similarity_psi', 'score_similarity_ss8', 'n_residues_acc', 'n_residues_phi', 'n_residues_psi', 'n_residues_ss8', 'correl_coeffi_acc', 'correl_p_value_acc', 'correl_coeffi_phi', 'correl_p_value_phi', 'correl_coeffi_psi', 'correl_p_value_psi', 'prop_pdb_evidence_query', 'prop_pdb_evidence_target', 'structure_id_query', 'count_structure_id_query', 'structure_id_target', 'count_structure_id_target', 'most_frequent_ss8_query', 'count_most_frequent_ss8_query', 'most_frequent_ss8_target', 'count_most_frequent_ss8_target' ] 
    l_l_value = list( )
    for window_size, id_alignment, source, query_accession, target_accession, e_value, identity, alignment_start, alignment_end, query_start, query_end, target_start, target_end, query_subsequence, target_subsequence, score_blosum, score_blosum_weighted, sum_of_weights, score_similarity_acc, score_similarity_phi, score_similarity_psi, score_similarity_ss8, n_residues_acc, n_residues_phi, n_residues_psi, n_residues_ss8, correl_coeffi_acc, correl_p_value_acc, correl_coeffi_phi, correl_p_value_phi, correl_coeffi_psi, correl_p_value_psi, prop_pdb_evidence_query, prop_pdb_evidence_target, structure_id_query, count_structure_id_query, structure_id_target, count_structure_id_target, most_frequent_ss8_query, count_most_frequent_ss8_query, most_frequent_ss8_target, count_most_frequent_ss8_target,  in df[ l_col ].values :
        max_overlap_query, id_alignment_structure_query, query_structure_start, query_structure_end, structure_id_query_start, structure_id_query_end, alignment_structure_query_start, alignment_structure_query_end = Retrieve_Alignment( arr_data_blastp_pdb_query[ dict_index_blastp_pdb_query[ query_accession, structure_id_query ] ], query_start, query_end ) if isinstance( structure_id_query, str ) else ( np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan )
        max_overlap_target, id_alignment_structure_target, target_structure_start, target_structure_end, structure_id_target_start, structure_id_target_end, alignment_structure_target_start, alignment_structure_target_end = Retrieve_Alignment( arr_data_blastp_pdb_target[ dict_index_blastp_pdb_target[ target_accession, structure_id_target ] ], target_start, target_end ) if isinstance( structure_id_target, str ) else ( np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan )
        l_l_value.append( [ query_structure_start, query_structure_end, structure_id_query_start, structure_id_query_end, max_overlap_query, id_alignment_structure_query, alignment_structure_query_start, alignment_structure_query_end, target_structure_start, target_structure_end, structure_id_target_start, structure_id_target_end, max_overlap_target, id_alignment_structure_target, alignment_structure_target_start, alignment_structure_target_end ] )
    df = df.join( pd.DataFrame( l_l_value, columns = [ 'query_structure_start', 'query_structure_end', 'structure_id_query_start', 'structure_id_query_end', 'max_overlap_query', 'id_alignment_structure_query', 'alignment_structure_query_start', 'alignment_structure_query_end', 'target_structure_start', 'target_structure_end', 'structure_id_target_start', 'structure_id_target_end', 'max_overlap_target', 'id_alignment_structure_target', 'alignment_structure_target_start', 'alignment_structure_target_end' ] ) )
    df.to_csv( f"{dir_file_input}_output.tsv.gz", sep = '\t', index = False )
    
def Retrieve_Overlapping_Structures__PostProcessing( str_uuid, dir_temp, name_file, dir_folder_pipeline ) :
    dir_folder_pipeline_web = f'{dir_folder_pipeline}web_application/' # a working directory for exporting data for web applications
    
    l_dir_file = glob.glob( f"{dir_temp}{str_uuid}*_output.tsv.gz" ) # retrieve list of output files
    dir_file_combined = f'{dir_folder_pipeline_web}b_cell.{name_file}.1.position_aligned_structure_added.tsv.gz'
    OS_FILE_Combine_Files_in_order( l_dir_file, dir_file_combined, overwrite_existing_file = True, flag_use_header_from_first_file = True, remove_n_lines = 1 )

def Prepare_data_for_web_application( dir_file, dir_folder_pipeline, dir_folder_pipeline_temp, n_threads ) :
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp2/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package
    
    # define folder directories for exporting data for web applications
    dir_folder_pipeline_web = f'{dir_folder_pipeline}web_application/' # a working directory for exporting data for web applications
    os.makedirs( dir_folder_pipeline_web, exist_ok = True )
    
    """ Prepare_data_for_web_application """
    name_file = dir_file.rsplit( '/', 1 )[ 1 ].rsplit( '.tsv', 1 )[ 0 ] # retrieve name of the file
    
    """ add full-length alignments to RCSB_PDB structures """
    df = pd.read_csv( dir_file, sep = '\t', low_memory = False )
    Multiprocessing( df, Function = Retrieve_Overlapping_Structures, n_threads = n_threads, Function_PostProcessing = Retrieve_Overlapping_Structures__PostProcessing, dir_temp = dir_folder_pipeline_temp, global_arguments = [ name_file, dir_folder_pipeline ] )
    del df

    # modify aligned coordinates of PDB structures so that it accurately match that in 'label_seq_id' (PDB sequence is often fragment of its parent proteins) 
    df_subsequence_pdb_web = pd.read_csv( f"{dir_folder_pipeline_web}b_cell.{name_file}.1.position_aligned_structure_added.tsv.gz", sep = '\t' )
    Map = MAP.Map( pd.read_csv( f"{dir_folder_cressp}data/pdb/rcsb_pdb.label_seq_id.start_end.tsv.gz", sep = '\t' ).set_index( 'structure_id' ).int_index_residue_start.to_dict( ) ) # read structure_id -> 'label_seq_id' start position mapping
    df_subsequence_pdb_web.structure_id_query_start = df_subsequence_pdb_web.structure_id_query_start + df_subsequence_pdb_web.structure_id_query.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_query_end = df_subsequence_pdb_web.structure_id_query_end + df_subsequence_pdb_web.structure_id_query.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_target_start = df_subsequence_pdb_web.structure_id_target_start + df_subsequence_pdb_web.structure_id_target.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_target_end = df_subsequence_pdb_web.structure_id_target_end + df_subsequence_pdb_web.structure_id_target.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.to_csv( f'{dir_folder_pipeline_web}b_cell.{name_file}.2.residue_pos_corrected.tsv.gz', sep = '\t', index = False ) # save an intermediate version
    
    