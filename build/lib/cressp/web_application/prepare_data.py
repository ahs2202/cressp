from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

# define read-only global variables during multiprocessing
dict_index_blastp_pdb_query, arr_data_blastp_pdb_query, dict_index_blastp_pdb_target, arr_data_blastp_pdb_target = dict( ), dict( ), dict( ), dict( )

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
def Retrieve_Overlapping_Structures( dir_file_input, dir_file_b_cell, name_file_b_cell, dir_folder_pipeline, dir_folder_pipeline_temp ) :
    """ retrieve id_structure of most widely overlapping structures """
    global dict_index_blastp_pdb_query, arr_data_blastp_pdb_query, dict_index_blastp_pdb_target, arr_data_blastp_pdb_target # use global variables
    
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # a working directory of estimating structural properties
    
    # reetrieve start and end indices of lines
    int_index_line_start, int_index_line_end = pd.read_csv( dir_file_input, sep = '\t', header = None ).values[ 0 ]
    
    str_uuid = UUID( ) # set uuid of the current process
    
    ''' initialize the output file '''
    l_col = [ 'window_size', 'id_alignment', 'source', 'query_accession', 'target_accession', 'e_value', 'identity', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'score_blosum', 'score_blosum_weighted', 'sum_of_weights', 'score_similarity_acc', 'score_similarity_phi', 'score_similarity_psi', 'score_similarity_ss8', 'n_residues_acc', 'n_residues_phi', 'n_residues_psi', 'n_residues_ss8', 'correl_coeffi_acc', 'correl_p_value_acc', 'correl_coeffi_phi', 'correl_p_value_phi', 'correl_coeffi_psi', 'correl_p_value_psi', 'prop_pdb_evidence_query', 'prop_pdb_evidence_target', 'structure_id_query', 'count_structure_id_query', 'structure_id_target', 'count_structure_id_target', 'most_frequent_ss8_query', 'count_most_frequent_ss8_query', 'most_frequent_ss8_target', 'count_most_frequent_ss8_target', 'query_structure_start', 'query_structure_end', 'structure_id_query_start', 'structure_id_query_end', 'max_overlap_query', 'id_alignment_structure_query', 'alignment_structure_query_start', 'alignment_structure_query_end', 'target_structure_start', 'target_structure_end', 'structure_id_target_start', 'structure_id_target_end', 'max_overlap_target', 'id_alignment_structure_target', 'alignment_structure_target_start', 'alignment_structure_target_end' ]
    newfile = gzip.open( f"{dir_folder_pipeline_temp}{str_uuid}.{name_file_b_cell}.1.position_aligned_structure_added.tsv.gz", 'wb' ) # open an output file
    newfile.write( ( '\t'.join( l_col ) + '\n' ).encode( ) )
    
    ''' search maximum overlap with RCSB_PDB structures for each record in the input data '''
    index_line = -1
    with gzip.open( dir_file_b_cell, 'rb' ) as file :
        file.readline( ) # read header
        while True :
            line = file.readline( )
            if len( line ) == 0 :
                break
            index_line += 1 # increase line index number by 1
            if index_line < int_index_line_start : # read until the target region is reached
                continue
            elif index_line >= int_index_line_end : # stop reading if processing of the target region is completed
                break
            ''' parse a single record '''
            line_without_newline = line.decode( ).strip( ) # retrieve line without a newline character
            window_size, id_alignment, source, query_accession, target_accession, e_value, identity, alignment_start, alignment_end, query_start, query_end, target_start, target_end, query_subsequence, target_subsequence, score_blosum, score_blosum_weighted, sum_of_weights, score_similarity_acc, score_similarity_phi, score_similarity_psi, score_similarity_ss8, n_residues_acc, n_residues_phi, n_residues_psi, n_residues_ss8, correl_coeffi_acc, correl_p_value_acc, correl_coeffi_phi, correl_p_value_phi, correl_coeffi_psi, correl_p_value_psi, prop_pdb_evidence_query, prop_pdb_evidence_target, structure_id_query, count_structure_id_query, structure_id_target, count_structure_id_target, most_frequent_ss8_query, count_most_frequent_ss8_query, most_frequent_ss8_target, count_most_frequent_ss8_target = Parse_Line( line_without_newline, [ int, int, str, str, str, float, float, int, int, int, int, int, int, str, str, int, float, float, float, float, float, float, int, int, int, int, float, float, float, float, float, float, float, float, str, int, str, int, str, int, str, int ], delimiter = '\t' )
            ''' search maximum overlap with RCSB_PDB structures for each record in the input data '''
            max_overlap_query, id_alignment_structure_query, query_structure_start, query_structure_end, structure_id_query_start, structure_id_query_end, alignment_structure_query_start, alignment_structure_query_end = Retrieve_Alignment( arr_data_blastp_pdb_query[ dict_index_blastp_pdb_query[ query_accession, structure_id_query ] ], query_start, query_end ) if isinstance( structure_id_query, str ) else ( np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan )
            max_overlap_target, id_alignment_structure_target, target_structure_start, target_structure_end, structure_id_target_start, structure_id_target_end, alignment_structure_target_start, alignment_structure_target_end = Retrieve_Alignment( arr_data_blastp_pdb_target[ dict_index_blastp_pdb_target[ target_accession, structure_id_target ] ], target_start, target_end ) if isinstance( structure_id_target, str ) else ( np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan )
            ''' write processed result '''
            newfile.write( ( line_without_newline + '\t' + '\t'.join( list( map( str, [ query_structure_start, query_structure_end, structure_id_query_start, structure_id_query_end, max_overlap_query, id_alignment_structure_query, alignment_structure_query_start, alignment_structure_query_end, target_structure_start, target_structure_end, structure_id_target_start, structure_id_target_end, max_overlap_target, id_alignment_structure_target, alignment_structure_target_start, alignment_structure_target_end ] ) ) ) + '\n' ).encode( ) )
def Prepare_data_for_web_application( dir_file_b_cell, dir_file_t_cell, dict_cressp_setting ) :
    """ Prepare_data_for_web_application """
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package
    
    dir_folder_output = dict_cressp_setting[ 'dir_folder_output' ]
    dir_folder_pipeline = dict_cressp_setting[ 'dir_folder_pipeline' ]
    dir_folder_pipeline_temp = dict_cressp_setting[ 'dir_folder_pipeline_temp' ]
    n_threads = dict_cressp_setting[ 'n_threads' ]
    dir_folder_pipeline_web = dict_cressp_setting[ 'dir_folder_pipeline_web' ]
    dir_folder_pipeline_struc = dict_cressp_setting[ 'dir_folder_pipeline_struc' ]
    dir_folder_web = dict_cressp_setting[ 'dir_folder_web' ]
    flag_deduplicate_based_on_aligned_subsequences_for_visualization = dict_cressp_setting[ 'flag_deduplicate_based_on_aligned_subsequences_for_visualization' ]
    
    # create folders if they do not exist
    for dir_folder in [ dir_folder_pipeline_web, dir_folder_pipeline_struc, dir_folder_web ] :
        os.makedirs( dir_folder_pipeline_web, exist_ok = True )
    
    ''' deduplicate records based on aligned subsequences '''
    if flag_deduplicate_based_on_aligned_subsequences_for_visualization :
        ''' deduplicate b-cell records '''
        dir_file_b_cell_deduplicated = f"{dir_file_b_cell.rsplit( '.tsv.gz', 1 )[ 0 ]}.deduplicated_by_subsequences.tsv.gz" # define an output file
        DF_Deduplicate_without_loading_in_memory( dir_file_b_cell, dir_file_b_cell_deduplicated, [ 'query_subsequence', 'target_subsequence' ], flag_header_is_present = True, str_delimiter = '\t' )
        dir_file_b_cell = dir_file_b_cell_deduplicated # use the deduplicated b-cell record file 
    
    """ Retrieve file names """
    name_file_b_cell = dir_file_b_cell.rsplit( '/', 1 )[ 1 ].rsplit( '.tsv', 1 )[ 0 ] # retrieve name of the input file of b-cell
    name_file_t_cell = dir_file_t_cell.rsplit( '/', 1 )[ 1 ].rsplit( '.tsv', 1 )[ 0 ] # retrieve name of the input file of t-cell
    # add to the setting dictionary
    dict_cressp_setting[ 'name_file_b_cell' ] = name_file_b_cell
    dict_cressp_setting[ 'name_file_t_cell' ] = name_file_t_cell
    
    """ 
    Process B-cell data 
    """
    """ add full-length alignments to RCSB_PDB structures """
    # update global read-only variables
    global dict_index_blastp_pdb_query, arr_data_blastp_pdb_query, dict_index_blastp_pdb_target, arr_data_blastp_pdb_target # use global variables
    df = pd.read_csv( dir_file_b_cell, sep = '\t', usecols = [ 'target_accession', 'query_accession' ] ) # retrieve target & query accessions
    
    # load pdb aligned to target proteins
    df_blastp_pdb_target = pd.read_csv( f'{dir_folder_pipeline_struc}protein_target.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz', sep = '\t' )
    df_blastp_pdb_target[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_target ) ) # retrieve integer index (id) of each alignment
    df_blastp_pdb_target = PD_Select( df_blastp_pdb_target, qaccver = set( df.target_accession.values ) ) # subset for the target accession in the current input
    
    # load pdb aligned to query proteins
    df_blastp_pdb_query = pd.read_csv( f'{dir_folder_pipeline_struc}protein_query.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz', sep = '\t' )
    df_blastp_pdb_query[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_query ) ) # retrieve integer index (id) of each alignment
    df_blastp_pdb_query = PD_Select( df_blastp_pdb_query, qaccver = set( df.query_accession.values ) ) # subset for the query accession in the current input
   
    # retrieve dictionary indices of the blastp alignment dataframes for more efficient accessing of rows
    dict_index_blastp_pdb_query = DF_Build_Index_Using_Dictionary( df_blastp_pdb_query, [ 'qaccver', 'saccver' ] )
    arr_data_blastp_pdb_query = df_blastp_pdb_query.values
    dict_index_blastp_pdb_target = DF_Build_Index_Using_Dictionary( df_blastp_pdb_target, [ 'qaccver', 'saccver' ] )
    arr_data_blastp_pdb_target = df_blastp_pdb_target.values
    
    # build bookmarks for simultaneous access of the file
    int_n_records = len( df )
    l = list( range( 0, int_n_records, int( np.ceil( int_n_records / n_threads ) ) ) ) + [ int_n_records ] # list of indices (for building bookmarks)
    l_bookmarks = list( [ l[ i ], l[ i + 1 ] ] for i in range( len( l ) - 1 ) ) # retrieve list of bookmarks
    del df
    
    Multiprocessing( l_bookmarks, Function = Retrieve_Overlapping_Structures, n_threads = len( l_bookmarks ), dir_temp = dir_folder_pipeline_temp, global_arguments = [ dir_file_b_cell, name_file_b_cell, dir_folder_pipeline, dir_folder_pipeline_temp ] )
    # combine output files
    OS_FILE_Combine_Files_in_order( glob.glob( f"{dir_folder_pipeline_temp}*.{name_file_b_cell}.1.position_aligned_structure_added.tsv.gz" ), f'{dir_folder_pipeline_web}{name_file_b_cell}.1.position_aligned_structure_added.tsv.gz', overwrite_existing_file = True, flag_use_header_from_first_file = True, remove_n_lines = 1 )


    """ Modify coordinates of RCSB_PDB structures """
    # modify aligned coordinates of PDB structures so that it accurately match that in 'label_seq_id' (PDB sequence is often fragment of its parent proteins) 
    df_subsequence_pdb_web = pd.read_csv( f"{dir_folder_pipeline_web}{name_file_b_cell}.1.position_aligned_structure_added.tsv.gz", sep = '\t' )
    PKG.Download_Data( "data/pdb/rcsb_pdb.label_seq_id.start_end.tsv.gz", dir_remote, name_package ) # download data
    Map = MAP.Map( pd.read_csv( f"{dir_folder_cressp}data/pdb/rcsb_pdb.label_seq_id.start_end.tsv.gz", sep = '\t' ).set_index( 'structure_id' ).int_index_residue_start.to_dict( ) ) # read structure_id -> 'label_seq_id' start position mapping
    df_subsequence_pdb_web.structure_id_query_start = df_subsequence_pdb_web.structure_id_query_start + df_subsequence_pdb_web.structure_id_query.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_query_end = df_subsequence_pdb_web.structure_id_query_end + df_subsequence_pdb_web.structure_id_query.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_target_start = df_subsequence_pdb_web.structure_id_target_start + df_subsequence_pdb_web.structure_id_target.apply( Map.a2b ) - 1
    df_subsequence_pdb_web.structure_id_target_end = df_subsequence_pdb_web.structure_id_target_end + df_subsequence_pdb_web.structure_id_target.apply( Map.a2b ) - 1
    # calculate score for sorting
    df_subsequence_pdb_web[ 'score_for_sorting' ] = ( df_subsequence_pdb_web.correl_coeffi_acc * df_subsequence_pdb_web.score_blosum_weighted ).fillna( -1 ) 
    # save an intermediate version
    df_subsequence_pdb_web.to_csv( f'{dir_folder_pipeline_web}{name_file_b_cell}.2.residue_pos_corrected.tsv.gz', sep = '\t', index = False ) 

    """ 
    Process T-cell data 
    """

    """ 
    Export skeletons of B-Tell and T-cell CrossReactivity data for Web application 
    """
    # prepare B-Tell and T-cell CrossReactivity data for Web-based visualization application (step2: encode bulky data with integer indices and write compact data for the web application) # 2021-01-05 14:44:49 
    str_col_id = 'id_protein'
    n_record_max = 500000

    df_subsequence_pdb_web = pd.read_csv( f'{dir_folder_pipeline_web}{name_file_b_cell}.2.residue_pos_corrected.tsv.gz', sep = '\t', low_memory = False )
    df_mhc_web = pd.read_csv( dir_file_t_cell, sep = '\t' )
    df_matched = pd.read_csv( f"{dir_folder_pipeline}matched.tsv.gz", sep = '\t' )
    df_blastp_pdb_query = pd.read_csv( f"{dir_folder_pipeline_struc}protein_query.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz", sep = '\t' )
    df_blastp_pdb_target = pd.read_csv( f"{dir_folder_pipeline_struc}protein_target.blastp_rcsb_pdb.with_aligned_seq.filtered.tsv.gz", sep = '\t' )
    df_blastp_pdb_target[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_target ) ) # retrieve integer index (id) of each alignment
    df_blastp_pdb_query[ 'id_alignment' ] = np.arange( len( df_blastp_pdb_query ) ) # retrieve integer index (id) of each alignment


    # subset alignments for records in B-cell CrossReactivity data
    # count alignment_id and sort alignment_id based on counts, so that the id_alignment with the largest counts are assigned with smaller integer index, which can reduce the file size
    df_subsequence_pdb_web = df_subsequence_pdb_web[ df_subsequence_pdb_web.correl_coeffi_acc > 0 ] # retrieve only significant records

    s_count_id_alignment_structure_target = LIST_COUNT( df_subsequence_pdb_web.id_alignment_structure_target.dropna( ).astype( int ), duplicate_filter = None ) 
    s_count_id_alignment_structure_query = LIST_COUNT( df_subsequence_pdb_web.id_alignment_structure_query.dropna( ).astype( int ), duplicate_filter = None )
    s_count_id_alignment = LIST_COUNT( df_subsequence_pdb_web.id_alignment.dropna( ).astype( int ), duplicate_filter = None )

    # 'df_blastp_pdb_target' should have been loaded from the previous cells
    df_blastp_pdb_target_subset = PD_Select( df_blastp_pdb_target, id_alignment = s_count_id_alignment_structure_target.index.values ).set_index( 'id_alignment' ).loc[ s_count_id_alignment_structure_target.index.values ].reset_index( )
    df_blastp_pdb_query_subset = PD_Select( df_blastp_pdb_query, id_alignment = s_count_id_alignment_structure_query.index.values ).set_index( 'id_alignment' ).loc[ s_count_id_alignment_structure_query.index.values ].reset_index( )
    df_matched_subset = df_matched.loc[ s_count_id_alignment.index.values ]
    df_matched_subset.index.name = 'id_alignment'
    df_matched_subset.reset_index( inplace = True )
    # save intermediate results
    df_blastp_pdb_target_subset.to_csv( f'{dir_folder_pipeline_web}alignment_target_pdb.source.tsv.gz', sep = '\t', index = False )
    df_blastp_pdb_query_subset.to_csv( f'{dir_folder_pipeline_web}alignment_query_pdb.source.tsv.gz', sep = '\t', index = False )
    df_matched_subset.to_csv( f'{dir_folder_pipeline_web}alignment_query_target.source.tsv.gz', sep = '\t', index = False )
    # rename columns and reset index to set new id_alignment
    df_blastp_pdb_target_subset = df_blastp_pdb_target_subset.rename( columns = { 'query_seq_aligned' : 'query_alignment', 'subject_seq_aligned' : 'target_alignment' } ).reset_index( drop = True )
    df_blastp_pdb_query_subset = df_blastp_pdb_query_subset.rename( columns = { 'query_seq_aligned' : 'query_alignment', 'subject_seq_aligned' : 'target_alignment' } ).reset_index( drop = True )
    # save compact dataframe for web application (id_alignment in the B-Cell CrossReactivity data contain 0-based row-index of each dataframe)
    df_blastp_pdb_target_subset[ [ 'query_alignment', 'target_alignment' ] ].to_csv( f'{dir_folder_pipeline_web}alignment_target_pdb.tsv', sep = '\t', index = False )
    df_blastp_pdb_query_subset[ [ 'query_alignment', 'target_alignment' ] ].to_csv( f'{dir_folder_pipeline_web}alignment_query_pdb.tsv', sep = '\t', index = False )
    df_matched_subset[ [ 'query_alignment', 'target_alignment' ] ].to_csv( f'{dir_folder_pipeline_web}alignment_query_target.tsv', sep = '\t', index = False )
    # replace id_alignment columns in the B-cell CrossReactivity data with 0-based row-index of each dataframe
    mo = MAP.Map( dict( ( id_alignment, index ) for index, id_alignment in enumerate( df_matched_subset.id_alignment.values ) ) )
    df_subsequence_pdb_web.id_alignment = df_subsequence_pdb_web.id_alignment.apply( mo.a2b )
    mo = MAP.Map( dict( ( id_alignment, index ) for index, id_alignment in enumerate( df_blastp_pdb_query_subset.id_alignment.values ) ) )
    df_subsequence_pdb_web.id_alignment_structure_query = df_subsequence_pdb_web.id_alignment_structure_query.apply( mo.a2b )
    mo = MAP.Map( dict( ( id_alignment, index ) for index, id_alignment in enumerate( df_blastp_pdb_target_subset.id_alignment.values ) ) )
    df_subsequence_pdb_web.id_alignment_structure_target = df_subsequence_pdb_web.id_alignment_structure_target.apply( mo.a2b )

    # encode accessions to integers ((1) make dataframe more compact and (2) to aid efficient visualization)
    # read input fasta files (with key=accession)
    dict_fasta_query = FASTA_Read( f'{dir_folder_pipeline}protein_query.fasta', header_split_at_space = True )
    dict_fasta_target = FASTA_Read( f'{dir_folder_pipeline}protein_target.fasta', header_split_at_space = True )
    # retrieve valid target and query accession (if accession has been changed to something else due to handling error throughput the pipeline, they will be filtered out)
    set_query_accession = set( dict_fasta_query )
    set_target_accession = set( dict_fasta_target )
    # retrieve records of valid target and query accessions
    df_subsequence_pdb_web = PD_Select( df_subsequence_pdb_web, query_accession = set_query_accession, target_accession = set_target_accession )
    df_mhc_web = PD_Select( df_mhc_web, query_accession = set_query_accession, target_accession = set_target_accession )

    def __Encode_to_Integer__( l, spread_based_on_weight = True ) :
        ''' encode given list of values to integer.
        'spread_based_on_weight': count the values in the given list to assign integer more effectively (spread integers for efficient visualization on 'rainbow' colormap) '''
        df = LIST_COUNT( l, duplicate_filter = None ).reset_index( ) # count values and sort so that value with largest count is located at the first row
        df.columns = [ 'value', 'n_counts' ]
        df[ 'encoding' ] = INTEGER_Spread( len( df ) ) if spread_based_on_weight else np.arange( len( df ) )
        return df

    df_acc_query = __Encode_to_Integer__( np.concatenate( [ df_subsequence_pdb_web.query_accession.values, df_mhc_web.query_accession.values ] ), spread_based_on_weight = True ).sort_values( 'encoding' ) # encode accession in BCellCrossReactivity data and TCellCrossReactivity into integers
    df_acc_target = __Encode_to_Integer__( np.concatenate( [ df_subsequence_pdb_web.target_accession.values, df_mhc_web.target_accession.values ] ), spread_based_on_weight = True ).sort_values( 'encoding' ) # sort in order to use integer index to access records in the web application
    df_acc_pdb = __Encode_to_Integer__( np.concatenate( [ df_subsequence_pdb_web.structure_id_query.values, df_subsequence_pdb_web.structure_id_target.values ] ), spread_based_on_weight = True ).sort_values( 'encoding' )

    # retrieve fasta header and sequences for query and target proteins in the exported data
    def __Read_Fasta_Header__( dir_file_fasta ) :
        ''' read fasta headers of an unzipped fasta file '''
        l = [ ]
        with open( dir_file_fasta, 'r' ) as file :
            while True :
                line = file.readline( )
                if len( line ) == 0 :
                    break
                if line[ 0 ] == '>' : # identify line containing fasta header
                    l.append( line.strip( )[ 1 : ] )
        return l

    for df_acc, dir_file_fasta in zip( [ df_acc_query, df_acc_target ], [  f'{dir_folder_pipeline}protein_query.fasta',  f'{dir_folder_pipeline}protein_target.fasta' ] ) :
        dict_fasta = FASTA_Read( dir_file_fasta )
        mo = MAP.Map( dict( ( h.split( ' ', 1 )[ 0 ], h ) for h in dict_fasta ) )
        df_acc[ 'fasta_header' ] = df_acc.value.apply( mo.a2b )
        mo = MAP.Map( dict( ( h.split( ' ', 1 )[ 0 ], dict_fasta[ h ] ) for h in dict_fasta ) )
        df_acc[ 'sequence' ] = df_acc.value.apply( mo.a2b )

    # save intermediate results
    df_acc_query.to_csv( f'{dir_folder_pipeline_web}acc_query.source.tsv.gz', sep = '\t', index = False )
    df_acc_target.to_csv( f'{dir_folder_pipeline_web}acc_target.source.tsv.gz', sep = '\t', index = False )
    df_acc_pdb.to_csv( f'{dir_folder_pipeline_web}acc_pdb.source.tsv.gz', sep = '\t', index = False )
    # save accessions for web application
    df_acc_query.to_csv( f'{dir_folder_pipeline_web}acc_query.tsv', columns = [ 'value', 'n_counts', 'fasta_header', 'sequence' ], sep = '\t', index = False )
    df_acc_target.to_csv( f'{dir_folder_pipeline_web}acc_target.tsv', columns = [ 'value', 'n_counts', 'fasta_header', 'sequence' ], sep = '\t', index = False )

    # replace accession_id columns in the B-cell CrossReactivity data with 0-based row-index of each dataframe
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_acc_query.value.values ) ) )
    df_subsequence_pdb_web.query_accession = df_subsequence_pdb_web.query_accession.apply( mo.a2b )
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_acc_target.value.values ) ) )
    df_subsequence_pdb_web.target_accession = df_subsequence_pdb_web.target_accession.apply( mo.a2b )
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_acc_pdb.value.values ) ) )
    df_subsequence_pdb_web.structure_id_query = df_subsequence_pdb_web.structure_id_query.apply( mo.a2b )
    df_subsequence_pdb_web.structure_id_target = df_subsequence_pdb_web.structure_id_target.apply( mo.a2b )


    df_subsequence_pdb_web[ 'source_is_blastp' ] = ( df_subsequence_pdb_web.source == 'blastp' ).astype( int )
    l_col_positive_integer_datatype = [ 'structure_id_query_start', 'structure_id_query_end', 'structure_id_target_start', 'structure_id_target_end', 'source_is_blastp', 'id_alignment', 'query_accession', 'target_accession', 'alignment_start', 'alignment_end', 'query_start', 'query_end', 'target_start', 'target_end', 'n_residues_acc', 'n_residues_phi', 'n_residues_psi', 'n_residues_ss8', 'structure_id_query', 'count_structure_id_query', 'structure_id_target', 'count_structure_id_target', 'window_size', 'query_structure_start', 'query_structure_end', 'max_overlap_query', 'id_alignment_structure_query', 'alignment_structure_query_start', 'alignment_structure_query_end', 'target_structure_start', 'target_structure_end', 'max_overlap_target', 'id_alignment_structure_target', 'alignment_structure_target_start', 'alignment_structure_target_end' ]
    for col in l_col_positive_integer_datatype : df_subsequence_pdb_web[ col ] = df_subsequence_pdb_web[ col ].fillna( -1 ).astype( int ) # nan values in columns containing positive integers will be represented by -1 value
    l_col_integer_datatype = [ 'score_blosum' ]
    for col in l_col_integer_datatype : df_subsequence_pdb_web[ col ] = df_subsequence_pdb_web[ col ].fillna( 0 ).astype( int ) # nan values will be represented by 0 value
    l_col_float_scientific_notation = [ 'e_value', 'correl_p_value_acc', 'correl_p_value_phi', 'correl_p_value_psi' ] # convert very small float values to scientific notation to reduce the size of the file while retainint the accuracy
    for col in l_col_float_scientific_notation : df_subsequence_pdb_web[ col ] = list( '' if np.isnan( value ) else "{:.2e}".format( value ) for value in df_subsequence_pdb_web[ col ].values )
    l_col_optional_data = [ 'score_similarity_acc', 'score_similarity_phi', 'score_similarity_psi', 'n_residues_acc', 'n_residues_phi', 'n_residues_psi', 'n_residues_ss8', 'max_overlap_query', 'max_overlap_target', 'count_structure_id_query', 'count_structure_id_target' ] # list of columns containing non-essential information, which will be dropped to reduce the file size and make the delivery of the dataframe in the web application more efficient # similarity scores of datatypes where correlation coefficients are available were considered optional
    df_subsequence_pdb_web.drop( columns = [ 'source', 'query_subsequence', 'target_subsequence' ] + l_col_optional_data, inplace = True ) # drop entirely unnecessary columns
    # write compact dataframes of the BCellCrossReactivity data for web application
    float_thres_margin = 1.5 # a threshold for skipping writing a file if the number of records for subset is too close to the maximum number
    df_subsequence_pdb_web.to_csv( f'{dir_folder_pipeline_web}BCellCrossReactivity.source.tsv.gz', sep = '\t', index = False ) # save an intermediate result
    for df, name in zip( [ df_subsequence_pdb_web ], [ 'nonredundant_sequence' ] ) : # for each redundancy level
        for int_window_size in df.window_size.unique( ) : # write the BCellCrossReactivity data for each window_size
            df_for_a_window_size = df[ df.window_size == int_window_size ]
            df_for_a_window_size = df_for_a_window_size.sample( n = min( n_record_max, len( df_for_a_window_size ) ) ) # subsample epitopes before sorting with a composite score if current dataset has too many number of records
            df_for_a_window_size.sort_values( 'score_for_sorting', ascending = False, inplace = True ) # sort records by 'score_for_sorting' so that the most significant record is situated at the the first row.
            df_for_a_window_size.drop( columns = [ 'score_for_sorting' ], inplace = True )
            for n_records in [ 1000, 2000, 5000, 10000, 20000, 40000, 80000, 160000 ] : # write the TCellCrossReactivity data for multiple numbers of records
                if len( df_for_a_window_size ) < n_records * float_thres_margin :
                    n_records = len( df_for_a_window_size )
                    df_for_a_window_size.to_csv( f'{dir_folder_pipeline_web}BCellCrossReactivity.{name}.windowSize_{int_window_size}.top_{n_records}.tsv', sep = '\t', index = False, float_format = '%.3f' ) # write only 3 digits below floating point to reduce file size (except for floats in columns in 'l_col_float_scientific_notation')
                    break  # stop subsetting after writing all records as a file
                else : df_for_a_window_size.iloc[ : n_records ].to_csv( f'{dir_folder_pipeline_web}BCellCrossReactivity.{name}.windowSize_{int_window_size}.top_{n_records}.tsv', sep = '\t', index = False, float_format = '%.3f' ) # write only 3 digits below floating point to reduce file size (except for floats in columns in 'l_col_float_scientific_notation')


    """ export prediction result of cross-reactivity of mhc ligands """
    # replace accession_id columns in the T-cell CrossReactivity data with 0-based row-index of each dataframe
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_acc_query.value.values ) ) )
    df_mhc_web.query_accession = df_mhc_web.query_accession.apply( mo.a2b )
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_acc_target.value.values ) ) )
    df_mhc_web.target_accession = df_mhc_web.target_accession.apply( mo.a2b )
    # for a couple of columns, replace repeating values with 0-based integer index of a list
    df_mhc_allele = __Encode_to_Integer__( df_mhc_web.mhc_allele, spread_based_on_weight = False )
    df_mhc_allele.to_csv( f'{dir_folder_pipeline_web}mhc_allele.tsv', columns = [ 'value' ], sep = '\t', index = False )
    mo = MAP.Map( dict( ( value, index ) for index, value in enumerate( df_mhc_allele.value.values ) ) )
    df_mhc_web.mhc_allele = df_mhc_web.mhc_allele.apply( mo.a2b )

    # convert the TCellCrossReactivity data values into more compact form
    df_mhc_web[ 'source_is_blastp' ] = ( df_mhc_web.source == 'blastp' ).astype( int )
    l_col_positive_integer_datatype = [ 'mhc_class', 'mhc_allele', 'source_is_blastp', 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'window_size' ]
    for col in l_col_positive_integer_datatype : df_mhc_web[ col ] = df_mhc_web[ col ].fillna( -1 ).astype( int ) # nan values in columns containing positive integers will be represented by -1 value
    l_col_integer_datatype = [ 'score_blosum' ]
    for col in l_col_integer_datatype : df_mhc_web[ col ] = df_mhc_web[ col ].fillna( 0 ).astype( int ) # nan values will be represented by 0 value
    l_col_float_scientific_notation = [ 'e_value' ] # convert very small float values to scientific notation to reduce the size of the file while retainint the accuracy
    for col in l_col_float_scientific_notation : df_mhc_web[ col ] = list( '' if np.isnan( value ) else "{:.2e}".format( value ) for value in df_mhc_web[ col ].values )
    l_col_optional_data = [ 'source', 'average_score_blosum', 'n_identical_residues', 'proportion_identical', 'n_gaps_in_alignment', 'query_subsequence_without_gap', 'target_subsequence_without_gap', 'length_query_subsequence_without_gap', 'length_target_subsequence_without_gap' ]
    df_mhc_web.drop( columns = list( set( df_mhc_web.columns.values ).intersection( l_col_optional_data ) ), inplace = True ) # drop unnecessary columns
    # write compact dataframes of the TCellCrossReactivity data for web application
    df_mhc_web.to_csv( f'{dir_folder_pipeline_web}TCellCrossReactivity.source.tsv.gz', sep = '\t', index = False ) # save an intermediate result
    # reduce redundancy of the alignments in two subsequent levels
    for df, name in zip( [ df_mhc_web ], [ 'all' ] ) :
        df.sort_values( 'score_for_sorting', ascending = False, inplace = True )
        for n_records in [ 1000, 2000, 5000, 10000, 20000, 40000, 80000 ] : # write the TCellCrossReactivity data for multiple numbers of records
            if len( df ) < n_records * float_thres_margin :
                n_records = len( df )
                df.drop( columns = [ 'score_for_sorting' ] ).iloc[ : n_records ].to_csv( f'{dir_folder_pipeline_web}TCellCrossReactivity.{name}.top_{n_records}.tsv', sep = '\t', index = False, float_format = '%.3f' ) # write only 3 digits below floating point to reduce file size (except for floats in columns in 'l_col_float_scientific_notation')
                break  # stop subsetting after writing all records as a file
            else : df.drop( columns = [ 'score_for_sorting' ] ).iloc[ : n_records ].to_csv( f'{dir_folder_pipeline_web}TCellCrossReactivity.{name}.top_{n_records}.tsv', sep = '\t', index = False, float_format = '%.3f' ) # write only 3 digits below floating point to reduce file size (except for floats in columns in 'l_col_float_scientific_notation')

    ''' add additional information for using the CoordinateServer for displaying structures '''
    # read pdb accessions and parse data
    df_acc_pdb = pd.read_csv( f'{dir_folder_pipeline_web}acc_pdb.source.tsv.gz', sep = '\t' )
    arr_value = df_acc_pdb.value.values # retrieve the encoded values
    df_acc_pdb = df_acc_pdb.reset_index( drop = True ).join( pd.DataFrame( list( e.split( '_' ) for e in df_acc_pdb.value.values ), columns = [ 'id_pdb', 'id_chain', 'id_model' ] ) )
    df_acc_pdb.id_pdb = df_acc_pdb.id_pdb.str.lower( )

    # read entity records
    PKG.Download_Data( "data/pdb/rcsb_pdb.id_model_and_label_entity_id.tsv.gz", dir_remote, name_package ) # download data
    df_pdb_entity = pd.read_csv( f"{dir_folder_cressp}data/pdb/rcsb_pdb.id_model_and_label_entity_id.tsv.gz", sep = '\t', keep_default_na = False )
    df_pdb_entity = PD_Select( df_pdb_entity, id_pdb = df_acc_pdb.id_pdb.values ) # subset the entity records to accelerate the downstream operations
    # build index and value for faster accession
    dict_index = DF_Build_Index_Using_Dictionary( df_pdb_entity, [ 'id_pdb', 'id_chain' ] )
    arr_value = df_pdb_entity[ [ 'id_model', 'id_entity' ] ].values
    l_l_from_resource = [ ]
    for id_pdb, id_chain, id_model in df_acc_pdb[ [ 'id_pdb', 'id_chain', 'id_model' ] ].values :
        # retrieve list of id_pdb and id_chain pairs with given id_pdb and id_chain
        l_l = arr_value[ dict_index[ id_pdb, id_chain ] ]

        # when id_model is -1 (not available for a given id_pdb and id_chain pair)
        if len( l_l ) == 1 and l_l[ 0 ][ 0 ] == -1 :
            l_l_from_resource.append( l_l[ 0 ] )
            continue

        # for each available id_model, search matched id_model and retrieve matched id_entity
        for id_model_from_resource, id_entity in l_l :
            if id_model_from_resource == id_model :
                l_l_from_resource.append( [ id_model_from_resource, id_entity ] )
                continue

        # when matched id_entity is not found, print error message and put a value indicating id_entity is invalid
        print( f"no matching id_entity for {id_pdb} and {id_chain} pair" )
        l_l_from_resource.append( [ np.nan, np.nan ] )
    # add id_entity information to df_acc_pdb
    df_acc_pdb = df_acc_pdb.reset_index( drop = True ).join( pd.DataFrame( l_l_from_resource, columns = [ 'id_model_from_resource', 'id_entity' ] ) )
    df_acc_pdb.to_csv( f'{dir_folder_pipeline_web}acc_pdb.source.tsv.gz', sep = '\t', index = False ) # update the data for the web_application
    df_acc_pdb.id_model = df_acc_pdb.id_model_from_resource.replace( -1, np.nan ) # replace -1 to np.nan, representing an invalid value
    df_acc_pdb.to_csv( f'{dir_folder_pipeline_web}acc_pdb.tsv', columns = [ 'value', 'id_pdb', 'id_chain', 'id_model', 'id_entity' ], sep = '\t', index = False )

    ''' save BLOSUM62 scores for web-application '''
    # read dict_blosum62 from the tsv file
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score
    # save blosum62 score matrix for the usage in the web application 
    l_l_value = list( )
    for key in dict_blosum62 :
        l_l_value.append( [ ''.join( key ), dict_blosum62[ key ] ] )
    df_blosum62 = pd.DataFrame( l_l_value, columns = [ 'pair_of_amino_acids', 'score_blosum62' ] )
    df_blosum62.to_csv( f'{dir_folder_pipeline_web}blosum62.tsv', sep = '\t', index = False )

    """ 
    Compress structural properties and alignments data
    """
    # retrieve RSA, Phi, Psi, SS8, and datatype mask for accession linked to BCellCrossReactivity data # 2020-08-12 22:19:43 
    str_col_id = 'id_protein' # a column name representing str_col_id

    # convert strings encoding RSA, Phi, Psi values using two ASCII characters to those using one ASCII character for web application usage (to reduce complexity of code and reduce file size) 
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )
    dict_kw_rsa_for_web_application = dict( ascii_min = 33, ascii_max = 126, n_char = 1, value_min = 0, value_max = 1 ) # since string encoding structure information for web application will not be delivered in FASTA format, '>' character can be used in the encoding
    dict_kw_torsion_angle_for_web_application = dict( ascii_min = 33, ascii_max = 126, n_char = 1, value_min = -180, value_max = 180 )

    df_fasta_acc_query = pd.read_csv( f'{dir_folder_pipeline}protein_query.tsv.gz', sep = '\t' ) # load structural property data for target sequences 
    df_fasta_acc_target = pd.read_csv( f'{dir_folder_pipeline}protein_target.tsv.gz', sep = '\t' ) # load structural property data for target sequences 
    
    # subset structural data and order records in the same order as 'df_acc_query' or 'df_acc_target'
    df_fasta_acc_query_subset = df_fasta_acc_query.set_index( str_col_id ).loc[ df_acc_query.value.values ].reset_index( )
    df_fasta_acc_target_subset = df_fasta_acc_target.set_index( str_col_id ).loc[ df_acc_target.value.values ].reset_index( )
    
    # parse structural data of query sequences
    df_fasta_acc_query_subset.set_index( str_col_id, inplace = True )
    df_fasta_acc_query_subset[ 'rsa___ascii_encoding_1_character_from_33_to_126__from_0_to_1__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_query_subset[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ].dropna( ).to_dict( ), ** dict_kw_rsa ), ** dict_kw_rsa_for_web_application ) )
    df_fasta_acc_query_subset[ 'phi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_query_subset[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle ), ** dict_kw_torsion_angle_for_web_application ) )
    df_fasta_acc_query_subset[ 'psi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_query_subset[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle ), ** dict_kw_torsion_angle_for_web_application ) )
    dict_fasta = df_fasta_acc_query_subset[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
    df_fasta_acc_query_subset[ 'structure_id' ] = pd.Series( dict( ( acc, Encode_List_of_Strings( Decode_List_of_Strings( dict_fasta[ acc ] ), chr_representing_repeated_string = None ) ) for acc in dict_fasta ) )
    df_fasta_acc_query_subset.reset_index( inplace = True )
    l_col_acc_query = [ 'rsa___ascii_encoding_1_character_from_33_to_126__from_0_to_1__for_web_application', 'phi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application', 'psi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB', 'structure_id' ]

    # parse structural data of target sequences
    df_fasta_acc_target_subset.set_index( str_col_id, inplace = True )
    df_fasta_acc_target_subset[ 'rsa___ascii_encoding_1_character_from_33_to_126__from_0_to_1__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_target_subset[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ].dropna( ).to_dict( ), ** dict_kw_rsa ), ** dict_kw_rsa_for_web_application ) )
    df_fasta_acc_target_subset[ 'phi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_target_subset[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle ), ** dict_kw_torsion_angle_for_web_application ) )
    df_fasta_acc_target_subset[ 'psi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application' ] = pd.Series( ASCII_Encode( ASCII_Decode( df_fasta_acc_target_subset[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].dropna( ).to_dict( ), ** dict_kw_torsion_angle ), ** dict_kw_torsion_angle_for_web_application ) )
    dict_fasta = df_fasta_acc_target_subset[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
    df_fasta_acc_target_subset[ 'structure_id' ] = pd.Series( dict( ( acc, Encode_List_of_Strings( Decode_List_of_Strings( dict_fasta[ acc ] ), chr_representing_repeated_string = None ) ) for acc in dict_fasta ) )
    df_fasta_acc_target_subset.reset_index( inplace = True )
    l_col_acc_target = [ 'rsa___ascii_encoding_1_character_from_33_to_126__from_0_to_1__for_web_application', 'phi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application', 'psi___ascii_encoding_1_character_from_33_to_126__from_-180_to_180__for_web_application', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB', 'structure_id' ]

    # 2020-08-09 03:18:37 
    # find intervals of alignments and protein sequences linked to BCellCrossReactivity data
    dict_dict_it = dict( ) 
    # retrieve intervals from BCellCrossReactivity data
    for filename, l_col in zip( [ 'alignment_query_target', 'alignment_query_pdb', 'alignment_target_pdb', 'structural_property_query', 'structural_property_target' ], [ [ 'id_alignment', 'alignment_start', 'alignment_end' ], [ 'id_alignment_structure_query', 'alignment_structure_query_start', 'alignment_structure_query_end' ], [ 'id_alignment_structure_target', 'alignment_structure_target_start', 'alignment_structure_target_end' ], [ 'query_accession', 'query_start', 'query_end' ], [ 'target_accession', 'target_start', 'target_end' ] ] ) :
        dict_it = dict( )
        for index, start, end in df_subsequence_pdb_web[ l_col ].values :
            if index == -1 : continue # ignore invalid 'index'
            if index not in dict_it : dict_it[ index ] = intervaltree.IntervalTree( )
            dict_it[ index ].addi( start - 1, end ) # 1-based to 0-based system
        for index in dict_it : dict_it[ index ].merge_overlaps( )
        dict_dict_it[ filename ] = dict_it
    # write compact alignment and structural property records by discarding the portion of data values that were not linked to the records in the BCellCrossReactivity data
    for df, filename, l_col in zip( [ df_matched_subset, df_blastp_pdb_query_subset, df_blastp_pdb_target_subset, df_fasta_acc_query_subset, df_fasta_acc_target_subset ], [ 'alignment_query_target', 'alignment_query_pdb', 'alignment_target_pdb', 'structural_property_query', 'structural_property_target' ], [ [ 'query_alignment', 'target_alignment' ], [ 'query_alignment', 'target_alignment' ], [ 'query_alignment', 'target_alignment' ], l_col_acc_query, l_col_acc_target ] ) : # dict_dict_it :
        arr_data = df[ l_col ].values # retrieve data values from the dataframe
        # use pre-initialized list of list to compose a dataframe since iterating a dictionary does not guarantee the order of keys
        l_l_value = list( list( ) for index in range( len( arr_data ) ) ) # initialize a list that will store compact data values after discarding regions that are not linked to BCellCrossReactivity data
        dict_it = dict_dict_it[ filename ]
        for index in dict_it :
            it = dict_it[ index ]
            arr = np.array( list( list( interval )[ : 2 ] for interval in it ) ) # ignore the data field in the interval object
            arr = arr[ arr[ :, 0 ].argsort( ) ] # sort intervals with the start positions of intervals
            if len( arr ) == 1 :
                str_discarded_regions = "0:{}".format( arr[ 0 ][ 0 ] )
            else :
                n_intervals = len( arr )
                arr_raveled = arr.ravel( )
                arr_interval_end = np.zeros( n_intervals, dtype = int )
                arr_interval_end[ 1 : ] = arr_raveled[ 1 : - 1 : 2 ]
                arr_interval_start = arr_raveled[ : : 2 ]
                str_discarded_regions = ';'.join( list( "{}:{}".format( int_discarded_start, len_discarded ) for int_discarded_start, len_discarded in zip( arr_interval_end, arr_interval_start - arr_interval_end ) ) )
            l_l_value[ index ].append( str_discarded_regions )
            for col, value in zip( l_col, arr_data[ index ] ) :
                str_compact = ''
                if isinstance( value, float ) and np.isnan( value ) : pass
                elif col == 'structure_id' :
                    l_structure_id = Decode_List_of_Strings( value, chr_separator = ';', chr_representing_repeated_string = None )
                    l_structure_id_compact = list( )
                    for start, end in arr : l_structure_id_compact.extend( l_structure_id[ start : end ] )
                    str_compact = '' if len( list( value for value in l_structure_id_compact if isinstance( value, float ) ) ) == len( l_structure_id_compact ) else Encode_List_of_Strings( l_structure_id_compact, chr_separator = ';', chr_representing_repeated_string = '=' ) # if 'l_structure_id_compact' only contains np.nan, set 'str_compact' to an empty string
                else :
                    for start, end in arr : str_compact += value[ start : end ]
                l_l_value[ index ].append( str_compact )
        l_empty_line = list( np.nan for index in range( 1 + len( l_col ) ) )
        l_l_value = list( l_empty_line if len( l_value ) == 0 else l_value for l_value in l_l_value )
        df_compact = pd.DataFrame( l_l_value, columns = [ 'discarded_regions' ] + list( col + '__compact__for_web_application' for col in l_col ) )
        df_compact.to_csv( f'{dir_folder_pipeline_web}{filename}__compact__for_web_application.tsv'.format( filename = filename ), sep = '\t', index = False )
        
        
    """ 
    Compress output files with gzip and base64 encoding
    """
    # compress tsv files into base64-encoded gzipped files to reduce the file size 
    # make directories for gzipped files
    dir_folder_pipeline_web_tmp = f"{dir_folder_pipeline_web}tmp/"
    dir_folder_pipeline_web_base64 = f"{dir_folder_pipeline_web}base64/"

    os.makedirs( dir_folder_pipeline_web_tmp, exist_ok = True )
    os.makedirs( dir_folder_pipeline_web_base64, exist_ok = True )
    for dir_file in glob.glob( f"{dir_folder_pipeline_web}*.tsv" ) : # export tsv files in the 'dir_folder_pipeline_web' folder
        print( f"compressing {dir_file} and export data for web application" )
        l = dir_file.rsplit( '/', 1 )
        l.insert( 1, 'tmp' )
        dir_file_binary = '/'.join( l ) + '.gz'
        with open( dir_file_binary, 'wb' ) as newfile : # create gzipped file using command line
            result = subprocess.run( [ 'gzip', '-c', dir_file ], stdin = PIPE, stdout = PIPE, stderr = PIPE )
            newfile.write( result.stdout )
        l[ 1 ] = 'base64'
        dir_file_binary_base64 = '/'.join( l ) + '.gz.base64.txt'
        Base64_Encode( dir_file_binary, dir_file_binary_base64 ) # convert binary file into text using base64 encoding

    """ 
    Filter out files with large number of records 
    """
    # identify and delete files with large number of records
    int_thres_large_files = 41000 # number of records for filtering out files with too many number of records (so that it is impractical to download, parse, and plot the data)
    df = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_pipeline_web_base64}*.top_*.tsv.gz.base64.txt", retrieve_file_size = True )
    df.wildcard_1 = df.wildcard_1.astype( int )
    df_large_file = PD_Threshold( df, wildcard_1a = int_thres_large_files )
    print( f"{len( df_large_file )}/{len( df )} files ({round( float( df_large_file.size_in_bytes.sum( ) / 1e6 ), 2 )}MB/{round( float( df.size_in_bytes.sum( ) / 1e6 ), 2 )}MB) with more than {int_thres_large_files} records will be filtered out for web deployment" )
    for dir_file in df_large_file.dir.values :
        os.remove( dir_file )

    """ 
    Record number of record for each window size and data type
    """
    # compose a dictionary showing the number of maximum records of BCellCrossReactivityData # 2021-01-07 07:46:10 
    df_file = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_pipeline_web_base64}*.*.windowSize_*.top_*.tsv.gz.base64.txt", retrieve_file_size = True )
    for col in [ 'wildcard_2', 'wildcard_3' ] :
        df_file[ col ] = df_file[ col ].astype( int )
    df_file_largest_n_records = df_file.sort_values( [ 'wildcard_3' ], ascending = False ).drop_duplicates( subset = [ "wildcard_0", "wildcard_1", "wildcard_2" ], keep = 'first' ) # retrieve files with largest number of records

    dict_data = dict( )
    for wc_1 in df_file_largest_n_records[ 'wildcard_1' ].unique( ) :
        df_wc_1 = PD_Select( df_file_largest_n_records, wildcard_1 = wc_1 )
        df_wc_1 = df_wc_1.sort_values( 'wildcard_2' )
        df_wc_1.wildcard_2 = df_wc_1.wildcard_2.astype( str ).astype( object )
        dict_data[ wc_1 ] = df_wc_1.set_index( 'wildcard_2' ).wildcard_3.to_dict( )
    # record number of record for each window size and data type (B-cell data)
    dict_cressp_setting[ 'dict_bcell_crossreactivity_data_to_n_records' ] = dict_data 


    # compose a dictionary showing the number of maximum records of TCellCrossReactivityData # 2021-01-07 07:46:10 
    df_file = GLOB_Retrive_Strings_in_Wildcards( f"{dir_folder_pipeline_web_base64}*.*.top_*.tsv.gz.base64.txt", retrieve_file_size = True )
    df_file = PD_Select( df_file, wildcard_0 = 'TCellCrossReactivity' )
    df_file.wildcard_2 = df_file.wildcard_2.astype( int )
    df_file_largest_n_records = df_file.sort_values( [ 'wildcard_2' ], ascending = False ).drop_duplicates( subset = [ "wildcard_0", "wildcard_1" ], keep = 'first' ) # retrieve files with largest number of records

    dict_data = df_file_largest_n_records.set_index( 'wildcard_1' ).wildcard_2.to_dict( )
    # record number of record for each window size and data type (T-cell data)
    dict_cressp_setting[ 'dict_tcell_crossreactivity_data_to_n_records' ] = dict_data        

    """ 
    export data, setting, and metadata for visualization on a web application 
    """
    name_file_cressp_web_viewer_setting = 'cressp_web_viewer_setting.json'
    # collect a list of file names exported for visualization on a web application
    l_name_file_web = [ name_file_cressp_web_viewer_setting ] # put the name of json cressp setting file before writing the setting as an json file
    # collect file names of base64-encoded files
    for dir_file in glob.glob( f"{dir_folder_pipeline_web_base64}*.base64.txt" ) :
        name_file = dir_file.rsplit( '/', 1 )[ 1 ]
        l_name_file_web.append( name_file )
        shutil.copyfile( dir_file, f"{dir_folder_web}{name_file}" )
    dict_cressp_setting[ 'l_name_file_web' ] = l_name_file_web

    # record completed time
    dict_cressp_setting[ 'str_time_completed' ] = datetime.datetime.now( ).strftime( "%Y.%m.%d (%H:%M)" )
    
    with open( f"{dir_folder_web}{name_file_cressp_web_viewer_setting}", 'w' ) as newfile :
        json.dump( dict_cressp_setting, newfile, indent = 6 )

    # copy web viewer to 
    shutil.copyfile( f"{dir_folder_cressp}web_application/cressp_web_viewer.html", f"{dir_folder_output}CRESSP_Web_Viewer.html" )
    # write read-me file
    with open( f'{dir_folder_output}readme.txt', 'w' ) as newfile :
        newfile.write( '\n'.join( [ 
            'CRESSP_Web_Viewer.html\t(web application)\t(1) please open CRESSP web viewer with a compatible web browser (tested on Google Chrome and MicroSoft Edge)',
            '\t\t(2) drag & drop files in web_application/ folder (>16 files, depending on the input file size) and click "start analysis" button',
            'web_application/\t(folder)\ta folder containing input files for CRESSP web viewer. gzipped base64 encoded text files were used to reduce the file size. unzipped original files can be found at pipeline/web_application/ folder',
            'pipeline/\t(folder)\ta folder containing intermediate data files, including alignment results and structural property estimation and prediction results.' ] ) + '\n' )