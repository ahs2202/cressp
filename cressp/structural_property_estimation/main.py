from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

import tensorflow as tf
from tensorflow.keras import layers

# define read-only global variables during multiprocessing
dict_index_df_blastp, arr_data_df_blastp, dict_acc_to_arr_acc_dssp, dict_acc_to_arr_phi_dssp, dict_acc_to_arr_psi_dssp, dict_acc_to_arr_ss8_dssp, dict_kw_rsa, dict_kw_torsion_angle, dict_kw_ss8, dict_kw_datatype, dict_fasta_protein = [ dict( ) ] * 11

def __Transfer_DSSP_Structural_Property_Through_BLAST__( dir_file_input, dir_folder_cressp, dir_folder_pipeline_temp ) : # 2020-07-29 01:09:01 
    """
    # 2021-05-31 15:15:54 
    transfer structural properties according to the BLAST alignment results
    
    'dir_file_input' : the directory to the output of 'BLAST_Parse_BTOP_String' function
    """
    ''' retrieve UUID of the current process '''
    str_uuid = UUID( )
    
    """ input BLASTP results """
    l_col_blastp = [ 'qaccver', 'saccver', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'btop', 'query_seq_aligned', 'subject_seq_aligned' ]
    global dict_index_df_blastp, arr_data_df_blastp, dict_acc_to_arr_acc_dssp, dict_acc_to_arr_phi_dssp, dict_acc_to_arr_psi_dssp, dict_acc_to_arr_ss8_dssp, dict_kw_rsa, dict_kw_torsion_angle, dict_kw_ss8, dict_kw_datatype, dict_fasta_protein # use global variable
    l_qaccver = pd.read_csv( dir_file_input, sep = '\t', keep_default_na = False, header = None ).iloc[ :, 0 ].values.ravel( ) # retrieve list of qaccver for the current process
    
    """ read BLOSUM62 score matrix """
    # read dict_blosum62 from the tsv file
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score
    
    """ settings """
    int_large_window_size = 25
    int_middle_window_size = 9
    int_small_window_size = 3
    float_thres_weight = 10000

    dict_saccver_to_arr_acc = dict_acc_to_arr_acc_dssp
    dict_saccver_to_arr_phi = dict_acc_to_arr_phi_dssp
    dict_saccver_to_arr_psi = dict_acc_to_arr_psi_dssp
    dict_saccver_to_arr_ss8 = dict_acc_to_arr_ss8_dssp
    dict_qaccver_to_seq = dict_fasta_protein
    
    """ main part """
    int_large_window_radius, int_middle_window_radius, int_small_window_radius = int( ( int_large_window_size - 1 ) / 2 ), int( ( int_middle_window_size - 1 ) / 2 ), int( ( int_small_window_size - 1 ) / 2 ) # retrieve radii of the different window sizes
    newfile_df_acc = gzip.open( f"{dir_folder_pipeline_temp}{str_uuid}.transferred.tsv.gz", 'wb' ) # open output file
    
    ''' process each qaccver '''
    for qaccver in l_qaccver :
        
        ''' initialize '''
        if qaccver not in dict_qaccver_to_seq : continue # if subject sequence of the current alignment is not included in the given 'dict_qaccver_to_seq', ignore the record.
        seq_query = dict_qaccver_to_seq[ qaccver ] # retrieve sequence of the query protein
        len_seq_query = len( dict_qaccver_to_seq[ qaccver ] ) # retrieve total length of the subject sequence
        l_l_acc = [ ]
        l_l_phi = [ ]
        l_l_psi = [ ]
        l_saccver = [ ]
        l_l_weight = [ ]
        # initialize dictionary that will contain structural property values of query sequences aligned to subject sequence # initialize dictionary of lists that will store query_id for each transfer so that sugject_id that contributed the most can be retrieved for each position later 
        arr_ss8_weight = np.zeros( ( 8, len_seq_query ) ) # initialize dictionary of arrays that will store weights for each of 8-state secondary structure classification results
        
        ''' iterate alignments '''
        for saccver, qaccver, pident, sstart, send, qstart, qend, subject_seq_aligned, query_seq_aligned in arr_data_df_blastp[ dict_index_df_blastp[ qaccver ] ][ :, [ 1, 0, 2, 8, 9, 6, 7, 14, 13 ] ] : # for each subject sequence, collect structural property of query sequences aligned to the subject sequence
            subject_seq_aligned = subject_seq_aligned.replace( 'U', 'C' ).replace( 'O', 'K' ) # consider selenocysteine (Sec; U) as cysteine # consider pyrrolysine (Pyl; O) as lysine (even though lysine and pyrrolysine are quite different, it seems lysine is the most similar amino acid to pyrrolysine in structure.) # consider all ambiguous amino acids ('BZJ') as 'X'
            query_seq_aligned = query_seq_aligned.replace( 'U', 'C' ).replace( 'O', 'K' )
            len_alignment = len( subject_seq_aligned )
            arr_acc_subject, arr_phi_subject, arr_psi_subject, arr_ss8_subject = dict_saccver_to_arr_acc[ saccver ], dict_saccver_to_arr_phi[ saccver ], dict_saccver_to_arr_psi[ saccver ], dict_saccver_to_arr_ss8[ saccver ] # retrieve structural property values of the query sequence
            arr_acc_query, arr_phi_query, arr_psi_query, arr_weight_query = np.full( ( 3 + 1, len_seq_query ), np.nan, dtype = np.float32 )
            arr_ss8_weight_query = np.zeros( ( 8, len_seq_query ) )

            arr_score = np.full( len_alignment, np.nan, dtype = float ) # initialize an array for saving a BLOSUM62 score for each position in the alignment
            for index_pos_in_alignment, residue_subject, residue_query in zip( np.arange( len_alignment ), subject_seq_aligned, query_seq_aligned ) : arr_score[ index_pos_in_alignment ] = dict_blosum62[ residue_subject, residue_query ] # retrieve a BLOSUM62 score for each position in the alignment
            int_pos_in_subject = sstart - 1 # 0-based coordinate
            int_pos_in_query = qstart - 1 # 0-based coordinate
            for index_pos_in_alignment, residue_subject, residue_query in zip( np.arange( len_alignment ), subject_seq_aligned, query_seq_aligned ) : # for each position in the alignment, calculate weight and transfer structural property of the query sequence to the position on the subject sequence.
                if residue_query == '-' : int_pos_in_subject += 1
                elif residue_subject == '-' : int_pos_in_query += 1
                else :
                    slice_large_window = slice( max( 0, index_pos_in_alignment - int_large_window_radius ), min( len_alignment, index_pos_in_alignment + 1 + int_large_window_radius ) )
                    slice_middle_window = slice( max( 0, index_pos_in_alignment - int_middle_window_radius ), min( len_alignment, index_pos_in_alignment + 1 + int_middle_window_radius ) )
                    slice_small_window = slice( max( 0, index_pos_in_alignment - int_small_window_radius ), min( len_alignment, index_pos_in_alignment + 1 + int_small_window_radius ) )
                    weight = max( 0, arr_score[ slice_large_window ].sum( ) ) * max( 0, arr_score[ slice_middle_window ].sum( ) ) * max( 0, arr_score[ slice_small_window ].sum( ) ) * ( ( pident / 100 ) ** 3 ) * ( qend - qstart + 1 ) # calculate weight of the position in the alignment based on local alignment scores, pident **3, and the alignment length of subject sequence

                    arr_acc_query[ int_pos_in_query ] = arr_acc_subject[ int_pos_in_subject ] # transfer structural property values of a position on the query to the corresponding position on the subject
                    arr_phi_query[ int_pos_in_query ] = arr_phi_subject[ int_pos_in_subject ] 
                    arr_psi_query[ int_pos_in_query ] = arr_psi_subject[ int_pos_in_subject ] 
                    ss8_residue = arr_ss8_subject[ int_pos_in_subject ]
                    if not np.isnan( ss8_residue ) : arr_ss8_weight_query[ int( ss8_residue ), int_pos_in_query ] = weight # if retrieved ss8 is valid, add ss8 of the residue to the array multiplied with the weight
                    arr_weight_query[ int_pos_in_query ] = weight # record a calculated weight for each position on the subject
                    int_pos_in_query += 1
                    int_pos_in_subject += 1

            l_l_acc.append( arr_acc_query )
            l_l_phi.append( arr_phi_query )
            l_l_psi.append( arr_psi_query )
            arr_ss8_weight += arr_ss8_weight_query # accumulate weight for each ss8 classification
            l_saccver.append( saccver )
            l_l_weight.append( arr_weight_query )
        
        ''' combine information from alignments '''
        arr_acc = np.vstack( l_l_acc ).astype( np.float64 )  # build numpy array for efficient calculation # convert datatype to float64 to increase precision
        arr_phi = np.vstack( l_l_phi ).astype( np.float64 ) 
        arr_psi = np.vstack( l_l_psi ).astype( np.float64 ) 
        arr_weight = np.vstack( l_l_weight ).astype( np.float64 ) 
        mask_invalid_weight = arr_weight < float_thres_weight # discard accessibility scores and weights where weights are too low for the transfer of structural properties from protein structures
        # mark residues with invalid weights
        arr_weight[ mask_invalid_weight ] = np.nan
        arr_acc[ mask_invalid_weight ] = np.nan
        arr_phi[ mask_invalid_weight ] = np.nan
        arr_psi[ mask_invalid_weight ] = np.nan
        marr_acc = np.ma.masked_array( data = arr_acc, mask = np.isnan( arr_acc ) ) # mask np.nan values to ignore positions containing invalid values during calculations
        marr_phi = np.ma.masked_array( data = arr_phi, mask = np.isnan( arr_phi ) )
        marr_psi = np.ma.masked_array( data = arr_psi, mask = np.isnan( arr_psi ) )
        marr_weight = np.ma.masked_array( data = arr_weight, mask = np.isnan( arr_weight ) )
        marr_weight_normalized = marr_weight / marr_weight.sum( axis = 0 ) # calculate normalized weights (sum of weights in each column = 1)
        arr_acc_query = ( marr_acc * marr_weight_normalized ).sum( axis = 0 ).filled( np.nan )
        arr_phi_query = ( marr_phi * marr_weight_normalized ).sum( axis = 0 ).filled( np.nan )
        arr_psi_query = ( marr_psi * marr_weight_normalized ).sum( axis = 0 ).filled( np.nan )

        # retrieved accumulated weight for each ss8 classification # zero-filled array
        arr_ss8_weight[ arr_ss8_weight < float_thres_weight ] = 0 # ss8 with accumulated weights less than 'float_thres_weight' will be ignored
        arr_ss8 = arr_ss8_weight.argmax( axis = 0 ).astype( float ) # Retrieving consensus of ss8 classification by weighted voting (ss8 classification with the largest weight will be chosen as consensus ss8 classification)
        arr_ss8[ arr_ss8_weight.sum( axis = 0 ) == 0 ] = np.nan # retrieve a mask for residues without transferred ss8 classification, and put np.nan representing an invalid value on the positions of these residues.

        arr_weight_normalized = marr_weight_normalized.filled( 0 ) # zero-filled array containing normalized weights # argmax does not work with np.nan, and even with masked array, residues without any weights should be identified by different means other than checking the output of np.argmax
        arr_saccver_aligned_to_qaccver = np.array( l_saccver, dtype = object ) # retrieve list of saccver that transferred structural property to qaccver
        arr_saccver = np.array( list( np.nan if bool_invalid_residue else arr_saccver_aligned_to_qaccver[ index_saccver_largest_weight ] for index_saccver_largest_weight, bool_invalid_residue in zip( arr_weight_normalized.argmax( axis = 0 ), arr_weight_normalized.sum( axis = 0 ) == 0 ) ), dtype = object ) # retrieve saccver of largest weight (saccver that was most significant in the transfer of structural property to subject sequence) for each residue position of subject sequence # ignore the sequence if it contains invalid residue

        # encode structural properties using ASCII characters and save it as a tabular data
        # l_col = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'structure_id___redundancy_reduced' ]
        
        __Encode_and_Write_Structural_Properties__( newfile_df_acc, { 'id_protein' : qaccver, 'seq' : seq_query, 'acc' : arr_acc_query, 'phi' : arr_phi_query, 'psi' : arr_psi_query, 'ss8' : arr_ss8, 'structure_id' : arr_saccver } ) # encode transferred and combined structural properties into ASCII strings using ASCII encoding
    
    newfile_df_acc.close( ) # close file
def __Iterate_and_Parse_Structural_Properties__( dir_file_sp, int_datatype = None ) :
    """
    # 2022-01-09 23:19:14 
    iterate through 'dir_file_sp' structural property database for each record (protein sequence)
    
    parse structural properties with typical parameters for parsing ascii-encoded structural properties and typical column names  
    Also, initializes a array containing the identifier of the current dataset of 'acc' datatype
    
    'dir_file_sp' : a TSV file 
    'int_datatype' : integer identifying the current dataset. Used for initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype. if 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' column already exists, 'int_datatype' argument will not be used.
    """
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    ''' define internal functions '''
    def __remove_quotation__( e, str_quotation = '"' ) :
        """ remove quotation from value in the CSV file """
        if len( e ) < 2 :
            return e
        ''' if a given string is wrapped in a quotation and an escape character of quotation character is found, remove quotation and replace the escape character with normal quotation character '''
        if e[ 0 ] == str_quotation and e[ - 1 ] == str_quotation and str_quotation * 2 in e[ 1 : -1 ] :
            return e[ 1 : -1 ].replace( str_quotation * 2, str_quotation )
        return e 

    with gzip.open( dir_file_sp, 'rb' ) as file :
        ''' read the header line '''
        l_name_col = file.readline( ).decode( )[ : -1 ].split( '\t' )
        flag_datatype_acc_does_not_exists = len( l_name_col ) == 7 # retrieve a flag indicating datatype_acc does not exist and should be initialized
        ''' iterate contents header '''
        while True :
            line = file.readline( ).decode( )
            if len( line ) == 0 :
                break

            """ parse line """
            l_data = list( __remove_quotation__( e ) for e in line[ : -1 ].split( '\t' ) )
            ''' handle when 'datatype_acc' column does not exist '''
            if flag_datatype_acc_does_not_exists : # parse with datatype does not exist
                id_protein, seq, rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1, phi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180, psi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180, ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C, structure_id___redundancy_reduced = list( __remove_quotation__( e ) for e in line.strip( ).split( '\t' ) ) # parse line
            else :
                id_protein, seq, rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1, phi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180, psi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180, ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C, rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB, structure_id___redundancy_reduced = list( __remove_quotation__( e ) for e in line.strip( ).split( '\t' ) ) # parse line
                ''' parse arr_datatype '''
                arr_datatype_acc = ASCII_Decode( [ rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB ], ** dict_kw_datatype )[ 0 ]

            ''' parse structural properties '''
            # decode structural properties
            arr_acc = ASCII_Decode( [ rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1 ], ** dict_kw_rsa )[ 0 ]
            arr_phi = ASCII_Decode( [ phi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180 ], ** dict_kw_torsion_angle )[ 0 ]
            arr_psi = ASCII_Decode( [ psi___ascii_encoding_2_characters_from_33_to_126__from_minus180_to_180 ], ** dict_kw_torsion_angle )[ 0 ]
            arr_ss8 = ASCII_Decode( [ ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C ], ** dict_kw_ss8 )[ 0 ]
            arr_structure_id = Decode_List_of_Strings( structure_id___redundancy_reduced ) if len( structure_id___redundancy_reduced ) > 0 else np.nan
            
            ''' initialize arr_datatype '''
            if flag_datatype_acc_does_not_exists : # parse with datatype does not exist
                arr_datatype_acc = np.zeros_like( arr_acc ) # default = 0 (predicted structural property)
                arr_datatype_acc[ ~ np.isnan( arr_acc ) ] = int_datatype
            
            ''' yield parsed record '''
            yield { 'id_protein' : id_protein, 'seq' : seq, 'acc' : arr_acc, 'phi' : arr_phi, 'psi' : arr_psi, 'ss8' : arr_ss8, 'structure_id' : arr_structure_id, 'datatype_acc' : arr_datatype_acc }
def __Encode_and_Write_Structural_Properties__( file_handle, dict_record ) :
    """
    # 2021-05-03 16:18:16 
    Encode structural properties of a given record to the given file_handle (should be gzipped file)
    
    This function corresponds to the writing version of '__Iterate_and_Parse_Structural_Properties__'
    
    
    'file_handle' : a gzipped TSV file to be written
    'dict_record' : structural properties and sequence information of a single record
    """
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    ''' define internal functions '''
    def __add_quotation__( e, str_quotation = '"' ) :
        """ add quotation from value to write a entry of a CSV file """
        ''' if the quotation character is found, change replace it with escape character and wrap the string with quotations '''
        if str_quotation in e :
            return str_quotation + e.replace( str_quotation, str_quotation * 2 ) + str_quotation
        else : 
            return e
    
    # l_col = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB', 'structure_id___redundancy_reduced' ] # when datatype_acc is present
    # l_col = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'structure_id___redundancy_reduced' ] # when datatype_acc is not present

    ''' encode structural properties '''
    l_data = [ dict_record[ 'id_protein' ], dict_record[ 'seq' ], ASCII_Encode( [ dict_record[ 'acc' ] ], ** dict_kw_rsa )[ 0 ], ASCII_Encode( [ dict_record[ 'phi' ] ], ** dict_kw_torsion_angle )[ 0 ], ASCII_Encode( [ dict_record[ 'psi' ] ], ** dict_kw_torsion_angle )[ 0 ], ASCII_Encode( [ dict_record[ 'ss8' ] ], ** dict_kw_ss8 )[ 0 ] ] + ( [ ASCII_Encode( [ dict_record[ 'datatype_acc' ] ], ** dict_kw_datatype )[ 0 ] ] if 'datatype_acc' in dict_record else [ ] ) + [ Encode_List_of_Strings( dict_record[ 'structure_id' ] if not isinstance( dict_record[ 'structure_id' ], float ) else [ '' ] * len( dict_record[ 'seq' ] ) ) ] # if 'datatype_acc' column does not exist, does not include it in the output
    file_handle.write( ( '\t'.join( list( map( __add_quotation__, l_data ) ) ) + '\n' ).encode( ) ) # write a record
def __Parse_Structural_Properties__( dir_file_sp, int_datatype ) :
    """
    # 2021-05-03 16:18:16 
    parse structural properties with typical parameters for parsing ascii-encoded structural properties and typical column names  
    Also, initializes a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
    'int_datatype' : integer identifying the current dataset. Used for initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
    """
    # read structural property database
    df_sp = pd.read_csv( dir_file_sp, sep = '\t' )
    df_sp.set_index( 'id_protein', inplace = True )

    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    dict_sp = dict( )
    dict_sp[ 'seq' ] = df_sp.seq.to_dict( )
    dict_sp[ 'acc' ] = ASCII_Decode( df_sp.rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1.to_dict( ), ** dict_kw_rsa ) # decode mkdssp outputs of RCSB PDB data
    dict_sp[ 'phi' ] = ASCII_Decode( df_sp[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
    dict_sp[ 'psi' ] = ASCII_Decode( df_sp[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
    dict_sp[ 'ss8' ] = ASCII_Decode( df_sp[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].to_dict( ), ** dict_kw_ss8 )
    dict_fasta = df_sp[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
    dict_sp[ 'structure_id' ] = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )

    # initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
    dict_datatype_acc = dict( )
    for h in dict_sp[ 'acc' ] :
        arr = np.zeros_like( dict_sp[ 'acc' ][ h ] ) # default = 0 (predicted structural property)
        arr[ ~ np.isnan( dict_sp[ 'acc' ][ h ] ) ] = int_datatype
        dict_datatype_acc[ h ] = arr
    dict_sp[ 'datatype_acc' ] = dict_datatype_acc
    return dict_sp
def __Iterate_and_Combine_Structural_Properties__( dir_file_sp_1, int_datatype_1, dir_file_sp_2, int_datatype_2, dir_file_sp_combined, flag_debugging = False ) :
    """
    # 2021-05-29 16:12:09 
    Combine structural properties from the two different sources. 
    Structural properties from the source with the lower priority will be overwritten with that of the higher priority
    the smaller database will be loaded on the memory and the larger database will be iterated without loading the database entirely on the memory.
    Also, the output file will be written line by line without saving results on the memory.

    'dir_file_sp_1' : the file directory of the structural property database 1 (larger)
    'int_datatype_1' : priority of 'dir_file_sp_1' (higher number means higher priority)
    'dir_file_sp_2' : the file directory of the structural property database 2 (smaller)
    'int_datatype_2' : priority of 'dir_file_sp_2' (higher number means higher priority)
    'dir_file_sp_combined' : the file directory of the structural property database containing combined information of structural property database 1 and 2
    
    """
    ''' open an output file '''
    newfile = gzip.open( dir_file_sp_combined, 'wb' )
    l_col = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB', 'structure_id___redundancy_reduced' ]
    ''' write a header '''
    newfile.write( ( '\t'.join( l_col ) + '\n' ).encode( ) )
    
    
    ''' load a smaller structural property database on the memory '''
    dict_sp = __Parse_Structural_Properties__( dir_file_sp_2, int_datatype = int_datatype_2 )
    for dict_record_1 in __Iterate_and_Parse_Structural_Properties__( dir_file_sp_1, int_datatype_1 ) : # iterate each record
        id_protein = dict_record_1[ 'id_protein' ] # retrieve accession
        seq = dict_record_1[ 'seq' ] 
        len_seq = len( seq ) # retrieve the length of the sequence
        
        ''' if the current protein does not exist in the smaller structural property database, write the record to the combined database as-is '''
        if id_protein not in dict_sp[ 'acc' ] :
            if flag_debugging :
                print( f'{id_protein} record from {dir_file_sp_1} was written' )
            __Encode_and_Write_Structural_Properties__( newfile, dict_record_1 )
        else :
            if flag_debugging :
                print( f'{id_protein} records from {dir_file_sp_1} and {dir_file_sp_2} were combined and written' )
            """ if the current protein exists in the smaller structural property database, pop the record out of the smaller structural property database and combine structural properties, and write the record to the combined database """
            ''' pop the record out of the smaller structural property database and combine structural properties, and write the record to the combined database '''
            dict_record_2 = { 'id_protein' : id_protein, 'seq' : dict_sp[ 'seq' ].pop( id_protein ), 'acc' : dict_sp[ 'acc' ].pop( id_protein ), 'phi' : dict_sp[ 'phi' ].pop( id_protein ), 'psi' : dict_sp[ 'psi' ].pop( id_protein ), 'ss8' : dict_sp[ 'ss8' ].pop( id_protein ), 'structure_id' : dict_sp[ 'structure_id' ].pop( id_protein ), 'datatype_acc' : dict_sp[ 'datatype_acc' ].pop( id_protein ) }
            
            # initialize combined dict_record
            dict_record_combined = { 'id_protein' : id_protein, 'seq' : seq }
            ''' for each datatype '''
            for str_datatype in [ 'acc', 'phi', 'psi', 'ss8', 'structure_id', 'datatype_acc' ] :
                # initialize array containing combined data
                arr_combined = np.full( len_seq, np.nan, dtype = object ) if str_datatype == 'structure_id' else np.full( len_seq, np.nan ) 
            
                # update data in the order of increasing priorities (so that data with the highest priority is updated last and thus not overwritten by any other data)                
                l_dict_record_for_structure_id = [ dict_record_2 ] if int_datatype_1 < int_datatype_2 else [ dict_record_1 ]
                l_dict_record_for_others = [ dict_record_1, dict_record_2 ] if int_datatype_1 < int_datatype_2 else [ dict_record_2, dict_record_1 ]
                for dict_record in l_dict_record_for_structure_id if str_datatype == 'structure_id' else l_dict_record_for_others : # for 'structure_id', only use structure_id of the database of the higher priority ('RCSB_PDB' to be exact) 
                    arr_current = dict_record[ str_datatype ]
                    mask_valid = ~ pd.isnull( arr_current ) if str_datatype == 'structure_id' else ~ np.isnan( arr_current )
                    arr_combined[ mask_valid ] = arr_current[ mask_valid ] # update data
                dict_record_combined[ str_datatype ] = arr_combined
            ''' write combined record '''
            __Encode_and_Write_Structural_Properties__( newfile, dict_record_combined )
            
    """ write the structural property records from the smaller structural property database that were absent in the larger structural property database and thus not merged with records from the larger database """
    for id_protein in list( dict_sp[ 'acc' ] ) :
        if flag_debugging :
            print( f'{id_protein} record from {dir_file_sp_2} was written' )
        ''' pop the record out of the smaller structural property database and write them to the combined database '''
        dict_sp[ 'structure_id' ].pop( id_protein )
        __Encode_and_Write_Structural_Properties__( newfile, { 'id_protein' : id_protein, 'seq' : dict_sp[ 'seq' ].pop( id_protein ), 'acc' : dict_sp[ 'acc' ].pop( id_protein ), 'phi' : dict_sp[ 'phi' ].pop( id_protein ), 'psi' : dict_sp[ 'psi' ].pop( id_protein ), 'ss8' : dict_sp[ 'ss8' ].pop( id_protein ), 'structure_id' : np.nan, 'datatype_acc' : dict_sp[ 'datatype_acc' ].pop( id_protein ) } )
    newfile.close( )
def __Combine_Structural_Properties__( dict_sp_1, dict_sp_2 ) :
    """
    # 2021-05-29 16:12:09 
    Combine structural properties from the two different sources. 
    Structural properties from the source with the lower priority will be overwritten with that of the higher priority

    'df_sp_1' : structural properties of the lower priority
    'df_sp_2' : structural properties with the higher priority
    """
    dict_sp_combined = dict( )

    for str_datatype in [ 'acc', 'phi', 'psi', 'ss8', 'structure_id', 'datatype_acc' ] :
        dict_arr_1 = dict_sp_1[ str_datatype ]
        dict_arr_2 = dict_sp_2[ str_datatype ]
        dict_arr_combined = dict( )
        for h in set( dict_arr_1 ).union( set( dict_arr_2 ) ) :
            # retrieve length of the protein sequence
            len_seq = 0
            for dict_arr in [ dict_arr_1, dict_arr_2 ] : 
                if h in dict_arr :
                    len_seq = len( dict_arr[ h ] )

            arr_combined = np.full( len_seq, np.nan, dtype = object ) if str_datatype == 'structure_id' else np.full( len_seq, np.nan ) # initialize array containing combined data
            # update data in the order of increasing priorities (so that data with the highest priority is updated last and thus not overwritten by any other data)
            for dict_arr in [ dict_arr_2 ] if str_datatype == 'structure_id' else [ dict_arr_1, dict_arr_2 ] : # for 'structure_id', only use structure_id of the database of the higher priority ('RCSB_PDB' to be exact) 
                if h in dict_arr :
                    arr_current = dict_arr[ h ]
                    mask_valid = ~ pd.isnull( arr_current ) if str_datatype == 'structure_id' else ~ np.isnan( arr_current )
                    arr_combined[ mask_valid ] = arr_current[ mask_valid ] # update data
            dict_arr_combined[ h ] = arr_combined
        dict_sp_combined[ str_datatype ] = dict_arr_combined
    return dict_sp_combined
def __Encode_Structural_Properties__( dict_sp, dict_fasta_protein ) :
    """
    # 2021-05-29 23:56:59 
    Compose a dataframe containing encoded structural properties using given 'dict_sp' and 'dict_fasta_protein'

    'dict_sp' : dictionary containing structural properties
    'dict_fasta_protein' : dictionary containing protein sequences

    """
    dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
    dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
    dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
    dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

    ''' initialize dataframe with protein sequences '''
    df_sp = pd.Series( dict_fasta_protein, name = 'seq' ).reset_index( ).rename( columns = { 'index' : 'id_protein' } ).set_index( 'id_protein' )

    # encode combined structural properties into ASCII strings using ASCII encoding
    df_sp[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ] = pd.Series( ASCII_Encode( dict_sp[ 'acc' ], ** dict_kw_rsa ) )
    df_sp[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'phi' ], ** dict_kw_torsion_angle ) )
    df_sp[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'psi' ], ** dict_kw_torsion_angle ) )
    df_sp[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ] = pd.Series( ASCII_Encode( dict_sp[ 'ss8' ], ** dict_kw_ss8 ) )
    df_sp[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ] = pd.Series( ASCII_Encode( dict_sp[ 'datatype_acc' ], ** dict_kw_datatype ) )
    df_sp[ 'structure_id___redundancy_reduced' ] = pd.Series( dict( ( h, Encode_List_of_Strings( dict_sp[ 'structure_id' ][ h ] ) ) for h in dict_sp[ 'structure_id' ] ) )
    df_sp.reset_index( drop = False, inplace = True )
    return df_sp
def __Predict_Structural_Properties_of_Remaining_Residues__( dir_file_protein, dir_folder_pipeline = None, dir_folder_pipeline_temp = '/tmp/', int_number_of_proteins_in_a_batch_during_dnn_prediction = 100, flag_use_all_gpu_devices = False ) :
    """
    Parse arguments
    """
    # directories of file and folders
    dir_file_protein = os.path.abspath( dir_file_protein )
    dir_folder_pipeline = os.path.abspath( dir_folder_pipeline )
    dir_folder_pipeline_temp = os.path.abspath( dir_folder_pipeline_temp )
    if dir_folder_pipeline[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline += '/'
    if dir_folder_pipeline_temp[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline_temp += '/'
    
    # fixed arguments
    float_search_thres_e_value = 30 # set e-value threshold for search (fixed for maximum sensitivity)

    # By default, set 'dir_folder_pipeline' as the folder where the given protein fasta file is located
    if dir_folder_pipeline is None :
        dir_folder_pipeline = dir_file_protein.rsplit( '/', 1 )[ 0 ] + '/'
    name_file = dir_file_protein.rsplit( '/', 1 )[ 1 ].rsplit( '.', 1 )[ 0 ] # retrieve name of file excluding file extension
    dir_file_protein_property = f"{dir_folder_pipeline}{name_file}.tsv.gz" # define output file (a file containing protein sequences and their estimated structural properties)
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # create a working directory of estimating structural properties
    os.makedirs( dir_folder_pipeline_struc, exist_ok = True )
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package
    
    """ read input proteins """
    dict_fasta_protein = FASTA_Read( dir_file_protein )
    dict_header = dict( ( h.split( ' ', 1 )[ 0 ], h ) for h in dict_fasta_protein ) # retrieve a dictionary of header (key = id_protein, value = header)
    dict_fasta_protein = dict( ( h.split( ' ', 1 )[ 0 ], dict_fasta_protein[ h ] ) for h in dict_fasta_protein ) # retrieve sequence_id by spliting the header at the first space (to makes sequence_id consistent with that used with blastp)
    
    """
    Prediction of Structural Properties Using Machine-Learning Models
    """
    # 2021-05-29 23:41:46 
    # use machine-learning model to predict structural properties
    ''' settings '''
    int_length_flanking = 74 # number of protein residues flanking the residue whose structural property values are predicted

    ''' download models '''
    file_name_model_rsa = "RSA.Conv1D.5_large_and_wide_dense.epoch_082.val_loss_0.1185.hdf5" # name of machine learning model used for structural property prediction (RSA)
    file_name_model_ss8 = "SS8.Conv1D.5_large_and_wide_dense.epoch_031.val_loss_0.7378.hdf5" # name of machine learning model used for structural property prediction (SS8)

    PKG.Download_Data( f"structural_property_estimation/{file_name_model_rsa}", dir_remote, name_package ) # download data
    PKG.Download_Data( f"structural_property_estimation/{file_name_model_ss8}", dir_remote, name_package ) # download data

    ''' prepare encoding '''
    def Index_List( l, start = 0 ) :
        ''' # 2021-05-06 22:32:53 
        index a given list containing hashable element with integer index starting with 'start' (default: 0)
        can be used for 'one-hot encoding'
        '''
        return dict( ( e, i ) for e, i in zip( l, np.arange( start, len( l ) + start, 1, dtype = int ) ) )

    # input structure : amino-acids (one-hot encoding) + ss8 (one-hot encoding) + rsa (0~1)
    int_index_start_amino_acid = 0
    l_amino_acids = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
    dict_amino_acid_to_index = Index_List( l_amino_acids, int_index_start_amino_acid ) # one-hot encoding dictionary for 20 amino acids
    int_index_start_ss8 = len( l_amino_acids )
    l_ss8 = [ 'G', 'H', 'I', 'E', 'B', 'T', 'S', 'C' ]
    dict_ss8_to_index = Index_List( l_ss8, int_index_start_ss8 ) # one-hot encoding dictionary for 20 amino acids
    index_rsa = len( l_amino_acids ) + len( l_ss8 ) # index of rsa data in the input array

    ''' predict and update structural properties using ML models for each batch '''
    import tensorflow as tf
    from tensorflow.keras import layers

    ''' initialize '''
    int_index_batch = 0 # record batch index for recording
    l_dict_record_for_a_batch = [ ]
    ''' open an iterator for iterating each record in the database '''
    iter_dict_record = __Iterate_and_Parse_Structural_Properties__( f"{dir_folder_pipeline_struc}{name_file}_transferred_combined.tsv.gz", int_datatype = 2 ) # default datatype is that of RCSB_PDB (int_datatype = 2)
    flag_predicting_sp_of_proteins_solely_based_on_aa = False # this flag will be set to True when 'iter_dict_record' iterator is exhausted
    flag_all_records_were_retrieved = False # this flag will be set to True when all protein records are loaded onto the last batch queue

    while True :
        if not flag_predicting_sp_of_proteins_solely_based_on_aa :
            try :
                ''' iterate each record in the database '''
                dict_record = next( iter_dict_record )
                dict_fasta_protein.pop( dict_record[ 'id_protein' ] ) # remove 'id_protein' from the fasta records to mark the 'id_protein' has been processed
            except StopIteration :
                ''' when all records in the database were retrieved, retrieve records from fasta protein sequences. perform initialization of the process as below '''
                flag_predicting_sp_of_proteins_solely_based_on_aa = True
        if flag_predicting_sp_of_proteins_solely_based_on_aa :
            l_id_protein = list( dict_fasta_protein )
            ''' when 'dict_fasta_protein' is empty, set a flag indicating all records has been retrieved '''
            if len( l_id_protein ) == 0 :
                flag_all_records_were_retrieved = True
            else :
                ''' initialize record of a protein with no aligned structures '''
                id_protein = l_id_protein[ 0 ]
                seq = dict_fasta_protein.pop( id_protein ) # remove 'id_protein' from the fasta records to mark the 'id_protein' has been processed
                len_seq = len( seq ) # retrieve length of protein sequence
                dict_record = { 'id_protein' : id_protein, 'seq' : seq }
                # { 'id_protein' : id_protein, 'seq' : seq, 'acc' : arr_acc, 'phi' : arr_phi, 'psi' : arr_psi, 'ss8' : arr_ss8, 'structure_id' : arr_structure_id, 'datatype_acc' : arr_datatype_acc }              
                for str_datatype in [ 'acc', 'phi', 'psi', 'ss8', 'structure_id', 'datatype_acc' ] :
                    ''' set default array for each datatype '''
                    if str_datatype == 'datatype_acc' :
                        arr = np.full( len_seq, 0, dtype = float )
                    elif str_datatype == 'structure_id' :
                        arr = np.full( len_seq, np.nan, dtype = object )
                    else :
                        arr = np.full( len_seq, np.nan )
                    dict_record[ str_datatype ] = arr
                
        ''' collect records for a batch '''
        if not flag_all_records_were_retrieved :
            l_dict_record_for_a_batch.append( dict_record )
        ''' when collection is complete, operate on the batch '''
        if len( l_dict_record_for_a_batch ) >= int_number_of_proteins_in_a_batch_during_dnn_prediction or flag_all_records_were_retrieved :
            ''' defint an output file directory '''
            int_n_records_prev_processed = int_index_batch * int_number_of_proteins_in_a_batch_during_dnn_prediction # retrieve the number of processed records
            int_n_records_processed = int_n_records_prev_processed + len( l_dict_record_for_a_batch )
            dir_file_output_for_a_batch_being_written = f"{dir_folder_pipeline_temp}{name_file}.batch{int_index_batch}_{int_n_records_prev_processed + 1}-{int_n_records_processed}.being_written.tsv.gz"
            dir_file_output_for_a_batch = f"{dir_folder_pipeline_temp}{name_file}.batch{int_index_batch}_{int_n_records_prev_processed + 1}-{int_n_records_processed}.tsv.gz"
            
            ''' perform encoding for the current batch '''
            if not os.path.exists( dir_file_output_for_a_batch ) : # if an output file does not exist
                # 'l_dict_record_for_a_batch': id_protein -> integer encoding to reduce memory usage
                # perform encoding and record metadata of each position
                x_data = [ ]
                meta_data = [ ]

                for index_protein, dict_record in enumerate( l_dict_record_for_a_batch ) :

                    id_protein = dict_record[ 'id_protein' ]
                    seq = dict_record[ 'seq' ]
                    arr_acc = dict_record[ 'acc' ]
                    arr_ss8 = dict_record[ 'ss8' ]
                    length = len( seq ) # retrieve the length of the protein

                    ''' one-hot encoding (+ RSA) of the current protein '''
                    arr_protein = np.zeros( ( length, len( l_amino_acids ) + len( l_ss8 ) + 1 ), dtype = np.float32 ) # initialize one-hot encoding of the proteinarray
                    for i in range( length ) :
                        # one-hot encoding of amino acid
                        aa = seq[ i ] # retrieve amino acid
                        if aa in dict_amino_acid_to_index : # invalid amino acid will be considered as a missing value (empty row)
                            arr_protein[ i, dict_amino_acid_to_index[ aa ] ] = 1 
                        # one-hot encoding of ss8
                        if not np.isnan( arr_ss8[ i ] ) : # missing value = empty row
                            arr_protein[ i, dict_ss8_to_index[ l_ss8[ int( arr_ss8[ i ] ) ] ] ] = 1
                        # encode rsa information
                        if not np.isnan( arr_acc[ i ] ) : # missing value = empty row
                            arr_protein[ i, index_rsa ] = arr_acc[ i ]

                    ''' prepare one-hot encoding input array for each missing position '''
                    ''' assumes positions of missing values in 'arr_acc' are identical to those of 'arr_ss8' '''
                    for int_pos in np.where( np.isnan( arr_acc ) )[ 0 ] : 
                        start = max( 0, int_pos - int_length_flanking )
                        end = min( length, int_pos + int_length_flanking + 1 )

                        arr = np.zeros( ( int_length_flanking + 1 + int_length_flanking, len( l_amino_acids ) + len( l_ss8 ) + 1 ), dtype = np.float32 ) # initialize input array
                        start_input = max( 0, int_length_flanking - int_pos ) # index of start position of the input array
                        arr[ start_input : start_input + end - start ] = arr_protein[ start : end ]
                        x_data.append( arr )
                        meta_data.append( [ index_protein, int_pos ] )
                x_data = np.array( x_data, dtype = np.float32 )

                ''' predict and update structural properties using ML models for the current batch '''
                for dir_file_model, str_datatype in zip( [ f"{dir_folder_cressp}structural_property_estimation/{file_name_model_rsa}", f"{dir_folder_cressp}structural_property_estimation/{file_name_model_ss8}" ], [ 'acc', 'ss8' ] ) :
                    model = tf.keras.models.load_model( dir_file_model )
                    try :
                        dir_file_during_gpu_usage = f"{dir_folder_cressp}structural_property_estimation/{UUID( )}.flag"
                        with tf.distribute.MirroredStrategy( ).scope( ) if flag_use_all_gpu_devices else open( dir_file_during_gpu_usage, 'w' ) : # distributed calculation if 'flag_use_all_gpus' is True
                            y_pred = model.predict( x_data ) # predict structural properties using a ML model
                            y_pred = y_pred.argmax( axis = 1 ) if str_datatype == 'ss8' else y_pred # multiclass classification for 'ss8' and regression for 'acc'
                            os.remove( dir_file_during_gpu_usage )
                    except tf.errors.InternalError : # when out of memory error occurs during training, use CPU instead
                        os.environ[ 'CUDA_VISIBLE_DEVICES' ] = '-1' # forcing tensorflow to use CPU
                        y_pred = model.predict( x_data ) # predict structural properties using a ML model
                        y_pred = y_pred.argmax( axis = 1 ) if str_datatype == 'ss8' else y_pred # multiclass classification for 'ss8' and regression for 'acc'

                    for i in range( len( y_pred ) ) :
                        index_protein, int_pos = meta_data[ i ] # parse metadata
                        l_dict_record_for_a_batch[ index_protein ][ str_datatype ][ int_pos ] = y_pred[ i ] # update structural property data with predicted data

                ''' write an output file of the current batch '''
                with gzip.open( dir_file_output_for_a_batch_being_written, 'wb' ) as newfile :
                    for dict_record in l_dict_record_for_a_batch :
                        __Encode_and_Write_Structural_Properties__( newfile, dict_record )
                os.rename( dir_file_output_for_a_batch_being_written, dir_file_output_for_a_batch ) # rename the file once the output file has been written completely

            ''' initialize the next batch'''
            l_dict_record_for_a_batch = [ ] 
            int_index_batch += 1
            print( f"{int_n_records_processed} proteins has been processed at " + TIME_GET_timestamp( True ) ) # report progress

            if flag_all_records_were_retrieved :
                print( f"all {int_n_records_processed} proteins has been processed at " + TIME_GET_timestamp( True ) ) # report progress
                break
                
    ''' open output files of individual batches into a single output file '''
    print( f"start combining output files at " + TIME_GET_timestamp( True ) ) # report progress
    l_col = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB', 'structure_id___redundancy_reduced' ]
    OS_FILE_Combine_Files_in_order( glob.glob( f"{dir_folder_pipeline_temp}{name_file}.batch*_*-*.tsv.gz" ), f"{dir_folder_pipeline}{name_file}.tsv.gz", overwrite_existing_file = True, header = '\t'.join( l_col ) + '\n', delete_input_files = True ) # combine output files
def __Combine_Structural_Properties_and_Predict_Remaining_Residues__( dir_file_protein, dir_folder_pipeline = None, dir_folder_pipeline_temp = '/tmp/', flag_use_rcsb_pdb_only = False, int_number_of_proteins_in_a_batch_during_dnn_prediction = 1000, flag_use_all_gpu_devices = False ) :
    
    """
    Parse arguments
    """
    # directories of file and folders
    dir_file_protein = os.path.abspath( dir_file_protein )
    dir_folder_pipeline = os.path.abspath( dir_folder_pipeline )
    dir_folder_pipeline_temp = os.path.abspath( dir_folder_pipeline_temp )
    if dir_folder_pipeline[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline += '/'
    if dir_folder_pipeline_temp[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline_temp += '/'
    
    # fixed arguments
    float_search_thres_e_value = 30 # set e-value threshold for search (fixed for maximum sensitivity)

    # By default, set 'dir_folder_pipeline' as the folder where the given protein fasta file is located
    if dir_folder_pipeline is None :
        dir_folder_pipeline = dir_file_protein.rsplit( '/', 1 )[ 0 ] + '/'
    name_file = dir_file_protein.rsplit( '/', 1 )[ 1 ].rsplit( '.', 1 )[ 0 ] # retrieve name of file excluding file extension
    dir_file_protein_property = f"{dir_folder_pipeline}{name_file}.tsv.gz" # define output file (a file containing protein sequences and their estimated structural properties)
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # create a working directory of estimating structural properties
    os.makedirs( dir_folder_pipeline_struc, exist_ok = True )
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package
    
    """ read input proteins """
    dict_fasta_protein = FASTA_Read( dir_file_protein )
    dict_header = dict( ( h.split( ' ', 1 )[ 0 ], h ) for h in dict_fasta_protein ) # retrieve a dictionary of header (key = id_protein, value = header)
    dict_fasta_protein = dict( ( h.split( ' ', 1 )[ 0 ], dict_fasta_protein[ h ] ) for h in dict_fasta_protein ) # retrieve sequence_id by spliting the header at the first space (to makes sequence_id consistent with that used with blastp)
    
    """ 
    Load and Combine Structural Properties
    """
    # 2021-05-29 17:25:14 

    def __Parse_Structural_Properties__( df_sp, int_datatype ) :
        """
        # 2021-05-03 16:18:16 
        parse structural properties with typical parameters for parsing ascii-encoded structural properties and typical column names  
        Also, initializes a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
        'int_datatype' : integer identifying the current dataset. Used for initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
        """
        dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
        dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
        dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
        dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

        dict_sp = dict( )
        dict_sp[ 'acc' ] = ASCII_Decode( df_sp.set_index( 'header' ).rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1.to_dict( ), ** dict_kw_rsa ) # decode mkdssp outputs of RCSB PDB data
        dict_sp[ 'phi' ] = ASCII_Decode( df_sp.set_index( 'header' )[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
        dict_sp[ 'psi' ] = ASCII_Decode( df_sp.set_index( 'header' )[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
        dict_sp[ 'ss8' ] = ASCII_Decode( df_sp.set_index( 'header' )[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].to_dict( ), ** dict_kw_ss8 )
        dict_fasta = df_sp.set_index( 'header' )[ 'structure_id___redundancy_reduced' ].dropna( ).to_dict( )
        dict_sp[ 'structure_id' ] = dict( ( acc, Decode_List_of_Strings( dict_fasta[ acc ] ) ) for acc in dict_fasta )

        # initialization of a dictionary of arrays containing the identifier of the current dataset of 'acc' datatype
        dict_datatype_acc = dict( )
        for h in dict_sp[ 'acc' ] :
            arr = np.zeros_like( dict_sp[ 'acc' ][ h ] ) # default = 0 (predicted structural property)
            arr[ ~ np.isnan( dict_sp[ 'acc' ][ h ] ) ] = int_datatype
            dict_datatype_acc[ h ] = arr
        dict_sp[ 'datatype_acc' ] = dict_datatype_acc
        return dict_sp

    def __Combine_Structural_Properties__( dict_sp_1, dict_sp_2 ) :
        """
        # 2021-05-29 16:12:09 
        Combine structural properties from the two different sources. 
        Structural properties from the source with the lower priority will be overwritten with that of the higher priority

        'df_sp_1' : structural properties of the lower priority
        'df_sp_2' : structural properties with the higher priority
        """
        dict_sp_combined = dict( )

        for str_datatype in [ 'acc', 'phi', 'psi', 'ss8', 'structure_id', 'datatype_acc' ] :
            dict_arr_1 = dict_sp_1[ str_datatype ]
            dict_arr_2 = dict_sp_2[ str_datatype ]
            dict_arr_combined = dict( )
            for h in set( dict_arr_1 ).union( set( dict_arr_2 ) ) :
                # retrieve length of the protein sequence
                len_seq = 0
                for dict_arr in [ dict_arr_1, dict_arr_2 ] : 
                    if h in dict_arr :
                        len_seq = len( dict_arr[ h ] )

                arr_combined = np.full( len_seq, np.nan, dtype = object ) if str_datatype == 'structure_id' else np.full( len_seq, np.nan ) # initialize array containing combined data
                # update data in the order of increasing priorities (so that data with the highest priority is updated last and thus not overwritten by any other data)
                for dict_arr in [ dict_arr_2 ] if str_datatype == 'structure_id' else [ dict_arr_1, dict_arr_2 ] : # for 'structure_id', only use structure_id of the database of the higher priority ('RCSB_PDB' to be exact) 
                    if h in dict_arr :
                        arr_current = dict_arr[ h ]
                        mask_valid = ~ pd.isnull( arr_current ) if str_datatype == 'structure_id' else ~ np.isnan( arr_current )
                        arr_combined[ mask_valid ] = arr_current[ mask_valid ] # update data
                dict_arr_combined[ h ] = arr_combined
            dict_sp_combined[ str_datatype ] = dict_arr_combined
        return dict_sp_combined

    df_sp = pd.read_csv( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_rcsb_pdb.tsv.gz", sep = '\t' )
    dict_sp = __Parse_Structural_Properties__( df_sp, 2 ) # parse structural data from rcsb_pdb

    if not flag_use_rcsb_pdb_only :
        df_sp_swiss = pd.read_csv( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_swiss_model.tsv.gz", sep = '\t' )
        dict_sp_swiss = __Parse_Structural_Properties__( df_sp_swiss, 1 ) # parse structural data from swiss-model

        dict_sp = __Combine_Structural_Properties__( dict_sp_swiss, dict_sp ) # combine structural property data from RCSB_PDB and SWISS-MODEL

    ''' initialize arrays of proteins with no aligned structures '''
    for id_protein in dict_fasta_protein :
        len_seq = len( dict_fasta_protein[ id_protein ] ) # retrieve length of protein
        for str_datatype in [ 'acc', 'phi', 'psi', 'ss8', 'structure_id', 'datatype_acc' ] :
            if id_protein not in dict_sp[ str_datatype ] :
                ''' set default array for each datatype '''
                if str_datatype == 'datatype_acc' :
                    arr = np.full( len_seq, 0, dtype = float )
                elif str_datatype == 'structure_id' :
                    arr = np.full( len_seq, np.nan, dtype = object )
                else :
                    arr = np.full( len_seq, np.nan )
                dict_sp[ str_datatype ][ id_protein ] = arr

    """
    Prediction of Structural Properties Using Machine-Learning Models
    """
    # 2021-05-29 23:41:46 
    # use machine-learning model to predict structural properties
    ''' settings '''
    int_length_flanking = 74 # number of protein residues flanking the residue whose structural property values are predicted

    ''' download models '''
    file_name_model_rsa = "RSA.Conv1D.5_large_and_wide_dense.epoch_082.val_loss_0.1185.hdf5" # name of machine learning model used for structural property prediction (RSA)
    file_name_model_ss8 = "SS8.Conv1D.5_large_and_wide_dense.epoch_031.val_loss_0.7378.hdf5" # name of machine learning model used for structural property prediction (SS8)

    PKG.Download_Data( f"structural_property_estimation/{file_name_model_rsa}", dir_remote, name_package ) # download data
    PKG.Download_Data( f"structural_property_estimation/{file_name_model_ss8}", dir_remote, name_package ) # download data

    ''' prepare encoding '''
    def Index_List( l, start = 0 ) :
        ''' # 2021-05-06 22:32:53 
        index a given list containing hashable element with integer index starting with 'start' (default: 0)
        can be used for 'one-hot encoding'
        '''
        return dict( ( e, i ) for e, i in zip( l, np.arange( start, len( l ) + start, 1, dtype = int ) ) )

    # input structure : amino-acids (one-hot encoding) + ss8 (one-hot encoding) + rsa (0~1)
    int_index_start_amino_acid = 0
    l_amino_acids = [ 'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y' ]
    dict_amino_acid_to_index = Index_List( l_amino_acids, int_index_start_amino_acid ) # one-hot encoding dictionary for 20 amino acids
    int_index_start_ss8 = len( l_amino_acids )
    l_ss8 = [ 'G', 'H', 'I', 'E', 'B', 'T', 'S', 'C' ]
    dict_ss8_to_index = Index_List( l_ss8, int_index_start_ss8 ) # one-hot encoding dictionary for 20 amino acids
    index_rsa = len( l_amino_acids ) + len( l_ss8 ) # index of rsa data in the input array


    ''' predict and update structural properties using ML models for each batch '''
    import tensorflow as tf
    from tensorflow.keras import layers

    l_id_protein_all = list( dict_sp[ 'acc' ] ) # retrieve a list containing all 'id_protein'
    n_batches = int( np.ceil( len( l_id_protein_all ) / int_number_of_proteins_in_a_batch_during_dnn_prediction ) ) # number of batches
    ''' for each batch '''
    for l_id_protein in LIST_Split( l_id_protein_all, n_split = n_batches ) :
        ''' perform encoding for the current batch '''
        # 'l_id_protein': id_protein -> integer encoding to reduce memory usage
        # perform encoding and record metadata of each position
        x_data = [ ]
        meta_data = [ ]

        for index_protein, id_protein in enumerate( l_id_protein ) :
            seq = dict_fasta_protein[ id_protein ]
            length = len( seq ) # retrieve the length of the protein

            arr_acc = dict_sp[ 'acc' ][ id_protein ]
            arr_ss8 = dict_sp[ 'ss8' ][ id_protein ]

            ''' one-hot encoding (+ RSA) of the current protein '''
            arr_protein = np.zeros( ( length, len( l_amino_acids ) + len( l_ss8 ) + 1 ), dtype = np.float32 ) # initialize one-hot encoding of the proteinarray
            for i in range( length ) :
                # one-hot encoding of amino acid
                aa = seq[ i ] # retrieve amino acid
                if aa in dict_amino_acid_to_index : # invalid amino acid will be considered as a missing value (empty row)
                    arr_protein[ i, dict_amino_acid_to_index[ aa ] ] = 1 
                # one-hot encoding of ss8
                if not np.isnan( arr_ss8[ i ] ) : # missing value = empty row
                    arr_protein[ i, dict_ss8_to_index[ l_ss8[ int( arr_ss8[ i ] ) ] ] ] = 1
                # encode rsa information
                if not np.isnan( arr_acc[ i ] ) : # missing value = empty row
                    arr_protein[ i, index_rsa ] = arr_acc[ i ]

            ''' prepare one-hot encoding input array for each missing position '''
            ''' assumes positions of missing values in 'arr_acc' are identical to those of 'arr_ss8' '''
            for int_pos in np.where( np.isnan( arr_acc ) )[ 0 ] : 
                start = max( 0, int_pos - int_length_flanking )
                end = min( length, int_pos + int_length_flanking + 1 )

                arr = np.zeros( ( int_length_flanking + 1 + int_length_flanking, len( l_amino_acids ) + len( l_ss8 ) + 1 ), dtype = np.float32 ) # initialize input array
                start_input = max( 0, int_length_flanking - int_pos ) # index of start position of the input array
                arr[ start_input : start_input + end - start ] = arr_protein[ start : end ]
                x_data.append( arr )
                meta_data.append( [ index_protein, int_pos ] )
        x_data = np.array( x_data, dtype = np.float32 )

        ''' predict and update structural properties using ML models for the current batch '''
        for dir_file_model, str_datatype in zip( [ f"{dir_folder_cressp}structural_property_estimation/{file_name_model_rsa}", f"{dir_folder_cressp}structural_property_estimation/{file_name_model_ss8}" ], [ 'acc', 'ss8' ] ) :

            model = tf.keras.models.load_model( dir_file_model )

            try :
                dir_file_during_gpu_usage = f"{dir_folder_cressp}structural_property_estimation/{UUID( )}.flag"
                with tf.distribute.MirroredStrategy( ).scope( ) if flag_use_all_gpu_devices else open( dir_file_during_gpu_usage, 'w' ) : # distributed calculation if 'flag_use_all_gpus' is True
                    y_pred = model.predict( x_data ) # predict structural properties using a ML model
                    y_pred = y_pred.argmax( axis = 1 ) if str_datatype == 'ss8' else y_pred # multiclass classification for 'ss8' and regression for 'acc'
                    os.remove( dir_file_during_gpu_usage )
            except tf.errors.InternalError : # when out of memory error occurs during training, use CPU instead
                os.environ[ 'CUDA_VISIBLE_DEVICES' ] = '-1' # forcing tensorflow to use CPU
                y_pred = model.predict( x_data ) # predict structural properties using a ML model
                y_pred = y_pred.argmax( axis = 1 ) if str_datatype == 'ss8' else y_pred # multiclass classification for 'ss8' and regression for 'acc'

            for i in range( len( y_pred ) ) :
                index_protein, int_pos = meta_data[ i ] # parse metadata
                id_protein = l_id_protein[ index_protein ]
                dict_sp[ str_datatype ][ id_protein ][ int_pos ] = y_pred[ i ] # update structural property data with predicted data

    """
    Save Result
    """
    def __Encode_Structural_Properties__( dict_sp, dict_fasta_protein ) :
        """
        # 2021-05-29 23:56:59 
        Compose a dataframe containing encoded structural properties using given 'dict_sp' and 'dict_fasta_protein'

        'dict_sp' : dictionary containing structural properties
        'dict_fasta_protein' : dictionary containing protein sequences

        """
        dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
        dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
        dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
        dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

        ''' initialize dataframe with protein sequences '''
        df_sp = pd.Series( dict_fasta_protein, name = 'seq' ).reset_index( ).rename( columns = { 'index' : 'id_protein' } ).set_index( 'id_protein' )

        # encode combined structural properties into ASCII strings using ASCII encoding
        df_sp[ 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1' ] = pd.Series( ASCII_Encode( dict_sp[ 'acc' ], ** dict_kw_rsa ) )
        df_sp[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'phi' ], ** dict_kw_torsion_angle ) )
        df_sp[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ] = pd.Series( ASCII_Encode( dict_sp[ 'psi' ], ** dict_kw_torsion_angle ) )
        df_sp[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ] = pd.Series( ASCII_Encode( dict_sp[ 'ss8' ], ** dict_kw_ss8 ) )
        df_sp[ 'rsa_datatype___ascii_encoding_1_character_from_33_to_36__states_Pred_Model_PDB' ] = pd.Series( ASCII_Encode( dict_sp[ 'datatype_acc' ], ** dict_kw_datatype ) )
        df_sp[ 'structure_id___redundancy_reduced' ] = pd.Series( dict( ( h, Encode_List_of_Strings( dict_sp[ 'structure_id' ][ h ] ) ) for h in dict_sp[ 'structure_id' ] ) )
        df_sp.reset_index( drop = False, inplace = True )
        return df_sp

    df_sp = __Encode_Structural_Properties__( dict_sp, dict_fasta_protein ) # encode structural properties into a dataframe
    df_sp[ 'fasta_header' ] = list( dict_header[ i ] for i in df_sp.id_protein.values ) # retrieve headers
    df_sp.to_csv( f"{dir_folder_pipeline}{name_file}.tsv.gz", sep = '\t', index = False ) # save structural properties of input proteins as a tabular file
def Estimate_structural_property( dir_file_protein, n_threads, dir_folder_pipeline = None, dir_folder_pipeline_temp = '/tmp/', flag_use_rcsb_pdb_only = False, int_number_of_proteins_in_a_batch_during_dnn_prediction = 1000, flag_use_all_gpu_devices = False, flag_predict_structural_properties_of_remaining_residues_using_dnn = True ) :
    """
    # 2021-05-31 15:13:13 
    Estimate structural properties of given proteins, and write a tsv file containing structural properties of the proteins.
    When CUDA-enabled GPU is available to TensorFlow, GPUs will be used to predict structural properties of protein residues not covered by known and/or predicted protein structures through homology-modeling.
    
    'dir_file_protein' : unzipped fasta file
    'flag_predict_structural_properties_of_remaining_residues_using_dnn' : (default: True) predict structural properties of the remaining residues. if set to False, output database will only contains transferred structural properties
    """
    
    """
    Parse arguments
    """
    # By default, set 'dir_folder_pipeline' as the folder where the given protein fasta file is located
    if dir_folder_pipeline is None :
        dir_folder_pipeline = dir_file_protein.rsplit( '/', 1 )[ 0 ] + '/'
    
    # directories of file and folders
    dir_file_protein = os.path.abspath( dir_file_protein )
    dir_folder_pipeline = os.path.abspath( dir_folder_pipeline )
    dir_folder_pipeline_temp = os.path.abspath( dir_folder_pipeline_temp )
    if dir_folder_pipeline[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline += '/'
    if dir_folder_pipeline_temp[ -1 ] != '/' : # last character of a directory should be '/'
        dir_folder_pipeline_temp += '/'
    
    # fixed arguments
    float_search_thres_e_value = 30 # set e-value threshold for search (fixed for maximum sensitivity)


    name_file = dir_file_protein.rsplit( '/', 1 )[ 1 ].rsplit( '.', 1 )[ 0 ] # retrieve name of file excluding file extension
    dir_file_protein_property = f"{dir_folder_pipeline}{name_file}.tsv.gz" # define output file (a file containing protein sequences and their estimated structural properties)
    dir_folder_pipeline_struc = f'{dir_folder_pipeline}struc/' # create a working directory of estimating structural properties
    os.makedirs( dir_folder_pipeline_struc, exist_ok = True )

    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package


    """
    RCSB_PDB 
    """
    dir_file_db_rcsb_pdb = f'{dir_folder_cressp}data/pdb/rcsb_pdb.tsv.gz'
    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline_struc}{name_file}.blastp_rcsb_pdb.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :
        # download and parse associated data (RCSB_PDB)
        PKG.Download_Data( "data/pdb/rcsb_pdb.tsv.gz", dir_remote, name_package ) # download data
        PKG.Download_Data( "data/pdb/rcsb_pdb.label_seq_id.start_end.tsv.gz", dir_remote, name_package ) # download data

        df_protein_rcsb_pdb = pd.read_csv( dir_file_db_rcsb_pdb, sep = '\t' )

        dir_file_protein_rcsb_pdb = f'{dir_folder_cressp}data/pdb/rcsb_pdb.fa' 
        if not os.path.exists( dir_file_protein_rcsb_pdb ) : # if rcsb_pdb fasta file is not available, extract the sequence from the dataframe
            dict_fasta_protein_rcsb_pdb = df_protein_rcsb_pdb.set_index( 'header' ).seq.to_dict( )
            FASTA_Write( dir_file_protein_rcsb_pdb, dict_fasta = dict_fasta_protein_rcsb_pdb )

        # create blastp_db using rcsb_pdb protein sequences
        dir_prefix_blastdb_protein_rcsb_pdb = f"{dir_folder_cressp}data/pdb/makeblastdb_out/protein_rcsb_pdb"
        if not os.path.exists( f"{dir_prefix_blastdb_protein_rcsb_pdb}.psq" ) : # check whether blastdb exists
            os.makedirs( f"{dir_folder_cressp}data/pdb/makeblastdb_out/", exist_ok = True )  
            OS_Run( [ "makeblastdb", "-in", dir_file_protein_rcsb_pdb, '-dbtype', 'prot', '-parse_seqids', '-max_file_sz', '1GB', '-out', dir_prefix_blastdb_protein_rcsb_pdb ], dir_file_stdout = f"{dir_prefix_blastdb_protein_rcsb_pdb}.makeblastdb.stdout.txt", dir_file_stderr = f"{dir_prefix_blastdb_protein_rcsb_pdb}.makeblastdb.stderr.txt", return_output = False ) # make blast db for protein_query

        # run blastp
        dir_file_blastp_output = f'{dir_folder_pipeline_struc}{name_file}.blastp_rcsb_pdb.tsv'
        OS_Run( [ 'blastp', '-query', dir_file_protein, '-db', dir_prefix_blastdb_protein_rcsb_pdb, '-out', dir_file_blastp_output, '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop', '-num_threads', f'{n_threads}', '-evalue', f'{float_search_thres_e_value}' ], dir_file_stdout = f"{dir_file_blastp_output}.blastp.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.blastp.stderr.txt", return_output = False ) # run blastp
        OS_Run( [ 'gzip', dir_file_blastp_output ], dir_file_stdout = f"{dir_file_blastp_output}.gzip.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.gzip.stderr.txt", return_output = False ) # compress blastp output
        dir_file_blastp_output += '.gz'
        
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )


    """
    SWISS-MODEL
    """
    dir_file_db_swiss_model = f'{dir_folder_cressp}data/pdb/swiss_model.tsv.gz'
    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline_struc}{name_file}.blastp_swiss_model.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :
        if not flag_use_rcsb_pdb_only :
            # download and parse associated data (RCSB_PDB)
            PKG.Download_Data( "data/pdb/swiss_model.tsv.gz", dir_remote, name_package ) # download data

            df_protein_swiss_model = pd.read_csv( dir_file_db_swiss_model, sep = '\t' )

            dir_file_protein_swiss_model = f'{dir_folder_cressp}data/pdb/swiss_model.fa' 
            if not os.path.exists( dir_file_protein_swiss_model ) : # if rcsb_pdb fasta file is not available, extract the sequence from the dataframe
                dict_fasta_protein_swiss_model = df_protein_swiss_model.set_index( 'header' ).seq.to_dict( )
                FASTA_Write( dir_file_protein_swiss_model, dict_fasta = dict_fasta_protein_swiss_model )

            # create blastp_db using rcsb_pdb protein sequences
            dir_prefix_blastdb_protein_swiss_model = f"{dir_folder_cressp}data/pdb/makeblastdb_out/protein_swiss_model"
            if not os.path.exists( f"{dir_prefix_blastdb_protein_swiss_model}.psq" ) : # check whether blastdb exists
                os.makedirs( f"{dir_folder_cressp}data/pdb/makeblastdb_out/", exist_ok = True )  
                OS_Run( [ "makeblastdb", "-in", dir_file_protein_swiss_model, '-dbtype', 'prot', '-parse_seqids', '-max_file_sz', '1GB', '-out', dir_prefix_blastdb_protein_swiss_model ], dir_file_stdout = f"{dir_prefix_blastdb_protein_swiss_model}.makeblastdb.stdout.txt", dir_file_stderr = f"{dir_prefix_blastdb_protein_swiss_model}.makeblastdb.stderr.txt", return_output = False ) # make blast db for protein_query

            # run blastp
            dir_file_blastp_output = f'{dir_folder_pipeline_struc}{name_file}.blastp_swiss_model.tsv'
            OS_Run( [ 'blastp', '-query', dir_file_protein, '-db', dir_prefix_blastdb_protein_swiss_model, '-out', dir_file_blastp_output, '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop', '-num_threads', f'{n_threads}', '-evalue', f'{float_search_thres_e_value}' ], dir_file_stdout = f"{dir_file_blastp_output}.blastp.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.blastp.stderr.txt", return_output = False ) # run blastp
            OS_Run( [ 'gzip', dir_file_blastp_output ], dir_file_stdout = f"{dir_file_blastp_output}.gzip.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.gzip.stderr.txt", return_output = False ) # compress blastp output
            dir_file_blastp_output += '.gz'
            
            """ set flag """
            with open( dir_file_flag, 'w' ) as newfile :
                newfile.write( 'completed\n' )
                
    """
    AlphaFold Structure Database
    """
    dir_file_db_alphafold_structure_db = f'{dir_folder_cressp}data/pdb/AlphaFoldStructureDB.tsv.gz'
    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline_struc}{name_file}.blastp_alphafold_structure_db.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :
        if not flag_use_rcsb_pdb_only :
            # download and parse associated data (RCSB_PDB)
            PKG.Download_Data( "data/pdb/AlphaFoldStructureDB.tsv.gz", dir_remote, name_package ) # download data

            df_protein_alphafold_structure_db = pd.read_csv( dir_file_db_alphafold_structure_db, sep = '\t' )

            dir_file_protein_alphafold_structure_db = f'{dir_folder_cressp}data/pdb/alphafold_structure_db.fa' 
            if not os.path.exists( dir_file_protein_alphafold_structure_db ) : # if rcsb_pdb fasta file is not available, extract the sequence from the dataframe
                dict_fasta_protein_alphafold_structure_db = df_protein_alphafold_structure_db.set_index( 'header' ).seq.to_dict( )
                FASTA_Write( dir_file_protein_alphafold_structure_db, dict_fasta = dict_fasta_protein_alphafold_structure_db )

            # create blastp_db using rcsb_pdb protein sequences
            dir_prefix_blastdb_protein_alphafold_structure_db = f"{dir_folder_cressp}data/pdb/makeblastdb_out/protein_alphafold_structure_db"
            if not os.path.exists( f"{dir_prefix_blastdb_protein_alphafold_structure_db}.psq" ) : # check whether blastdb exists
                os.makedirs( f"{dir_folder_cressp}data/pdb/makeblastdb_out/", exist_ok = True )  
                OS_Run( [ "makeblastdb", "-in", dir_file_protein_alphafold_structure_db, '-dbtype', 'prot', '-parse_seqids', '-max_file_sz', '1GB', '-out', dir_prefix_blastdb_protein_alphafold_structure_db ], dir_file_stdout = f"{dir_prefix_blastdb_protein_alphafold_structure_db}.makeblastdb.stdout.txt", dir_file_stderr = f"{dir_prefix_blastdb_protein_alphafold_structure_db}.makeblastdb.stderr.txt", return_output = False ) # make blast db for protein_query

            # run blastp
            dir_file_blastp_output = f'{dir_folder_pipeline_struc}{name_file}.blastp_alphafold_structure_db.tsv'
            OS_Run( [ 'blastp', '-query', dir_file_protein, '-db', dir_prefix_blastdb_protein_alphafold_structure_db, '-out', dir_file_blastp_output, '-outfmt', '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop', '-num_threads', f'{n_threads}', '-evalue', f'{float_search_thres_e_value}' ], dir_file_stdout = f"{dir_file_blastp_output}.blastp.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.blastp.stderr.txt", return_output = False ) # run blastp
            OS_Run( [ 'gzip', dir_file_blastp_output ], dir_file_stdout = f"{dir_file_blastp_output}.gzip.stdout.txt", dir_file_stderr = f"{dir_file_blastp_output}.gzip.stderr.txt", return_output = False ) # compress blastp output
            dir_file_blastp_output += '.gz'
            
            """ set flag """
            with open( dir_file_flag, 'w' ) as newfile :
                newfile.write( 'completed\n' )

    """ 
    Transfer Structural Properties from Experimental & Modeled Protein Structures and Combine Structural Properties
    """
    """ check flag """
    dir_file_flag2 = f"{dir_folder_pipeline_struc}{name_file}_transferred_combined.tsv.gz.completed.flag"
    if not os.path.exists( dir_file_flag2 ) :
        # transfer surface accessibility data from structures to proteins
        """ initialize """
        global dict_index_df_blastp, arr_data_df_blastp, dict_acc_to_arr_acc_dssp, dict_acc_to_arr_phi_dssp, dict_acc_to_arr_psi_dssp, dict_acc_to_arr_ss8_dssp, dict_kw_rsa, dict_kw_torsion_angle, dict_kw_ss8, dict_kw_datatype, dict_fasta_protein # use global variables
        # set global setting for encoding & decoding structural info.
        dict_kw_rsa = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = 0, value_max = 1 )
        dict_kw_torsion_angle = dict( ascii_min = 33, ascii_max = 126, l_ascii_to_exclude = [ 62 ], n_char = 2, value_min = -180, value_max = 180 )
        dict_kw_ss8 = dict( ascii_min = 33, ascii_max = 41, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 8 )
        dict_kw_datatype = dict( ascii_min = 33, ascii_max = 36, l_ascii_to_exclude = [ 62 ], n_char = 1, value_min = 0, value_max = 3 )

        """ read input proteins """
        dict_fasta_protein = FASTA_Read( dir_file_protein )
        dict_header = dict( ( h.split( ' ', 1 )[ 0 ], h ) for h in dict_fasta_protein ) # retrieve a dictionary of header (key = id_protein, value = header)
        dict_fasta_protein = dict( ( h.split( ' ', 1 )[ 0 ], dict_fasta_protein[ h ] ) for h in dict_fasta_protein ) # retrieve sequence_id by spliting the header at the first space (to makes sequence_id consistent with that used with blastp)

        for dir_file_db, name_dataset in zip( [ dir_file_db_rcsb_pdb ] if flag_use_rcsb_pdb_only else [ dir_file_db_rcsb_pdb, dir_file_db_swiss_model, dir_file_db_alphafold_structure_db ], [ 'rcsb_pdb' ] if flag_use_rcsb_pdb_only else [ 'rcsb_pdb', 'swiss_model', 'alphafold_structure_db' ] ) :
            print( dir_file_db, name_dataset )
            """ check flag """
            dir_file_flag = f"{dir_folder_pipeline_struc}{name_file}_transferred_from_{name_dataset}.tsv.gz.completed.flag"
            if not os.path.exists( dir_file_flag ) :
                """ load blastp output """
                dir_file_blastp_output = f'{dir_folder_pipeline_struc}{name_file}.blastp_{name_dataset}.tsv.gz'
                dir_prefix_blastp_output = dir_file_blastp_output.rsplit( '.tsv.gz', 1 )[ 0 ]


                if not os.path.exists( f"{dir_prefix_blastp_output}.with_aligned_seq.filtered.tsv.gz" ) :
                    """ process blastp output """
                    """ filter alignment  """
                    # set threshold values for transferring structural properties from structures to protein sequences
                    float_transfer_pidenta = 70 # percent identity
                    float_transfer_evalueb = 1e-8 # significance of the alignment 
                    df_blastp = BLAST_Read( dir_file_blastp_output, dict_qaccver_to_seq = dict_fasta_protein, dir_file_output = f"{dir_prefix_blastp_output}.with_aligned_seq.filtered.saving.tsv.gz", float_transfer_pidenta = float_transfer_pidenta, float_transfer_evalueb = float_transfer_evalueb ) # save filtered BLASTP result
                    os.rename( f"{dir_prefix_blastp_output}.with_aligned_seq.filtered.saving.tsv.gz", f"{dir_prefix_blastp_output}.with_aligned_seq.filtered.tsv.gz" ) # rename the saved file once write operation is completed
                df_blastp = pd.read_csv( f"{dir_prefix_blastp_output}.with_aligned_seq.filtered.tsv.gz", sep = '\t' )

                ''' check whether the blastp result is empty '''
                if len( df_blastp ) == 0 :
                    pd.DataFrame( [ ], columns = l_col_df_transferred ).to_csv( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_{name_dataset}.tsv.gz", sep = '\t', index = False ) # save an empty output file
                    continue
                l_col_df_transferred = [ 'id_protein', 'seq', 'rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1', 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180', 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C', 'structure_id___redundancy_reduced' ]
                ''' load blastp result globally '''
                dict_index_df_blastp = DF_Build_Index_Using_Dictionary( df_blastp, l_col_for_index = 'qaccver' )
                arr_data_df_blastp = df_blastp.values

                """ load dssp db globally """
                df_dssp = PD_Select( pd.read_csv( dir_file_db, sep = '\t' ), header = df_blastp.saccver.unique( ) ) # subset structural data so that only structures aligned to the input proteins are included
                dict_acc_to_arr_acc_dssp = ASCII_Decode( df_dssp.set_index( 'header' ).rsa___ascii_encoding_2_characters_from_33_to_126__from_0_to_1.to_dict( ), ** dict_kw_rsa ) # decode mkdssp outputs of RCSB PDB data
                dict_acc_to_arr_phi_dssp = ASCII_Decode( df_dssp.set_index( 'header' )[ 'phi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
                dict_acc_to_arr_psi_dssp = ASCII_Decode( df_dssp.set_index( 'header' )[ 'psi___ascii_encoding_2_characters_from_33_to_126__from_-180_to_180' ].to_dict( ), ** dict_kw_torsion_angle )
                dict_acc_to_arr_ss8_dssp = ASCII_Decode( df_dssp.set_index( 'header' )[ 'ss8___ascii_encoding_1_character_from_33_to_41__states_G_H_I_E_B_T_S_C' ].to_dict( ), ** dict_kw_ss8 )

                Multiprocessing( df_blastp.qaccver.unique( ), __Transfer_DSSP_Structural_Property_Through_BLAST__, n_threads, dir_temp = dir_folder_pipeline_temp, global_arguments = [ dir_folder_cressp, dir_folder_pipeline_temp ] ) # transfer structural properties with multiple processing
                OS_FILE_Combine_Files_in_order( glob.glob( f"{dir_folder_pipeline_temp}*.transferred.tsv.gz" ), f"{dir_folder_pipeline_struc}{name_file}_transferred_from_{name_dataset}.tsv.gz", overwrite_existing_file = True, header = '\t'.join( l_col_df_transferred ) + '\n', delete_input_files = True ) # combine output files

                """ set flag """
                with open( dir_file_flag, 'w' ) as newfile :
                    newfile.write( 'completed\n' )

        ''' combine structural properties transferred from SWISS-MODEL and RCSB_PDB '''
        if flag_use_rcsb_pdb_only :
            ''' when only RCSB_PDB was used, simply rename the output file '''
            os.rename( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_rcsb_pdb.tsv.gz", f"{dir_folder_pipeline_struc}{name_file}_transferred_combined.tsv.gz" )
        else :
            ''' combine structural properties from RCSB with those from SWISS-MODELs '''
            __Iterate_and_Combine_Structural_Properties__( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_rcsb_pdb.tsv.gz", 2, f"{dir_folder_pipeline_struc}{name_file}_transferred_from_swiss_model.tsv.gz", 1, f"{dir_folder_pipeline_struc}{name_file}_transferred_from_rcsb_pdb_and_swiss_model.tsv.gz", flag_debugging = False )
            ''' combine structural properties from RCSB + SWISS-MODELS with those from AlphaFold structure DB '''
            __Iterate_and_Combine_Structural_Properties__( f"{dir_folder_pipeline_struc}{name_file}_transferred_from_rcsb_pdb_and_swiss_model.tsv.gz", 2, f"{dir_folder_pipeline_struc}{name_file}_transferred_from_alphafold_structure_db.tsv.gz", 1, f"{dir_folder_pipeline_struc}{name_file}_transferred_combined.tsv.gz", flag_debugging = False )

        """ set flag """
        with open( dir_file_flag2, 'w' ) as newfile :
            newfile.write( 'completed\n' )

    """ 
    Predict Structural Properties of Remaining Residues
    """  
    """ check flag """
    dir_file_flag = f"{dir_folder_pipeline}{name_file}.tsv.gz.completed.flag"
    if not os.path.exists( dir_file_flag ) :
        if flag_predict_structural_properties_of_remaining_residues_using_dnn :
            __Predict_Structural_Properties_of_Remaining_Residues__( dir_file_protein, dir_folder_pipeline, dir_folder_pipeline_temp, int_number_of_proteins_in_a_batch_during_dnn_prediction, flag_use_all_gpu_devices ) # predict RSA
        else :
            shutil.copyfile( f"{dir_folder_pipeline_struc}{name_file}_transferred_combined.tsv.gz", f"{dir_folder_pipeline}{name_file}.tsv.gz" ) # simply copy the the file containing combined transferred structural properties as the output file
            
        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )


