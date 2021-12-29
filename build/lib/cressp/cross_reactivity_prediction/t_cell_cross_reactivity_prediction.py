from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

def Predict_T_cell_cross_reactivity( dir_folder_pipeline, float_thres_avg_blosum62_score_for_mhc, float_thres_min_mhc_allele_frequency, float_thres_binding_affinities_in_nM, flag_replace_unconventional_acid_code ) :
    """
    Predict_T_cell_cross_reactivity
    """

    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package

    ''' read dict_blosum62 from the tsv file '''
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score

    """ 
    Extract subsequences from the alignments, filter, predict binding affinity, and prepare an output file containing individual records
    """

    # setting
    int_thres_min_length_peptide_mhc_i = 8
    int_thres_max_length_peptide_mhc_i = 15

    ''' define and open an output folder '''
    dir_file_output = f"{dir_folder_pipeline}t_cell.mhc_binding.tsv.gz"
    newfile = gzip.open( dir_file_output, 'wb' )
    l_col = [ 'mhc_allele', 'affinity_query', 'affinity_target', 'source', 'query_accession', 'target_accession', 'e_value', 'score_blosum', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'window_size', 'average_score_blosum', 'query_processing_score', 'target_processing_score', 'mhc_class', 'score_for_sorting' ]
    newfile.write( ( '\t'.join( l_col ) + '\n' ).encode( ) ) # write a header line for the output file

    ''' load MHCFlurry binding affinity predictor '''
    from mhcflurry import Class1PresentationPredictor # load MHCflurry predictor
    predictor = Class1PresentationPredictor.load( ) # load MHCflurry predictor for benchmarking

    ''' retrieve MHC-I alleles with allele frequency larger then the given threshold in at least one population '''
    PKG.Download_Data( "data/mhc_population_allele_frequency.tsv.gz", dir_remote, name_package ) # download data
    PKG.Download_Data( "data/autoimmune_disease_associated_mhc_alleles.tsv.gz", dir_remote, name_package ) # download data
    df_mhc_af = pd.read_csv( f"{dir_folder_cressp}data/mhc_population_allele_frequency.tsv.gz", sep = '\t', index_col = [ 0, 1 ] )
    l_mhc_i_allele = list( PD_Threshold( ( df_mhc_af.loc[ 'I' ] > float_thres_min_mhc_allele_frequency ).sum( axis = 1 ), a = 0 ).index.values ) + list( pd.read_csv( f"{dir_folder_cressp}data/autoimmune_disease_associated_mhc_alleles.tsv.gz", sep = '\t' ).mhc_allele.unique( ) ) # retrieve mhc_i_alleles with allele frequency > 'float_thres_min_mhc_allele_frequency' in at least one population # also retrieve autoimmune-associated alleles
    l_mhc_i_allele = list( set( l_mhc_i_allele ).intersection( predictor.supported_alleles ) ) # retrive valid alleles for MHCflurry

    ''' prepare df_matched before iteration '''
    if 'df_matched' not in locals( ) :
        df_matched = pd.read_csv( f'{dir_folder_pipeline}matched.tsv.gz', sep = '\t' ) # read alignments between query and target protein sequences
    df_matched_for_mhc_i = df_matched[ df_matched.query_alignment.apply( len ) >= int_thres_max_length_peptide_mhc_i ]
    iter_arr_df_matched_for_mhc_i = ( arr for arr in df_matched_for_mhc_i[ [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'source' ] ].values )

    ''' initialize a bucket '''
    flag_file_end_reached = False
    l_l_value = [ ] # initialize a bucket
    int_max_num_line_before_start_prediction = 50000 # since predicting multiple records in a batch is much faster, records will be collected in a bucket before they are processed together

    while True :
        ''' iterate 'df_matched_for_mhc_i' '''
        try :
            query_accession, target_accession, query_start, query_end, target_start, target_end, query_alignment, target_alignment, e_value, source = next( iter_arr_df_matched_for_mhc_i )
        except StopIteration :
            flag_file_end_reached = True

        ''' extract and filter subsequences '''
        if not flag_file_end_reached :
            try : l_similarity_score = list( dict_blosum62[ ( amino_acid_query, amino_acid_target ) ] for amino_acid_query, amino_acid_target in zip( query_alignment.upper( ), target_alignment.upper( ) ) ) # calculate similarity score for each residue
            except : 
                continue # in case of ambiguous amino acids (e.g. J = I or L, B = D or N), skip the protein itself
            ''' iterate different int_window_size '''
            for int_window_size in np.arange( int_thres_min_length_peptide_mhc_i, int_thres_max_length_peptide_mhc_i + 1 ) :
                n_windows = len( l_similarity_score ) + 1 - int_window_size
                for int_start in range( n_windows ) : 
                    int_similarity_score_for_window = sum( l_similarity_score[ int_start : int_start + int_window_size ] ) 
                    query_alignment_subsequence = query_alignment[ int_start : int_start + int_window_size ]
                    target_alignment_subsequence = target_alignment[ int_start : int_start + int_window_size ]
                    int_gap_count_subsequence_query, int_gap_count_subsequence_target = query_alignment_subsequence.count( '-' ), target_alignment_subsequence.count( '-' )
                    int_gap_count_before_subsequence_query, int_gap_count_before_subsequence_target = query_alignment[ : int_start ].count( '-' ), target_alignment[ : int_start ].count( '-' )
                    # [ 'source', 'query_accession', 'target_accession', 'e_value', 'score_blosum', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'window_size', 'average_score_blosum' ]
                    ''' filter subsequences based on average blosum62 scores '''
                    average_score_blosum = int_similarity_score_for_window / int_window_size # retrieve an average blosum score
                    if average_score_blosum >= float_thres_avg_blosum62_score_for_mhc :
                        l_l_value.append( [ source, query_accession, target_accession, e_value, int_similarity_score_for_window, query_start + int_start - int_gap_count_before_subsequence_query, query_start + int_start + int_window_size - 1 - int_gap_count_subsequence_query - int_gap_count_before_subsequence_query, target_start + int_start - int_gap_count_before_subsequence_target, target_start + int_start + int_window_size - 1 - int_gap_count_subsequence_target - int_gap_count_before_subsequence_target, query_alignment_subsequence, target_alignment_subsequence, int_window_size, average_score_blosum ] )        

        ''' check whether the bucket has been filled, and predict binding affinity of the peptides '''
        if len( l_l_value ) >= int_max_num_line_before_start_prediction or flag_file_end_reached : # detect whether the bucket has been filled
            ''' compose a dataframe from the collected records '''
            df_subsequence_mhc_i = pd.DataFrame( l_l_value, columns = [ 'source', 'query_accession', 'target_accession', 'e_value', 'score_blosum', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence', 'window_size', 'average_score_blosum' ] )
            ''' initialize the next bucket ''' 
            l_l_value = [ ]

            ''' remove subsequences containing invalid amino acids, and retrieve subsequences without gaps '''
            l_aa_invalid = list( "XBJZ" ) if flag_replace_unconventional_acid_code else list( "OBJUZX" ) # define invalid amino acids based on settings 'flag_replace_unconventional_acid_code'
            df_subsequence_mhc_i = PD_Search( df_subsequence_mhc_i, query_subsequence = l_aa_invalid, target_subsequence = l_aa_invalid, is_negative_query = True ) # retrive subset of df_subsequence for MHC_I allele binding prediction # invalid residues 'B' and 'Z' should be removed 
            df_subsequence_mhc_i[ 'query_subsequence_without_gap' ] = list( seq.replace( '-', '' ).replace( 'U', 'C' ).replace( 'O', 'Y' ) for seq in df_subsequence_mhc_i.query_subsequence.values ) if flag_replace_unconventional_acid_code else list( seq.replace( '-', '' ) for seq in df_subsequence_mhc_i.query_subsequence.values ) # remove gap in the aligment and use the sequence as inputs
            df_subsequence_mhc_i[ 'target_subsequence_without_gap' ] = list( seq.replace( '-', '' ).replace( 'U', 'C' ).replace( 'O', 'Y' ) for seq in df_subsequence_mhc_i.target_subsequence.values ) if flag_replace_unconventional_acid_code else list( seq.replace( '-', '' ) for seq in df_subsequence_mhc_i.target_subsequence.values )

            int_min_len_peptide = 5 # lower limit of the length of the peptide for binding affinity prediction
            df_subsequence_mhc_i = df_subsequence_mhc_i[ ( df_subsequence_mhc_i.query_subsequence_without_gap.apply( len ) >= int_min_len_peptide ) & ( df_subsequence_mhc_i.target_subsequence_without_gap.apply( len ) >= int_min_len_peptide ) ] # drop aligned peptides with peptide shorter than the minimum peptide length of prediction algorithm (MHCflurry)

            ''' predict binding affinity of peptides pairs with MHCflurry '''
            df_subsequence_mhc_i.reset_index( drop = True, inplace = True ) # reset index before join operation
            dict_data = dict( )
            for str_allele in l_mhc_i_allele :
                df_predicted_query_subsequence = predictor.predict( peptides = list( df_subsequence_mhc_i.query_subsequence_without_gap.values ), alleles = [ str_allele ], verbose = 0 ) 
                df_predicted_target_subsequence = predictor.predict( peptides = list( df_subsequence_mhc_i.target_subsequence_without_gap.values ), alleles = [ str_allele ], verbose = 0 ) # remove gap in the aligment and use the sequence as inputs
                if 'query_processing_score' not in dict_data : dict_data[ 'query_processing_score' ] = df_predicted_query_subsequence.processing_score.values # retreve processing score only once for query and target peptides
                if 'target_processing_score' not in dict_data : dict_data[ 'target_processing_score' ] = df_predicted_target_subsequence.processing_score.values
                dict_data[ f'query_affinity_{str_allele}' ] = df_predicted_query_subsequence.affinity.values
                dict_data[ f'target_affinity_{str_allele}' ] = df_predicted_target_subsequence.affinity.values
            df_subsequence_mhc_i = df_subsequence_mhc_i.join( pd.DataFrame( dict_data ) ) # combine predicted result

            ''' prepare a dataframe containing individual records '''
            float_thres_product_binding_affinities = float_thres_binding_affinities_in_nM ** 2 # calculate a threshold for the product of binding affinities for cross-reactivity prediction
             # prepare an array of indices
            df_subsequence_mhc_i.reset_index( drop = True, inplace = True ) # reset index before join operation
            arr_index = df_subsequence_mhc_i.index.values
            l = [ ]
            for mhc_i_allele in l_mhc_i_allele :
                for index, affinity_query, affinity_target in zip( arr_index, dict_data[ f'query_affinity_{mhc_i_allele}' ], dict_data[ f'target_affinity_{mhc_i_allele}' ] ) :
                    ''' perform filtering '''
                    if affinity_query * affinity_target <= float_thres_product_binding_affinities : # drop a pair of peptides if the geometric average of predicted binding affinity values (IC50) is above the threshold
                        l.append( [ index, mhc_i_allele, affinity_query, affinity_target ] )    
            df_subsequence_mhc_i.drop( columns = list( c for c in df_subsequence_mhc_i.columns.values if '_affinity_' in c ) + [ 'query_subsequence_without_gap', 'target_subsequence_without_gap' ], inplace = True ) # drop unnecessary columns
            df_mhc_web = pd.DataFrame( l, columns = [ 'index', 'mhc_allele', 'affinity_query', 'affinity_target' ] ).set_index( 'index' ).join( df_subsequence_mhc_i, how = 'left' ) 
            df_mhc_web[ 'mhc_class' ] = 1 # MHC-Class is I
            df_mhc_web[ 'score_for_sorting' ] = df_mhc_web.score_blosum / ( df_mhc_web.affinity_query * df_mhc_web.affinity_target ) ** 0.5 # calculate score for sorting
            df_mhc_web.to_csv( newfile, sep = '\t', index = False, header = False, mode = 'wb' )

            ''' complete the operation if the end of the file has been reached '''
            if flag_file_end_reached :
                break
            break

    newfile.close( )
