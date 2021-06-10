from biobookshelf.main import *
from biobookshelf import *

pd.options.mode.chained_assignment = None  # default='warn' # to disable worining

def Predict_T_cell_cross_reactivity( dir_folder_pipeline, float_thres_avg_blosum62_score_for_mhc, float_thres_min_mhc_allele_frequency ) :
    """
    Predict_T_cell_cross_reactivity
    """
    
    """
    Package settings
    """
    name_package = 'cressp'
    dir_remote = 'https://github.com/ahs2202/cressp2/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package

    ''' read dict_blosum62 from the tsv file '''
    df_blosum62 = pd.read_csv( f'{dir_folder_cressp}data/blosum62.tsv.gz', sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score
    
    """ 
    Calculate T-cell epitope similarity
    """
    
    """ retrieve subsequences from aligned query and target protein sequences """
    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :
        # setting
        int_thres_min_length_peptide_mhc_i = 8
        int_thres_max_length_peptide_mhc_i = 15

        if 'df_matched' not in locals( ) :
            df_matched = pd.read_csv( f'{dir_folder_pipeline}matched.tsv.gz', sep = '\t' ) # read alignments between query and target protein sequences

        s_query_alignment_len = df_matched.query_alignment.apply( len )
        df_matched_for_mhc_i = df_matched[ ( s_query_alignment_len >= int_thres_min_length_peptide_mhc_i ) & ( s_query_alignment_len <= int_thres_max_length_peptide_mhc_i ) ]
        s_query_alignment_len = df_matched_for_mhc_i.query_alignment.apply( len )
        l_df = list( )
        for int_window_size in np.arange( int_thres_min_length_peptide_mhc_i, int_thres_max_length_peptide_mhc_i + 1 ) :
            df = df_matched_for_mhc_i[ s_query_alignment_len >= int_window_size ] # filter entries with length below window size
            l_l_value = list( )
            for query_accession, target_accession, query_start, query_end, target_start, target_end, query_alignment, target_alignment, e_value, source in df[ [ 'query_accession', 'target_accession', 'query_start', 'query_end', 'target_start', 'target_end', 'query_alignment', 'target_alignment', 'e_value', 'source' ] ].values :
                try : l_similarity_score = list( dict_blosum62[ ( amino_acid_query, amino_acid_target ) ] for amino_acid_query, amino_acid_target in zip( query_alignment.upper( ), target_alignment.upper( ) ) ) # calculate similarity score for each residue
                except : 
                    continue # in case of ambiguous amino acids (e.g. J = I or L, B = D or N), skip the protein itself
                n_windows = len( l_similarity_score ) + 1 - int_window_size
                for int_start in range( n_windows ) : 
                    int_similarity_score_for_window = sum( l_similarity_score[ int_start : int_start + int_window_size ] ) 
                    query_alignment_subsequence = query_alignment[ int_start : int_start + int_window_size ]
                    target_alignment_subsequence = target_alignment[ int_start : int_start + int_window_size ]
                    int_gap_count_subsequence_query, int_gap_count_subsequence_target = query_alignment_subsequence.count( '-' ), target_alignment_subsequence.count( '-' )
                    int_gap_count_before_subsequence_query, int_gap_count_before_subsequence_target = query_alignment[ : int_start ].count( '-' ), target_alignment[ : int_start ].count( '-' )
                    l_l_value.append( [ source, query_accession, target_accession, e_value, int_similarity_score_for_window, query_start + int_start - int_gap_count_before_subsequence_query, query_start + int_start + int_window_size - 1 - int_gap_count_subsequence_query - int_gap_count_before_subsequence_query, target_start + int_start - int_gap_count_before_subsequence_target, target_start + int_start + int_window_size - 1 - int_gap_count_subsequence_target - int_gap_count_before_subsequence_target, query_alignment_subsequence, target_alignment_subsequence ] )        
            if len( l_l_value ) == 0 : # if no records were retrieved, continue
                continue 
            df_subsequence = pd.DataFrame( l_l_value, columns = [ 'source', 'query_accession', 'target_accession', 'e_value', 'score_blosum', 'query_start', 'query_end', 'target_start', 'target_end', 'query_subsequence', 'target_subsequence' ] ).sort_values( 'e_value' )
            df_subsequence[ 'window_size' ] = int_window_size
            l_df.append( df_subsequence )
        df_subsequence_mhc_i = pd.concat( l_df )

        df_subsequence_mhc_i.to_csv( f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.tsv.gz', sep = '\t', index = False )

        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )

    """ Filter predicted cross-reactive T-cell epitope pairs """
    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.filtered.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :
        # read 'df_subsequence_mhc_i' if it has not been loaded yet.
        if 'df_subsequence_mhc_i' not in locals( ) :
            df_subsequence_mhc_i = pd.read_csv( f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.tsv.gz', sep = '\t' )
        # calculate average blosum62 score of the alignment between query and target peptides
        df_subsequence_mhc_i[ 'average_score_blosum' ] = df_subsequence_mhc_i.score_blosum / df_subsequence_mhc_i.window_size 
        df_subsequence_mhc_i = df_subsequence_mhc_i[ df_subsequence_mhc_i.average_score_blosum >= float_thres_avg_blosum62_score_for_mhc ] # perform filtering
        df_subsequence_mhc_i.to_csv( f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.filtered.tsv.gz', sep = '\t', index = False )

        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )


    """ 
    Calculate T-cell epitope similarity
    """

    """ check flag """
    dir_file_flag = f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.filtered.binding_affinity.tsv.gz.completed.flag'
    if not os.path.exists( dir_file_flag ) :

        if 'df_subsequence_mhc_i' not in locals( ) :
            df_subsequence_mhc_i = pd.read_csv( f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.filtered.tsv.gz', sep = '\t' )

        # retrieve list of common MHC allele names
        from mhcflurry import Class1PresentationPredictor # load MHCflurry predictor

        ''' load MHCFlurry binding affinity predictor '''
        predictor = Class1PresentationPredictor.load( ) # load MHCflurry predictor for benchmarking

        ''' retrieve MHC-I alleles with allele frequency larger then the given threshold in at least one population '''
        df_mhc_af = pd.read_csv( f"{dir_folder_cressp}data/mhc_population_allele_frequency.tsv.gz", sep = '\t', index_col = [ 0, 1 ] )
        l_mhc_i_allele = PD_Threshold( ( df_mhc_af.loc[ 'I' ] > float_thres_min_mhc_allele_frequency ).sum( axis = 1 ), a = 0 ).index.values # retrieve mhc_i_alleles with allele frequency > 'float_thres_min_mhc_allele_frequency' in at least one population
        l_mhc_i_allele = list( set( l_mhc_i_allele ).intersection( predictor.supported_alleles ) ) # retrive valid alleles for MHCflurry

        ''' remove subsequences containing invalid amino acids, and retrieve subsequences without gaps '''
        df_subsequence_mhc_i = PD_Search( df_subsequence_mhc_i, query_subsequence = [ 'B' ], target_subsequence = [ 'B' ], is_negative_query = True ) # retrive subset of df_subsequence for MHC_I allele binding prediction # invalid residues 'B' should be removed 
        df_subsequence_mhc_i[ 'query_subsequence_without_gap' ] = list( seq.replace( '-', '' ) for seq in df_subsequence_mhc_i.query_subsequence.values ) # remove gap in the aligment and use the sequence as inputs
        df_subsequence_mhc_i[ 'target_subsequence_without_gap' ] = list( seq.replace( '-', '' ) for seq in df_subsequence_mhc_i.target_subsequence.values )

        ''' predict binding affinity of peptides pairs with MHCflurry '''
        dict_data = dict( )
        for str_allele in l_mhc_i_allele :
            df_predicted_query_subsequence = predictor.predict( peptides = list( df_subsequence_mhc_i.query_subsequence_without_gap.values ), alleles = [ str_allele ], verbose = 0 ) 
            df_predicted_target_subsequence = predictor.predict( peptides = list( df_subsequence_mhc_i.target_subsequence_without_gap.values ), alleles = [ str_allele ], verbose = 0 ) # remove gap in the aligment and use the sequence as inputs
            if 'query_processing_score' not in dict_data : dict_data[ 'query_processing_score' ] = df_predicted_query_subsequence.processing_score.values # retreve processing score only once for query and target peptides
            if 'target_processing_score' not in dict_data : dict_data[ 'target_processing_score' ] = df_predicted_target_subsequence.processing_score.values
            dict_data[ f'query_affinity_{str_allele}' ] = df_predicted_query_subsequence.affinity.values
            dict_data[ f'target_affinity_{str_allele}' ] = df_predicted_target_subsequence.affinity.values
        df_subsequence_mhc_i = df_subsequence_mhc_i.join( pd.DataFrame( dict_data ) ) # combine predicted result
        df_subsequence_mhc_i.to_csv( f'{dir_folder_pipeline}t_cell.mhc_i.subsequence.filtered.binding_affinity.tsv.gz', sep = '\t', index = False )

        """ set flag """
        with open( dir_file_flag, 'w' ) as newfile :
            newfile.write( 'completed\n' )