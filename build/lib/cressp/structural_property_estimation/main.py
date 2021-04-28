from biobookshelf.main import *
from biobookshelf import *

def Estimate_structural_property( dir_file_protein, n_threads, dir_folder_output, dir_folder_pipeline, dir_folder_pipeline_temp ) :
    """
    # 2021-04-21 20:33:45 
    Estimate structural properties of given proteins, and write a tsv file containing structural properties of the proteins
    
    'dir_file_protein' : unzipped fasta file
    """
    
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
    dir_remote = 'https://github.com/ahs2202/cressp2/raw/main/cressp/' # remote directory from which datafiles will be downloaded
    dir_folder_cressp = f"{pkg_resources.resource_filename( name_package, '' )}/" # directory of the current installed package
    
    
    """
    RCSB_PDB 
    """
    # download and parse associated data (RCSB_PDB)
    PKG.Download_Data( "data/pdb/rcsb_pdb.tsv.gz", dir_remote, name_package ) # download data
    
    df_protein_rcsb_pdb = pd.read_csv( f'{dir_folder_cressp}data/pdb/rcsb_pdb.tsv.gz', sep = '\t' )
        
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
    
    