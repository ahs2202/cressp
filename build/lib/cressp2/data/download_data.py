from biobookshelf.main import *
from biobookshelf import *

import pkg_resources

def Download_Data( ) :
    # download data files from Git-LFS
    dir_file_download_flag = pkg_resources.resource_filename( "cressp2", 'data/pdb/rcsb_pdb.tsv.gz.download_completed.flag' )
    if not os.path.exists( dir_file_download_flag ) :
        OS_Download( 'https://github.com/ahs2202/cressp2/raw/main/cressp2/data/pdb/rcsb_pdb.tsv.gz', pkg_resources.resource_filename( "cressp2", 'data/pdb/rcsb_pdb.tsv.gz' ) )
        with open( dir_file_download_flag, 'w' ) as file :
            file.write( 'download completed at ' + TIME_GET_timestamp( ) )
    else :
        print( "data already exists" )