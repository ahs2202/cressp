#!/usr/bin/env python
# 20210419

from biobookshelf.main import *
from biobookshelf import *

from cressp2.data import Download_Data

import argparse
import os, sys, getopt
from io import StringIO
import time
import math


def main( ) :
    parser = argparse.ArgumentParser( description = "a program to find cross-reactive epitopes with structural information from known protein structures." )
#     parser.add_argument( "-r", "--dir_file_fasta_ref", help = "(Required) directory of a fasta file containing the reference sequence. (e.g. a fasta sequence from AddGene)" )
#     parser.add_argument( "-i", "--dir_file_fastq", help = "(Required) directory of a fastq file from a nanopore sequencing" )
#     parser.add_argument( "-o", "--dir_folder_output", help = "(Default: subdirectory of the folder containing the given fastq file) directory of output folder", default = 'default' )
#     parser.add_argument( "-t", "--threads", help = "(Default: 10) Number of threads to use in the current compute node.", default = '10' )

    args = parser.parse_args( )
#     if args.dir_file_fasta_ref is None or args.dir_file_fastq is None  :
#         print( "required arguments are not given, exiting" )
#         sys.exit( )

    # [input] parse arguments
#     dir_file_fasta_ref = Program__Get_Absolute_Path_of_a_File( args.dir_file_fasta_ref )
#     dir_file_fastq = Program__Get_Absolute_Path_of_a_File( args.dir_file_fastq )
#     dir_folder_output = args.dir_folder_output
#     n_threads = int( args.threads )

    # read dict_blosum62 from the tsv file
    df_blosum62 = pd.read_csv( pkg_resources.resource_filename( "cressp2", 'data/blosum62.tsv.gz' ), sep = '\t' )
    dict_blosum62 = dict( )
    for aa_0, aa_1, score in df_blosum62.values : # sould be in [ 'aa_0', 'aa_1', 'BLOSUM62_score' ] order
        dict_blosum62[ aa_0, aa_1 ] = score
        
    Download_Data( ) # download data
    
    

    
    
if __name__ == "__main__" :
    main( )
    