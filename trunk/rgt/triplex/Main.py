# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#sys.path.append(lib_path)

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import GenomeData, OverlapType, Html

dir = os.getcwd()
"""
Statistical tests and plotting tools for triplex binding site analysis

Author: Joseph Kuo

"""

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################
    

######### Universal functions



def main():
    #################################################################################################
    ##### PARAMETERS ################################################################################
    #################################################################################################
    
    
    
    parser = argparse.ArgumentParser(description='Provides various Statistical tests for triplex \
                                     binding site analysis \
                                     \nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Random test ##########################################
    parser_randomtest = subparsers.add_parser('randomtest', help='Test the TTS are due to chance \
                                              or not by randomization')
    parser_randomtest.add_argument('-i', help="Input file name (.txp) ")
    parser_randomtest.add_argument('-o', help="Output directory name")
    parser_randomtest.add_argument('-n', type=int, default=10000, 
                                   help="Number of times for randomization (Default: 10000)")
    parser_randomtest.add_argument('-bg', help="Define a BED file as background. If not defined, the background is whole genome according to the given organism.")
    parser_randomtest.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_randomtest.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
