# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
import argparse 
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
    
    ################### Triplex search ##########################################

    parser_search = subparsers.add_parser('search', help='Search the possible triplex loci between \
                                                          single strand (RNA) and double strand (DNA)')
    parser_search.add_argument('-rna', help="Input file name for RNA (in fasta or bed format)")
    parser_search.add_argument('-dna', help="Input file name for DNA (in fasta or bed format)")
    
    parser_search.add_argument('-it', choices= ['fasta', 'bed'], default='fasta', help="Input file type (fasta or bed)")
    parser_search.add_argument('-o', help="Output directory name")
    parser_search.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')
    
    parser_search.add_argument('-min',type=int, default=4, help="Minimum length of TFO (Default: 4)")
    parser_search.add_argument('-max',type=int, default=0, help="Maxmum length of TFO (Default is infinite)")
    

    
    parser_search.add_argument('-bg', help="Define a BED file as background. If not defined, \
                                                the background is whole genome according to the given organism.")

    ################### Random test ##########################################
    parser_randomtest = subparsers.add_parser('randomtest', help='Test the TTS are due to chance \
                                              or not by randomization')
    parser_randomtest.add_argument('-i', help="Input file name (.txp) ")
    parser_randomtest.add_argument('-o', help="Output directory name")
    parser_randomtest.add_argument('-ns', type=int, default=1000, 
                                   help="Number of sequences for each randomization (Default: 1000)")
    parser_randomtest.add_argument('-l', type=int, default=1000, 
                                   help="Length of random sequences (Default: 1000 bp)")
    parser_randomtest.add_argument('-n', type=int, default=10000, 
                                   help="Number of times for randomization (Default: 10000)")
    parser_randomtest.add_argument('-bg', help="Define a BED file as background. If not defined, \
                                                the background is whole genome according to the given organism.")

    

    ################### Parsing the arguments ################################
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2: 
        # retrieve subparsers from parser
        subparsers_actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
        # there will probably only be one subparser_action,but better save than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                if choice == sys.argv[1]:
                    print("\nYou need more arguments.")
                    print("\nSubparser '{}'".format(choice))        
                    subparser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()
        
    #################################################################################################
    ##### Main #####################################################################################
    #################################################################################################

    if args.mode == 'search':
        if args.rna:
            if args.dna:
                os.system("triplexator -ss "+args.rna+" -ds "+args.dna+" > "+args.o)
            else:
                os.system("triplexator -ss "+args.rna+" > "+args.o)
        

    if args.mode == 'randomtest':
        pass
    #if args.mode == 'search':
