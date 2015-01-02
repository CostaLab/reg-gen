# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
import argparse 

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import GenomicRegionSet
from triplexTools import TriplexSearch, FischerTest, RandomTest
from SequenceSet import Sequence, SequenceSet
from Util import SequenceType


dir = os.getcwd()
"""
Statistical tests and plotting tools for triplex binding site analysis

Author: Joseph Kuo

"""

def main():
    ##########################################################################
    ##### PARAMETERS #########################################################
    ##########################################################################
    
    
    
    parser = argparse.ArgumentParser(description='Provides \
                                     triplex binding sites searching tool and \
                                     various Statistical tests for analysis. \
                                     \nAuthor: Joseph Kuo', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Triplex search ##########################################

    parser_search = subparsers.add_parser('search', help='Search the possible triplex binding sites \
                                                          between single strand (RNA) and \
                                                          double strand (DNA)')
    parser_search.add_argument('-r', '-RNA', type=str, help="Input file name for RNA (in fasta or bed format)")
    parser_search.add_argument('-d', '-DNA', type=str, help="Input file name for DNA (in fasta or bed format)")
    
    parser_search.add_argument('-rt', choices= ['fasta', 'bed'], default='fasta', 
                               help="Input file type (fasta or bed)")
    parser_search.add_argument('-dt', choices= ['fasta', 'bed'], default='fasta', 
                               help="Input file type (fasta or bed)")

    parser_search.add_argument('-o', type=str, help="Output directory name")
    parser_search.add_argument('-genome',type=str, help='Define the directory where the genome FASTA files locate.')
    
    parser_search.add_argument('-min',type=int, default=4, help="Minimum length of binding site (Default: 4)")
    parser_search.add_argument('-max',type=int, default=None, help="Maxmum length of binding site (Default is infinite)")
    parser_search.add_argument('-m',type=str, default="RYMPA", help="Define the motif for binding site searching (Default is RYMPA)")

    ################### Fischer exact test ##########################################
    parser_fischertest = subparsers.add_parser('fischer', help='Test the TTS are due to chance \
                                              or not by randomization')
    parser_fischertest.add_argument('-r', '-RNA', type=str, help="Input file name for RNA (in fasta or bed format)")
    parser_fischertest.add_argument('-de', help="Input file for defferentially expression gene list ")
    parser_fischertest.add_argument('-pl', type=int, default=0, 
                                   help="Define the promotor length (Default: 0)")
    parser_fischertest.add_argument('-o', help="Output directory name for all the results and temporary files")
    parser_fischertest.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')

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
        if not args.o: 
            print("Please define the output diractory name. \n")
            sys.exit(1)
        
        # Both RNA and DNA input
        if args.r and args.d: 
            pass


        # Only RNA input
        elif args.r and not args.d:
            rnaname = os.path.basename(args.r).split(".")[0]
            print("Read files")
            # Input is FASTA
            if args.rt == 'fasta':
                rnas = SequenceSet(name=rnaname, seq_type=SequenceType.RNA)
                rnas.read_fasta(args.r)

            # Input is BED
            else:
                if not args.genome: 
                    print("Please add the directory where the genome FASTA files locate.\n")
                    sys.exit(1)

                rnas = SequenceSet(name=rnaname, seq_type=SequenceType.RNA)
                rnas.read_bed(os.path.join(dir, args.r), args.genome)
            
            print("Search files")
            triplex = TriplexSearch()
            bs = triplex.search_bindingsites(sequence_set=rnas, seq_type=SequenceType.RNA, 
                                             motif=args.m, min_len=args.min, max_len=args.max)

            print("Write files")
            bs.write_rbs(os.path.join(args.o, rnaname+".rbs"))


        # Only DNA input
        elif args.d and not args.r:
            pass
        # No input
        else:
            print("Please define either RNA strand or DNA strand (or both) as input\n")
        
            """
            if args.d:
                os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -ss "+args.r+" -ds "+args.d+" > "+args.o)
            else:
                os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -ss "+args.r+" > "+args.o)
            """

    if args.mode == 'fischer':
        #print(args.de)
        #print(args.organism)
        #print(args.pl)
        fischer = FischerTest(gene_list_file=args.de, organism=args.organism, promoterLength=args.pl)
        fischer.search_triplex(rna=args.r, temp=os.path.join(args.o,"fischer_test"))
    #if args.mode == 'search':
