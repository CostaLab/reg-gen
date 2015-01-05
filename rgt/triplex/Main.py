# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
import argparse 
import time, datetime, getpass, fnmatch
# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import GenomicRegionSet
from triplexTools import TriplexSearch, FischerTest, RandomTest
from SequenceSet import Sequence, SequenceSet
from Util import SequenceType

dir = os.getcwd()

# To do: merge absolute path and relative path

"""
Statistical tests and plotting tools for triplex binding site analysis

Author: Joseph Kuo
"""

##########################################################################
##### UNIVERSAL FUNCTIONS ################################################
##########################################################################

def print2(summary, string):
    """ Show the message on the console and also save in summary. """
    print(string)
    summary.append(string)
    
def output_summary(summary, directory, filename):
    """Save the summary log file into the defined directory"""
    pd = os.path.join(dir,directory)
    try: os.stat(pd)
    except: os.mkdir(pd)    
    if summary:
        with open(os.path.join(pd,"parameters.txt"),'w') as f:
            print("#########  RGT Triplex: Summary information #########", file=f)
            for s in parameter:
                print(s, file=f)
        
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
    
    ################### Triplex search #######################################

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
    
    parser_search.add_argument('-min',type=int, default=8, help="Minimum length of binding site (Default: 4)")
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
        if not args.o: 
            print("Please define the output diractory name. \n")
            sys.exit(1)
            
        args = parser.parse_args()
        
        t0 = time.time()
        # Normalised output path
        args.o = os.path.normpath(os.path.join(dir,args.o))
        
        # Input parameters dictionary
        summary = []
        summary.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        summary.append("User: " + getpass.getuser())
        summary.append("\nCommand:\n   $ " + " ".join(sys.argv))

    ################################################################################
    ##### Search ###################################################################
    ################################################################################

    if args.mode == 'search':
        
        ############################################################################
        ##### Both RNA and DNA input ###############################################
        if args.r and args.d: 
            print2(summary, "\nSearch potential triplex forming binding sites between RNA and DNA")
            # Normalised paths
            args.r = os.path.normpath(os.path.join(dir,args.r)) 
            args.d = os.path.normpath(os.path.join(dir,args.d))
            triplex = TriplexSearch()
            
            ##### Read RNA sequences ###############################################
            print2(summary, "Step 1: Searching potential binding sites on RNA")
            rnaname = os.path.basename(args.r).split(".")[0]
            rnas = SequenceSet(name=rnaname, seq_type=SequenceType.RNA)
            
            if args.rt == 'fasta': # Input is FASTA
                print2(summary, "\tRead RNA in FASTA: "+args.r)
                rnas.read_fasta(args.r)
            else: # Input is BED
                if not args.genome: 
                    print("Please add the directory where the genome FASTA files locate.\n")
                    sys.exit(1)
                
                args.genome = os.path.normpath(os.path.join(dir,args.genome))    # Normalised paths
                print2(summary, "\tRead RNA in BED: "+args.r)
                print2(summary, "\tRefer to genome in FASTA: "+args.genome)
                rnas.read_bed(args.r, args.genome)
            
            ##### Search RNA potential binding sites ###############################
            # Save parameters
            print2(summary, "\tMotif: "+args.m)
            print2(summary, "\tMinimum length: "+args.min+" bp")
            if args.max: print2(summary, "\tMaximum length: "+args.max+" bp")
            else: print2(summary, "\tMaximum length: infinite")
            
            rbs = triplex.search_bindingsites(sequence_set=rnas, seq_type=SequenceType.RNA, 
                                              motif=args.m, min_len=args.min, max_len=args.max)
            rbs.write_bs(os.path.join(args.o, rnaname+".rbs"))
            t1 = time.time()
            print2(summary, "\tRunning time : " + str(datetime.timedelta(seconds=round(t1-t0))))
            print2(summary, "\tRNA binding sites are saved in: "+os.path.join(args.o, rnaname+".rbs"))
            
            
            ##### Read DNA sequences ###############################################
            print2(summary, "Step 2: Searching potential binding sites on DNA")
            dnaname = os.path.basename(args.d).split(".")[0]
            dnas = SequenceSet(name=dnaname, seq_type=SequenceType.DNA)
            
            if args.dt == 'fasta': # Input is FASTA
                print2(summary, "\tRead DNA in FASTA: "+args.d)
                dnas.read_fasta(args.d)
            else: # Input is BED
                if not args.genome: 
                    print("Please add the directory where the genome FASTA files locate.\n")
                    sys.exit(1)
                
                args.genome = os.path.normpath(os.path.join(dir,args.genome))    # Normalised paths
                print2(summary, "\tRead DNA in BED: "+args.d)
                print2(summary, "\tRefer to genome in FASTA: "+args.genome)
                dnas.read_bed(args.d, args.genome)
            
            ##### Search DNA potential binding sites ###############################
            # Save parameters
            print2(summary, "\tMinimum length: "+args.min+" bp")
            if args.max: print2(summary, "\tMaximum length: "+args.max+" bp")
            else: print2(summary, "\tMaximum length: infinite")
            
            dbs = triplex.search_bindingsites(sequence_set=dnas, seq_type=SequenceType.DNA, 
                                              motif=args.m, min_len=args.min, max_len=args.max)
            dbs.write_bs(os.path.join(args.o, dnaname+".dbs"))
            t2 = time.time()
            print2(summary, "\tRunning time : " + str(datetime.timedelta(seconds=round(t2-t1))))
            print2(summary, "\tDNA binding sites are saved in: "+os.path.join(args.o, dnaname+".dbs"))
            
            ##### Compare the binding sites between RNA and DNA ####################
            rbs = 
            dbs
            
            
            
            output_summary(summary, args.o, "summary.log")
            
            
            
            
            
            
            
            
        ############################################################################
        ##### Only RNA input #######################################################
        elif args.r and not args.d:
            print2(summary, "\nSearch potential triplex forming binding sites on RNA")
            
            args.r = os.path.normpath(os.path.join(dir,args.r))   # Normalised paths
            rnaname = os.path.basename(args.r).split(".")[0]
            rnas = SequenceSet(name=rnaname, seq_type=SequenceType.RNA)
            
            # Input is FASTA
            if args.rt == 'fasta':
                print2(summary, "\tRead RNA in FASTA: "+args.r)
                rnas.read_fasta(args.r)

            # Input is BED
            else:
                if not args.genome: 
                    print("Please add the directory where the genome FASTA files locate.\n")
                    sys.exit(1)
                
                args.genome = os.path.normpath(os.path.join(dir,args.genome))    # Normalised paths
                print2(summary, "\tRead RNA in BED: "+args.r)
                print2(summary, "\tRefer to genome in FASTA: "+args.genome)
                rnas.read_bed(args.r, args.genome)
            
            triplex = TriplexSearch()
            print2(summary, "\tMotif: "+args.m)
            print2(summary, "\tMinimum length: "+args.min)
            print2(summary, "\tMaximum length: "+args.max)
            bs = triplex.search_bindingsites(sequence_set=rnas, seq_type=SequenceType.RNA, 
                                             motif=args.m, min_len=args.min, max_len=args.max)

            bs.write_bs(os.path.join(args.o, rnaname+".rbs"))
            t1 = time.time()
            print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            print2(summary, "Results are saved in: "+os.path.join(args.o, rnaname+".rbs"))
            output_summary(summary, args.o, "summary.log")
         
        ############################################################################
        ##### Only DNA input #######################################################
        elif args.d and not args.r:
            print2(summary, "\nSearch potential triplex forming binding sites on DNA")
            
            args.d = os.path.normpath(os.path.join(dir,args.d))   # Normalised paths
            dnaname = os.path.basename(args.r).split(".")[0]
            dnas = SequenceSet(name=dnaname, seq_type=SequenceType.DNA)
            
            # Input is FASTA
            if args.dt == 'fasta':
                print2(summary, "\tRead DNA in FASTA: "+args.d)
                dnas.read_fasta(args.d)

            # Input is BED
            else:
                if not args.genome: 
                    print("Please add the directory where the genome FASTA files locate.\n")
                    sys.exit(1)
                args.genome = os.path.normpath(os.path.join(dir,args.genome))   # Normalised paths
                print2(summary, "\tRead DNA in BED: "+args.d)
                print2(summary, "\tRefer to genome in FASTA: "+args.genome)
                dnas.read_bed(os.path.join(dir, args.d), args.genome)
            
            triplex = TriplexSearch()
            print2(summary, "\tMotif: "+args.m)
            print2(summary, "\tMinimum length: "+args.min)
            print2(summary, "\tMaximum length: "+args.max)
            bs = triplex.search_bindingsites(sequence_set=dnas, seq_type=SequenceType.DNA, 
                                             motif=args.m, min_len=args.min, max_len=args.max)

            bs.write_bs(os.path.join(args.o, dnaname+".dbs"))
            t1 = time.time()
            print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            print2(summary, "Results are saved in: "+os.path.join(args.o, dnaname+".dbs"))
            output_summary(summary, args.o, "summary.log")
            
            
        # No input
        else:
            print("Please define either RNA strand or DNA strand (or both) as inputs\n")
        
    ################################################################################
    ##### Fischer ##################################################################
    ################################################################################
    if args.mode == 'fischer':
        #print(args.de)
        #print(args.organism)
        #print(args.pl)
        fischer = FischerTest(gene_list_file=args.de, organism=args.organism, promoterLength=args.pl)
        fischer.search_triplex(rna=args.r, temp=os.path.join(args.o,"fischer_test"))
    
    ################################################################################
    ##### Random ###################################################################
    ################################################################################
    if args.mode == 'random':
        pass
