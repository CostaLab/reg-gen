# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import argparse 
import time, datetime, getpass, fnmatch
# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import GenomicRegionSet
from triplexTools import TriplexSearch, PromoterTest, RandomTest
from rgt.SequenceSet import Sequence, SequenceSet
from rgt.Util import SequenceType, Html

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
        with open(os.path.join(pd,filename),'w') as f:
            print("********* RGT Triplex: Summary information *********", file=f)
            for s in summary:
                print(s, file=f)
    
def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try: os.stat(path)
    except: os.mkdir(path)

def list_all_index(path):
    """Creat an 'index.html' in the defined directory """
    dirname = os.path.basename(path)
    link_d = {dirname:os.path.join(path,"index.html")}
    
    html = Html(name="Triplex Domain Finder", links_dict=link_d, 
                fig_dir=os.path.join(path,"style"), fig_rpath="./style", RGT_header=False)
    header_list = ["Experiments"]
    html.add_heading("All experiments in: "+dirname+"/")
    data_table = []
    type_list = 's'
    col_size_list = [10]
    for root, dirnames, filenames in os.walk(path):
        roots = root.split('/')
        for filename in fnmatch.filter(filenames, '*.html'):
            if filename == 'index.html' and root.split('/')[-1] != dirname:
                data_table.append(['<a href="'+os.path.join(root.split('/')[-1], filename)+'"><font size='+'"4"'+'>'+root.split('/')[-1]+"</a>"])
                #print(link_d[roots[-1]])
    html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=50, cell_align="left")
    html.write(os.path.join(path,"index.html"))

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
    parser_search.add_argument('-r', '-RNA', type=str, help="Input file name for RNA (in fasta format)")
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
    parser_search.add_argument('-mp', action="store_true", help="Perform multiprocessing for faster computation.")
    
    ################### Promoter test ##########################################

    h_promotor = "Promoter test evaluates the association between the given lncRNA to the target promoters."
    parser_promotertest = subparsers.add_parser('promotertest', help=h_promotor)
    parser_promotertest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_promotertest.add_argument('-rn', type=str, default=None, metavar='  ', help="Define the RNA name")
    parser_promotertest.add_argument('-de', default=False, metavar='  ', help="Input file for defferentially expression gene list (gene symbols or Ensembl ID)")
    parser_promotertest.add_argument('-bed', default=False, metavar='  ', help="Input BED file of the promoter regions of defferentially expression genes")
    parser_promotertest.add_argument('-bg', default=False, metavar='  ', help="Input BED file of the promoter regions of background genes")
    parser_promotertest.add_argument('-o', metavar='  ', help="Output directory name for all the results and temporary files")
    
    parser_promotertest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_promotertest.add_argument('-genome_path',type=str, metavar='  ', help='Define the path of genome FASTA file')

    parser_promotertest.add_argument('-pl', type=int, default=1000, metavar='  ', help="Define the promotor length (Default: 1000)")
    
    parser_promotertest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_promotertest.add_argument('-a', type=int, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (Default: 0.05)")
    parser_promotertest.add_argument('-ccf', type=int, default=20, metavar='  ', help="Define the cut off value for promoter counts (Default: 20)")
    parser_promotertest.add_argument('-rt', action="store_true", default=False, help="Remove temporary files (fa, txp...etc)")
    parser_promotertest.add_argument('-log', action="store_true", default=False, help="Set the plots in log scale")
    parser_promotertest.add_argument('-ac', type=str, default=False, metavar='  ', help="Input file for RNA accecibility ")
    parser_promotertest.add_argument('-accf', type=float, default=500, metavar='  ', help="Define the cut off value for RNA accecibility")
    parser_promotertest.add_argument('-obed', action="store_true", default=False, help="Output the BED files for DNA binding sites.")
    parser_promotertest.add_argument('-showpa', action="store_true", default=False, help="Show parallel and antiparallel bindings in the plot separately.")
    
    parser_promotertest.add_argument('-l', type=int, default=15, metavar='  ', help="[Triplexator] Define the minimum length of triplex (Default: 15)")
    parser_promotertest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (Default: 20)")
    parser_promotertest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (Default: 2)")
    parser_promotertest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (Default: off)")
    parser_promotertest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (Default 0).'0' = greedy approach; '1' = q-gram filtering.")
    parser_promotertest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (Default: 1)")
    parser_promotertest.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_promotertest.add_argument('-rm', action="store_true", default=False, help="[Triplexator] Set the multiprocessing")
    
    
    ################### Random test ##########################################
    h_random = "Region test evaluates the association between the given lncRNA to the target regions by randomization."
    parser_randomtest = subparsers.add_parser('regiontest', help=h_random)
    parser_randomtest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_randomtest.add_argument('-rn', type=str, default=False, metavar='  ', help="Define the RNA name")
    parser_randomtest.add_argument('-bed', metavar='  ', help="Input BED file for interested regions on DNA")
    parser_randomtest.add_argument('-o', metavar='  ', help="Output directory name for all the results and temporary files")
    
    parser_randomtest.add_argument('-n', type=int, default=10000, metavar='  ', 
                                   help="Number of times for randomization (Default: 10000)")

    parser_randomtest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_randomtest.add_argument('-genome_path',type=str, metavar='  ', help='Define the path of genome FASTA file.')
    
    parser_randomtest.add_argument('-a', type=int, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (Default: 0.05)")
    parser_randomtest.add_argument('-ccf', type=int, default=20, metavar='  ', help="Define the cut off value for DBS counts (Default: 20)")
    parser_randomtest.add_argument('-rt', action="store_true", default=False, help="Remove temporary files (fa, txp...etc)")
    parser_randomtest.add_argument('-log', action="store_true", default=False, help="Set the plots in log scale")
    parser_randomtest.add_argument('-f', type=str, default=False, metavar='  ', help="Input BED file as mask in randomization")
    parser_randomtest.add_argument('-ac', type=str, default=False, metavar='  ', help="Input file for RNA accecibility ")
    parser_randomtest.add_argument('-accf', type=float, default=500, metavar='  ', help="Define the cut off value for RNA accecibility")
    parser_randomtest.add_argument('-obed', action="store_true", default=False, help="Output the BED files for DNA binding sites.")
    parser_randomtest.add_argument('-showpa', action="store_true", default=False, help="Show parallel and antiparallel bindings in the plot separately.")
    
    parser_randomtest.add_argument('-l', type=int, default=15, metavar='  ', help="[Triplexator] Define the minimum length of triplex (Default: 15)")
    parser_randomtest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (Default: 20)")
    parser_randomtest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (Default: 2)")
    parser_randomtest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (Default: off)")
    parser_randomtest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (Default 0).'0' = greedy approach; '1' = q-gram filtering.")
    parser_randomtest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (Default: 1)")
    parser_randomtest.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_randomtest.add_argument('-rm', action="store_true", default=False, help="[Triplexator] Set the multiprocessing")
    

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
        if not args.o: 
            print("Please define the output diractory name. \n")
            sys.exit(1)
        if not args.organism: 
            print("Please define the organism. (hg19 or mm9)")
            sys.exit(1)
        if not args.rn: 
            print("Please define RNA sequence name.")
            sys.exit(1)
        if not args.genome_path: 
            print("Please define the path of genome FASTA file.")
            sys.exit(1)
        if not args.l: 
            print("Please define the minimum length and other parameters for Triplexator.")
            sys.exit(1)
        
        t0 = time.time()
        # Normalised output path
        args.o = os.path.normpath(os.path.join(dir,args.o))
        check_dir(args.o)
        # Input parameters dictionary
        summary = []
        summary.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        summary.append("User: " + getpass.getuser())
        summary.append("\nCommand:\n\t$ " + " ".join(sys.argv))

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
            output_summary(summary, args.o, "summary.txt")
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
            print2(summary, "\tMinimum length: "+str(args.min))
            print2(summary, "\tMaximum length: "+str(args.max))

            bs = triplex.search_bindingsites(sequence_set=rnas, seq_type=SequenceType.RNA, 
                                             motif=args.m, min_len=args.min, max_len=args.max, multiprocess=args.mp)

            bs.write_rbs(os.path.join(args.o, rnaname+".rbs"))
            t1 = time.time()
            print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            print2(summary, "Results are saved in: "+os.path.join(args.o, rnaname+".rbs"))
            output_summary(summary, args.o, "summary.txt")
         
        ############################################################################
        ##### Only DNA input #######################################################
        elif args.d and not args.r:
            print2(summary, "\nSearch potential triplex forming binding sites on DNA")
            
            args.d = os.path.normpath(os.path.join(dir,args.d))   # Normalised paths
            dnaname = os.path.basename(args.d).split(".")[0]
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
            print2(summary, "\tMinimum length: "+str(args.min))
            print2(summary, "\tMaximum length: "+str(args.max))
            bs = triplex.search_bindingsites(sequence_set=dnas, seq_type=SequenceType.DNA, 
                                             motif=args.m, min_len=args.min, max_len=args.max, multiprocess=args.mp)

            bs.write_dbs(os.path.join(args.o, dnaname+".dbs"))
            t1 = time.time()
            print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            print2(summary, "Results are saved in: "+os.path.join(args.o, dnaname+".dbs"))
            output_summary(summary, args.o, "summary.txt")
            
            
        # No input
        else:
            print("Please define either RNA strand or DNA strand (or both) as inputs\n")
        
    ################################################################################
    ##### Promoter Test ############################################################
    ################################################################################
    if args.mode == 'promotertest':
        print2(summary, "\n"+"*************** Promoter Test ****************")
        print2(summary, "*** Input RNA sequence: "+args.r)
        print2(summary, "*** Output directory: "+os.path.basename(args.o))

        args.r = os.path.normpath(os.path.join(dir,args.r))
        
        if args.de: args.de = os.path.normpath(os.path.join(dir,args.de))
        if args.bed: args.bed = os.path.normpath(os.path.join(dir,args.bed))
        if args.bg: args.bg = os.path.normpath(os.path.join(dir,args.bg))

        # Get GenomicRegionSet from the given genes
        print2(summary, "Step 1: Calculate the triplex forming sites on RNA and DNA.")
        promoter = PromoterTest(gene_list_file=args.de, rna_name=args.rn, bed=args.bed, bg=args.bg, organism=args.organism, 
                                promoterLength=args.pl, summary=summary, genome_path=args.genome_path, temp=dir,
                                showdbs=args.showdbs)
        promoter.search_triplex(rna=args.r, temp=args.o, l=args.l, e=args.e, remove_temp=args.rt,
                                c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf)
        t1 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t1-t0))))

        print2(summary, "Step 2: Calculate the frequency of DNA binding sites within the promotors.")
        if args.obed: obedp = os.path.basename(args.o)
        else: obedp = None
        promoter.count_frequency(temp=args.o, remove_temp=args.rt, obedp=obedp, cutoff=args.ccf)
        promoter.fisher_exact(alpha=args.a)
        t2 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t2-t1))))

        print2(summary, "Step 3: Establishing promoter profile.")
        promoter.promoter_profile()
        t3 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t3-t2))))

        print2(summary, "Step 4: Generate plot and output html files.")
        promoter.plot_lines(txp=promoter.txp_de, rna=args.r, dirp=args.o, ac=args.ac, 
                            cut_off=args.accf, log=args.log, showpa=args.showpa,
                            sig_region=promoter.sig_region_promoter,
                            ylabel="Number of target promoters with DBSs", 
                            linelabel="No. promoters", filename="plot_promoter.png")
        promoter.barplot(dirp=args.o, filename="bar_promoter.png", sig_region=promoter.sig_region_promoter)
        if args.showdbs:
            promoter.plot_lines(txp=promoter.txp_def, rna=args.r, dirp=args.o, ac=args.ac, 
                                cut_off=args.accf, log=args.log, showpa=args.showpa,
                                sig_region=promoter.sig_region_dbs,
                                ylabel="Number of DBSs on target promoters", 
                                linelabel="No. DBSs", filename="plot_dbss.png")
            promoter.barplot(dirp=args.o, filename="bar_dbss.png", sig_region=promoter.sig_region_dbs)
            
        promoter.gen_html(directory=args.o, parameters=args, ccf=args.ccf, align=50, alpha=args.a)
        promoter.gen_html_genes(directory=args.o, align=50, alpha=args.a, nonDE=False)
        t4 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t4-t3))))
        print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t4-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        list_all_index(path=os.path.dirname(args.o))
        
    ################################################################################
    ##### Random ###################################################################
    ################################################################################
    if args.mode == 'regiontest':
        print2(summary, "\n"+"*************** Region Test ***************")
        print2(summary, "*** Input RNA sequence: "+args.r)
        print2(summary, "*** Input regions in BED: "+os.path.basename(args.bed))
        print2(summary, "*** Number of randomization: "+str(args.n))
        print2(summary, "*** Output directoey: "+os.path.basename(args.o))
        
        args.r = os.path.normpath(os.path.join(dir,args.r))
        
        print2(summary, "\nStep 1: Calculate the triplex forming sites on RNA and the given regions")
        randomtest = RandomTest(rna_fasta=args.r, rna_name=args.rn, dna_region=args.bed, 
                                organism=args.organism, genome_path=args.genome_path)
        if args.obed: obed = os.path.basename(args.o)
        else: obed=False
        randomtest.target_dna(temp=args.o, remove_temp=args.rt, l=args.l, e=args.e, obed=obed,
                              c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, cutoff=args.ccf )
        t1 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
        

        print2(summary, "Step 2: Randomization and counting number of binding sites")
        randomtest.random_test(repeats=args.n, temp=args.o, remove_temp=args.rt, l=args.l, e=args.e,
                               c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, rm=args.rm,
                               filter_bed=args.f, alpha=args.a)
        t2 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t2-t1))))
        
        print2(summary, "Step 3: Generating plot and output HTML")
        randomtest.lineplot(txp=randomtest.txp, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
                            log=args.log, ylabel="Number of target regions with DBS", 
                            sig_region=randomtest.data["region"]["sig_region"],
                            linelabel="No. target regions", filename="lineplot_region.png")
        randomtest.lineplot(txp=randomtest.txpf, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
                            log=args.log, ylabel="Number of DBS on target regions",
                            sig_region=randomtest.data["dbs"]["sig_region"], 
                            linelabel="No. DBS", filename="lineplot_dbs.png")

        randomtest.boxplot(dir=args.o, matrix=randomtest.region_matrix, 
                           sig_region=randomtest.data["region"]["sig_region"], 
                           truecounts=[r[0] for r in randomtest.counts_tr.values()],
                           sig_boolean=randomtest.data["region"]["sig_boolean"], 
                           ylabel="Number of target regions with DBS",
                           filename="boxplot_regions" )
        
        randomtest.boxplot(dir=args.o, matrix=randomtest.dbss_matrix, 
                           sig_region=randomtest.data["dbs"]["sig_region"], 
                           truecounts=randomtest.counts_dbs.values(),
                           sig_boolean=randomtest.data["dbs"]["sig_boolean"], 
                           ylabel="Number of DBS on target regions",
                           filename="boxplot_dbs" )

        randomtest.gen_html(directory=args.o, parameters=args, align=50, alpha=args.a)
        t3 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t3-t2))))
        
        print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t3-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        list_all_index(path=os.path.dirname(args.o))
        
