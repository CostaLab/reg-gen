# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import argparse 
import time, datetime, getpass, fnmatch
import subprocess
# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import GenomicRegionSet
from triplexTools import PromoterTest, RandomTest
from rgt.SequenceSet import Sequence, SequenceSet
from rgt.Util import SequenceType, Html, ConfigurationFile

dir = os.getcwd()

"""
Triplex Domain Finder (TDF) provides statistical tests and plotting tools for triplex binding site analysis

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
    
    link_d = {"List":"index.html"}
    html = Html(name="Directory: "+dirname, links_dict=link_d, 
                fig_dir=os.path.join(path,"style"), fig_rpath="./style", RGT_header=False, other_logo="TDF")
    
    html.add_heading("All experiments in: "+dirname+"/")
    data_table = []
    type_list = 'sssss'
    col_size_list = [10, 10, 10, 10, 10]
    c = 0
    for root, dirnames, filenames in os.walk(path):
        #roots = root.split('/')
        for filename in fnmatch.filter(filenames, '*.html'):
            if filename == 'index.html' and root.split('/')[-1] != dirname:
                c += 1
                if "_" in root.split('/')[-1]:
                    tags = root.split('/')[-1].split("_")
                    p1 = tags[0]
                    p2 = tags[-1]
                    data_table.append([str(c), '<a href="'+os.path.join(root.split('/')[-1], filename)+'">'+root.split('/')[-1]+"</a>",
                                       p1, p2])
                    header_list = ["No.", "Experiments", "Tag1", "Tag2"]
                else:
                    data_table.append([str(c), '<a href="'+os.path.join(root.split('/')[-1], filename)+'">'+root.split('/')[-1]+"</a>"])
                    header_list = ["No.", "Experiments"]
                #print(link_d[roots[-1]])
    html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=50, cell_align="left", sortable=True)
    html.add_fixed_rank_sortable()
    html.write(os.path.join(path,"index.html"))

def main():
    ##########################################################################
    ##### PARAMETERS #########################################################
    ##########################################################################
    
    parser = argparse.ArgumentParser(description='Triplex Domain Finder is a statistical framework \
                                                  for detection of triple helix potential of \
                                                  lncRNAs from genome-wide functional data. \
                                                  Author: Chao-Chung Kuo\
                                                  \nVersion: 0.1.1', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
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

    parser_promotertest.add_argument('-pl', type=int, default=1000, metavar='  ', help="Define the promotor length (Default: 1000)")
    
    parser_promotertest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_promotertest.add_argument('-score', action="store_true", help="Load score column from input gene list of BED file for analysis.")
    parser_promotertest.add_argument('-scoreh', action="store_true", help="Use the header of scores from the given gene list or BED file.")
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
    
    
    ################### Genomic Region Test ##########################################
    h_region = "Genomic region test evaluates the association between the given lncRNA to the target regions by randomization."
    parser_randomtest = subparsers.add_parser('regiontest', help=h_region)
    parser_randomtest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_randomtest.add_argument('-rn', type=str, default=False, metavar='  ', help="Define the RNA name")
    parser_randomtest.add_argument('-bed', metavar='  ', help="Input BED file for interested regions on DNA")
    parser_randomtest.add_argument('-o', metavar='  ', help="Output directory name for all the results and temporary files")
    
    parser_randomtest.add_argument('-n', type=int, default=10000, metavar='  ', 
                                   help="Number of times for randomization (Default: 10000)")

    parser_randomtest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
 
    parser_randomtest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_randomtest.add_argument('-score', action="store_true", help="Load score column from input BED file")
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
    
    ##########################################################################
    parser_triplexator = subparsers.add_parser('triplexator', help="Setting Triplexator.")
    parser_triplexator.add_argument('-path',type=str, metavar='  ', help='Define the path of Triplexator.')
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
        cf = ConfigurationFile()
        
        process = subprocess.Popen(["triplexator", "--help"], shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait for the process to terminate
        out, err = process.communicate()
        errcode = process.returncode
        if errcode == 0:
            print("** Triplexator: OK")
        else:
            print("** Error: Triplexator cannot be found. Please export Triplexator path into $PATH")
            print(err)
            sys.exit(1)

        if not args.o: 
            print("Please define the output diractory name. \n")
            sys.exit(1)
        if not args.organism: 
            print("Please define the organism. (hg19 or mm9)")
            sys.exit(1)
        if not args.rn: 
            print("Please define RNA sequence name.")
            sys.exit(1)

        t0 = time.time()
        # Normalised output path
        args.o = os.path.normpath(os.path.join(dir,args.o))
        check_dir(os.path.dirname(os.path.dirname(args.o)))
        check_dir(os.path.dirname(args.o))
        check_dir(args.o)
        # Input parameters dictionary
        summary = []
        summary.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        summary.append("User: " + getpass.getuser())
        summary.append("\nCommand:\n\t$ " + " ".join(sys.argv))
           
    ################################################################################
    ##### Promoter Test ############################################################
    ################################################################################
    if args.mode == 'promotertest':
        if args.bed and not args.bg:
            print("Please add background promoters in BED format. (-bg)")
            sys.exit(1)
        if args.scoreh and not args.score:
            print("Score header (-scoreh) can only be used when scores (-score) are loaded.")
            print("Please add '-score'.")
            sys.exit(1)

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
                                promoterLength=args.pl, summary=summary, temp=dir,
                                showdbs=args.showdbs, score=args.score, scoreh=args.scoreh)
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
            promoter.barplot(dirp=args.o, filename="bar_dbss.png", sig_region=promoter.sig_region_dbs, dbs=True)
            
        promoter.gen_html(directory=args.o, parameters=args, ccf=args.ccf, align=50, alpha=args.a)
        promoter.gen_html_genes(directory=args.o, align=50, alpha=args.a, nonDE=False)
        t4 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t4-t3))))
        print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t4-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        list_all_index(path=os.path.dirname(args.o))
        
    ################################################################################
    ##### Genomic Region Test ######################################################
    ################################################################################
    if args.mode == 'regiontest':
        print2(summary, "\n"+"*************** Genomic Region Test ***************")
        print2(summary, "*** Input RNA sequence: "+args.r)
        print2(summary, "*** Input regions in BED: "+os.path.basename(args.bed))
        print2(summary, "*** Number of randomization: "+str(args.n))
        print2(summary, "*** Output directoey: "+os.path.basename(args.o))
        
        args.r = os.path.normpath(os.path.join(dir,args.r))
        
        print2(summary, "\nStep 1: Calculate the triplex forming sites on RNA and the given regions")
        randomtest = RandomTest(rna_fasta=args.r, rna_name=args.rn, dna_region=args.bed, 
                                organism=args.organism, showdbs=args.showdbs)
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
        
        randomtest.boxplot(dir=args.o, matrix=randomtest.region_matrix, 
                           sig_region=randomtest.data["region"]["sig_region"], 
                           truecounts=[r[0] for r in randomtest.counts_tr.values()],
                           sig_boolean=randomtest.data["region"]["sig_boolean"], 
                           ylabel="Number of target regions with DBS",
                           filename="boxplot_regions" )
        if args.showdbs:
            randomtest.lineplot(txp=randomtest.txpf, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
                                log=args.log, ylabel="Number of DBS on target regions",
                                sig_region=randomtest.data["dbs"]["sig_region"], 
                                linelabel="No. DBS", filename="lineplot_dbs.png")
            
            randomtest.boxplot(dir=args.o, matrix=randomtest.dbss_matrix, 
                               sig_region=randomtest.data["dbs"]["sig_region"], 
                               truecounts=randomtest.counts_dbs.values(),
                               sig_boolean=randomtest.data["dbs"]["sig_boolean"], 
                               ylabel="Number of DBS on target regions",
                               filename="boxplot_dbs" )

        randomtest.gen_html(directory=args.o, parameters=args, align=50, alpha=args.a, 
                            score=args.score)
        t3 = time.time()
        print2(summary, "\tRunning time is : " + str(datetime.timedelta(seconds=round(t3-t2))))
        
        print2(summary, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t3-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        list_all_index(path=os.path.dirname(args.o))
        
