# Python Libraries
from __future__ import print_function
from __future__ import division
import os
import sys
import time
import getpass
import argparse
import datetime
import subprocess
import matplotlib
matplotlib.use('Agg', warn=False)
from collections import OrderedDict

# Local Libraries
# Distal Libraries
from .. import __version__
# from rgt.Util import Html
from triplexTools import get_dbss, check_dir,generate_rna_exp_pv_table, revise_index, \
                         no_binding_response, integrate_stat, update_profile, integrate_html, \
                         merge_DBD_regions, silentremove, summerize_stat, shorten_dir, merge_DBSs

# from tdf_promotertest import PromoterTest
# from tdf_regiontest import RandomTest
from ..tdf.Input import Input
from ..tdf.Triplexes import Triplexes
from ..tdf.Statistics import Statistics
from ..tdf.Report import Report


dir = os.getcwd()

"""
Triplex Domain Finder (TDF) provides statistical tests and plotting tools for 
triplex binding site analysis

Author: Joseph C.C. Kuo

Triplexator
https://github.com/zbarni/triplexator
Author: Barna Zajzon
"""

def main():
    ##########################################################################
    ##### PARAMETERS #########################################################
    ##########################################################################
    version_message = "TDF - Regulatory Analysis Toolbox (RGT). Version: " + str(__version__)
    parser = argparse.ArgumentParser(description='Triplex Domain Finder offers a statistical framework \
                                                  for detection of triple helix potential of \
                                                  lncRNAs from genome-wide functional data. \
                                                  Author: Chao-Chung Kuo',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter, version=version_message)
    
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Promoter test ##########################################

    h_promotor = "Promoter test evaluates the association between the given lncRNA to the target promoters."
    parser_promotertest = subparsers.add_parser('promotertest', help=h_promotor)
    parser_promotertest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    # parser_promotertest.add_argument('-rl', type=str, default=None, metavar='  ', help="Input list for name and file path to all RNA sequences (in fasta format)")
    parser_promotertest.add_argument('-rn', type=str, default=None, metavar='  ', help="Define the RNA name")
    parser_promotertest.add_argument('-de', default=False, metavar='  ', help="Input file for target gene list (gene symbols or Ensembl ID)")
    parser_promotertest.add_argument('-bed', default=False, metavar='  ', help="Input BED file of the promoter regions of target genes")
    parser_promotertest.add_argument('-bg', default=False, metavar='  ', help="Input BED file of the promoter regions of background genes")
    parser_promotertest.add_argument('-o', metavar='  ', help="Output directory name for all the results")
    parser_promotertest.add_argument('-t', metavar='  ', default=False, help="Define the title name for the results under the Output name. (default: %(default)s)")
    
    parser_promotertest.add_argument('-organism', metavar='  ', help='Define the organism')
    parser_promotertest.add_argument('-gtf', metavar='  ', default=None, help='Define the GTF file for annotation (optional)')
    parser_promotertest.add_argument('-tss', type=int, default=0, metavar='  ', help="Define the distance between the promoter regions and TSS along gene body (default: %(default)s)")
    parser_promotertest.add_argument('-pl', type=int, default=1000, metavar='  ', help="Define the promotor length (default: %(default)s)")
    
    parser_promotertest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_promotertest.add_argument('-score', action="store_true", help="Load score column from input gene list or BED file for analysis.")
    parser_promotertest.add_argument('-scoreh', action="store_true", help="Use the header of scores from the given gene list or BED file.")
    parser_promotertest.add_argument('-a', type=float, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (default: %(default)s)")
    parser_promotertest.add_argument('-ccf', type=int, default=100, metavar='  ', help="Define the cut off value for promoter counts (default: %(default)s)")
    parser_promotertest.add_argument('-rt', action="store_true", default=False, help="Remove temporary files (fa, txp...etc)")
    parser_promotertest.add_argument('-log', action="store_true", default=False, help="Set the plots in log scale")
    parser_promotertest.add_argument('-ac', type=str, default=False, metavar='  ', help="Input file for RNA accecibility ")
    parser_promotertest.add_argument('-accf', type=float, default=500, metavar='  ', help="Define the cut off value for RNA accecibility")
    parser_promotertest.add_argument('-obed', action="store_true", default=True, help="Output the BED files for DNA binding sites.")
    parser_promotertest.add_argument('-showpa', action="store_true", default=False, help="Show parallel and antiparallel bindings in the plot separately.")
    # parser_promotertest.add_argument('-motif', action="store_true", default=False, help="Show motif of binding sites.")
    parser_promotertest.add_argument('-filter_havana', type=str, default="F", metavar='  ', help="Apply filtering to remove HAVANA entries.")
    parser_promotertest.add_argument('-protein_coding', type=str, default="F", metavar='  ', help="Apply filtering to get only protein coding genes.")
    parser_promotertest.add_argument('-known_only', type=str, default="F", metavar='  ', help="Apply filtering to get only known genes.")
    parser_promotertest.add_argument('-dump', action="store_true", default=False, help="Only dump the experimental file and leave the program.")
    parser_promotertest.add_argument('-rnaexp', type=str, default=None, metavar='  ',
                                     help="Given a file with RNA name and the expression value")

    parser_promotertest.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexes] Define the minimum length of triplex (default: %(default)s)")
    parser_promotertest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexes] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_promotertest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexes] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_promotertest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexes] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_promotertest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexes] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_promotertest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexes] Define output formats of Triplexator (default: %(default)s)")
    parser_promotertest.add_argument('-mf', action="store_true", default=False, help="[Triplexes] Merge overlapping features into a cluster and report the spanning region.")
    parser_promotertest.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexes] Set the multiprocessing")
    parser_promotertest.add_argument('-par', type=str, default="", metavar='  ', help="[Triplexes] Define other parameters for Triplexator")
    
    ################### Genomic Region Test ##########################################
    h_region = "Genomic region test evaluates the association between the given lncRNA to the target regions by randomization."
    parser_randomtest = subparsers.add_parser('regiontest', help=h_region)
    parser_randomtest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_randomtest.add_argument('-rl', type=str, default=None, metavar='  ', help="Input list for paths to all RNA sequences (in fasta format)")
    parser_randomtest.add_argument('-rn', type=str, default=False, metavar='  ', help="Define the RNA name")
    parser_randomtest.add_argument('-bed', metavar='  ', help="Input BED file for interested regions on DNA")
    parser_randomtest.add_argument('-o', metavar='  ', help="Output directory name for all the results and temporary files")
    parser_randomtest.add_argument('-t', metavar='  ', default=False, help="Define the title name for the results under the Output name. (default: %(default)s)")
    
    parser_randomtest.add_argument('-n', type=int, default=10000, metavar='  ', 
                                   help="Number of times for randomization (default: %(default)s)")

    parser_randomtest.add_argument('-organism', metavar='  ', help='Define the organism')
 
    parser_randomtest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_randomtest.add_argument('-score', action="store_true", help="Load score column from input BED file")
    parser_randomtest.add_argument('-a', type=float, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (default: %(default)s)")
    parser_randomtest.add_argument('-ccf', type=int, default=40, metavar='  ', help="Define the cut off value for DBS counts (default: %(default)s)")
    parser_randomtest.add_argument('-rt', action="store_true", default=False, help="Remove temporary files (fa, txp...etc)")
    parser_randomtest.add_argument('-log', action="store_true", default=False, help="Set the plots in log scale")
    parser_randomtest.add_argument('-f', type=str, default=False, metavar='  ', help="Input BED file as mask in randomization")
    parser_randomtest.add_argument('-ac', type=str, default=False, metavar='  ', help="Input file for RNA accecibility ")
    parser_randomtest.add_argument('-accf', type=float, default=500, metavar='  ', help="Define the cut off value for RNA accecibility")
    parser_randomtest.add_argument('-obed', action="store_true", default=True, help="Output the BED files for DNA binding sites.")
    parser_randomtest.add_argument('-showpa', action="store_true", default=False, help="Show parallel and antiparallel bindings in the plot separately.")
    parser_randomtest.add_argument('-mp', type=int, default=1, metavar='  ', help="Define the number of threads for multiprocessing")
    
    parser_randomtest.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexes] Define the minimum length of triplex (default: %(default)s)")
    parser_randomtest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexes] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_randomtest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexes] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_randomtest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexes] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_randomtest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexes] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_randomtest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexes] Define output formats of Triplexator (default: %(default)s)")
    parser_randomtest.add_argument('-mf', action="store_true", default=False, help="[Triplexes] Merge overlapping features into a cluster and report the spanning region.")
    parser_randomtest.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexes] Set the multiprocessing")
    parser_randomtest.add_argument('-par', type=str, default="", metavar='  ', help="[Triplexes] Define other parameters for Triplexator")

    ##########################################################################
    parser_bed2bed = subparsers.add_parser('get_dbss', help="Get DBSs in BED format from the single BED file")
    parser_bed2bed.add_argument('-i',type=str, metavar='  ', help='Input BED file of the target regions')
    parser_bed2bed.add_argument('-dbs',type=str, metavar='  ', help='Output BED file of the DBSs')
    parser_bed2bed.add_argument('-rbs',type=str, metavar='  ', help='Output BED file of the RBSs')
    parser_bed2bed.add_argument('-r',type=str, metavar='  ', help='Input FASTA file of the RNA')
    parser_bed2bed.add_argument('-organism', metavar='  ', help='Define the organism')
    parser_bed2bed.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexes] Define the minimum length of triplex (default: %(default)s)")
    parser_bed2bed.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexes] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_bed2bed.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexes] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_bed2bed.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexes] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_bed2bed.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexes] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_bed2bed.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexes] Define output formats of Triplexator (default: %(default)s)")
    parser_bed2bed.add_argument('-mf', action="store_true", default=False, help="[Triplexes] Merge overlapping features into a cluster and report the spanning region.")
    parser_bed2bed.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexes] Set the multiprocessing")
    
    ##########################################################################
    # rgt-TDF integrate -path 
    parser_integrate = subparsers.add_parser('integrate', help="Integrate the project's links and generate project-level statistics.")
    parser_integrate.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
    parser_integrate.add_argument('-exp', action="store_true", default=False, help='Include expression score for ranking.')
    ##########################################################################
    parser_updatehtml = subparsers.add_parser('updatehtml', help="Update the project's html.")
    parser_updatehtml.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
    parser_updatehtml.add_argument('-exp', type=str, metavar='  ', help='Define file with expression data.')

    ################### Parsing the arguments ################################
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "-v" or sys.argv[1] == "--version":
            print(version_message)
            sys.exit(0)
        else:
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

        ####################################################################################
        ######### Integration
        if args.mode == "integrate":

            targets = []
            for dirpath, dnames, fnames in os.walk(args.path):
                for f in fnames:
                    if f == "stat.txt":
                        targets.append(os.path.dirname(dirpath))
            targets = list(set(targets))
            # print(targets)
            # Build tabs for each condition
            link_d = OrderedDict()
            for tar in sorted([os.path.basename(t) for t in targets]):
                link_d[tar] = "../" + tar + "/index.html"

            # For each condition
            for target in targets:
                merge_DBD_regions(path=target)
                merge_DBSs(path=target)
                # stat
                integrate_stat(path=target)
                summerize_stat(target=target, link_d=link_d, score=args.exp)
            # Project level index file

            for item in os.listdir(args.path):
                print("\t"+item)
                if item == "style": continue
                elif os.path.isfile(os.path.join(args.path, item)):
                    continue
                elif os.path.isdir(os.path.join(args.path, item)):

                    for it in os.listdir(os.path.join(args.path, item)):
                        # print("\t" + item + "\t" + it)
                        if it == "style": continue
                        elif os.path.isfile(os.path.join(args.path, item, it)):
                            continue
                        elif os.path.isdir(os.path.join(args.path, item, it)):
                            integrate_html(os.path.join(args.path, item, it))
                    integrate_html(os.path.join(args.path, item))
            integrate_html(args.path)

            # gen_heatmap(path=args.path)
            # generate_rna_exp_pv_table(root=args.path, multi_corr=False)
            # merge_DBD_regions(path=args.path)

            sys.exit(0)

        ####################################################################################
        ######### updatehtml
        elif args.mode == "updatehtml":
            for item in os.listdir(args.path):
                pro = os.path.join(args.path, item, "profile.txt")
                if os.path.isfile(pro): update_profile(dirpath=os.path.join(args.path, item),
                                                       expression=args.exp)
            revise_index(root=args.path)
            generate_rna_exp_pv_table(root=args.path, multi_corr=True)
            sys.exit(0)
        
        ####################################################################################
        ######### get_dbss
        elif args.mode == "get_dbss":

            get_dbss(input_BED=args.i,output_BED=args.dbs,rna_fasta=args.r,output_rbss=args.rbs,
                     organism=args.organism,l=args.l,e=args.e,c=args.c,
                     fr=args.fr,fm=args.fm,of=args.of,mf=args.mf,rm=args.rm,temp=dir)
            os.remove("dna_targeted_region.fa")
            os.remove("dna_targeted_region.txp")
            os.remove("rna_temp.fa")
            sys.exit(0)


        #######################################################################
        #### Checking arguments
        if not args.o: 
            print("Please define the output directory name. \n")
            sys.exit(1)
        if not args.organism: 
            print("Please define the organism. (hg19 or mm9)")
            sys.exit(1)
        if not args.rn:
            print("Please define RNA sequence name.")
            sys.exit(1)
        # if args.r and args.rl:
        #     print("Both -r and -rl are given. TDF will skip -r and process -rl ")
        # if args.rl:
        #     with open(args.rl) as f:
        #         for line in f:
        #             line = line.strip()
        #             rn = os.path.basename(line).rpartition(".")[0]
        #             print("\tProcessing: "+rn)
        #             command = ["rgt-TDF", args.mode,
        #                        "-r", line, "-rn", rn,
        #                        "-o", os.path.join(args.o, rn),
        #                        "-organism", args.organism ]
        #             if args.de and not args.bed: command += ["-de", args.de]
        #             if args.bed and args.bg: command += ["-bed", args.bed, "-bg", args.bg]
        #
        #             if args.score: command += ["-score"]
        #             if args.rt: command += ["-rt" ]
        #             if args.pl != 1000: command += ["-pl", args.pl]
        #             if args.ccf != 40: command += ["-ccf", args.ccf]
        #             if args.obed: command += ["-obed"]
        #             if args.a != 0.05: command += ["-a", args.a]
        #             if args.filter_havana == 'F': command += ["-filter_havana", 'F']
        #             if args.protein_coding == 'T': command += ["-protein_coding", 'T']
        #             if args.known_only == 'F': command += ["-known_only", 'F']
        #
        #             if args.rm > 0: command += ["-rm", args.rm ]
        #             if args.fr != 'off': command += ["-fr", args.fr ]
        #             if args.c != 2: command += ["-c", args.c ]
        #             if args.e != 20: command += ["-e", args.e ]
        #             if args.of != 1: command += ["-of", args.of ]
        #             if args.l != 15: command += ["-l", args.l ]
        #             if args.fr != 'off': command += ["-fr", args.fr ]
        #             if args.fr != 'off': command += ["-fr", args.fr ]
        #             if args.fr != 'off': command += ["-fr", args.fr ]
        #             subprocess.call(command)
        #     sys.exit(0)

        t0 = time.time()
        # Normalised output path
        if not args.t: args.t = args.rn
        # else: title = args.t
        args.r = os.path.normpath(os.path.join(dir, args.r))
        args.o = os.path.normpath(os.path.join(dir, args.o, args.t))
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
        # Checking all arguments
        if args.bed and not args.bg:
            print("Please add background promoters in BED format. (-bg)")
            sys.exit(1)
        if args.scoreh and not args.score:
            print("Score header (-scoreh) can only be used when scores (-score) are loaded.")
            print("Please add '-score'.")
            sys.exit(1)
        if args.de: args.de = os.path.normpath(os.path.join(dir, args.de))
        if args.bed: args.bed = os.path.normpath(os.path.join(dir, args.bed))
        if args.bg: args.bg = os.path.normpath(os.path.join(dir, args.bg))

        if args.filter_havana == "T": args.filter_havana = True
        else: args.filter_havana = False
        if args.protein_coding == "T": args.protein_coding = True
        else: args.protein_coding = False
        if args.known_only == "T": args.known_only = True
        else: args.known_only = False

        print("\n*************** Promoter Test ****************")
        print("*** Input RNA sequence: " + args.r)
        print("*** Output directory: " + shorten_dir(args.o))
        print("Step 1: Calculate the triplex forming sites on RNA and DNA.")
        #######################################
        # Input
        tdf_input = Input(pars=args)
        if args.de:
            tdf_input.dna.degenes()
        elif args.bed:
            tdf_input.dna.de_bed_input()
        #######################################
        # Triplexes
        triplexes = Triplexes(organism=args.organism, pars=args)
        tpx_de = triplexes.search_triplex(target_regions=tdf_input.dna.target_regions,
                                          prefix="target_promoters", remove_temp=True)
        tpx_nde = triplexes.search_triplex(target_regions=tdf_input.dna.nontarget_regions,
                                           prefix="nontarget_promoters", remove_temp=True)
        t1 = time.time()
        print("\tRunning time: " + str(datetime.timedelta(seconds=round(t1-t0))))
        #######################################
        # Statistics
        print("Step 2: Calculate the frequency of DNA binding sites within the promoters.")
        stat = Statistics(pars=args)
        stat.count_frequency_promoters(target_regions=tdf_input.dna.target_regions,
                                       background=tdf_input.dna.nontarget_regions,
                                       file_tpx_de=tpx_de, file_tpx_nde=tpx_nde)
        triplexes.find_autobinding(rbss=stat.rbss)
        stat.fisher_exact_de()
        stat.dbs_motif(tpx=stat.tpx_def)
        stat.uniq_motif(tpx=stat.tpx_def, rnalen=tdf_input.rna.seq_length)
        stat.dbd_regions(rna_exons=tdf_input.rna.regions)
        stat.output_bed(input=tdf_input, tpx=stat.tpx_def)
        stat.summary_stat(input=tdf_input, triplexes=triplexes, mode="promotertest")
        stat.write_stat(filename=os.path.join(args.o, "stat.txt"))
        t2 = time.time()
        print("\tRunning time: " + str(datetime.timedelta(seconds=round(t2 - t1))))
        #######################################
        # Reports
        print("Step 3: Generate plot and output html files.")
        if len(stat.rbss) == 0:
            no_binding_response(args=args,  stat=stat.stat)
        else:
            reports = Report(pars=args, input=tdf_input, triplexes=triplexes, stat=stat)
            reports.plot_lines(tpx=stat.tpx_def, ylabel="Number of DBSs",
                               linelabel="No. DBSs", filename=args.rn + "_lineplot.png")
            reports.barplot(filename=args.rn+"_barplot.png")
            reports.gen_html_promotertest()
            reports.gen_html_genes()
        t3 = time.time()
        print("\tRunning time: " + str(datetime.timedelta(seconds=round(t3 - t2))))
        silentremove(os.path.join(args.o, "rna_temp.fa"))
        silentremove(os.path.join(args.o, "rna_temp.fa.fai"))
        silentremove(os.path.join(args.o, "de.fa"))
        silentremove(os.path.join(args.o, "nde.fa"))
        silentremove(os.path.join(args.o, "de.txp"))
        silentremove(os.path.join(args.o, "autobinding.tpx"))
        print("\nTotal running time: " + str(datetime.timedelta(seconds=round(t3 - t0))))


    ################################################################################
    ##### Genomic Region Test ######################################################
    ################################################################################
    if args.mode == 'regiontest':
        if args.bed:
            args.bed = os.path.normpath(os.path.join(dir, args.bed))

        print("\n"+"*************** Genomic Region Test ***************")
        print("*** Input RNA sequence: " + args.r)
        print("*** Input regions in BED: " + os.path.basename(args.bed))
        print("*** Number of randomization: " + str(args.n))
        print("*** Output directory: " + os.path.basename(args.o))
        print("Step 1: Calculate the triplex forming sites on RNA and DNA.")

        #######################################
        # Input
        tdf_input = Input(pars=args)
        tdf_input.dna.bed_input(bed=args.bed)

        #######################################
        # Triplexes
        triplexes = Triplexes(organism=args.organism, pars=args)
        stat = Statistics(pars=args)
        stat.tpx = triplexes.get_tpx(rna_fasta_file=os.path.join(args.o,"rna_temp.fa"),
                                     target_regions=tdf_input.dna.target_regions,
                                     prefix="target_regions", remove_temp=args.rt, dna_fine_posi=False)

        stat.tpxf = triplexes.get_tpx(rna_fasta_file=os.path.join(args.o,"rna_temp.fa"),
                                      target_regions=tdf_input.dna.target_regions,
                                      prefix="target_regions_fine", remove_temp=args.rt, dna_fine_posi=True)
        t1 = time.time()
        print("\tRunning time: " + str(datetime.timedelta(seconds=round(t1 - t0))))


        #######################################
        # Statistics
        print("Step 2: Permutation by randomization the target regions for "+str(args.n)+ " times.")
        stat.target_stat(target_regions=tdf_input.dna.target_regions, tpx=stat.tpx, tpxf=stat.tpxf)
        triplexes.find_autobinding(rbss=stat.rbss)

        if len(stat.rbss) == 0:
            stat.summary_stat(input=tdf_input, triplexes=triplexes, mode="regiontest", no_binding=True)
            no_binding_response(args=args,  stat=stat.stat)
        else:
            stat.random_test(repeats=args.n, target_regions=tdf_input.dna.target_regions,
                             filter_bed=args.f, mp=args.mp, genome_fasta=triplexes.genome.get_genome())

            stat.dbs_motif(tpx=stat.tpxf)
            stat.uniq_motif(tpx=stat.tpxf, rnalen=tdf_input.rna.seq_length)
            stat.dbd_regions(rna_exons=tdf_input.rna.regions)
            stat.output_bed(input=tdf_input, tpx=stat.tpxf)
            stat.summary_stat(input=tdf_input, triplexes=triplexes, mode="regiontest")
            stat.write_stat(filename=os.path.join(args.o, "stat.txt"))
            t2 = time.time()
            print("\tRunning time: " + str(datetime.timedelta(seconds=round(t2 - t1))))
            #######################################
            # Reports
            print("Step 3: Generate plot and output html files.")
            if len(stat.rbss) == 0:
                no_binding_response(args=args, stat=stat.stat)

            reports = Report(pars=args, input=tdf_input, triplexes=triplexes, stat=stat)
            reports.plot_lines(tpx=stat.tpx, ylabel="Number of DBSs",
                               linelabel="No. DBSs", filename=args.rn + "_lineplot.png")
            reports.boxplot(filename=args.rn + "_boxplot.png", matrix=stat.region_matrix, sig_region=stat.sig_DBD,
                            truecounts=stat.counts_dbs.values(), sig_boolean=stat.data["region"]["sig_boolean"],
                            ylabel="Number of DBS on target regions")
            reports.gen_html_regiontest()

            t3 = time.time()
            print("\tRunning time: " + str(datetime.timedelta(seconds=round(t3 - t2))))
            silentremove(os.path.join(args.o, "rna_temp.fa"))
            silentremove(os.path.join(args.o, "rna_temp.fa.fai"))
            silentremove(os.path.join(args.o, "de.fa"))
            # silentremove(os.path.join(args.o, "de.tpx"))
            silentremove(os.path.join(args.o, "autobinding.tpx"))
            print("\nTotal running time: " + str(datetime.timedelta(seconds=round(t3 - t0))))
