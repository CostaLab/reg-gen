# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import re
import time
import shutil
import getpass
import argparse
import datetime
import subprocess
import matplotlib
matplotlib.use('Agg', warn=False)
from collections import *
import natsort as natsort_ob

# Local Libraries
# Distal Libraries
from rgt import __version__
from rgt.Util import Html
from triplexTools import rna_associated_gene, get_dbss, check_dir,\
                         gen_heatmap, generate_rna_exp_pv_table, revise_index, print2, \
                         output_summary, list_all_index, no_binding_response, write_stat, \
                         integrate_stat, update_profile, integrate_html, \
                         merge_DBD_regions, rank_array, silentremove, summerize_stat

from tdf_promotertest import PromoterTest
from tdf_regiontest import RandomTest

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
    
    parser = argparse.ArgumentParser(description='Triplex Domain Finder is a statistical framework \
                                                  for detection of triple helix potential of \
                                                  lncRNAs from genome-wide functional data. \
                                                  Author: Chao-Chung Kuo\
                                                  \nVersion: ' + __version__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Promoter test ##########################################

    h_promotor = "Promoter test evaluates the association between the given lncRNA to the target promoters."
    parser_promotertest = subparsers.add_parser('promotertest', help=h_promotor)
    parser_promotertest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_promotertest.add_argument('-rl', type=str, default=None, metavar='  ', help="Input list for paths to all RNA sequences (in fasta format)")
    parser_promotertest.add_argument('-rn', type=str, default=None, metavar='  ', help="Define the RNA name")
    parser_promotertest.add_argument('-de', default=False, metavar='  ', help="Input file for target gene list (gene symbols or Ensembl ID)")
    parser_promotertest.add_argument('-bed', default=False, metavar='  ', help="Input BED file of the promoter regions of target genes")
    parser_promotertest.add_argument('-bg', default=False, metavar='  ', help="Input BED file of the promoter regions of background genes")
    parser_promotertest.add_argument('-o', metavar='  ', help="Output directory name for all the results")
    parser_promotertest.add_argument('-t', metavar='  ', default=False, help="Define the title name for the results under the Output name. (default: %(default)s)")
    
    parser_promotertest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_promotertest.add_argument('-gtf', metavar='  ', default=None, help='Define the GTF file for annotation (optional)')

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

    parser_promotertest.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexator] Define the minimum length of triplex (default: %(default)s)")
    parser_promotertest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_promotertest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_promotertest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_promotertest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_promotertest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (default: %(default)s)")
    parser_promotertest.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_promotertest.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexator] Set the multiprocessing")
    parser_promotertest.add_argument('-par', type=str, default="", metavar='  ', help="[Triplexator] Define other parameters for Triplexator")
    
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

    parser_randomtest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
 
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
    
    parser_randomtest.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexator] Define the minimum length of triplex (default: %(default)s)")
    parser_randomtest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_randomtest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_randomtest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_randomtest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_randomtest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (default: %(default)s)")
    parser_randomtest.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_randomtest.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexator] Set the multiprocessing")
    parser_randomtest.add_argument('-par', type=str, default="", metavar='  ', help="[Triplexator] Define other parameters for Triplexator")

    ##########################################################################
    parser_bed2bed = subparsers.add_parser('get_dbss', help="Get DBSs in BED format from the single BED file")
    parser_bed2bed.add_argument('-i',type=str, metavar='  ', help='Input BED file of the target regions')
    parser_bed2bed.add_argument('-dbs',type=str, metavar='  ', help='Output BED file of the DBSs')
    parser_bed2bed.add_argument('-rbs',type=str, metavar='  ', help='Output BED file of the RBSs')
    parser_bed2bed.add_argument('-r',type=str, metavar='  ', help='Input FASTA file of the RNA')
    parser_bed2bed.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_bed2bed.add_argument('-l', type=int, default=20, metavar='  ', help="[Triplexator] Define the minimum length of triplex (default: %(default)s)")
    parser_bed2bed.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (default: %(default)s)")
    parser_bed2bed.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (default: %(default)s)")
    parser_bed2bed.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (default: %(default)s)")
    parser_bed2bed.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (default: %(default)s).'0' = greedy approach; '1' = q-gram filtering.")
    parser_bed2bed.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (default: %(default)s)")
    parser_bed2bed.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_bed2bed.add_argument('-rm', type=int, default=0, metavar='  ', help="[Triplexator] Set the multiprocessing")
    
    ##########################################################################
    # rgt-TDF integrate -path 
    parser_integrate = subparsers.add_parser('integrate', help="Integrate the project's links and generate project-level statistics.")
    parser_integrate.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
    ##########################################################################
    parser_updatehtml = subparsers.add_parser('updatehtml', help="Update the project's html.")
    parser_updatehtml.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
    parser_updatehtml.add_argument('-exp', type=str, metavar='  ', help='Define file with expression data.')

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
                # stat
                integrate_stat(path=target)
                summerize_stat(target=target, link_d=link_d)
            # Project level index file
            integrate_html(args.path)
            for item in os.listdir(args.path):
                # print("\t"+item)
                if item == "style": continue
                elif os.path.isfile(os.path.join(args.path, item)):
                    continue
                elif os.path.isdir(os.path.join(args.path, item)):
                    integrate_html(os.path.join(args.path, item))
                    for it in os.listdir(os.path.join(args.path, item)):
                        # print("\t" + item + "\t" + it)
                        if it == "style": continue
                        elif os.path.isfile(os.path.join(args.path, item, it)):
                            continue
                        elif os.path.isdir(os.path.join(args.path, item, it)):
                            integrate_html(os.path.join(args.path, item, it))


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
        if not args.rn and not args.rl: 
            print("Please define RNA sequence name.")
            sys.exit(1)
        if args.r and args.rl:
            print("Both -r and -rl are given. TDF will skip -r and process -rl ")
        if args.rl:
            with open(args.rl) as f:
                for line in f:
                    line = line.strip()
                    rn = os.path.basename(line).rpartition(".")[0]
                    print("\tProcessing: "+rn)
                    command = ["rgt-TDF", args.mode, 
                               "-r", line, "-rn", rn,
                               "-o", os.path.join(args.o, rn),
                               "-organism", args.organism ]
                    if args.de and not args.bed: command += ["-de", args.de]
                    if args.bed and args.bg: command += ["-bed", args.bed, "-bg", args.bg]

                    if args.score: command += ["-score"]
                    if args.rt: command += ["-rt" ]
                    if args.pl != 1000: command += ["-pl", args.pl]
                    if args.ccf != 40: command += ["-ccf", args.ccf]
                    if args.obed: command += ["-obed"]
                    if args.a != 0.05: command += ["-a", args.a]
                    if args.filter_havana == 'F': command += ["-filter_havana", 'F']
                    if args.protein_coding == 'T': command += ["-protein_coding", 'T']
                    if args.known_only == 'F': command += ["-known_only", 'F']
                    
                    if args.rm > 0: command += ["-rm", args.rm ]
                    if args.fr != 'off': command += ["-fr", args.fr ]
                    if args.c != 2: command += ["-c", args.c ]
                    if args.e != 20: command += ["-e", args.e ]
                    if args.of != 1: command += ["-of", args.of ]
                    if args.l != 15: command += ["-l", args.l ]
                    if args.fr != 'off': command += ["-fr", args.fr ]
                    if args.fr != 'off': command += ["-fr", args.fr ]
                    if args.fr != 'off': command += ["-fr", args.fr ]
                    subprocess.call(command)
            sys.exit(0)

        t0 = time.time()
        # Normalised output path
        if not args.t: title = args.rn
        else: title = args.t
        
        args.o = os.path.normpath(os.path.join(dir,args.o,title))
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


################################################################################################3

        if args.bed and not args.bg:
            print("Please add background promoters in BED format. (-bg)")
            sys.exit(1)
        if args.scoreh and not args.score:
            print("Score header (-scoreh) can only be used when scores (-score) are loaded.")
            print("Please add '-score'.")
            sys.exit(1)

        print2(summary, "\n"+"*************** Promoter Test ****************")
        print2(summary, "*** Input RNA sequence: "+args.r)
        
        if args.o.count("/") < 3:
            print2(summary, "*** Output directory: "+ args.o)
        else:
            n = args.o.count("/") - 3 + 1
            print2(summary, "*** Output directory: "+ args.o.split("/",n)[-1] )

        args.r = os.path.normpath(os.path.join(dir, args.r))
        
        if args.de: args.de = os.path.normpath(os.path.join(dir,args.de))
        if args.bed: args.bed = os.path.normpath(os.path.join(dir,args.bed))
        if args.bg: args.bg = os.path.normpath(os.path.join(dir,args.bg))

        # Get GenomicRegionSet from the given genes
        print2(summary, "Step 1: Calculate the triplex forming sites on RNA and DNA.")
        promoter = PromoterTest(gene_list_file=args.de, gtf=args.gtf, rna_name=args.rn, bed=args.bed, bg=args.bg, 
                                organism=args.organism, promoterLength=args.pl, summary=summary, 
                                temp=dir, output=args.o, showdbs=args.showdbs, score=args.score, 
                                scoreh=args.scoreh, filter_havana=args.filter_havana, 
                                protein_coding=args.protein_coding, known_only=args.known_only)
        if args.dump: sys.exit(0)
        promoter.get_rna_region_str(rna=args.r, expfile=args.rnaexp)
        promoter.connect_rna(rna=args.r, temp=args.o)
        promoter.search_triplex(temp=args.o, l=args.l, e=args.e, remove_temp=args.rt, 
                                c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, par=args.par)

        t1 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t1-t0))))

        print2(summary, "Step 2: Calculate the frequency of DNA binding sites within the promotors.")
        if args.obed: obedp = os.path.basename(args.o)
        else: obedp = None
        promoter.count_frequency(temp=args.o, remove_temp=args.rt, obedp=obedp, cutoff=args.ccf, l=args.l)
        promoter.fisher_exact(alpha=args.a)
        t2 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t2-t1))))
        promoter.dbd_regions(output=args.o)
        promoter.autobinding(output=args.o, l=args.l, e=args.e,
                             c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, par=args.par)
        # promoter.dbs_motif(outdir=args.o)
        if len(promoter.rbss) == 0:
            no_binding_response(args=args, rna_regions=promoter.rna_regions,
                                rna_name=promoter.rna_name, organism=promoter.organism,
                                stat=promoter.stat, expression=promoter.rna_expression)


        print2(summary, "Step 3: Establishing promoter profile.")
        t3 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t3-t2))))

        print2(summary, "Step 4: Generate plot and output html files.")
        promoter.plot_lines(txp=promoter.txp_def, rna=args.r, dirp=args.o, ac=args.ac, 
                            cut_off=args.accf, log=args.log, showpa=args.showpa,
                            sig_region=promoter.sig_DBD,
                            ylabel="Number of DBSs", 
                            linelabel="No. DBSs", filename=promoter.rna_name+"_lineplot.png")

        promoter.barplot(dirp=args.o, filename=promoter.rna_name+"_barplot.png", sig_region=promoter.sig_DBD)
        promoter.uniq_motif()
        #if args.showdbs:
        #    promoter.plot_lines(txp=promoter.txp_def, rna=args.r, dirp=args.o, ac=args.ac, 
        #                        cut_off=args.accf, log=args.log, showpa=args.showpa,
        #                        sig_region=promoter.sig_region_dbs,
        #                        ylabel="Number of DBSs on target promoters", 
        #                        linelabel="No. DBSs", filename="plot_dbss.png")
        #    promoter.barplot(dirp=args.o, filename="bar_dbss.png", sig_region=promoter.sig_region_dbs, dbs=True)
        # if args.motif: promoter.gen_motifs(temp=args.o)

        promoter.gen_html(directory=args.o, parameters=args, ccf=args.ccf, align=50, alpha=args.a)
        promoter.gen_html_genes(directory=args.o, align=50, alpha=args.a, nonDE=False)
        # promoter.save_table(path=os.path.dirname(args.o), table=promoter.ranktable,
        #                         filename="lncRNA_target_ranktable.txt")
        # promoter.save_table(path=os.path.dirname(args.o), table=promoter.dbstable,
        #                         filename="lncRNA_target_dbstable.txt")

        #promoter.heatmap(table="ranktable.txt", temp=os.path.dirname(args.o))

        t4 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t4-t3))))
        print2(summary, "\nTotal running time is: " + str(datetime.timedelta(seconds=round(t4-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        # save_profile(rna_regions=promoter.rna_regions, rna_name=promoter.rna_name,
        #              organism=promoter.organism, output=args.o, bed=args.bed,
        #              geneset=args.de, stat=promoter.stat, topDBD=promoter.topDBD,
        #              sig_DBD=promoter.sig_DBD, expression=promoter.rna_expression)
        # revise_index(root=os.path.dirname(os.path.dirname(args.o)))


        silentremove(os.path.join(args.o,"rna_temp.fa"))
        silentremove(os.path.join(args.o,"rna_temp.fa.fai"))
        silentremove(os.path.join(args.o, "de.fa"))
        silentremove(os.path.join(args.o, "nde.fa"))
        silentremove(os.path.join(args.o, "de.txp"))
        silentremove(os.path.join(args.o, "autobinding.txp"))


        write_stat(stat=promoter.stat, filename=os.path.join(args.o, "stat.txt"))


    ################################################################################
    ##### Genomic Region Test ######################################################
    ################################################################################
    if args.mode == 'regiontest':
        def no_binding_code():
            print("*** Find no triple helices binding on the given RNA")

            pro_path = os.path.join(os.path.dirname(args.o), "profile.txt")
            exp = os.path.basename(args.o)
            tar_reg = os.path.basename(args.bed)
            r_genes = rna_associated_gene(rna_regions=randomtest.rna_regions, name=randomtest.rna_name, organism=randomtest.organism)
            newlines = []
            if os.path.isfile(pro_path):
                with open(pro_path,'r') as f:
                    new_exp = True
                    for line in f:
                        line = line.strip()
                        line = line.split("\t")
                        if line[0] == exp:
                            newlines.append([exp, args.rn, args.o.split("_")[-1],
                                             args.organism, tar_reg, "0",
                                             "-", "1.0", r_genes, "No triplex found" ])
                            new_exp = False
                        else:
                            newlines.append(line)
                    if new_exp:
                        newlines.append([exp, args.rn, args.o.split("_")[-1],
                                         args.organism, tar_reg,"0",
                                         "-", "1.0", r_genes, "No triplex found" ])
            else:
                newlines.append(["Experiment","RNA_names","Tag","Organism","Target_region","No_sig_DBDs",
                                 "Top_DBD", "p-value","closest_genes"])
                newlines.append([exp, args.rn, args.o.split("_")[-1],
                                 args.organism, tar_reg, "0",
                                 "-", "1.0", r_genes, "No triplex found" ])
            with open(pro_path,'w') as f:
                for lines in newlines:
                    print("\t".join(lines), file=f)

            #shutil.rmtree(args.o)
            list_all_index(path=os.path.dirname(args.o), show_RNA_ass_gene=randomtest.rna_regions)
            shutil.rmtree(args.o)
            sys.exit(1)

            #########################################################
        print2(summary, "\n"+"*************** Genomic Region Test ***************")
        print2(summary, "*** Input RNA sequence: "+args.r)
        print2(summary, "*** Input regions in BED: "+os.path.basename(args.bed))
        print2(summary, "*** Number of randomization: "+str(args.n))
        print2(summary, "*** Output directoey: "+os.path.basename(args.o))

        args.r = os.path.normpath(os.path.join(dir,args.r))

        print2(summary, "\nStep 1: Calculate the triplex forming sites on RNA and the given regions")
        randomtest = RandomTest(rna_fasta=args.r, rna_name=args.rn, dna_region=args.bed,
                                organism=args.organism, showdbs=args.showdbs)
        randomtest.get_rna_region_str(rna=args.r)
        obed = os.path.basename(args.o)
        randomtest.connect_rna(rna=args.r, temp=args.o)

        randomtest.target_dna(temp=args.o, remove_temp=args.rt, l=args.l, e=args.e, obed=obed,
                              c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, par=args.par, cutoff=args.ccf )
        t1 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t1-t0))))
        # print(args.par)
        randomtest.autobinding(output=args.o, l=args.l, e=args.e,
                               c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, par=args.par)
        # randomtest.dbs_motif(outdir=args.o)
        if len(randomtest.rbss) == 0:
            # no_binding_code()
            no_binding_response(args=args, rna_regions=randomtest.rna_regions,
                                rna_name=randomtest.rna_name, organism=randomtest.organism,
                                stat=randomtest.stat, expression=randomtest.rna_expression)

        print2(summary, "Step 2: Randomization and counting number of binding sites")
        randomtest.random_test(repeats=args.n, temp=args.o, remove_temp=args.rt, l=args.l, e=args.e,
                               c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, par=args.par, rm=args.rm,
                               filter_bed=args.f, alpha=args.a, mp=args.mp)
        randomtest.dbd_regions(sig_region=randomtest.data["region"]["sig_region"], output=args.o)

        t2 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t2-t1))))

        print2(summary, "Step 3: Generating plot and output HTML")

        randomtest.lineplot(txp=randomtest.txpf, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
                            log=args.log, ylabel="Number of DBS",
                            sig_region=randomtest.data["region"]["sig_region"],
                            linelabel="No. DBS", filename=randomtest.rna_name+"_lineplot.png")

        #randomtest.lineplot(txp=randomtest.txp, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
        #                    log=args.log, ylabel="Number of target regions with DBS", 
        #                    sig_region=randomtest.data["region"]["sig_region"],
        #                    linelabel="No. target regions", filename="lineplot_region.png")

        randomtest.boxplot(dir=args.o, matrix=randomtest.region_matrix,
                           sig_region=randomtest.data["region"]["sig_region"],
                           truecounts=[r[0] for r in randomtest.counts_tr.values()],
                           sig_boolean=randomtest.data["region"]["sig_boolean"],
                           ylabel="Number of target regions",
                           filename=randomtest.rna_name+"_boxplot" )
        #if args.showdbs:
        #    randomtest.lineplot(txp=randomtest.txpf, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
        #                        log=args.log, ylabel="Number of DBS on target regions",
        #                        sig_region=randomtest.data["dbs"]["sig_region"], 
        #                        linelabel="No. DBS", filename="lineplot_dbs.png")

        #    randomtest.boxplot(dir=args.o, matrix=randomtest.dbss_matrix, 
        #                       sig_region=randomtest.data["dbs"]["sig_region"], 
        #                       truecounts=randomtest.counts_dbs.values(),
        #                       sig_boolean=randomtest.data["dbs"]["sig_boolean"], 
        #                       ylabel="Number of DBS on target regions",
        #                       filename="boxplot_dbs" )
        randomtest.uniq_motif()
        randomtest.gen_html(directory=args.o, parameters=args, align=50, alpha=args.a,
                            score=args.score, obed=obed)

        t3 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t3-t2))))

        print2(summary, "\nTotal running time is: " + str(datetime.timedelta(seconds=round(t3-t0))))

        output_summary(summary, args.o, "summary.txt")

        for f in os.listdir(args.o):
            if re.search("dna*.fa", f) or re.search("dna*.txp", f):
                silentremove(os.path.join(args.o, f))
        silentremove(os.path.join(args.o, "rna_temp.fa"))
        silentremove(os.path.join(args.o, "rna_temp.fa.fai"))
        silentremove(os.path.join(args.o, randomtest.rna_name+"_autobinding.txp"))

        write_stat(stat=randomtest.stat, filename=os.path.join(args.o, "stat.txt"))
