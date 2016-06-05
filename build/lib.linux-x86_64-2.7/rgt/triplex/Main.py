# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import re
import argparse
import shutil
import time, datetime, getpass, fnmatch
import subprocess
import pickle
from collections import OrderedDict, defaultdict
import numpy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import colors
import scipy.cluster.hierarchy as sch

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from triplexTools import PromoterTest, RandomTest, value2str,\
                         split_gene_name, rna_associated_gene, get_dbss
from rgt.SequenceSet import Sequence, SequenceSet
from rgt.Util import SequenceType, Html, ConfigurationFile
from rgt.motifanalysis.Statistics import multiple_test_correction

dir = os.getcwd()

"""
Triplex Domain Finder (TDF) provides statistical tests and plotting tools for 
triplex binding site analysis

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
            print("********* RGT Triplex: Summary information *********", 
                  file=f)
            for s in summary:
                print(s, file=f)
    
def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try: os.stat(path)
    except: 
        try: os.mkdir(path)
        except: pass

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp

def list_all_index(path, link_d=None, show_RNA_ass_gene=False):
    """Creat an 'index.html' in the defined directory """

    dirname = os.path.basename(path)
    
    if link_d: pass
    else: link_d = {"List":"index.html"}

    html = Html(name="Directory: "+dirname, links_dict=link_d, 
                fig_rpath="./style", fig_dir=os.path.join(path,"style"), 
                RGT_header=False, other_logo="TDF", homepage="../index.html")
    
    html.add_heading("All experiments in: "+dirname+"/")

    
    data_table = []
    type_list = 'sssssssssssss'
    col_size_list = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
    c = 0
    if show_RNA_ass_gene:
        header_list = ["No.", "Experiments", "RNA", "Closest genes", 
                       "No sig. DBD",
                       "Top DBD", "p-value", "Organism", "Target region"]
    else:
        header_list = ["No.", "Experiments", "RNA", "No sig. DBD",
                       "Top DBD", "p-value", "Organism", #"Condition", 
                       "Target region" ]

    profile_f = open(os.path.join(path, "profile.txt"),'r')
    profile = {}
    for line in profile_f:
        line = line.strip()
        line = line.split("\t")
        profile[line[0]] = line[1:]
    #profile = pickle.load(profile_f)
    
    #for root, dirnames, filenames in os.walk(path):
        #roots = root.split('/')
        #for filename in fnmatch.filter(filenames, '*.html'):
        #    if filename == 'index.html' and root.split('/')[-1] != dirname:
    #    for i, dirname in enumerate(dirnames):
            
            #if dirname in profile.keys():
    for i, exp in enumerate(profile.keys()):
        #print(exp)
        c += 1
        #exp = root.split('/')[-1]
        #exp = dirn
        #if profile[exp][5] == "-":
        #    new_line = [ str(c), exp, profile[exp][0] ]
        #else:
        try:
            if profile[exp][5] == "-":
                new_line = [ str(c), exp,
                             profile[exp][0] ]
            else:
                new_line = [ str(c), 
                             '<a href="'+os.path.join(exp, "index.html")+\
                             '">'+exp+"</a>",
                             profile[exp][0] ]

            if show_RNA_ass_gene: 
                new_line.append( 
                    split_gene_name(gene_name=profile[exp][7], 
                                    org=profile[exp][2]) 
                    )

            if profile[exp][6] == "-":
                new_line += [ profile[exp][4], 
                              profile[exp][5], profile[exp][6],
                              profile[exp][2], profile[exp][3] ]

            elif float(profile[exp][6]) < 0.05:
                new_line += [ profile[exp][4], 
                              profile[exp][5], 
                              "<font color=\"red\">"+\
                              profile[exp][6]+"</font>",
                              profile[exp][2], profile[exp][3] ]
            else:
                new_line += [ profile[exp][4], 
                              profile[exp][5], profile[exp][6],
                              profile[exp][2], profile[exp][3] ]
            data_table.append(new_line)
        except:
            if exp != "Experiment":
                print("Error in loading profile: "+exp)
            continue

    html.add_zebra_table( header_list, col_size_list, type_list, data_table, 
                          align=10, cell_align="left", sortable=True)
    #if os.path.isfile(os.path.join(path,"lncRNA_target_dendrogram.png")):
    #    html.add_figure("lncRNA_target_dendrogram.png", align="center", width="90%")

    html.add_fixed_rank_sortable()
    html.write(os.path.join(path,"index.html"))

def revise_index(root, show_RNA_ass_gene=False):
    "Revise other index.html in the same project"
    
    dirlist = {}
    plist = {}
    for item in os.listdir(root):
        h = os.path.join(root, item, "index.html")
        pro = os.path.join(root, item, "profile.txt")
        if os.path.isfile(pro):
            dirlist[os.path.basename(item)] = "../"+item+"/index.html"
            plist[os.path.basename(item)] = h
    dirlist = OrderedDict(sorted(dirlist.items()))
    #print(dirlist)
    for d, p in plist.iteritems():
        list_all_index(path=os.path.dirname(p), 
                       link_d=dirlist, show_RNA_ass_gene=show_RNA_ass_gene)

def gen_heatmap(path):
    """Generate the heatmap to show the sig RNA in among conditions"""
    def fmt(x, pos):
        # a, b = '{:.0e}'.format(x).split('e')
        # b = int(b)
        # return r'${} \times 10^{{{}}}$'.format(a, b)
        a = -numpy.log10(x)
        return '{:.0f}'.format(a)


    matrix = OrderedDict()
    rnas = []
    for item in os.listdir(path):
        # print(item)
        # print(os.path.isdir(os.path.join(path,item)))
        if not os.path.isdir(os.path.join(path,item)): continue
        if item == "style": continue
        #if item == "index.html": continue
        matrix[item] = {}
        pro = os.path.join(path, item, "profile.txt")
        with open(pro) as f:
            for line in f:
                line = line.strip().split("\t")
                if line[0] == "Experiment": continue
                if line[6] == "-": continue
                matrix[item][line[0]] = float(line[7])
                rnas.append(line[0])
    rnas = list(set(rnas))
    # rnas.sort()

    # Convert into array
    ar = []
    exps = natsorted(matrix.keys())
    # print(exps)
    for exp in exps:
        row = []
        for rna in rnas:
            try: row.append(matrix[exp][rna])
            except: row.append(1)
        ar.append(row)
    ar = numpy.array(ar)
    ar = numpy.transpose(ar)

    # print(ar.shape)
    data = ar[~numpy.all(ar == 1, axis=1)]
    # print(data.shape)
    

    fig = plt.figure(figsize=(len(matrix.keys())*1.5,len(rnas)*2.5))
    # fig = plt.figure()
    # ax1 = fig.add_axes([0.09,0.2,0.2,0.6])
    # Y = sch.linkage(data, method='single')
    # Z1 = sch.dendrogram(Y, orientation='right')
    # ax1.set_xticks([])
    # ax1.set_yticks([])

    # Compute and plot second dendrogram.
    # ax2 = fig.add_axes([0.3, 1.1, 0.55,0.04])
    # Y = sch.linkage(data.T, method='single')
    # Z2 = sch.dendrogram(Y)
    # ax2.set_xticks([])
    # ax2.set_yticks([])
    bounds = []
    for n in range(-8,0):
        # print(10**n)
        bounds.append(10**n)
    norm = colors.BoundaryNorm(bounds, plt.cm.YlOrRd_r.N)

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.1,0.2,0.8,0.6])
    #axmatrix = fig.add_axes()
    # idx1 = Z1['leaves']
    #idx2 = Z2['leaves']
    # data = data[idx1,:]
    #data = data[:,idx2]
    im = axmatrix.matshow(data, aspect='auto', origin='lower', cmap=plt.cm.YlOrRd_r, norm=norm)

    axmatrix.set_xticks(range(data.shape[1]))
    axmatrix.set_xticklabels( exps, minor=False, ha="left")
    axmatrix.xaxis.set_label_position('top')
    axmatrix.xaxis.tick_top()
    plt.xticks(rotation=40, fontsize=10)

    axmatrix.set_yticks(range(data.shape[0]))
    axmatrix.set_yticklabels( rnas, minor=False)
    # axmatrix.set_yticklabels( [ rnas[i] for i in idx1 ], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(rotation=0, fontsize=10)
    # axmatrix.tight_layout()

    # Plot colorbar.
    axcolor = fig.add_axes([0.1,0.1,0.8,0.02])

    plt.colorbar(im, cax=axcolor, orientation='horizontal', norm=norm, 
                 boundaries=bounds, ticks=bounds, format=matplotlib.ticker.FuncFormatter(fmt))
    axcolor.set_xlabel('p value (-log10)')
    
    numpy.savetxt(os.path.join(path,'matrix_p.txt'), ar, delimiter='\t')
    fig.savefig(os.path.join(path,'condition_lncRNA_dendrogram.png'))
    fig.savefig(os.path.join(path,'condition_lncRNA_dendrogram.pdf'), format="pdf")

def generate_rna_exp_pv_table(root, multi_corr=True):
    "Generate p value table for Experiments vs RNA in the same project"
    
    nested_dict = lambda: defaultdict(nested_dict)
    #nested_dict = lambda: defaultdict(lambda: 'n.a.')

    data = nested_dict()
    rnas = []
    
    for item in os.listdir(root):
        pro = os.path.join(root, item, "profile.txt")
        if os.path.isfile(pro):
            with open(pro) as f:
                for line in f:
                    if line.startswith("Experiment"): continue
                    else:
                        line = line.strip().split("\t")
                        data[item][line[0]] = float(line[7])
                        rnas.append(line[0])
                        
    
    exp_list = sorted(data.keys())
    rnas = sorted(list(set(rnas)))
    
    pvs = []
    for rna in rnas:
        for exp in exp_list:
            if data[exp][rna]: pvs.append(data[exp][rna])
    reject, pvals_corrected = multiple_test_correction(pvs, alpha=0.05, method='indep')

    with open(os.path.join(root, "table_exp_rna_pv.txt"), "w") as t:
        print("\t".join(["RNA_ID"] + exp_list), file=t)
        i = 0
        for rna in rnas:
            newline = [rna]
            for exp in exp_list:
                if data[exp][rna]:
                    newline.append(str(pvals_corrected[i]))
                    i += 1
                else:
                    newline.append("n.a.")
            print("\t".join(newline), file=t)

    
    for d, p in plist.iteritems():
        list_all_index(path=os.path.dirname(p), 
                       link_d=dirlist, show_RNA_ass_gene=show_RNA_ass_gene)

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
    parser_promotertest.add_argument('-rl', type=str, default=None, metavar='  ', help="Input list for paths to all RNA sequences (in fasta format)")
    parser_promotertest.add_argument('-rn', type=str, default=None, metavar='  ', help="Define the RNA name")
    parser_promotertest.add_argument('-de', default=False, metavar='  ', help="Input file for target gene list (gene symbols or Ensembl ID)")
    parser_promotertest.add_argument('-bed', default=False, metavar='  ', help="Input BED file of the promoter regions of target genes")
    parser_promotertest.add_argument('-bg', default=False, metavar='  ', help="Input BED file of the promoter regions of background genes")
    parser_promotertest.add_argument('-o', metavar='  ', help="Output directory name for all the results")
    parser_promotertest.add_argument('-t', metavar='  ', default=False, help="Define the title name for the results under the Output name. Default is -rn.")
    
    parser_promotertest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_promotertest.add_argument('-gtf', metavar='  ', default=None, help='Define the GTF file for annotation (optional)')

    parser_promotertest.add_argument('-pl', type=int, default=1000, metavar='  ', help="Define the promotor length (Default: 1000)")
    
    parser_promotertest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_promotertest.add_argument('-score', action="store_true", help="Load score column from input gene list or BED file for analysis.")
    parser_promotertest.add_argument('-scoreh', action="store_true", help="Use the header of scores from the given gene list or BED file.")
    parser_promotertest.add_argument('-a', type=float, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (Default: 0.05)")
    parser_promotertest.add_argument('-ccf', type=int, default=40, metavar='  ', help="Define the cut off value for promoter counts (Default: 40)")
    parser_promotertest.add_argument('-rt', action="store_true", default=False, help="Remove temporary files (fa, txp...etc)")
    parser_promotertest.add_argument('-log', action="store_true", default=False, help="Set the plots in log scale")
    parser_promotertest.add_argument('-ac', type=str, default=False, metavar='  ', help="Input file for RNA accecibility ")
    parser_promotertest.add_argument('-accf', type=float, default=500, metavar='  ', help="Define the cut off value for RNA accecibility")
    parser_promotertest.add_argument('-obed', action="store_true", default=False, help="Output the BED files for DNA binding sites.")
    parser_promotertest.add_argument('-showpa', action="store_true", default=False, help="Show parallel and antiparallel bindings in the plot separately.")
    parser_promotertest.add_argument('-motif', action="store_true", default=False, help="Show motif of binding sites.")
    parser_promotertest.add_argument('-filter_havana', type=str, default="F", metavar='  ', help="Apply filtering to remove HAVANA entries.")
    parser_promotertest.add_argument('-protein_coding', type=str, default="F", metavar='  ', help="Apply filtering to get only protein coding genes.")
    parser_promotertest.add_argument('-known_only', type=str, default="F", metavar='  ', help="Apply filtering to get only known genes.")
    

    parser_promotertest.add_argument('-l', type=int, default=15, metavar='  ', help="[Triplexator] Define the minimum length of triplex (Default: 15)")
    parser_promotertest.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (Default: 20)")
    parser_promotertest.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (Default: 2)")
    parser_promotertest.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (Default: off)")
    parser_promotertest.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (Default 0).'0' = greedy approach; '1' = q-gram filtering.")
    parser_promotertest.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (Default: 1)")
    parser_promotertest.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_promotertest.add_argument('-rm', action="store_true", default=True, help="[Triplexator] Set the multiprocessing")
    
    
    ################### Genomic Region Test ##########################################
    h_region = "Genomic region test evaluates the association between the given lncRNA to the target regions by randomization."
    parser_randomtest = subparsers.add_parser('regiontest', help=h_region)
    parser_randomtest.add_argument('-r', type=str, metavar='  ', help="Input file name for RNA sequence (in fasta format)")
    parser_randomtest.add_argument('-rl', type=str, default=None, metavar='  ', help="Input list for paths to all RNA sequences (in fasta format)")
    parser_randomtest.add_argument('-rn', type=str, default=False, metavar='  ', help="Define the RNA name")
    parser_randomtest.add_argument('-bed', metavar='  ', help="Input BED file for interested regions on DNA")
    parser_randomtest.add_argument('-o', metavar='  ', help="Output directory name for all the results and temporary files")
    parser_randomtest.add_argument('-t', metavar='  ', default=False, help="Define the title name for the results under the Output name. Default is -rn.")
    
    parser_randomtest.add_argument('-n', type=int, default=10000, metavar='  ', 
                                   help="Number of times for randomization (Default: 10000)")

    parser_randomtest.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
 
    parser_randomtest.add_argument('-showdbs', action="store_true", help="Show the plots and statistics of DBS (DNA Binding sites)")
    parser_randomtest.add_argument('-score', action="store_true", help="Load score column from input BED file")
    parser_randomtest.add_argument('-a', type=float, default=0.05, metavar='  ', help="Define significance level for rejection null hypothesis (Default: 0.05)")
    parser_randomtest.add_argument('-ccf', type=int, default=40, metavar='  ', help="Define the cut off value for DBS counts (Default: 20)")
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
    ##########################################################################
    parser_bed2bed = subparsers.add_parser('get_dbss', help="Get DBSs in BED format from the single BED file")
    parser_bed2bed.add_argument('-i',type=str, metavar='  ', help='Input BED file of the target regions')
    parser_bed2bed.add_argument('-dbs',type=str, metavar='  ', help='Output BED file of the DBSs')
    parser_bed2bed.add_argument('-rbs',type=str, metavar='  ', help='Output BED file of the RBSs')
    parser_bed2bed.add_argument('-r',type=str, metavar='  ', help='Input FASTA file of the RNA')
    parser_bed2bed.add_argument('-organism', metavar='  ', help='Define the organism (hg19 or mm9)')
    parser_bed2bed.add_argument('-l', type=int, default=15, metavar='  ', help="[Triplexator] Define the minimum length of triplex (Default: 15)")
    parser_bed2bed.add_argument('-e', type=int, default=20, metavar='  ', help="[Triplexator] Set the maximal error-rate in %% tolerated (Default: 20)")
    parser_bed2bed.add_argument('-c', type=int, default=2, metavar='  ', help="[Triplexator] Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro (Default: 2)")
    parser_bed2bed.add_argument('-fr', type=str, default="off", metavar='  ', help="[Triplexator] Activates the filtering of low complexity regions and repeats in the sequence data (Default: off)")
    parser_bed2bed.add_argument('-fm', type=int, default=0, metavar='  ', help="[Triplexator] Method to quickly discard non-hits (Default 0).'0' = greedy approach; '1' = q-gram filtering.")
    parser_bed2bed.add_argument('-of', type=int, default=1, metavar='  ', help="[Triplexator] Define output formats of Triplexator (Default: 1)")
    parser_bed2bed.add_argument('-mf', action="store_true", default=False, help="[Triplexator] Merge overlapping features into a cluster and report the spanning region.")
    parser_bed2bed.add_argument('-rm', action="store_true", default=True, help="[Triplexator] Set the multiprocessing")
    
    ##########################################################################
    parser_integrate = subparsers.add_parser('integrate', help="Integrate the project's links and generate project-level statistics.")
    parser_integrate.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
    ##########################################################################
    parser_updatehtml = subparsers.add_parser('updatehtml', help="Update the project's html.")
    parser_updatehtml.add_argument('-path',type=str, metavar='  ', help='Define the path of the project.')
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
        #cf = ConfigurationFile()
        ####################################################################################
        ######### Integration
        if args.mode == "integrate":
            condition_list = [] # name, link, no. tests, no. sig.
            for item in os.listdir(args.path):
                if item == "style": continue
                if os.path.isfile(os.path.join(args.path,item)): continue
                elif os.path.isdir(os.path.join(args.path,item)):
                    h = os.path.join(item, "index.html")
                    pro = os.path.join(args.path, item, "profile.txt")
                    if os.path.isfile(pro):
                        nt = 0
                        ns = 0
                        with open(pro) as f:
                            for line in f:
                                line = line.strip().split("\t")
                                if line[0] == "Experiment": continue
                                nt += 1
                                if float(line[7]) < 0.05: ns += 1
                        # print([item, h, str(nt), str(ns)])
                        condition_list.append( [item, h, str(nt), str(ns)] )
            # print(condition_list)
            link_d = {"List":"index.html"}
            fp = condition_list[0][0] + "/style"
            html = Html(name="Directory: "+args.path, links_dict=link_d, 
                        fig_rpath=fp, #fig_dir=fp, 
                        RGT_header=False, other_logo="TDF")
            html.add_heading("All conditions in: "+args.path+"/")
            data_table = []
            type_list = 'sssssssssssss'
            col_size_list = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
            c = 0
            header_list = ["No.", "Conditions", "No. tests", "No. sig. tests" ]
            for i, exp in enumerate(condition_list):
                c += 1
                data_table.append([str(c), 
                                   '<a href="'+exp[1]+'">'+exp[0]+"</a>",
                                   exp[2], exp[3] ])
            html.add_zebra_table( header_list, col_size_list, type_list, data_table, 
                                  align=10, cell_align="left", sortable=True)
            html.add_fixed_rank_sortable()
            html.write(os.path.join(args.path,"index.html"))
            #revise_index(root=args.path, show_RNA_ass_gene=True)
            gen_heatmap(path=args.path)
            sys.exit(0)
        ####################################################################################
        ######### updatehtml
        elif args.mode == "updatehtml":
            revise_index(root=args.path, show_RNA_ass_gene=True)
            generate_rna_exp_pv_table(root=args.path, multi_corr=True)
            sys.exit(0)


        #######################################################################################
        ######## Test triplexator
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
        
        ####################################################################################
        ######### get_dbss
        if args.mode == "get_dbss":
            get_dbss(input_BED=args.i,output_BED=args.dbs,rna_fasta=args.r,output_rbss=args.rbs,
                     organism=args.organism,l=args.l,e=args.e,c=args.c,
                     fr=args.fr,fm=args.fm,of=args.of,mf=args.mf,rm=args.rm,temp=dir)
            sys.exit(0)


        #######################################################################
        #### Checking arguments
        if not args.o: 
            print("Please define the output diractory name. \n")
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
                    #print(ps)
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
                    
                    if args.rm: command += ["-rm" ]
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
        def no_binding_code():
            print("*** Find no triple helices binding on the given RNA")

            pro_path = os.path.join(os.path.dirname(args.o), "profile.txt")
            exp = os.path.basename(args.o)
            if args.de: tar_reg = os.path.basename(args.de)
            else: tar_reg = os.path.basename(args.bed)
            r_genes = rna_associated_gene(rna_regions=promoter.rna_regions, name=promoter.rna_name, organism=promoter.organism)
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
            #list_all_index(path=os.path.dirname(args.o), show_RNA_ass_gene=promoter.rna_regions)
            revise_index(root=os.path.dirname(os.path.dirname(args.o)), show_RNA_ass_gene=promoter.rna_regions)
            shutil.rmtree(args.o)
            sys.exit(1)

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

        args.r = os.path.normpath(os.path.join(dir,args.r))
        
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
        promoter.get_rna_region_str(rna=args.r)
        promoter.connect_rna(rna=args.r, temp=args.o)
        promoter.search_triplex(temp=args.o, l=args.l, e=args.e, remove_temp=args.rt, 
                                c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf)
        
        t1 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t1-t0))))

        print2(summary, "Step 2: Calculate the frequency of DNA binding sites within the promotors.")
        if args.obed: obedp = os.path.basename(args.o)
        else: obedp = None
        promoter.count_frequency(temp=args.o, remove_temp=args.rt, obedp=obedp, cutoff=args.ccf, l=args.l)
        promoter.fisher_exact(alpha=args.a)
        t2 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t2-t1))))
        
        if len(promoter.rbss) == 0:  no_binding_code()

        promoter.dbd_regions(sig_region=promoter.sig_region_promoter, output=args.o)
        os.remove(os.path.join(args.o,"rna_temp.fa"))
        try: os.remove(os.path.join(args.o,"rna_temp.fa.fai"))
        except: pass
        print2(summary, "Step 3: Establishing promoter profile.")
        #promoter.promoter_profile()
        t3 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t3-t2))))

        print2(summary, "Step 4: Generate plot and output html files.")
        promoter.plot_lines(txp=promoter.txp_def, rna=args.r, dirp=args.o, ac=args.ac, 
                            cut_off=args.accf, log=args.log, showpa=args.showpa,
                            sig_region=promoter.sig_region_promoter,
                            ylabel="Number of DBSs", 
                            linelabel="No. DBSs", filename="plot_promoter.png")

        promoter.barplot(dirp=args.o, filename="bar_promoter.png", sig_region=promoter.sig_region_promoter)
        #if args.showdbs:
        #    promoter.plot_lines(txp=promoter.txp_def, rna=args.r, dirp=args.o, ac=args.ac, 
        #                        cut_off=args.accf, log=args.log, showpa=args.showpa,
        #                        sig_region=promoter.sig_region_dbs,
        #                        ylabel="Number of DBSs on target promoters", 
        #                        linelabel="No. DBSs", filename="plot_dbss.png")
        #    promoter.barplot(dirp=args.o, filename="bar_dbss.png", sig_region=promoter.sig_region_dbs, dbs=True)
        if args.motif: promoter.gen_motifs(temp=args.o)

        
        promoter.gen_html(directory=args.o, parameters=args, ccf=args.ccf, align=50, alpha=args.a)
        promoter.gen_html_genes(directory=args.o, align=50, alpha=args.a, nonDE=False)
        promoter.save_table(path=os.path.dirname(args.o), table=promoter.ranktable, 
                                filename="lncRNA_target_ranktable.txt")
        promoter.save_table(path=os.path.dirname(args.o), table=promoter.dbstable, 
                                filename="lncRNA_target_dbstable.txt")

        #promoter.heatmap(table="ranktable.txt", temp=os.path.dirname(args.o))

        t4 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t4-t3))))
        print2(summary, "\nTotal running time is: " + str(datetime.timedelta(seconds=round(t4-t0))))
    
        output_summary(summary, args.o, "summary.txt")
        promoter.save_profile(output=args.o, bed=args.bed, geneset=args.de)
        #list_all_index(path=os.path.dirname(args.o), show_RNA_ass_gene=promoter.rna_regions)
        revise_index(root=os.path.dirname(os.path.dirname(args.o)), show_RNA_ass_gene=promoter.rna_regions)
        try: os.remove(os.path.join(args.o, "de.fa"))
        except OSError: pass
        try: os.remove(os.path.join(args.o, "nde.fa"))
        except OSError: pass


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
                              c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, cutoff=args.ccf )
        t1 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t1-t0))))

        print2(summary, "Step 2: Randomization and counting number of binding sites")
        randomtest.random_test(repeats=args.n, temp=args.o, remove_temp=args.rt, l=args.l, e=args.e,
                               c=args.c, fr=args.fr, fm=args.fm, of=args.of, mf=args.mf, rm=args.rm,
                               filter_bed=args.f, alpha=args.a)
        
        if len(randomtest.rbss) == 0: 
            no_binding_code()

        t2 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t2-t1))))
        
        print2(summary, "Step 3: Generating plot and output HTML")
        randomtest.dbd_regions(sig_region=randomtest.data["region"]["sig_region"], output=args.o)
        
        os.remove(os.path.join(args.o,"rna_temp.fa"))
        try: os.remove(os.path.join(args.o,"rna_temp.fa.fai"))
        except: pass
        

        randomtest.lineplot(txp=randomtest.txpf, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
                            log=args.log, ylabel="Number of DBS",
                            sig_region=randomtest.data["region"]["sig_region"], 
                            linelabel="No. DBS", filename="lineplot_region.png")

        #randomtest.lineplot(txp=randomtest.txp, dirp=args.o, ac=args.ac, cut_off=args.accf, showpa=args.showpa,
        #                    log=args.log, ylabel="Number of target regions with DBS", 
        #                    sig_region=randomtest.data["region"]["sig_region"],
        #                    linelabel="No. target regions", filename="lineplot_region.png")
        
        randomtest.boxplot(dir=args.o, matrix=randomtest.region_matrix, 
                           sig_region=randomtest.data["region"]["sig_region"], 
                           truecounts=[r[0] for r in randomtest.counts_tr.values()],
                           sig_boolean=randomtest.data["region"]["sig_boolean"], 
                           ylabel="Number of target regions",
                           filename="boxplot_regions" )
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

        randomtest.gen_html(directory=args.o, parameters=args, align=50, alpha=args.a, 
                            score=args.score, obed=obed)

        t3 = time.time()
        print2(summary, "\tRunning time is: " + str(datetime.timedelta(seconds=round(t3-t2))))

        print2(summary, "\nTotal running time is: " + str(datetime.timedelta(seconds=round(t3-t0))))

        output_summary(summary, args.o, "summary.txt")
        randomtest.save_profile(output=args.o, bed=args.bed)
        list_all_index(path=os.path.dirname(args.o), show_RNA_ass_gene=False)
        for f in os.listdir(args.o):
            if re.search("dna*.fa", f) or re.search("dna*.txp", f):
                os.remove(os.path.join(args.o, f))
