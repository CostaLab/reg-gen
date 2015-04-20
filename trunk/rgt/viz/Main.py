# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import argparse
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import time, datetime, getpass, fnmatch
from shutil import copyfile

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.Util import GenomeData, OverlapType, Html
from plotTools import Projection, Jaccard, Intersect, Boxplot, Lineplot

dir = os.getcwd()
"""
Statistical analysis methods and plotting tools for ExperimentalMatrix

Author: Joseph Kuo
"""

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################
    

######### Universal functions
def print2(parameter, string):
    """ Show the message on the console and also save in a list for future backup. """
    print(string)
    parameter.append(string)
    
######### Projection test

def output(f, directory, folder, filename, extra=None, pdf=False, show=None):
    """Output the file in the defined folder """
    pd = os.path.normpath(os.path.join(dir,directory,folder))
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)    
    
    # Saving 
    if not extra:
        f.savefig(os.path.join(pd,filename), facecolor='w', edgecolor='w', bbox_inches='tight', dpi=300)
    else:
        f.savefig(os.path.join(pd,filename), facecolor='w', edgecolor='w', bbox_extra_artists=(extra), bbox_inches='tight',dpi=300)
    
    if pdf:
        try:
            pp = PdfPages(os.path.join(pd,filename) + '.pdf')
            pp.savefig(f, bbox_extra_artists=(extra),bbox_inches='tight') 
            pp.close()
        except:
            print("ERROR: Problem in PDF conversion. Skipped.")
    if show:
        plt.show()
        
def output_parameters(parameter, directory, folder, filename):
    pd = os.path.join(dir,directory,folder)
    try: os.stat(os.path.dirname(pd))
    except: os.mkdir(os.path.dirname(pd))
    try: os.stat(pd)
    except: os.mkdir(pd)    
    if parameter:
        with open(os.path.join(pd,"parameters.txt"),'w') as f:
            for s in parameter:
                print(s, file=f)

def copy_em(em, directory, folder, filename="experimental_matrix.txt"):
    copyfile(em, os.path.join(dir,directory,folder,filename))

def list_all_index(path):
    """Creat an 'index.html' in the defined directory """
    dirname = os.path.basename(path)
    
    link_d = {"List":"index.html"}
    html = Html(name="Directory: "+dirname, links_dict=link_d, 
                fig_dir=os.path.join(path,"style"), fig_rpath="./style", RGT_header=False, other_logo="viz")
    header_list = ["No.", "Experiments"]
    html.add_heading("All experiments in: "+dirname+"/")
    data_table = []
    type_list = 'ssss'
    col_size_list = [10, 10, 10]
    c = 0
    for root, dirnames, filenames in os.walk(path):
        #roots = root.split('/')
        for filename in fnmatch.filter(filenames, '*.html'):
            if filename == 'index.html' and root.split('/')[-1] != dirname:
                c += 1
                data_table.append([str(c), '<a href="'+os.path.join(root.split('/')[-1], filename)+'"><font size='+'"4"'+'>'+root.split('/')[-1]+"</a>"])
                #print(link_d[roots[-1]])
    html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=50, cell_align="left", sortable=True)
    html.add_fixed_rank_sortable()
    html.write(os.path.join(path,"index.html"))

def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try: os.stat(path)
    except: os.mkdir(path)

def main():
    ###############################################################################
    ##### PARAMETERS ##############################################################
    ###############################################################################
    
    # Some general help descriptions
    ######### Some general plotting arguments descriptions ###############
    helpinput = 'The file name of the input Experimental Matrix file. Recommended to add more columns for more information for ploting. For example, cell type or factors.'
    helpoutput = 'The directory name for the output files. For example, project name.'
    helptitle = 'The title shown on the top of the plot and also the folder name.'
    helpgroup = "Group the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."
    helpgroupbb = "Group the data by any optional column (for example, 'cell') of experimental matrix, or None."
    helpsort = "Sort the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."
    helpcolor = "Color the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."
    helpcolorbb = "Color the data by any optional column (for example, 'cell') of experimental matrix, or None."
    helpDefinedColot = 'Define the specific colors with the given column "color" in experimental matrix. The color should be in the format of matplotlib.colors. For example, "r" for red, "b" for blue, or "(100, 35, 138)" for RGB.'
    helpreference = 'The file name of the reference Experimental Matrix. Multiple references are acceptable.'
    helpquery = 'The file name of the query Experimental Matrix. Multiple queries are acceptable.'
    
    parser = argparse.ArgumentParser(description='Provides various Statistical analysis methods and plotting tools for ExperimentalMatrix.\
    \nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Projection test ##########################################
    parser_projection = subparsers.add_parser('projection',help='Projection test evaluates the association level by comparing to the random binomial model.')
    parser_projection.add_argument('-r', metavar='  ', help=helpreference)
    parser_projection.add_argument('-q', metavar='  ', help=helpquery)
    parser_projection.add_argument('-o', metavar='  ', help=helpoutput) 
    parser_projection.add_argument('-t', metavar='  ', default='projection_test', help=helptitle)
    parser_projection.add_argument('-g', metavar='  ', default=None, help=helpgroupbb +" (Default:None)")
    parser_projection.add_argument('-c', metavar='  ', default="regions", help=helpcolorbb +' (Default: regions)')
    parser_projection.add_argument('-bg', metavar='  ', default=None, help="Define a BED file as background. If not defined, the background is whole genome according to the given organism.")
    parser_projection.add_argument('-union', action="store_true", help='Take the union of references as background for binominal test.')
    parser_projection.add_argument('-organism', metavar='  ', default='hg19', help='Define the organism. (Default: hg19)')
    parser_projection.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
    parser_projection.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_projection.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_projection.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
    
    ################### Intersect Test ##########################################
    parser_intersect = subparsers.add_parser('intersect',help='Intersection test provides various modes of intersection to test the association between references and queries.')
    parser_intersect.add_argument('-r', metavar='  ', help=helpreference)
    parser_intersect.add_argument('-q', metavar='  ', help=helpquery)
    parser_intersect.add_argument('-o', help=helpoutput)
    parser_intersect.add_argument('-t', metavar='  ', default='intersection_test', help=helptitle)
    parser_intersect.add_argument('-g', metavar='  ', default=None, help=helpgroupbb +" (Default:None)")
    parser_intersect.add_argument('-c', metavar='  ', default="regions", help=helpcolorbb +' (Default: regions)')
    parser_intersect.add_argument('-organism', metavar='  ', default='hg19', help='Define the organism. (Default: hg19)')
    parser_intersect.add_argument('-bg', metavar='  ', help="Define a BED file as background. If not defined, the background is whole genome according to the given organism.")
    parser_intersect.add_argument('-m', metavar='  ', default="count", choices=['count','bp'],
                                  help="Define the mode of calculating intersection. \
                                  'count' outputs the number of overlapped regions.\
                                  'bp' outputs the coverage(basepair) of intersection.")
    parser_intersect.add_argument('-tc', metavar='  ', type=int, default=False, help="Define the threshold(in percentage) of reference length for intersection counting. For example, '20' means that the query which overlaps more than 20%% of reference is counted as intersection.")
    parser_intersect.add_argument('-ex', metavar='  ', type=int, default=0, help="Define the extension(in percentage) of reference length for intersection counting. For example, '20' means that each region of reference is extended by 20%% in order to include proximal queries.")
    parser_intersect.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
    parser_intersect.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_intersect.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_intersect.add_argument('-stest', metavar='  ', type=int, default= 0, help='Define the repetition time of random subregion test between reference and query.')
    
    ################### Jaccard test ##########################################
    
    parser_jaccard = subparsers.add_parser('jaccard',help='Jaccard test evaluates the association level by comparing with jaccard index from repeating randomization.')
    
    parser_jaccard.add_argument('-o', help=helpoutput) 
    parser_jaccard.add_argument('-r', '--reference',help=helpreference)
    parser_jaccard.add_argument('-q', '--query', help=helpquery)
    parser_jaccard.add_argument('-t','--title', default='jaccard_test', help=helptitle)
    parser_jaccard.add_argument('-rt','--runtime', type=int, default=500, help='Define how many times to run the randomization. (Default:500)')
    parser_jaccard.add_argument('-g', default=None, help=helpgroupbb +" (Default:None)")
    parser_jaccard.add_argument('-c', default="regions", help=helpcolorbb +' (Default: regions)')
    parser_jaccard.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')
    parser_jaccard.add_argument('-nlog', action="store_false", help='Set y axis of the plot not in log scale.')
    parser_jaccard.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_jaccard.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_jaccard.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')

    ################### Combinatorial Test ##########################################
    parser_combinatorial = subparsers.add_parser('combinatorial',help='Combinatorial test compare all combinatorial possibilities from reference to test the association between references and queries.')
    
    parser_combinatorial.add_argument('-o', help=helpoutput)
    parser_combinatorial.add_argument('-r', '--reference',help=helpreference)
    parser_combinatorial.add_argument('-q', '--query', help=helpquery)
    parser_combinatorial.add_argument('-t','--title', default='combinatorial_test', help=helptitle)
    parser_combinatorial.add_argument('-g', default=None, help=helpgroupbb +" (Default:None)")
    parser_combinatorial.add_argument('-c', default="regions", help=helpcolorbb +' (Default: regions)')
    parser_combinatorial.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')
    parser_combinatorial.add_argument('-bg', help="Define a BED file as background. If not defined, the background is whole genome according to the given organism.")
    parser_combinatorial.add_argument('-m', default="count", choices=['count','bp'],
                                      help="Define the mode of calculating intersection. \
                                      'count' outputs the number of overlapped regions.\
                                      'bp' outputs the coverage(basepair) of intersection.")
    parser_combinatorial.add_argument('-tc', type=int, default=False, help="Define the threshold(in percentage) of reference length for intersection counting. For example, '20' means that the query which overlaps more than 20%% of reference is counted as intersection.")
    parser_combinatorial.add_argument('-ex', type=int, default=0, help="Define the extension(in percentage) of reference length for intersection counting. For example, '20' means that each region of reference is extended by 20%% in order to include proximal queries.")
    parser_combinatorial.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
    parser_combinatorial.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_combinatorial.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_combinatorial.add_argument('-stest', type=int, default= 0, help='Define the repetition time of random subregion test between reference and query.')
    
    ################### Boxplot ##########################################
    
    parser_boxplot = subparsers.add_parser('boxplot',help='Boxplot based on the BAM and BED files for gene association analysis.')
    parser_boxplot.add_argument('input',help=helpinput)
    parser_boxplot.add_argument('-o', metavar='  ', help=helpoutput)
    parser_boxplot.add_argument('-t', metavar='  ', default='boxplot', help=helptitle)
    parser_boxplot.add_argument('-g', metavar='  ', default='reads', help=helpgroup + " (Default:reads)")
    parser_boxplot.add_argument('-c', metavar='  ', default='regions', help=helpcolor + " (Default:regions)")
    parser_boxplot.add_argument('-s', metavar='  ', default='None', help=helpsort + " (Default:None)")
    parser_boxplot.add_argument('-sy', action="store_true", help="Share y axis for convenience of comparison.")
    parser_boxplot.add_argument('-nlog', action="store_false", help='Set y axis of the plot not in log scale.')
    parser_boxplot.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_boxplot.add_argument('-nqn', action="store_true", help='No quantile normalization in calculation.')
    parser_boxplot.add_argument('-df', action="store_true", help="Show the difference of the two signals which share the same labels.The result is the subtraction of the first to the second.")
    parser_boxplot.add_argument('-ylim', metavar='  ', type=int, default=None, help="Define the limit of y axis.")
    parser_boxplot.add_argument('-p', metavar='  ', type=float, default=0.05, help='Define the significance level for multiple test. Default: 0.01')
    parser_boxplot.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_boxplot.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
    
    ################### Lineplot ##########################################
    parser_lineplot = subparsers.add_parser('lineplot', help='Generate lineplot with various modes.')
    
    choice_center = ['midpoint','leftend','rightend','bothends'] 
    # Be consist as the arguments of GenomicRegionSet.relocate_regions
    
    parser_lineplot.add_argument('input', help=helpinput)
    parser_lineplot.add_argument('-o', help=helpoutput)
    parser_lineplot.add_argument('-ga', action="store_true", help="Use genetic annotation data as input regions (e.g. TSS, TTS, exons and introns) instead of the BED files in the input matrix.")
    parser_lineplot.add_argument('-t', metavar='  ', default='lineplot', help=helptitle)
    parser_lineplot.add_argument('-center', metavar='  ', choices=choice_center, default='midpoint', 
                                 help='Define the center to calculate coverage on the regions. Options are: '+', '.join(choice_center) + 
                                 '.(Default:midpoint) The bothend mode will flap the right end region for calculation.')
    parser_lineplot.add_argument('-g', metavar='  ', default='reads', help=helpgroup + " (Default:reads)")
    parser_lineplot.add_argument('-c', metavar='  ', default='regions', help=helpcolor + " (Default:regions)")
    parser_lineplot.add_argument('-s', metavar='  ', default='None', help=helpsort + " (Default:None)")
    parser_lineplot.add_argument('-e', metavar='  ', type=int, default=2000, help='Define the extend length of interested region for plotting.(Default:2000)')
    parser_lineplot.add_argument('-rs', metavar='  ', type=int, default=200, help='Define the readsize for calculating coverage.(Default:200)')
    parser_lineplot.add_argument('-ss', metavar='  ', type=int, default=50, help='Define the stepsize for calculating coverage.(Default:50)')
    parser_lineplot.add_argument('-bs', metavar='  ', type=int, default=100, help='Define the binsize for calculating coverage.(Default:100)')
    parser_lineplot.add_argument('-sy', action="store_true", help="Share y axis for convenience of comparison.")
    parser_lineplot.add_argument('-sx', action="store_true", help="Share x axis for convenience of comparison.")
    parser_lineplot.add_argument('-organism', metavar='  ', default='hg19', help='Define the organism. (Default: hg19)')
    parser_lineplot.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_lineplot.add_argument('-mp', action="store_true", help="Perform multiprocessing for faster computation.")
    parser_lineplot.add_argument('-df', action="store_true", help="Show the difference of the two signals which share the same labels.The result is the subtraction of the first to the second.")
    parser_lineplot.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_lineplot.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
    
    ################### Heatmap ##########################################
    parser_heatmap = subparsers.add_parser('heatmap', help='Generate heatmap with various modes.')
    
    choice_center = ['midpoint','leftend','rightend','bothends'] 
    # Be consist as the arguments of GenomicRegionSet.relocate_regions
    
    parser_heatmap.add_argument('input', help=helpinput)
    parser_heatmap.add_argument('-o', metavar='  ', help=helpoutput)
    parser_heatmap.add_argument('-ga', action="store_true", help="Use genetic annotation data as input regions (e.g. TSS, TTS, exons and introns) instead of the BED files in the input matrix.")
    parser_heatmap.add_argument('-t', metavar='  ', default='heatmap', help=helptitle)
    parser_heatmap.add_argument('-center', metavar='  ', choices=choice_center, default='midpoint', 
                                 help='Define the center to calculate coverage on the regions. Options are: '+', '.join(choice_center) + 
                                 '.(Default:midpoint) The bothend mode will flap the right end region for calculation.')
    parser_heatmap.add_argument('-sort', metavar='  ', type=int, default=None, help='Define the way to sort the signals.'+
                                'Default is no sorting at all, the signals arrange in the order of their position; '+
                                '"0" is sorting by the average ranking of all signals; '+
                                '"1" is sorting by the ranking of 1st column; "2" is 2nd and so on... ')
    parser_heatmap.add_argument('-s', metavar='  ', default='None', help=helpsort + " (Default:None)")
    parser_heatmap.add_argument('-g', metavar='  ', default='regions', help=helpgroup + " (Default:regions)")
    parser_heatmap.add_argument('-c', metavar='  ', default='reads', help=helpcolor + " (Default:reads)")
    parser_heatmap.add_argument('-e', metavar='  ', type=int, default=2000, help='Define the extend length of interested region for plotting.(Default:2000)')
    parser_heatmap.add_argument('-rs', metavar='  ', type=int, default=200, help='Define the readsize for calculating coverage.(Default:200)')
    parser_heatmap.add_argument('-ss', metavar='  ', type=int, default=50, help='Define the stepsize for calculating coverage.(Default:50)')
    parser_heatmap.add_argument('-bs', metavar='  ', type=int, default=100, help='Define the binsize for calculating coverage.(Default:100)')
    parser_heatmap.add_argument('-organism', metavar='  ', default='hg19', help='Define the organism. (Default: hg19)')
    parser_heatmap.add_argument('-color', action="store_true", help=helpDefinedColot)
    parser_heatmap.add_argument('-log', action="store_true", help='Set colorbar in log scale.')
    parser_heatmap.add_argument('-mp', action="store_true", help="Perform multiprocessing for faster computation.")
    parser_heatmap.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_heatmap.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
    
    ################### Integration ##########################################
    parser_integration = subparsers.add_parser('integration', help='Provides some tools to deal with experimental matrix or other purposes.')
    parser_integration.add_argument('-ihtml', action="store_true", help='Integrate all the html files within the given directory and generate index.html for all plots.')
    parser_integration.add_argument('-l2m', help='Convert a given file list in txt format into a experimental matrix.')
    parser_integration.add_argument('-o', help='Define the folder of the output file.') 
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
            print("** Error: Please define the output directory (-o).")
            sys.exit(1)
        
        t0 = time.time()
        # Normalised output path
        args.o = os.path.normpath(os.path.join(dir,args.o))
        check_dir(args.o)
        check_dir(os.path.join(args.o, args.t))
        
        # Input parameters dictionary
        parameter = []
        parameter.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        parameter.append("User: " + getpass.getuser())
        parameter.append("\nCommand:\n\t$ " + " ".join(sys.argv))

        #################################################################################################
        ##### Main #####################################################################################
        #################################################################################################

        ################### Projection test ##########################################
        if args.mode == 'projection':
            # Fetching reference and query EM
            print2(parameter, "\n############# Projection Test #############")
            print2(parameter, "\tReference:        "+args.r)
            print2(parameter, "\tQuery:            "+args.q)
            print2(parameter, "\tOutput directory: "+os.path.basename(args.o))
            print2(parameter, "\tExperiment title: "+args.t)

            projection = Projection( args.r, args.q )
            projection.group_refque(args.g)
            projection.colors( args.c, args.color )
            if args.bg: projection.background(args.bg)
            if args.union: 
                projection.ref_union()
                projection.projection_test(organism = args.organism)
                print2(parameter, "\tTaking intersect of references as the background. ")
            else:
                projection.projection_test(organism = args.organism)
            
            # generate pdf
            projection.plot(args.log)
            output(f=projection.fig, directory = args.o, folder = args.t, filename="projection_test",
                   extra=plt.gci(),pdf=True,show=args.show)
            
            # generate html 
            projection.gen_html(args.o, args.t, args=args)
            
            if args.table:
                projection.table(directory = args.o, folder = args.t)
                
            print("\nAll related files are saved in:  "+ os.path.join(os.path.basename(args.o),args.t))
            t1 = time.time()
            print2(parameter,"\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            output_parameters(parameter, directory = args.o, folder = args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="query_experimental_matrix.txt")
            list_all_index(path=args.o)
            
        ################### Intersect Test ##########################################
        if args.mode == 'intersect':
            print2(parameter, "\n############ Intersection Test ############")
            print2(parameter, "\tReference:        "+args.r)
            print2(parameter, "\tQuery:            "+args.q)
            print2(parameter, "\tOutput directory: "+os.path.basename(args.o))
            print2(parameter, "\tExperiment title: "+args.t)
            # Fetching reference and query EM
            inter = Intersect(args.r,args.q, mode_count=args.m, organism=args.organism)
            # Setting background
            inter.background(args.bg)
            # Grouping
            inter.group_refque(args.g)
            
            # Extension
            if args.ex == 0: pass
            elif args.ex > 0: inter.extend_ref(args.ex)
            elif args.ex < 0: 
                print("\n**** extension percentage(-ex) should be positive value, not negative.\n")
                sys.exit(1)
            
            inter.colors(args.c, args.color)
            inter.count_intersect(threshold=args.tc)
            
            # generate pdf
            inter.barplot(logt=args.log)
            output(f=inter.bar, directory = args.o, folder = args.t, filename="intersection_bar",
                   extra=plt.gci(), pdf=True,show=args.show)
            inter.stackedbar()
            output(f=inter.sbar, directory = args.o, folder = args.t, filename="intersection_stackedbar",
                   extra=plt.gci(), pdf=True,show=args.show)
            inter.barplot(logt=args.log, percentage=True)
            output(f=inter.bar, directory = args.o, folder = args.t, filename="intersection_barp",
                   extra=plt.gci(), pdf=True,show=args.show)
            
            if args.stest > 0:
                inter.stest(repeat=args.stest,threshold=args.tc)
            
            # generate html
            inter.gen_html(args.o, args.t, align=50, args=args)
            
            t1 = time.time()
            print2(parameter, "\nAll related files are saved in:  "+ os.path.join(os.path.basename(args.o),args.t))
            print2(parameter,"\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            output_parameters(parameter, directory = args.o, folder = args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="query_experimental_matrix.txt")
            list_all_index(path=args.o)
            ################### Jaccard test ##########################################
        if args.mode == "jaccard":
            """Return the jaccard test of every possible comparisons between two ExperimentalMatrix. 
            
            Method:
            The distribution of random jaccard index is calculated by randomizing query for given times. 
            Then, we compare the real jaccard index to the distribution and formulate p-value as 
            p-value = (# random jaccard > real jaccard)/(# random jaccard)
            
            """
            print("\n############## Jaccard Test ###############")
            jaccard = Jaccard(args.reference,args.query)
            jaccard.group_refque(args.g)
            jaccard.colors(args.c, args.color)
            
            # jaccard test
            jaccard.jaccard_test(args.runtime, args.organism)
            parameter = parameter + jaccard.parameter
            t1 = time.time()
            # ploting and generate pdf
            jaccard.plot(logT=args.nlog)
            for i,f in enumerate(jaccard.fig):
                output(f=f, directory = args.o, folder = args.title, filename="jaccard_test"+str(i+1),extra=plt.gci(),pdf=True,show=args.show)
            # generate html
            jaccard.gen_html(args.o, args.title)
            
            if args.table:
                jaccard.table(directory = args.o, folder = args.title)
            
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.title))
            print2(parameter,"\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            output_parameters(parameter, directory = args.o, folder = args.title, filename="parameters.txt")
            copy_em(em=args.reference, directory=args.o, folder=args.title, filename="Reference_experimental_matrix.txt")
            copy_em(em=args.query, directory=args.o, folder=args.title, filename="Query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ################### Combinatorial Test ##########################################
        if args.mode == 'combinatorial':
            print("\n############ Combinatorial Test ############")
            # Fetching reference and query EM
            inter = Intersect(args.reference,args.query, mode_count=args.m, organism=args.organism)
            # Setting background
            inter.background(args.bg)
            # Grouping
            inter.group_refque(args.g)
            
            # Extension
            if args.ex == 0: pass
            elif args.ex > 0: inter.extend_ref(args.ex)
            elif args.ex < 0: 
                print("\n**** extension percentage(-ex) should be positive value, not negative.\n")
                sys.exit(1)
            # Combinatorial 
            print2(parameter, "Generating all combinatorial regions for further analysis...")
            inter.combinatorial()
            inter.count_intersect(threshold=args.tc, frequency=True)
            
            # generate pdf
            inter.colors_comb()
            #inter.barplot(args.log)
            #output(f=inter.bar, directory = args.output, folder = args.title, filename="intersection_bar",extra=plt.gci(),pdf=True,show=args.show)
            #if args.stackedbar:
            #inter.colors(args.c, args.color,ref_que = "ref")
            inter.comb_stacked_plot()
            output(f=inter.sbar, directory = args.o, folder = args.title, filename="intersection_stackedbar",extra=plt.gci(),pdf=True,show=args.show)
            #if args.lineplot:
            #    inter.comb_lineplot()
            if args.stest > 0:
                inter.stest(repeat=args.stest,threshold=args.tc)
            # generate html
            inter.gen_html_comb(args.o, args.title, align=50)
            
            parameter = parameter + inter.parameter
            t1 = time.time()
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.title))
            print2(parameter,"\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            output_parameters(parameter, directory = args.o, folder = args.title, filename="parameters.txt")
            copy_em(em=args.reference, directory=args.o, folder=args.title, filename="Reference_experimental_matrix.txt")
            copy_em(em=args.query, directory=args.o, folder=args.title, filename="Query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ################### Proportion Plot ##########################################
        if args.mode == 'proportionplot':
            print("\n############ Proportion Plot ############")
            # Fetching reference and query EM
            
            
            parameter = parameter + inter.parameter
            t1 = time.time()
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.title))
            print2(parameter,"\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1-t0))))
            output_parameters(parameter, directory = args.o, folder = args.title, filename="parameters.txt")
            copy_em(em=args.reference, directory=args.o, folder=args.title, filename="Reference_experimental_matrix.txt")
            copy_em(em=args.query, directory=args.o, folder=args.title, filename="Query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ################### Boxplot ##########################################
        if args.mode == 'boxplot':
            print("\n################# Boxplot #################")
            boxplot = Boxplot(args.input, title=args.t, df=args.df)
            
            print2(parameter,"\nStep 1/5: Combining all regions")
            boxplot.combine_allregions()
            print2(parameter,"    " + str(len(boxplot.all_bed)) + " regions from all bed files are combined.")
            t1 = time.time()
            print2(parameter,"    --- finished in {0} secs\n".format(round(t1-t0)))
            
            # Coverage of reads on all_bed
            print2(parameter,"Step 2/5: Calculating coverage of each bam file on all regions")
            boxplot.bedCoverage() 
            t2 = time.time()
            print2(parameter,"    --- finished in {0} (H:M:S)\n".format(datetime.timedelta(seconds=round(t2-t1))))
            
            # Quantile normalization
            print2(parameter,"Step 3/5: Quantile normalization of all coverage table")
            if args.nqn:
                print2(parameter,"    No quantile normalization.")
                boxplot.norm_table = boxplot.all_table
            else: boxplot.quantile_normalization()
            t3 = time.time()
            print2(parameter,"    --- finished in {0} secs\n".format(round(t3-t2)))
            
            # Generate individual table for each bed
            print2(parameter,"Step 4/5: Constructing different tables for box plot")
            boxplot.tables_for_plot()
            if args.table: boxplot.print_plot_table(directory = args.o, folder = args.t)
            t4 = time.time()
            print2(parameter,"    --- finished in {0} secs\n".format(round(t4-t3)))
            
            # Plotting
            print2(parameter,"Step 5/5: Plotting")
            boxplot.group_tags(groupby=args.g, sortby=args.s, colorby=args.c)
            
            boxplot.group_data(directory = args.o, folder = args.t, log=args.nlog)
            boxplot.color_map(colorby=args.c, definedinEM=args.color)
            boxplot.plot(title=args.t, logT=args.nlog, sy=args.sy, ylim=args.ylim)
            #if args.table: boxplot.print_table(directory=args.output, folder=args.title)
            output(f=boxplot.fig, directory = args.o, folder = args.t, filename="boxplot",extra=plt.gci(),pdf=True,show=args.show)
            # HTML
            boxplot.gen_html(args.o, args.t, align = 50)
            t5 = time.time()
            print2(parameter,"    --- finished in {0} secs\n".format(round(t5-t4)))
            print2(parameter,"Total running time is: " + str(datetime.timedelta(seconds=round(t5-t0))) + " (H:M:S)\n")
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.t))
            output_parameters(parameter, directory = args.o, folder = args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)

        ################### Lineplot #########################################
        if args.mode == 'lineplot':
            if args.sy and args.sx:
                print("** Err: -sy and -sx cannot be used simutaneously.")
                sys.exit(1)

            print("\n################ Lineplot #################")
            # Read experimental matrix
            t0 = time.time()
            if "reads" not in (args.g, args.c, args.s):
                print("Please add 'reads' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            if "regions" not in (args.g, args.c, args.s):
                print("Please add 'regions' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            print2(parameter, "Parameters:\tExtend length:\t"+str(args.e))
            print2(parameter, "\t\tRead size:\t"+str(args.rs))
            print2(parameter, "\t\tBin size:\t"+str(args.bs))
            print2(parameter, "\t\tStep size:\t"+str(args.ss))
            print2(parameter, "\t\tCenter mode:\t"+str(args.center+"\n"))
            
            lineplot = Lineplot(EMpath=args.input, title=args.t, annotation=args.ga, 
                                organism=args.organism, center=args.center, extend=args.e, rs=args.rs, 
                                bs=args.bs, ss=args.ss, df=args.df)
            # Processing the regions by given parameters
            print2(parameter, "Step 1/3: Processing regions by given parameters")
            lineplot.relocate_bed()
            t1 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t1-t0))))
            
            if args.mp: print2(parameter, "\nStep 2/3: Calculating the coverage to all reads and averaging with multiprocessing ")
            else: print2(parameter, "\nStep 2/3: Calculating the coverage to all reads and averaging")
            lineplot.group_tags(groupby=args.g, sortby=args.s, colorby=args.c)
            lineplot.gen_cues()
            lineplot.coverage(sortby=args.s, mp=args.mp)
            t2 = time.time()
            print2(parameter, "    --- finished in {0} (H:M:S)".format(str(datetime.timedelta(seconds=round(t2-t1)))))
            
            # Plotting
            print2(parameter, "\nStep 3/3: Plotting the lineplots")
            lineplot.colormap(colorby = args.c, definedinEM = args.color)
            lineplot.plot(groupby=args.g, colorby=args.c, output=args.o, printtable=args.table, sy=args.sy, sx=args.sx)
            output(f=lineplot.fig, directory = args.o, folder = args.t, filename="lineplot",extra=plt.gci(),pdf=True,show=args.show)
            lineplot.gen_html(args.o, args.t)
            t3 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t3-t2))))
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t3-t0))) + "(H:M:S)\n")
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.t))
            output_parameters(parameter, directory = args.o, folder = args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)

        ################### Heatmap ##########################################
        if args.mode=='heatmap':
            print("\n################# Heatmap #################")
            # Most part of heat map are the same as lineplot, so it share the same class as lineplot
            # Read experimental matrix
            t0 = time.time()
            if "reads" not in (args.g, args.c, args.s):
                print("Please add 'reads' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            if "regions" not in (args.g, args.c, args.s):
                print("Please add 'regions' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            print2(parameter, "Parameters:\tExtend length:\t"+str(args.e))
            print2(parameter, "\t\tRead size:\t"+str(args.rs))
            print2(parameter, "\t\tBin size:\t"+str(args.bs))
            print2(parameter, "\t\tStep size:\t"+str(args.ss))
            print2(parameter, "\t\tCenter mode:\t"+str(args.center+"\n"))
        
            lineplot = Lineplot(EMpath=args.input, title=args.t, annotation=args.ga, 
                                organism=args.organism, center=args.center, extend=args.e, rs=args.rs, 
                                bs=args.bs, ss=args.ss, df=False)
            # Processing the regions by given parameters
            print2(parameter, "Step 1/4: Processing regions by given parameters")
            lineplot.relocate_bed()
            t1 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t1-t0))))
            
            if args.mp: print2(parameter, "\nStep 2/4: Calculating the coverage to all reads and averaging with multiprocessing ")
            else: print2(parameter, "\nStep 2/4: Calculating the coverage to all reads and averaging")
            lineplot.group_tags(groupby=args.g, sortby=args.s, colorby=args.c)
            lineplot.gen_cues()
            lineplot.coverage(sortby=args.s, heatmap=True, logt=args.log, mp=args.mp)
            t2 = time.time()
            print2(parameter, "    --- finished in {0} (h:m:s)".format(str(datetime.timedelta(seconds=round(t2-t1)))))
            
            # Sorting 
            print2(parameter, "\nStep 3/4: Sorting the data for heatmap")
            lineplot.hmsort(sort=args.sort)
            t3 = time.time()
            print2(parameter, "    --- finished in {0} (h:m:s)".format(str(datetime.timedelta(seconds=round(t3-t2)))))
            
            # Plotting
            print2(parameter, "\nStep 4/4: Plotting the heatmap")
            lineplot.hmcmlist(colorby = args.c, definedinEM = args.color)
            lineplot.heatmap(args.log)
            for i, name in enumerate(lineplot.hmfiles):
                output(f=lineplot.figs[i], directory = args.o, folder = args.t, filename=name,pdf=True,show=args.show)
            lineplot.gen_htmlhm(args.o, args.t)
            t4 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t4-t3))))
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t4-t0))) + "(H:M:S)\n")
            print("\nAll related files are saved in:  "+ os.path.join(dir,args.o,args.t))
            output_parameters(parameter, directory = args.o, folder = args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)