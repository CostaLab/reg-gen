# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import argparse
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy
import time, datetime, getpass, fnmatch, HTML

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import GenomeData
from plottools import *

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

def output(f, directory, folder, filename, extra=None):
    """Output the file in the defined folder """
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)    
    
    # Saving 
    if extra == None:
        f.savefig(os.path.join(pd,filename), bbox_inches='tight',dpi=300)
    else:
        f.savefig(os.path.join(pd,filename), bbox_extra_artists=(extra), bbox_inches='tight',dpi=300)
    
    if args.pdf:
        pp = PdfPages(os.path.join(pd,filename) + '.pdf')
        pp.savefig(f, bbox_extra_artists=(extra),bbox_inches='tight') 
    
        # Set the file's metadata via the PdfPages object:
        d = pp.infodict()
        d['Title'] = args.title
        d['CreationDate'] = datetime.datetime.today()
        pp.close()
    
    if args.show:
        plt.show()
        
def output_parameters(parameter, directory, folder, filename):
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)    
    if parameter:
        logp = open(os.path.join(pd,"parameters.log"),'w')
        for s in parameter:
            logp.write("{}\n".format(s))
        logp.close()

def fetch_refandque(referenceEM, queryEM):
    rEM, qEM = ExperimentalMatrix(), ExperimentalMatrix()
    rEM.read(referenceEM)
    qEM.read(queryEM)
    references = rEM.get_regionsets()
    referencenames = rEM.get_regionsnames()
    query = qEM.get_regionsets()
    querynames = qEM.get_regionsnames()
    return [rEM, references, referencenames, qEM, query, querynames]

def group_refque(rEM, references, referencenames, qEM, query, querynames, groupby):
    groupedreference = OrderedDict()  # Store all bed names according to their types
    for r in references:
        ty = rEM.get_type(r.name,groupby)
        try: groupedreference[ty].append(r)
        except: groupedreference[ty] =[r]
    groupedquery = OrderedDict()  # Store all bed names according to their types
    for q in query:
        ty = qEM.get_type(q.name,groupby)
        try: groupedquery[ty].append(q)
        except: groupedquery[ty] =[q]
    return groupedreference, groupedquery

def projection_test(groupedreference, groupedquery, organism):
    qlist = OrderedDict()
    for ty in groupedquery.keys():
        qlist[ty] = OrderedDict()
        for i, r in enumerate(groupedreference[ty]):
            qlist[ty][r.name] = OrderedDict()
            for j, q in enumerate(groupedquery[ty]):
                background, ratio, p = r.projection_test(q, organism, extra=True)
                qlist[ty][r.name][q.name] = ratio
                if p < 0.025: 
                    if len(q) == 0:
                        print("    {0:25s}{1:25s}{2:.2e}\tEmpty query!".format(r.name,q.name,p))
                    else:
                        print("    {0:25s}{1:25s}{2:.2e}\tSignificantly unassociated!".format(r.name,q.name,p))
                elif p > 0.975:
                    if len(q) == 0:
                        print("    {0:25s}{1:25s}{2:.2e}\tEmpty query!".format(r.name,q.name,p))
                    else:
                        print("    {0:25s}{1:25s}{2:.2e}\tSignificantly associated!".format(r.name,q.name,p))
                else: print("    {0:25s}{1:25s}{2:.2e}".format(r.name,q.name,p))
        qlist[ty][r.name]['Background'] = background
    return qlist


######### 
#########
######### 
######### 
#########     
    

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

# Some general help descriptions
######### Some general plotting arguments descriptions ###############
helpinput = 'The file name of the input Experimental Matrix file. Recommended to add more columns for more information for ploting. For example, cell type or factors.'
helpoutput = 'The directory name for the output files. For example, project name.'
helptitle = 'The title shown on the top of the plot and also the folder name.'
helpgroup = "Group the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."
helpsort = "Sort the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."
helpcolor = "Color the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None."

parser = argparse.ArgumentParser(description='Provides various Statistical analysis methods and plotting tools for ExperimentalMatrix.\
\nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

################### Projection test ##########################################
parser_projection = subparsers.add_parser('projection',help='Projection test evaluates the association level by comparing to the random binomial model. \
The null hypothesis is that no association between reference and query and their distribution is random.')
parser_projection.add_argument('reference',help='The file name of the reference Experimental Matrix. Multiple references are acceptable.')
parser_projection.add_argument('query', help='The file name of the query Experimental Matrix. Multiple queries are acceptable.')
parser_projection.add_argument('-g', default=None, help=helpgroup +" (Default:None)")
parser_projection.add_argument('-c', default=None, help=helpcolor +' (Default: None)')
parser_projection.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')
parser_projection.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
parser_projection.add_argument('-output', default='projection_test', help='Define the filename of the output plot.(Default: projection_test)') 
parser_projection.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
parser_projection.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_projection.add_argument('-show', action="store_true", help='Show the figure in the screen.')

################### Jaccard test ##########################################

parser_jaccard = subparsers.add_parser('jaccard',help='Jaccard test evaluate the association level by comparing with jaccard index from repeating randomization.')
parser_jaccard.add_argument('reference',help='The file name of the reference Experimental Matrix file.')
parser_jaccard.add_argument('query', help='The file name of the query Experimental Matrix file.')
parser_jaccard.add_argument('-r', type=int, default=500, help='Repetition times of randomization.')
parser_jaccard.add_argument('-organism', default='hg19', help='Define the organism')
#parser_jaccard.add_argument('-plot', action="store_true", help='Generate the plot.') 
#parser_jaccard.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
#parser_jaccard.add_argument('-html', action="store_true", help='Save the figure in html format.')
#parser_jaccard.add_argument('-show', action="store_true", help='Show the figure in the screen.')

################### Boxplot ##########################################

parser_boxplot = subparsers.add_parser('boxplot',help='Boxplot based on the BAM and BED files for gene association analysis.')
parser_boxplot.add_argument('input',help=helpinput)
parser_boxplot.add_argument('output', help=helpoutput)
parser_boxplot.add_argument('-t','--title', default='Boxplot', help=helptitle)
parser_boxplot.add_argument('-g', default='reads', help=helpgroup + " (Default:reads)")
parser_boxplot.add_argument('-c', default='regions', help=helpcolor + " (Default:regions)")
parser_boxplot.add_argument('-s', default='cell', help=helpsort + " (Default:cell)")
parser_boxplot.add_argument('-color', action="store_true", help='Define the specific colors with the given column "color" in experimental matrix. The color should be in the format of matplotlib.colors. For example, "r" for red, "b" for blue, or "(100/255, 35/255, 138/255)" for RGB.')
parser_boxplot.add_argument('-nqn', action="store_true", help='No quantile normalization in calculation.')
parser_boxplot.add_argument('-pdf', action="store_true", help='Save the figure in pdf format.')
parser_boxplot.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_boxplot.add_argument('-p','--pvalue', type=float, default=0.05, help='Define the significance level for multiple test. Default: 0.01')
parser_boxplot.add_argument('-show', action="store_true", help='Show the figure in the screen.')


################### Lineplot ##########################################

################### Heatmap ##########################################

################### Integration ##########################################


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

t0 = time.time()
# Input parameters dictionary
parameter = [] 

parameter.append("Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
parameter.append("User: " + getpass.getuser())
parameter.append("Project: " + args.output)
parameter.append("\nCommand:\n   $ python " + " ".join(sys.argv))
parameter.append("")


if args.mode == 'projection':
    # Fetching reference and query EM
    [rEM, references, referencenames, qEM, query, querynames] = fetch_refandque(args.reference,args.query)
    
    if args.g:
        groupedreference, groupedquery = group_refque(rEM, references, referencenames, qEM, query, querynames)
    else:
        groupedreference = OrderedDict()
        groupedreference["All"] = references
        groupedquery = OrderedDict()
        groupedquery["All"] = query
    
    ############# Color #####################################
    #color_list = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan']
    colors = plt.cm.Paired(numpy.linspace(0, 1, 12))
    #colors = [(0, 35/255, 138/255),(132/255, 29/255, 20/255)]
    if args.c:
        color_tags = {}
        for q in query:
            c = qEM.get_type(q.name,qEM.fields[int(args.c[3])-1])
            color_tags[q.name] = c
    else:
        color_tags = {}
        for q in query:
            color_tags[q.name] = q.name
    
    #print(color_tags.values())
    color_list = {}
    for i, c in enumerate(set(color_tags.values())):
        for q in color_tags.keys():
            if color_tags[q] == c:
                color_list[q] = colors[i]
    color_tags['Background'] = 'Background'
    color_list['Background'] = '0.70'
    
    ################### Projection test ##########################################
    print("\nProjection test")
    print("    {0:25s}{1:25s}{2:s}".format("Reference","Query","p value"))
    
    qlist = projection_test(groupedreference, groupedquery, organism = args.organism)
    
## TO DO
    if args.pdf:
        f = projection_plot(args.log, qlist,color_list,groupedreference,color_tags)
        f.savefig(filename = args.output, bbox_extra_artists=(plt.gci()), bbox_inches='tight',dpi=300)

################### Jaccard test ##########################################
if args.mode == "jaccard":
    """Return the jaccard test of every possible comparisons between two ExperimentalMatrix. 
    
    Method:
    The distribution of random jaccard index is calculated by randomizing query for given times. 
    Then, we compare the real jaccard index to the distribution and formulate p-value as 
    p-value = (# random jaccard > real jaccard)/(# random jaccard)
    
    """
    print("\nJaccard Test")
    print("    {0:25s}{1:25s}{2:s}".format("Reference","Query","p-value"))
    
    for i, r in enumerate(referencenames):
        for j, q in enumerate(querynames):
            #t0 = time.clock()
            random_jaccards = [] # Store all the jaccard index from random regions
            for k in range(args.r):
                random = query[j].random_regions(organism=args.organism, multiply_factor=1, overlap_result=True, overlap_input=True, chrom_M=False)
                random_jaccards.append(references[i].jaccard(random))
            real_jaccard = query[j].jaccard(references[i]) # The real jaccard index from r and q
            # How many randomizations have higher jaccard index than the real index?
            p = len([x for x in random_jaccards if x > real_jaccard])/args.r
            print("    {0:25s}{1:25s}{2:.2e}".format(referencenames[i],querynames[j],p))
         
################### Boxplot ##########################################
if args.mode == 'boxplot':
    exps = ExperimentalMatrix()
    exps.read(args.input)
    boxplot = boxplot(exps, title="Boxplot")
    
    print2(parameter,"\nStep 1/5: Combining all regions")
    all_bed = boxplot.combine_allregions()
    print2(parameter,"    " + str(len(all_bed.sequences)) + " regions from all bed files are combined.")
    t1 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t1-t0))
    
    # Coverage of reads on all_bed
    print2(parameter,"Step 2/5: Calculating coverage of each bam file on all regions")
    all_table = boxplot.bedCoverage(all_bed) 
    t2 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t2-t1))
    
    # Quantile normalization
    print2(parameter,"Step 3/5: Quantile normalization of all coverage table")
    if args.nqn:
        print2(parameter,"    No quantile normalization.")
        norm_table = all_table
    else:
        norm_table = boxplot.quantile_normalization(all_table)
    t3 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t3-t2))
    
    # Generate individual table for each bed
    print2(parameter,"Step 4/5: Constructing different tables for box plot")
    tables = boxplot.tables_for_plot(norm_table,all_bed)
    t4 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t4-t3))
    
    # Plotting
    print2(parameter,"Step 5/5: Plotting")
    boxplot.group_tags(groupby=args.g, sortby=args.s, colorby=args.c)
    boxplot.group_data(tables=tables)
    boxplot.colormap(colorby=args.c, definedinEM=args.color)
    boxplot.plot(title=args.title,html=args.html)
    output(f=boxplot.fig, directory = args.output, folder = args.title, filename="boxplot",extra=plt.gci())
    # HTML
    if args.html: boxplot.gen_html(args.output, args.title, args.pvalue)
    t5 = time.time()
    print2(parameter,"    --- finished in {0:.3f} secs\n".format(t5-t4))
    print2(parameter,"Total running time is : " + str(datetime.timedelta(seconds=round(t5-t0))))
    output_parameters(parameter, directory = args.output, folder = args.title, filename="parameters.log")


################### Lineplot #########################################

################### Heatmap ##########################################

################### Integration ######################################
