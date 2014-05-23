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
import HTML

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import GenomeData


dir = os.getcwd()
"""
Statistical analysis methods for ExperimentalMatrix

Author: Joseph Kuo

"""

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

parser = argparse.ArgumentParser(description='Provides various statistical tools for association analysis.\nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho')

subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

# Projection test
choice_group=['col4','col5','col6']
parser_projection = subparsers.add_parser('projection',help='Projection test evaluates the association level by comparing to the random binomial model. \
The null hypothesis is that no association between reference and query and their distribution is random.')
parser_projection.add_argument('reference',help='The file name of the reference Experimental Matrix file. Multiple references are acceptable.')
parser_projection.add_argument('query', help='The file name of the query Experimental Matrix file. Multiple queries are acceptable.')
parser_projection.add_argument('-g', choices = choice_group, default=None, help='Group data of query by given column for plot. Default: None')
parser_projection.add_argument('-c', choices = choice_group, default=None, help='Define the column of query for various colors in plot. Default: None')
parser_projection.add_argument('-organism', default='hg19', help='Define the organism')
parser_projection.add_argument('-plot', action="store_true", help='Generate the plot.')
parser_projection.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
parser_projection.add_argument('-output', default='projection_test', help='Define the filename of the output plot.(Default: projection_test)') 
#parser_projection.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
#parser_projection.add_argument('-html', action="store_true", help='Save the figure in html format.')
#parser_projection.add_argument('-show', action="store_true", help='Show the figure in the screen.')

# Jaccard
parser_jaccard = subparsers.add_parser('jaccard',help='Jaccard test evaluate the association level by comparing with jaccard index from repeating randomization.')
parser_jaccard.add_argument('reference',help='The file name of the reference Experimental Matrix file.')
parser_jaccard.add_argument('query', help='The file name of the query Experimental Matrix file.')
parser_jaccard.add_argument('-r', type=int, default=500, help='Repetition times of randomization.')
parser_jaccard.add_argument('-organism', default='hg19', help='Define the organism')
#parser_jaccard.add_argument('-plot', action="store_true", help='Generate the plot.') 
#parser_jaccard.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
#parser_jaccard.add_argument('-html', action="store_true", help='Save the figure in html format.')
#parser_jaccard.add_argument('-show', action="store_true", help='Show the figure in the screen.')

args = parser.parse_args()

#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

rEM, qEM = ExperimentalMatrix(), ExperimentalMatrix()
rEM.read(args.reference)
qEM.read(args.query)
references = rEM.get_regionsets()
referencenames = rEM.get_regionsnames()
query = qEM.get_regionsets()
querynames = qEM.get_regionsnames()

if args.g:
    groupedreference = OrderedDict()  # Store all bed names according to their types
    for r in references:
        ty = rEM.get_type(r.name,rEM.fields[int(args.g[3])-1])
        try: groupedreference[ty].append(r)
        except: groupedreference[ty] =[r]
        
    groupedquery = OrderedDict()  # Store all bed names according to their types
    for q in query:
        ty = qEM.get_type(q.name,qEM.fields[int(args.g[3])-1])
        try: groupedquery[ty].append(q)
        except: groupedquery[ty] =[q]
        
else:
    groupedreference = OrderedDict()
    groupedreference["All"] = references
    groupedquery = OrderedDict()
    groupedquery["All"] = query

############# Color #####################################
#color_list = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan']
#colors = plt.cm.Paired(numpy.linspace(0, 1, 12))
colors = [(0, 35/255, 138/255),(132/255, 29/255, 20/255)]
if args.c:
    color_tags = {}
    for q in query:
        c = qEM.get_type(q.name,qEM.fields[int(args.c[3])-1])
        color_tags[q.name] = c
else:
    color_tags = {}
    for q in query:
        color_tags[q.name] = q.name

print(color_tags.values())
color_list = {}
for i, c in enumerate(set(color_tags.values())):
    for q in color_tags.keys():
        if color_tags[q] == c:
            color_list[q] = colors[i]
color_tags['Background'] = 'Background'
color_list['Background'] = '0.70'

if args.mode == "projection":
    print("\nProjection test")
    print("    {0:25s}{1:25s}{2:s}".format("Reference","Query","p value"))
    qlist = OrderedDict()
    for ty in groupedquery.keys():
        qlist[ty] = OrderedDict()
        for i, r in enumerate(groupedreference[ty]):
            qlist[ty][r.name] = OrderedDict()
            for j, q in enumerate(groupedquery[ty]):
                background, ratio, p = r.projection_test(q, args.organism, extra=True)
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

    if args.plot:
        #nr = len(referencenames)
        #nq = len(querynames)
        #nAll = nr * nq
        
        #x = numpy.arange(0.75,nr+0.75,1)
        f, ax = plt.subplots()
        if args.log:
            ax.set_yscale('log')
        else:
            ax.locator_params(axis = 'y', nbins = 2)
        
        for ind_ty, ty in enumerate(qlist.keys()):
                         
            for ind_r,r in enumerate(qlist[ty].keys()):
                width = 0.8/(len(qlist[ty].keys()) * len(qlist[ty][r].keys())+1) # Plus one background
                for ind_q, q in enumerate(qlist[ty][r].keys()):
                    ax.bar(ind_ty * len(qlist[ty].keys())+ ind_r * len(qlist[ty][r].keys()) + ind_q*width,
                           qlist[ty][r][q], width=width, color=color_list[q],align='edge')
            
        ax.set_ylabel("Percentage of associated regions",fontsize=12)
        ax.yaxis.tick_left()
        ax.set_xlim(-0.2,len(qlist.keys())-0.2)
        ax.set_xticks([i + 0.5 - width for i in range(len(groupedreference.keys()))])
        ax.set_xticklabels(groupedreference.keys())
        ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
        ax.legend(set(color_tags.values()), loc='center left', handlelength=1, handletextpad=1, 
                  columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
        for spine in ['top', 'right']:  # 'left', 'bottom'
            ax.spines[spine].set_visible(False)
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(filename = args.output, bbox_extra_artists=(plt.gci()), bbox_inches='tight',dpi=300)
        
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
            
            
    
