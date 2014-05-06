# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import argparse
import matplotlib.pyplot as plt
import numpy

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
def jaccard_test(query, reference, organism,replicates=500, ):
    """Return the jaccard test of every possible comparisons between two ExperimentalMatrix. 
    
    Method:
    The distribution of random jaccard index is calculated by randomizing query for given times. 
    Then, we compare the real jaccard index to the distribution and formulate p-value as 
    p-value = (# random jaccard > real jaccard)/(# random jaccard)
    
    """
    print("Jaccard test")
    print("query\treference\tp-value")
    
    result = []
    for s in query.objectsDict.keys():
        for ss in reference.objectsDict.keys():
            #t0 = time.clock()
            distribution = []
            for rep in range(replicates):
                random = query.objectsDict[s].random_regions(organism, multiply_factor=1, overlap_result=True, overlap_input=True, chrom_M=False)
                distribution.append(reference.objectsDict[ss].jaccard(random))
            real_jaccard = query.objectsDict[s].jaccard(reference.objectsDict[ss])
            p = sum(x for x in distribution if x > real_jaccard)/replicates
            print(s, ss, p, sep="\t")
            #t1 = time.clock()
            #print(t1 - t0, "randoming")
    

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

parser = argparse.ArgumentParser(description='Provides various statistical tools for association analysis.\nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho')

subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

# Projection test
parser_projection = subparsers.add_parser('projection',help='Projection test evaluate the association level by comparing to the random binomial model.')
parser_projection.add_argument('reference',help='The file name of the reference Experimental Matrix file.')
parser_projection.add_argument('query', help='The file name of the query Experimental Matrix file.')
parser_projection.add_argument('-organism', default='hg19', help='Define the organism')
parser_projection.add_argument('-plot', action="store_true", help='Generate the plot.')
parser_projection.add_argument('-output', default='projection_test', help='Define the filename of the output plot.(Default: projection_test)') 
#parser_projection.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
#parser_projection.add_argument('-html', action="store_true", help='Save the figure in html format.')
#parser_projection.add_argument('-show', action="store_true", help='Show the figure in the screen.')

# Jaccard
parser_jaccard = subparsers.add_parser('jaccard',help='Jaccard test evaluate the association level by comparing with jaccard index from repeating randomization.')
parser_jaccard.add_argument('reference',help='The file name of the reference Experimental Matrix file.')
parser_jaccard.add_argument('query', help='The file name of the query Experimental Matrix file.')
parser_jaccard.add_argument('-r', type=int, default=500, help='Repetition times of randomization.')
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
query = rEM.get_regionsets()
querynames = qEM.get_regionsnames()

if args.mode == "projection":
    print("\nProjection test")
    print("    {0:25s}{1:25s}{2:s}".format("Reference","Query","p value"))
    qlist = []
    for i, r in enumerate(references):
        for j, q in enumerate(query):
            q, p = r.projection_test(q, args.organism, proportion=True)
            qlist.append(q)
            if p < 0.025: print("    {0:25s}{1:25s}{2:.2e}\tSignificantly unassociated!".format(referencenames[i],querynames[j],p))
            elif p > 0.975: print("    {0:25s}{1:25s}{2:.2e}\tSignificantly associated!".format(referencenames[i],querynames[j],p))
    if args.plot:
        nr = len(referencenames)
        nq = len(querynames)
        nAll = nr * nq
        color_list = plt.cm.Set2(numpy.linspace(0, 1, 12))
        width = 0.5/nq
        x = numpy.arange(0.75,nr+0.75,1)
        f, ax = plt.subplots()
        for i,qy in enumerate(querynames):
            ax.bar(x+i*width,[qlist[j] for j in range(i,nAll,nq)],width=width,color=color_list[i],align='edge')
        ax.locator_params(axis = 'y', nbins = 2)
        ax.set_ylabel("Proportion of overlapped query",fontsize=10)
        ax.yaxis.tick_left()
        ax.set_xticks(range(1,nr+1))
        ax.set_xticklabels(referencenames)
        ax.legend(querynames,loc='center left', handlelength=1, handletextpad=1, 
                  columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
        for spine in ['top', 'right']:  # 'left', 'bottom'
            ax.spines[spine].set_visible(False)
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(filename = args.output, bbox_extra_artists=(plt.gci()), bbox_inches='tight',dpi=300)
        
if args.mode == "jaccard":
    #jaccard_test(args.reference, args.query, replicates=args.r, args.organism)
    print("done")