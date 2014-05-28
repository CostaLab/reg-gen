# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import numpy
from scipy.stats import mstats, wilcoxon, mannwhitneyu, rankdata
import time, datetime
import argparse
import HTML
from collections import *
import statsmodels.sandbox.stats.multicomp as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import getpass
import fnmatch

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import *
from rgt.CoverageSet import *

# Local test
dir = os.getcwd()


###########################################################################################
#                    Universal functions 
###########################################################################################

###########################################################################################
#                    Projection test
###########################################################################################
class projection_plot:
    def __init__(self):
        pass
    def plot(self,logt=None,qlist,color_list,groupedreference,color_tags):
        f, ax = plt.subplots()
        if logt:
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
        return f


###########################################################################################
#                    Jaccard test
###########################################################################################

###########################################################################################
#                    Boxplot 
###########################################################################################
class boxplot:
    """
    input:
        exps: input experimental matrix
        title: Default = boxplot
        groupby: Group the data by the given factor in the header of experimental matrix
        
        
        
    output:
        parameters: list of records
        figs: a list of figure(s)
        
        
    """
    def __init__(self,exps, title="Boxplot"):
        # Read the Experimental Matrix
        self.title = title
        self.beds = exps.get_regionsets() # A list of GenomicRegionSets
        self.bednames = exps.get_regionsnames()
        self.reads = exps.get_readsfiles()
        self.readsnames = exps.get_readsnames()
        self.fieldsDict = exps.fieldsDict
        self.parameter = []
        self.figs = []
        
    def 

###########################################################################################
#                    lineplot 
###########################################################################################
class lineplot:
###########################################################################################
#                    Heatmap 
###########################################################################################
class heatmap:







