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
    
    def plot(self,qlist,color_list,groupedreference,color_tags,logt=None):
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
        self.exps = exps
        self.beds = exps.get_regionsets() # A list of GenomicRegionSets
        self.bednames = exps.get_regionsnames()
        self.reads = exps.get_readsfiles()
        self.readsnames = exps.get_readsnames()
        self.fieldsDict = exps.fieldsDict
        self.parameter = []
    
    def combine_allregions(self):
        all_bed = GenomicRegionSet("All regions")
        for bed in self.beds:
            all_bed.combine(bed)
        all_bed.remove_duplicates() #all_bed is sorted!!
        return all_bed
    
    def bedCoverage(self,bed):
        """ Return coverage matrix of multiple reads on one bed. 
        bed --> GenomicRegionSet
        """
        c=[]
        for rp in self.reads:
            r = os.path.abspath(rp)   # Here change the relative path into absolute path
            cov = CoverageSet(r,bed)
            cov.coverage_from_genomicset(r)
            #cov.normRPM()
            c.append(cov.coverage)
            print("    processing: "+rp)
        return numpy.transpose(c)
      
    def quantile_normalization(self,matrix):
        """ Return the np.array which contains the normalized values
        """
        
        rank_matrix = []
        for c in range(matrix.shape[1]):
            col = matrix[:,c]
            rank_col = mstats.rankdata(col)
            rank_matrix.append(rank_col)
    
        ranks = numpy.array(rank_matrix)
        trans_rank = numpy.transpose(ranks)
        
        # Calculate for means of ranks
        print("    Calculating for the mean of ranked data...")
        sort_matrix = numpy.sort(matrix,axis=0)
        means = []
        for r in range(matrix.shape[0]):
            row = [x for x in sort_matrix[r,:]]
            means.append(numpy.mean(row))
    
        # Replace the value by new means
        print("    Replacing the data value by normalized mean...")
        normalized_table = numpy.around(trans_rank)
        for i, v  in enumerate(means):
            normalized_table[normalized_table == i+1] = v
        #print(rounded_rank)
        return normalized_table

###########################################################################################
#                    lineplot 
###########################################################################################

    def tables_for_plot(self,norm_table,all_bed):
        """ Return a Dict which stores all tables for each bed with file path(more unique) as its key. """
        tableDict = {} # Storage all tables for each bed with bedname as the key
        conList = []   # Store containers of beds
        iterList = []
        
        for i,bed in enumerate(self.beds):
            tableDict[bed.name] = []
            bed.sort()
            conList.append(bed.__iter__())
            iterList.append(conList[-1].next())
            
        for i, r in enumerate(all_bed.sequences):
            for j in range(len(self.beds)):
                while r > iterList[j]:
                    try:
                        iterList[j] = conList[j].next()
                    except:
                        break
                if r == iterList[j]:
                    tableDict[self.beds[j].name].append(norm_table[i])
                elif r < iterList[j]:
                    continue
        return tableDict

    def group_tags(self, groupby, sortby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        def unique(a):
            seen = set()
            return [seen.add(x) or x for x in a if x not in seen]
        
        def gen_tags(tag):
            if tag == "reads":
                try: l = [self.exps.get_type(i,"factor") for i in self.readsnames]
                except: 
                    print("You must define 'factor' column in experimental matrix for grouping.")
                    sys.exit(1)
            elif tag == "regions":
                try: l = [self.exps.get_type(i,"factor") for i in self.bednames]
                except: 
                    print("You must define 'factor' column in experimental matrix for grouping.")
                    sys.exit(1)
            else:
                try: l = self.exps.fieldsDict[tag]
                except: 
                    print('Cannot find the column "' + tag +'"')
                    sys.exit(1)
            return unique(l)
        
        if groupby == "None":
            self.group_tags = ["All"]
        else:
            self.group_tags = gen_tags(groupby)
        if sortby == "None":
            self.sort_tags = ["All"]
        else:
            self.sort_tags = gen_tags(sortby)
        if colorby == "None":
            self.color_tags = ["All"]
        else:
            self.color_tags = gen_tags(colorby)

    def group_data(self, tables):  
        plotDict = OrderedDict()  # Extracting the data from different bed_bams file
        cues = OrderedDict()   # Storing the cues for back tracking
        for bedname in tables.keys():
            plotDict[bedname] = {}
            mt = numpy.array(tables[bedname])
            for i,readname in enumerate(self.readsnames):
                plotDict[bedname][readname] = mt[:,i]
                #print(plotDict[bedname][readname])
                x = tuple(self.exps.get_types(readname) + self.exps.get_types(bedname))
                cues[x] = [bedname, readname]
        #print(cues.keys())
        sortDict = {}  # Storing the data by sorting tags
        for g in self.group_tags:
            #print("    "+g)
            sortDict[g] = {}
            for a in self.sort_tags:
                #print("        "+c)
                sortDict[g][a] = {}
                for c in self.color_tags:
                    #print("            "+a)
                    sortDict[g][a][c] = []
                    for k in cues.keys():
                        if set([g,a,c]).difference(set(['All'])) <= set(k):
                            sortDict[g][a][c] = plotDict[cues[k][0]][cues[k][1]]
        self.sortDict = sortDict

    def colormap(self, colorby, definedinEM=False):
        """Generate the self.colors in the format which compatible with matplotlib"""
        if definedinEM:
            if colorby == "reads":
                self.colors = []
                for i in self.readsnames:
                    c = self.exps.get_type(i,"color")
                    if c[0] == "(":
                        rgb = [ eval(j) for j in c.strip('()').split(',')]
                        self.colors.append(rgb)
                    else:
                        self.colors.append(self.exps.get_type(i,"color"))
            elif colorby == "regions":
                self.colors = []
                for i in self.bednames:
                    c = self.exps.get_type(i,"color")
                    if c[0] == "(":
                        rgb = [ eval(j) for j in c.strip('()').split(',')]
                        self.colors.append(rgb)
                    else:
                        self.colors.append(self.exps.get_type(i,"color"))
            else:
                self.colors = [self.exps.get_type(i,"color") for i in self.exps.fieldsDict[colorby]]
        if definedinEM == False:
            self.colors = plt.cm.Set2(numpy.linspace(0, 1, len(self.color_tags)))

    def plot(self, title, html=False):
        """ Return boxplot from the given tables.
        
        """
        f, axarr = plt.subplots(1, len(self.group_tags), dpi=300)
        canvas = FigureCanvas(f)
        canvas.set_window_title(title)
        try: axarr = axarr.reshape(-1)
        except: axarr = [axarr]
        plt.subplots_adjust(bottom=0.3)
        
        axarr[0].set_ylabel("Count number (log)")
        for i, g in enumerate(self.group_tags):
            axarr[i].set_title(g, y=0.94)
            axarr[i].set_yscale('log')
            axarr[i].tick_params(axis='y', direction='out')
            axarr[i].yaxis.tick_left()
            axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
            d = []  # Store data within group
            color_t = []  # Store tag for coloring boxes
            x_ticklabels = []  # Store ticklabels
            for k, a in enumerate(self.sort_tags):
                for j, c in enumerate(self.color_tags):
                    if self.sortDict[g][a][c] == []:  # When there is no matching data, skip it
                        continue
                    else:
                        d.append([x+1 for x in self.sortDict[g][a][c]])
                        color_t.append(self.colors[j])
                        x_ticklabels.append(a)  #  + "." + c
            # Fine tuning boxplot
            bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None, 
                                  patch_artist=True, bootstrap=None)
            z = 10 # zorder for bosplot
            plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8,zorder=z)
            plt.setp(bp['fliers'], markerfacecolor='gray',color='none',alpha=0.3,markersize=1.8,zorder=z)
            plt.setp(bp['caps'],color='none',zorder=z)
            plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
            legends = []
            for patch, color in zip(bp['boxes'], color_t):
                patch.set_facecolor(color) # When missing the data, the color patch will exceeds
                patch.set_zorder(z)
                legends.append(patch)
                
            # Fine tuning subplot
            axarr[i].set_xticks([len(self.color_tags)*n + 1 + (len(self.color_tags)-1)/2 for n,s in enumerate(self.sort_tags)])
            #plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
            axarr[i].set_xticklabels(self.sort_tags, rotation=0, fontsize=10)
            
            axarr[i].set_ylim(bottom=0.95)
            for spine in ['top', 'right', 'left', 'bottom']:
                axarr[i].spines[spine].set_visible(False)
            axarr[i].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
            
            if i > 0:
                plt.setp(axarr[i].get_yticklabels(),visible=False)
                #plt.setp(axarr[i].get_yticks(),visible=False)
                axarr[i].minorticks_off()
                axarr[i].tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
                    
        plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
        #plt.legend(colors, color_tags, loc=7)
        axarr[-1].legend(legends[0:len(self.color_tags)], self.color_tags, loc='center left', handlelength=1, 
                 handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10},
                 bbox_to_anchor=(1.05, 0.5))
        f.tight_layout(pad=2, h_pad=None, w_pad=None)
        self.fig = f
        
    def gen_html(self,outputname, title, pvalue):
        ########## HTML ###################
        pd = os.path.join(dir,outputname,title)
        try:
            os.stat(os.path.dirname(pd))
        except:
            os.mkdir(os.path.dirname(pd))
        try:
            os.stat(pd)
        except:
            os.mkdir(pd)    
        f = open(os.path.join(pd,'boxplot.html'),'w')
        table = []
        # Header
        table.append(['<font size="7">' + title + "</font>"])
        # Each row is a plot with its data
        table.append(["<img src='boxplot.png' width=800 >"])
        
        #### Calculate p value ####
            
        for g in self.group_tags:
            table.append(['<font size="5">' + g + "</font>"])
            indM = 0
            header = []
            data_p = []
            arr = []
            
            for s in self.sort_tags:
                for c in self.color_tags:
                    header.append("{0}.{1}".format(s,c))
                    data_p.append(self.sortDict[g][s][c])

                    for i, d in enumerate(data_p[:indM]):
                        u, p_value = mannwhitneyu(data_p[indM], d)
                        arr.append(p_value)
                    indM = indM + 1
            #print(len(arr))
            [h,pc,a,b] = sm.multipletests(arr, alpha=pvalue, returnsorted=False)
            ar = numpy.chararray([len(self.color_tags)*len(self.sort_tags),len(self.color_tags)*len(self.sort_tags)], itemsize=10)
            ar[:] = "-"
            k = 0
            for c in self.color_tags:
                for s in self.sort_tags:
                    for i, d in enumerate(header[:k]):
                        ar[k,i] = "{:3.1e}".format(pc[0.5*k*(k-1) + i])
                    k = k + 1
                        
            nrows, ncols = ar.shape
            subtable = '<style>table,th,td{border:1px solid black;border-collapse:collapse;text-align:center;table-layout: fixed;font-size:8pt;}\
            </style><table style="width:800px">'
            for r in range(nrows+1):
                subtable += '<tr>'
                for c in range(ncols+1):
                    if r == 0:
                        if c == 0: subtable += '<td></td>'
                        elif c > 0: subtable += '<td>'+header[c-1]+'</td>'
                    if r > 0:
                        if c == 0: subtable += '<td>'+header[r-1]+'</td>'
                        elif c > 0:
                            print(r,"  ",c,"   ", int(0.5*(r-1)*(r-2) + c -1))
                            if c < r and h[int(0.5*(r-1)*(r-2) + c -1)]:
                                subtable += '<td><font color="red">'+ar[r-1,c-1]+'</font></td>'
                            else: subtable += '<td>'+ar[r-1,c-1]+'</td>'
                subtable += '</tr>'
            subtable += '</table>'
            table.append([subtable])
        table.append(["<a href='"+os.path.join(dir, outputname,title,"parameters.log")+" '><font size="+'"5"'+">Parameters</a>"])
        htmlcode = HTML.table(table)
        for line in htmlcode: f.write(line)
        f.close()
        
class lineplot:
    def __init__(self,exps, title="Boxplot"):
        pass
    
###########################################################################################
#                    Heatmap 
###########################################################################################
class heatmap:
    
    def __init__(self,exps, title="Boxplot"):
        pass







