# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import numpy
from scipy.stats import mstats, wilcoxon, mannwhitneyu
import matplotlib.pyplot as plt
import time, datetime
import argparse

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import *
from rgt.CoverageSet import *

# Local test
tool_path = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/reggen/tools/plotTools.py"
input_path = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/box_plot/input_EM.txt"
dir = os.getcwd()
#sys.argv = ["plotTools","boxplot",input_path]

"""
Plot tool for generating analytic figures.

Author: Joseph Kuo

Method:


"""

def bedCoverage(bed,reads):
    """ Return coverage matrix of multiple reads on one bed. 
    bed --> GenomicRegionSet
    
    """
    c=[]
    for rp in reads:
        r = os.path.abspath(rp)   # Here change the relative path into absolute path
        cov = CoverageSet(r,bed)
        cov.coverage_from_genomicset(r)
        #cov.normRPM()
        c.append(cov.coverage)
        print("    processing: "+rp)
    return numpy.transpose(c)  

def printTable(namesCol,namesLines,table,fileName):
    """ Write a txt file from the given table. """
    d = os.path.dirname(fileName)
    try:
        os.stat(d)
    except:
        os.mkdir(d)
             
    f = open(fileName,"w")
    f.write("\t"+("\t".join(namesCol))+"\n")
    for i,line in enumerate(table):
        f.write(namesLines[i]+"\t"+("\t".join([str(j) for j in line]))+"\n")
    f.close()

# For boxplot

def quantile_normalization(matrix):
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

def tables_for_plot(norm_table,all_bed,beds):
    """ Return a Dict which stores all tables for each bed with bedname as its key. """
    tableDict = {} # Storage all tables for each bed with bedname as the key
    conList = []   # Store containers of beds
    iterList = []
    
    for i,bed in enumerate(beds):
        tableDict[bed.name] = []
        bed.sort()
        conList.append(bed.__iter__())
        iterList.append(conList[-1].next())
        
    for i, r in enumerate(all_bed.sequences):
        for j in range(len(beds)):
            while r > iterList[j]:
                try:
                    iterList[j] = conList[j].next()
                except:
                    break
            if r == iterList[j]:
                tableDict[beds[j].name].append(norm_table[i])
            elif r < iterList[j]:
                continue
    return tableDict

def group_tags(exps, groupby, colorby, arrangeby):
    """Generate the tags for the grouping of plot
    
    Parameters:
        groupby = 'bam','bed','cell',or 'factor'
        colorby = 'bam','bed','cell',or 'factor'
        arrangeby = 'bam','bed','cell',or 'factor'
    """
    
    if groupby == "bam":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        group_tags = [ i for i in set(l)]
    elif groupby == "bed":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        group_tags = [ i for i in set(l)]
    else:
        group_tags = exps.fieldsDict[groupby].keys()
        
    if colorby == "bam":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        color_tags = [ i for i in set(l)]
    elif colorby == "bed":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        color_tags = [ i for i in set(l)]
    else:
        color_tags = exps.fieldsDict[colorby].keys()
    
    if arrangeby == "bam":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        arrange_tags = [ i for i in set(l)]
    elif arrangeby == "bed":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        arrange_tags = [ i for i in set(l)]
    else:
        arrange_tags = exps.fieldsDict[arrangeby].keys()
    
    return group_tags, color_tags, arrange_tags

def group_data(tables, exps, group_tags, color_tags, arrange_tags):
    
    plotDict = {}  # Extracting the data from different bed_bams file
    cues = {}   # Storing the cues for back tracking
    for bedname in tables.keys():
        plotDict[bedname] = {}
        mt = numpy.array(tables[bedname])
        for i,readname in enumerate(exps.get_readsnames()):
            plotDict[bedname][readname] = mt[:,i]
            #print(plotDict[bedname][readname])
            x = tuple(exps.get_types(readname) + [exps.get_type(bedname,'factor')])
            cues[x] = [bedname, readname]
    #print(cues.keys())
    
    sortDict = {}  # Storing the data by sorting tags
    for g in group_tags:
        #print("    "+g)
        sortDict[g] = {}
        for c in color_tags:
            #print("        "+c)
            sortDict[g][c] = {}
            for a in arrange_tags:
                #print("            "+a)
                sortDict[g][c][a] = []
                for k in cues.keys():
                    if set([g,c,a]) == set(k):
                        sortDict[g][c][a] = plotDict[cues[k][0]][cues[k][1]]
                        
                    #else: sortDict[g][c][a] = []
    return sortDict 

def box_plot(sortDict, group_tags, color_tags, arrange_tags):
    """ Return boxplot from the given tables.
    
    """
    ## Prepare the data
    """        
    data = []
    for g in sortDict.keys():
        for c in sortDict[g].keys():
            try:
                for a in sortDict[g][c]:
                    data.append(sortDict[g][c][a])
            except:
                data.append(sortDict[g][c])
    """            
    
    # Calculate p value
    """
    p = numpy.ones((len(data),len(data)))  # Array for storing p values
    posi = numpy.array((len(data),len(data))) # A list storing position [x,y]
    for i, s1 in enumerate(data[1:]):
        for j, s2 in enumerate(data[i+1:]):
            t, p[i,j] = mannwhitneyu(s1,s2)
            posi[i,j] = ((i+1+j+1)*0.5,max(data[i,j])*1.1)
    """
    ## Initialize the figure
    #plt.title("Boxplot of Gene Association analysis")
    
    colors = [ 'pink', 'cyan', 'lightblue', 'lightgreen', 'tan']
    
    # Subplt by group_tags
    f, axarr = plt.subplots(1, len(group_tags))
    axarr = axarr.reshape(-1)
    plt.subplots_adjust(bottom=0.3)
    
    axarr[0].set_ylabel("Count number (log)")
    
    for i, g in enumerate(group_tags):
        axarr[i].set_title(g)
        axarr[i].set_yscale('log')
        #axarr[i].grid(color='w', linestyle='-', linewidth=2, dash_capstyle='round')
        d = []
        color_t = []
        for j, c in enumerate(color_tags):
            for a in arrange_tags:
                d.append([x for x in sortDict[g][c][a]])
                color_t.append(colors[j])
        # Fine tuning boxplot
        bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None, 
                              patch_artist=True, bootstrap=None)
        plt.setp(bp['whiskers'], color='black')
        plt.setp(bp['fliers'], markerfacecolor='gray',color='none',alpha=0.3,size=0.2)
        plt.setp(bp['medians'], color='black', linewidth=2)
        for patch, color in zip(bp['boxes'], color_t):
            patch.set_facecolor(color)
        # Fine tuning subplot
        axarr[i].set_xticklabels([ s  + "." + c for s in arrange_tags]*len(color_tags), rotation=90)
        for spine in ['top', 'right', 'left', 'bottom']:
            axarr[i].spines[spine].set_visible(False)
        axarr[i].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        #axarr[i].set_aspect(1)
#        for j, c in enumerate(color_tags):
#            bp['boxes'][j*len(arrange_tags):(j+1)*len(arrange_tags)].facecolor = colors[j]
                
    plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
    #plt.legend(colors, color_tags, loc=7)
    
    
    """
    gn = len(x_tick)/len(x_label)
    plt.xticks(range(0.5*(gn+1),len(x_tick),gn),x_label)
    # Label significance
    for i,r in enumerate(p):
        for j,pv in enumerate(r):
            if pv < 0.5:
                plt.annotate("*",xy=(posi[i,j][0],posi[i,j][1]+7),horizontalalignment='center',zorder=10)
                plt.annotate("",xy=(i+1,posi[i,j][1]),xytext=(j+1,posi[i,j][1]),arrowprops=props,zorder=10)
    """
    # Saving
    plt.savefig("boxplot.png")
    plt.show()

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################


parser = argparse.ArgumentParser(description='Provides various plotting tools.')

subparsers = parser.add_subparsers(help='sub-command help')

parser_boxplot = subparsers.add_parser('boxplot',help='Boxplot based on the bam and bed files for gene association analysis.')
parser_boxplot.add_argument('input',help='The name of the input Experimental Matrix file.')
parser_boxplot.add_argument('output',default='boxplot', help='The name of the output file.(pdf or html)')
parser_boxplot.add_argument('-o','--output', choices=["pdf","html"], default="pdf",help='Choose the output format, pdf or html.')
parser_boxplot.add_argument('-t','--title', default='The boxplot for gene association analysis', help='The title shown on the top of the plot.')
parser_boxplot.add_argument('-g','--grouped',choices=['bam','bed','cell','factor'], default='bam', help="The way to group data into different subplots.")
parser_boxplot.add_argument('-c','--colored',choices=['bam','bed','cell','factor'], default='bed', help="The way to color data within the same subplot.")
parser_boxplot.add_argument('-a','--arranged',choices=['bam','bed','cell','factor'], default='cell', help="The way to arranged data which shared the same color.")

parser_else = subparsers.add_parser('Else tool', help='others')
parser_else.add_argument('input',help='The name of the input Experimental Matrix file.')

args = parser.parse_args()
print(args)

#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

t0 = time.time()
# Input parameters dictionary

# Read the Experimental Matrix
exps = ExperimentalMatrix()
exps.read(args.input)
beds = exps.get_regionsets() # A list of GenomicRegionSets
bednames = exps.get_regionsnames()
reads = exps.get_readsfiles()
readsnames = exps.get_readsnames()
fieldsDict = exps.fieldsDict

if args.input:

    # Combine the beds
    print("\nStep 1/5: Combining all regions")
    all_bed = GenomicRegionSet("All regions")
    for bed in beds:
        all_bed.combine(bed)
    all_bed.remove_duplicates() #all_bed is sorted!!
    print("    " + str(len(all_bed.sequences)) + " regions from all bed files are combined.")
    t1 = time.time()
    
    print("    --- finished in {0:.3f} secs".format(t1-t0))
    
    # Coverage of reads on all_bed
    print("Step 2/5: Calculating coverage of each bam file on all regions")
    all_table = bedCoverage(all_bed,reads) 
    regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in all_bed]
    #printTable(readsnames,regionnames,all_table,os.path.join(dir,"cov/all_bed.txt"))
    t2 = time.time()
    print("    --- finished in {0:.3f} secs".format(t2-t1))
    
    # Quantile normalization
    print("Step 3/5: Quantile normalization of all coverage table")
    norm_table = quantile_normalization(all_table)
    #printTable(readsnames,regionnames,norm_table,os.path.join(dir,"cov/norm_bed.txt"))
    t3 = time.time()
    print("    --- finished in {0:.3f} secs".format(t3-t2))
    
    # Generate individual table for each bed
    print("Step 4/5: Constructing different tables for box plot")
    tables = tables_for_plot(norm_table,all_bed, beds)
    #print(tables.keys())
    #print(tables)
    
    #for bed in beds:
        #print(len(bed.sequences))
        #regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in bed]
        #printTable(readsnames,regionnames,tables[bed.name],dir+"/cov/"+bed.name+".txt")
    t4 = time.time()
    print("    --- finished in {0:.3f} secs".format(t4-t3))
    
    # Plotting
    print("Step 5/5: Plotting")
    group_tags, color_tags, arrange_tags = group_tags(exps, groupby=args.grouped, colorby=args.colored, arrangeby=args.arranged)
    #print(group_tags, "\n", color_tags, "\n", arrange_tags)
    sortDict = group_data(tables, exps, group_tags, color_tags, arrange_tags)
    #print(sortDict)
    box_plot(sortDict,group_tags, color_tags, arrange_tags)
    t5 = time.time()
    print("    --- finished in {0:.3f} secs".format(t5-t4))
    print("Total running time is : " + str(datetime.timedelta(seconds=t5-t0)))
    
    #print("Total running time is : {0:.3f} secs.\n".format(t5-t0))
    