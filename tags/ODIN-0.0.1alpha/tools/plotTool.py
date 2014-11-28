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
import HTML
from matplotlib.backends.backend_pdf import PdfPages

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import *
from rgt.CoverageSet import *

# Local test
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

def unique(a):
    seen = set()
    return [seen.add(x) or x for x in a if x not in seen]

def group_tags(exps, groupby, colorby, sortby):
    """Generate the tags for the grouping of plot
    
    Parameters:
        groupby = 'bam','bed','cell',or 'factor'
        colorby = 'bam','bed','cell',or 'factor'
        sortby = 'bam','bed','cell',or 'factor'
    """
    
    if groupby == "read":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        group_tags = unique(l)
    elif groupby == "region":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        group_tags = unique(l)
    else:
        group_tags = exps.fieldsDict[exps.fields[int(groupby[-1])-1]].keys()
        
    if colorby == "read":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        color_tags = unique(l)
    elif colorby == "region":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        color_tags = unique(l)
    else:
        color_tags = exps.fieldsDict[exps.fields[int(colorby[-1])-1]].keys()
    
    if sortby == "read":
        l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        sort_tags = unique(l)
    elif sortby == "region":
        l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        sort_tags = unique(l)
    else:
        sort_tags = exps.fieldsDict[str(exps.fields[int(sortby[-1])-1])].keys()
    return group_tags, color_tags, sort_tags

def group_data(tables, exps, group_tags, color_tags, sort_tags):
    
    plotDict = {}  # Extracting the data from different bed_bams file
    cues = {}   # Storing the cues for back tracking
    for bedname in tables.keys():
        plotDict[bedname] = {}
        mt = numpy.array(tables[bedname])
        for i,readname in enumerate(exps.get_readsnames()):
            plotDict[bedname][readname] = mt[:,i]
            #print(plotDict[bedname][readname])
            x = tuple(exps.get_types(readname) + exps.get_types(bedname))
            cues[x] = [bedname, readname]
    #print(cues.keys())
    
    sortDict = {}  # Storing the data by sorting tags
    for g in group_tags:
        #print("    "+g)
        sortDict[g] = {}
        for c in color_tags:
            #print("        "+c)
            sortDict[g][c] = {}
            for a in sort_tags:
                #print("            "+a)
                sortDict[g][c][a] = []
                for k in cues.keys():
                    if set([g,c,a]) == set(k):
                        sortDict[g][c][a] = plotDict[cues[k][0]][cues[k][1]]
                        #print(sortDict[g][c][a])
                    #else: sortDict[g][c][a] = []
    return sortDict 

def box_plot(sortDict, group_tags, color_tags, sort_tags):
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
    
    colors = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan']
    
    # Subplt by group_tags
    f, axarr = plt.subplots(1, len(group_tags), figsize=(10,7), dpi=100)
    f.canvas.set_window_title(args.t)
    f.suptitle(args.t, fontsize=20)
    axarr = axarr.reshape(-1)
    plt.subplots_adjust(bottom=0.3)
    
    axarr[0].set_ylabel("Count number")
    
    for i, g in enumerate(group_tags):
        axarr[i].set_title(g, y=0.94)
        axarr[i].set_yscale('log')
        axarr[i].tick_params(axis='y', direction='out')
        axarr[i].yaxis.tick_left()
        
        #axarr[i].grid(color='w', linestyle='-', linewidth=2, dash_capstyle='round')
        d = []  # Store data within group
        color_t = []  # Store tag for coloring boxes
        x_ticklabels = []  # Store ticklabels
        for j, c in enumerate(color_tags):
            for a in sort_tags:
                if sortDict[g][c][a] == []:  # When there is no matching data, skip it
                    continue
                else:
                    d.append([x+1 for x in sortDict[g][c][a]])
                    color_t.append(colors[j])
                    x_ticklabels.append(a + "." + c)

        # Fine tuning boxplot
        bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None, 
                              patch_artist=True, bootstrap=None)
        plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8)
        plt.setp(bp['fliers'], markerfacecolor='gray',color='none',alpha=0.3,markersize=1.8)
        plt.setp(bp['caps'],color='none')
        plt.setp(bp['medians'], color='black', linewidth=1.5)
        legends = []
        for patch, color in zip(bp['boxes'], color_t):
            patch.set_facecolor(color) # When missing the data, the color patch will exceeds
            legends.append(patch)
            
        # Fine tuning subplot
        axarr[i].set_xticklabels(x_ticklabels, rotation=90, fontsize=10)
        axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
        axarr[i].set_ylim(bottom=0.9)
        for spine in ['top', 'right', 'left', 'bottom']:
            axarr[i].spines[spine].set_visible(False)
        axarr[i].tick_params(axis='both', which='both', bottom='off', top='off', labelbottom='on')
        if i > 0:
            plt.setp(axarr[i].get_yticklabels(),visible=False)
            axarr[i].minorticks_off()
        
        #axarr[i].set_aspect(1)
#        for j, c in enumerate(color_tags):
#            bp['boxes'][j*len(sort_tags):(j+1)*len(sort_tags)].facecolor = colors[j]
                
    plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
    #plt.legend(colors, color_tags, loc=7)
    
    f.legend(legends[::len(sort_tags)], color_tags, loc='upper right', handlelength=1, 
             handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10},
             bbox_to_anchor=(0.96, 0.39))
             
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
    plt.savefig(args.output + ".eps", bbox_extra_artists=(legends), bbox_inches='tight')
    
    if args.pdf:
        with PdfPages(dir + args.output + '.pdf') as pp:
            #plt.savefig(pdf, format='pdf')
            pp.savefig(f)
            #pdf.savefig(f.suptitle)
            #pdf.savefig(axarr[:])
        
            # Set the file's metadata via the PdfPages object:
            d = pp.infodict()
            d['Title'] = args.t
            d['CreationDate'] = datetime.datetime.today()
            #pdf.close()
                
    if args.html:
        f = open(args.output+'.html','w')
        htmlcode = HTML.table(args.t, '<img src="' + args.output + '.png">')
        f.write(htmlcode)
        f.close()
    
    plt.show()

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################


parser = argparse.ArgumentParser(description='Provides various plotting tools.\nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho')

subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

parser_boxplot = subparsers.add_parser('boxplot',help='Boxplot based on the bam and bed files for gene association analysis.')
parser_boxplot.add_argument('input',help='The file name of the input Experimental Matrix file.')
parser_boxplot.add_argument('output',default='boxplot', help='The file name of the output file.(pdf or html)')
parser_boxplot.add_argument('-t', default='The boxplot', help='The title shown on the top of the plot.')
parser_boxplot.add_argument('-g',choices=['read','region','col4','col5','col6'], default='read', help="The way to group data into different subplots.(Default:read)")
parser_boxplot.add_argument('-c',choices=['read','region','col4','col5','col6'], default='region', help="The way to color data within the same subplot.(Default:region)")
parser_boxplot.add_argument('-s',choices=['read','region','col4','col5','col6'], default='col4', help="The way to sort data which shared the same color.(Default:col4)")
parser_boxplot.add_argument('-pdf', action="store_true", help='Save the figure in pdf format.')
parser_boxplot.add_argument('-html', action="store_true", help='Save the figure in html format.')
parser_else = subparsers.add_parser('Else tool', help='others')
parser_else.add_argument('input',help='The name of the input Experimental Matrix file.')

args = parser.parse_args()

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

if args.input and args.mode=='boxplot':
    
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
    group_tags, color_tags, sort_tags = group_tags(exps, groupby=args.g, colorby=args.c, sortby=args.s)
    #print(group_tags, "\n", color_tags, "\n", sort_tags)
    sortDict = group_data(tables, exps, group_tags, color_tags, sort_tags)
    #print(sortDict)
    box_plot(sortDict,group_tags, color_tags, sort_tags)
    t5 = time.time()
    print("    --- finished in {0:.3f} secs".format(t5-t4))
    print("Total running time is : " + str(datetime.timedelta(seconds=t5-t0)))
    
    #print("Total running time is : {0:.3f} secs.\n".format(t5-t0))
    