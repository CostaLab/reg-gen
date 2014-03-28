# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import numpy
from scipy.stats import mstats, wilcoxon
import matplotlib.pyplot as plt
import time

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
    """ Return coverage matrix of multiple reads on one bed. """
    c=[]
    for r in reads:
        cov = CoverageSet(r,bed)
        cov.coverage_from_genomicset(r)
        #cov.normRPM()
        c.append(cov.coverage)
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
    sort_matrix = numpy.sort(matrix,axis=0)
    means = []
    for r in range(matrix.shape[0]):
        row = [x for x in sort_matrix[r,:]]
        means.append(numpy.mean(row))

    # Replace the value by new means
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
    
def box_plot(tables,bednames,readnames,groupby = "bed"):
    """ Return boxplot from the given tables.
    
    Parameters:
        tables - a dictionary with bedname as key and normalized value as its value
        groupby - "bed" --> bed grouped; 
                  "bam" --> bam grouped
    
    """
    ## Prepare the data
    
    plotDict = {}
    for bedname in tables.keys():
        plotDict[bedname] = {}
        mt = numpy.array(tables[bedname])
        for i,readname in enumerate(readnames):
            plotDict[bedname][readname] = mt[:,i]
            
    data = []
    x_label = []
    x_tick = []
    nb = 8
    if groupby == "bed":
        for bed in plotDict.keys():
            x_label.append(bed)
            for bam in plotDict[bed].keys():
                data.append(plotDict[bed][bam])
                x_tick.append(bam)
    elif groupby == "bam":
        for bam in readnames:
            x_label.append(bam)
            for bed in plotDict.keys():
                data.append(plotDict[bed][bam])
                x_tick.append(bed)
    
    # Calculate p value
    p = numpy.ones(len(data))  # Array for storing p values
    posi = numpy.zero(len(data)) # A list storing position [x,y]
    for i, s1 in enumerate(data[1:]):
        for j, s2 in enumerate(data[i+1:]):
            t, p[i,j] = wilcoxon(s1,s2)
            posi[i,j] = ((i+1+j+1)*0.5,max(data[i,j])*1.1)
    
    ## Initialize the figure
    fig = plt.figure()
    fig.subplots_adjust(bottom=0.3)
    
    # Plotting
    bp = plt.boxplot(data, notch=False, sym='o', vert=True, whis=1.5,
                     positions=None, widths=None, patch_artist=False, bootstrap=None)
    
    # Editting for details
    plt.title("Boxplot of Gene Association analysis")
    plt.set_yscale('log')   
    plt.ylabel("Count (log)")

    plt.setp(bp['boxes'], color='black')
    plt.setp(bp['whiskers'], color='black')
    plt.setp(bp['fliers'], markerfacecolor='none',color='gray',alpha=0.6)
    plt.setp(bp['medians'], color='black', linewidth=2)
    #plt.xticks(range(1,len(x_tick)+1), x_tick, rotation=70)
    gn = len(x_tick)/len(x_label)
    plt.xticks(range(0.5*(gn+1),len(x_tick),gn),x_label)
    # Label significance
    for i,r in enumerate(p):
        for j,pv in enumerate(r):
            if pv < 0.5:
                plt.annotate("*",xy=(posi[i,j][0],posi[i,j][1]+7),horizontalalignment='center',zorder=10)
                plt.annotate("",xy=(i+1,posi[i,j][1]),xytext=(j+1,posi[i,j][1]),arrowprops=props,zorder=10)
    
    # Saving
    plt.savefig("boxplot.png")
    plt.show()

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

"""
Usage:  $ python plotTool.py <command> [options]
"""

# Parameters
params = []
params.append("\nPlotTool")
params.append("")
params.append("Usage:\t\tpython plotTool.py <command> [options]")
params.append("Command:\tboxplot\t\tBox plot demonstration")
params.append("")

# Box plot
params.append("\nUsage:\t\tpython plotTool.py boxplot <file>")
params.append("Options:\t<file>\tcoverage matrix which is generated from coverageFromGenomicSet.py")
params.append("")

if(len(sys.argv) == 1):
    for e in params[0:5]: print(e)
    sys.exit(0)
elif(len(sys.argv) == 2) and sys.argv[1] == "boxplot":
    for e in params[5:8]: print(e)
    sys.exit(0)

#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

t0 = time.time()
# Input parameters dictionary

# Read the Experimental Matrix
exps = ExperimentalMatrix()
exps.read(sys.argv[2])
beds = exps.get_regionsets() # A list of GenomicRegionSets
bednames = exps.get_regionsnames()
reads = exps.get_readsfiles()
readsnames = exps.get_readsnames()

if sys.argv[2] and sys.argv[1] == "boxplot":

    # Combine the beds
    all_bed = GenomicRegionSet("All regions")
    for bed in beds:
        all_bed.combine(bed)
    all_bed.remove_duplicates() #all_bed is sorted!!
    
    t1 = time.time()
    print("Step 1/5 finished: Combining " + str(len(all_bed.sequences)) + " regions from all bed files")
    print("    --- required time: "+ str(t1-t0)+ " secs")
    
    # Coverage of reads on all_bed
    all_table = bedCoverage(all_bed,reads) 
    regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in all_bed]
    printTable(readsnames,regionnames,all_table,os.path.join(dir,"cov/all_bed.txt"))
    print("Step 2/5 finished: Calculating coverage of each bam file on all regions (saved as 'cov/all_bed.txt')")
    t2 = time.time()
    print("    --- required time: "+ str(t2-t1)+ " secs")
    
    # Quantile normalization
    norm_table = quantile_normalization(all_table)
    printTable(readsnames,regionnames,norm_table,os.path.join(dir,"cov/norm_bed.txt"))
    print("Step 3/5 finished: Quantile normalization of all coverage table (saved as 'cov/norm_bed.txt')")
    t3 = time.time()
    print("    --- required time: "+ str(t3-t2)+ " secs")
    
    # Generate individual table for each bed
    tables = tables_for_plot(norm_table,all_bed, beds)
    for bed in beds:
        print(len(bed.sequences))
        regionnames = [r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in bed]
        printTable(readsnames,regionnames,tables[bed.name],dir+"/cov/plot/"+bed.name+".txt")
    print("Step 4/5 finished: Constructing different tables for box plot (saved in 'cov/plot/')")
    t4 = time.time()
    print("    --- required time: "+ str(t4-t3)+ " secs")
    
    # Plotting
    box_plot(tables,bednames,readsnames,groupby = "bed")
    print("Step 5/5 finished:")
    t5 = time.time()
    print("    --- required time: "+ str(t5-t4)+ " secs")
    