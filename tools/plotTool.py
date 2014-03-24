# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)
import numpy
from scipy.stats import mstats
import matplotlib.pyplot as plt

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import *
from rgt.CoverageSet import *

# Local test
tool_path = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/reg-gen/trunk/tools/plotTools.py"
input_path = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/box_plot/input_EM.txt"

sys.argv = ["plotTools","boxplot",input_path]

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
    f = open(fileName,"w")
    f.write("\t"+("\t".join(namesCol))+"\n")
    for i,line in enumerate(table):
        f.write(namesLines[i]+"\t"+("\t".join([str(j) for j in line]))+"\n")
    f.close()

# For boxplot
def quantile_normalization(matrix):
    """ Return the np.array which contains the normalized values
    """
    """
    f = open(fileName)
    data = []
    for i, l in enumerate(f.readlines()):
        if i == 0: 
            list_BAM = l.split()
            print("\tQuantile normalizing following datasets:")
            for b in list_BAM:
                print("\t\t",b)   
        else:
            data.append([int(j) for j in l.split()[1:]])
    f.close()
    #print(data)
    matrix = numpy.array(data)
    """
    #print(matrix)
    #matrix.astype(numpy.float)
    rank_matrix = []
    for c in range(matrix.shape[1]):
        col = matrix[:,c]
        #new_col = mstats.rankdata(numpy.ma.masked_invalid(col))
        rank_col = mstats.rankdata(col)
        rank_matrix.append(rank_col)

    ranks = numpy.array(rank_matrix)
    trans_rank = numpy.transpose(ranks)
    #print(trans_rank)
    #numpy.savetxt("rankfile.txt",ranks)
    
    # Calculate for means of ranks
    sort_matrix = numpy.sort(matrix,axis=0)
    means = []
    for r in range(matrix.shape[0]):
        row = [x for x in sort_matrix[r,:]]
        means.append(numpy.mean(row))
    #print(sort_matrix)
    #print(means)

    # Replace the value by new means
    normalized_table = numpy.around(trans_rank)
    for i, v  in enumerate(means):
        normalized_table[normalized_table == i+1] = v
    #print(rounded_rank)
    return normalized_table

def tables_for_plot(norm_table,all_bed,beds):
    """ Return a Dict which stores all tables for each bed with bedname as its key. """
    tableDict = {} # Storage all tables for each bed with bedname as the key
    for bedname in beds.keys():
        tableDict[bedname] = []
        for i,r in enumerate(all_bed.sequences):
            if r in beds[bedname].sequences:
                tableDict[bedname].append(norm_table[i])
    return tableDict
    
def box_plot(tables,groupby = "bed"):
    """ Return boxplot from the given tables.
    
    Parameters:
        tables - a dictionary with bedname as key and normalized value as its value
        groupby - Two ways to plot: 1. bed grouped; 2. bam grouped
        
    
    """
    plt.Figure()
    plt.boxplot(tables, notch=False, sym='o', vert=True, whis=1.5,
                positions=None, widths=None, patch_artist=False,
                bootstrap=None)
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

# Input parameters dictionary

# Read the Experimental Matrix
exps = ExperimentalMatrix()
exps.read(input_path)
beds = exps.get_regionsets()
bednames = beds.keys()
reads = exps.get_readsfiles()
readsnames = exps.get_readsnames()

if sys.argv[2] and sys.argv[1] == "boxplot":

    # Combine the beds
    all_bed = GenomicRegionSet("All regions")
    for bed in beds.values():
        all_bed.combine(bed)
    all_bed.remove_duplicates()
    
    # Coverage of reads on all_bed
    all_table = bedCoverage(all_bed,reads) 
    printTable(readsnames,bednames,all_table,"cov/all_bed.txt")
    # Quantile normalization
    norm_table = quantile_normalization(all_table)
    printTable(readsnames,bednames,norm_table,"cov/norm_bed.txt")
    
    # Generate individual table for each bed
    tables = tables_for_plot(norm_table,all_bed, beds)
    
    #box_plot(tables)
