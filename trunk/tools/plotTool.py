# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(lib_path)

# Local Libraries
# Distal Libraries
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.Util import GenomeData

"""
Plot tool for generating analytic figures.

Author: Joseph Kuo

"""
# For boxplot
def quantile_normalization(fileName):
    f = open(fileName)
    data = []
    for i, l in enumerate(f.readlines()):
        if i == 0: 
            list_BAM = l.split()   
        else:
            data.append([float(j) for j in l.split()[1:]])
    f.close()
    matrix = np.array(data)
    matrix.astype(numpy.float)
    rank_matrix = []
    for c in range(matrix.shape[1]):
        col = [x for x in matrix[:,c]]
        new_col = mstats.rankdata(np.ma.masked_invalid(col))
        rank_matrix.append(new_col)
    ranks = np.array(rank_matrix)
    np.transpose(ranks)
    
    # Calculate for means of ranks
    mdat = np.ma.masked_array(matrix,np.isnan(matrix))
    means = np.mean(mdat,axis=1)
    ff = open(path_input2,'w')
    ff.write(means)

    new_rank = numpy.copy(ranks)
    for i in range(len(means)):
        new_rank[ranks==i+1] = means[i]
    new_rank[ranks == 0] = np.nan
    
    for e in new_rank:
        print(e)
        ff.write(e)
    ff.close()

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

# Parameters
params = []
params.append("\nPlotTool")
params.append("")
params.append("Usage:\tpython plotTool.py <command> [options]\n")
params.append("Command:\n\tboxplot\t\tBox plot demonstration")

# Box plot
params.append("\nUsage:\tpython plotTool.py boxplot -cm <Path>\n")
params.append("Options:\t-cm\tcoverage matrix which is generated from coverageFromGenomicSet.py")

if(len(sys.argv) == 1):
    for e in params[0:4]: print(e)
    sys.exit(0)
elif(len(sys.argv) == 2) and sys.argv[1] == "boxplot":
    for e in params[4:6]: print(e)
    sys.exit(0)

#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

# Input parameters dictionary

# Input ExperimentalMatrix1
reference = ExperimentalMatrix()
if ("-r" in inputParameters.keys()):
    reference.read(inputParameters["-r"])
# Input ExperimentalMatrix2
query = ExperimentalMatrix()
if ("-q" in inputParameters.keys()):
    query.read(inputParameters["-q"])
# Input parameter
if inputParameters["command"] == "jaccard":
    try:
        repeats = inputParameters["-R"]
    except:
        repeats = 100
    jaccard_test(query, reference, replicates=repeats, organism=GenomeData.CHROMOSOME_SIZES)

if inputParameters["command"] == "projection":
    projection_test(query, reference)
    