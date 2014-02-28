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
Statistical analysis methods for ExperimentalMatrix

Author: Joseph Kuo

"""
def jaccard_test(query, reference, replicates=500, organism=GenomeData.CHROMOSOME_SIZES):
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
    
def projection_test(query, reference):
    """Return the projection test of each comparison of two ExperimentalMatrix.
    
    """
    print("Projection test")
    #print("query\treference\tp-value")

    for s in query.objectsDict.keys():
        for ss in reference.objectsDict.keys():
            inters = reference.objectsDict[ss].intersect(query.objectsDict[s])
            p = reference.objectsDict[ss].projection_test(query.objectsDict[s])
            print("%s\t%s\tp value: %.4f" % (s, ss, p))

#################################################################################################
##### PARAMETERS ################################################################################
#################################################################################################

# Parameters
params = []
params.append("\nAssociation Analysis")
params.append("")
params.append("Usage:\tpython AssociationAnalysis.py <command> [options]\n")
params.append("Command:\n\tjaccard\t\tJaccard test analysis")
params.append("\tprojection\tProjection test analysis\n")

# Jaccard
params.append("\nUsage:\tpython AssociationAnalysis.py jaccard -r <Path> -q <Path> [-R]\n")
params.append("Options:\t-r\treference input of ExperimetalMatrix")
params.append("\t\t-q\tquery input of ExperimentalMatrix")
params.append("\t\t-R\tnumber of randomizations for Jaccard test [100]\n")
# Projection
params.append("\nUsage:\tpython AssociationAnalysis.py projection -r <Path> -q <Path>\n")
params.append("Options:\t-r\treference input of ExperimetalMatrix")
params.append("\t\t-q\tquery input of ExperimentalMatrix\n")

if(len(sys.argv) == 1):
    for e in params[0:5]: print(e)
    sys.exit(0)
elif(len(sys.argv) == 2) and sys.argv[1] == "jaccard":
    for e in params[5:9]: print(e)
    sys.exit(0)
elif(len(sys.argv) == 2) and sys.argv[1] == "projection":
    for e in params[9:12]: print(e)
    sys.exit(0)
#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

# Input parameters dictionary
inputParameters={}
inputParameters["command"] = sys.argv[1]
i = 2
con_loop = True
while con_loop:
    try:
        inputParameters[sys.argv[i]] = sys.argv[i+1]
        i = i + 2
    except: con_loop = False

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
    