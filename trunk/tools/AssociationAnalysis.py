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
params.append("Analysis association between different experimental data by different statistical methods.")
params.append("\nAssociationAnalysis.py <method> <reference> <query> <others>\n")
params.append("    <method>")
params.append("        jaccard")
params.append("            Return the jaccard test of every possible comparisons between two ExperimentalMatrix.")
params.append("        projection")
params.append("            Return the projection test of each comparison of two ExperimentalMatrix.")
params.append("    <reference>")
params.append("        The text file contains ExperimentalMatrix. (As reference)")
params.append("        Format: txt defined by ExperimentalMatrix")
params.append("        Default: None.")
params.append("    <query>")
params.append("        The text file contains ExperimentalMatrix. (As query)")
params.append("        Format: txt defined by ExperimentalMatrix")
params.append("        Default: None.")
params.append("\nOther Parameters: ")
params.append("    <replicate>")
params.append("        For Jaccard test's replication times.")
params.append("        Default: 100")
params.append("")
if(len(sys.argv) == 1):
    for e in params: print(e)
    sys.exit(0)
elif(1 < len(sys.argv) < 4):
    print("Error: Please add",4 - len(sys.argv), "more parameters according to the following format.")
    for e in params[2:]: print(e)
    sys.exit(0)
#################################################################################################
##### INPUT #####################################################################################
#################################################################################################

# Input parameters dictionary
inputParameters={}
inputParameters["-method"] = sys.argv[1]
inputParameters["-reference"] = sys.argv[2]
inputParameters["-query"] = sys.argv[3]

# Input ExperimentalMatrix1
reference = ExperimentalMatrix()
if ("-reference" in inputParameters.keys()):
    reference.read(inputParameters["-reference"])
# Input ExperimentalMatrix2
query = ExperimentalMatrix()
if ("-query" in inputParameters.keys()):
    query.read(inputParameters["-query"])
# Input parameter
if ("-method" in inputParameters.keys()):
    if inputParameters["-method"] == "jaccard":
        try:
            repeats = sys.argv[4]
        except:
            repeats = 100
        jaccard_test(query, reference, replicates=repeats, organism=GenomeData.CHROMOSOME_SIZES)
    
    if inputParameters["-method"] == "projection":
        projection_test(query, reference)
    