from __future__ import print_function
from __future__ import division
import sys
import os.path
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
            t0 = time.clock()
            distribution = []
            for rep in range(replicates):
                random = query.objectsDict[s].random_regions(organism, multiply_factor=1, overlap_result=True, overlap_input=True, chrom_M=False)
                distribution.append(reference.objectsDict[ss].jaccard(random))
            real_jaccard = query.objectsDict[s].jaccard(reference.objectsDict[ss])
            p = sum(x for x in distribution if x > real_jaccard)/replicates
            print(s, ss, p, sep="\t")
            t1 = time.clock()
            print(t1 - t0, "randoming")
    
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
            print()

    
        