from __future__ import print_function
import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *

"""
Statistical analysis methods for ExperimentalMatrix

Author: Joseph Kuo

"""

class AssociationAnalysis:
    
    def jaccard_test(self, query, reference, replicates=500):
        """Return the jaccard indexes of every possible comparisons between two ExperimentalMatrix. 
        
        Method:
        The distribution of random jaccard index is calculated by randomizing query for certain times. 
        Then, we compare the real jaccard index to the distribution and formulate p-value as 
        p-value = (# random jaccard > real jaccard)/(# random jaccard)
        
        """
        print("Jaccard test")
        print("query\treference\tp-value")

        for s in query.get_regionsets():
            for ss in reference.get_regionsets():
                distribution = []
                for i in range(replicates):
                    random = s.random_regions(multiply_factor=1, overlap_result=False, overlap_input=True, chrom_M=False)
                    distribution.append(ss.jaccard(random))
                real_jaccard = s[1].jaccard(ss[1])
                p = sum(x > real_jaccard for x in distribution)/replicates
                print(s[0],"\t",ss[0],"\t",p)
        
    def projection_test(self, query, reference):
        """Return the projection test of each comparison of two ExperimentalMatrix.
        
        """
        print("Projection test")
        print("query\treference\tp-value")

        for s in query.get_regionsets():
            for ss in reference.get_regionsets():
                p = ss[1].projection_test(s[1])
                print(s[0],"\t",ss[0],"\t",p)
    
        
        