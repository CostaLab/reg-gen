from __future__ import print_function
import unittest
import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
import AssociationAnalysis

path_input1 = "/Users/Yu-ChingTsai/Documents/workspace/RG local test/data/20140120 ExpMat/input jaccard.txt"
path_input2 = "/Users/Yu-ChingTsai/Documents/workspace/RG local test/data/20140120 ExpMat/input.txt"

class Test(unittest.TestCase):
    
    def test_jaccard_test(self):
        matrix1 = ExperimentalMatrix()
        matrix2 = ExperimentalMatrix()
        matrix1.read(path_input2)
        matrix2.read(path_input2)
        #AssociationAnalysis.jaccard_test(matrix1, matrix2,replicates=2)

    def test_projection_test(self):
        matrix1 = ExperimentalMatrix()
        matrix2 = ExperimentalMatrix()
        matrix1.read(path_input2)
        matrix2.read(path_input2)
        AssociationAnalysis.projection_test(matrix1, matrix2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()