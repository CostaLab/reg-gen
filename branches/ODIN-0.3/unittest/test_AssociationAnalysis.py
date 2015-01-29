from __future__ import print_function
import unittest
import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
import AssociationAnalysis

path_input1 = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/20140120 ExpMat/input.txt"
path_input2 = "/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/projects/20140120 ExpMat/input2.txt"

class Test(unittest.TestCase):
    
    def test_jaccard_test(self):
        matrix1 = ExperimentalMatrix()
        matrix2 = ExperimentalMatrix()
        matrix1.read(path_input2)
        matrix2.read(path_input2)
        #AssociationAnalysis.jaccard_test(matrix1, matrix2,replicates=3)

    def test_projection_test(self):
        matrix1 = ExperimentalMatrix()
        matrix2 = ExperimentalMatrix()
        matrix1.read(path_input2)
        matrix2.read(path_input2)
        AssociationAnalysis.projection_test(matrix1, matrix2)

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()