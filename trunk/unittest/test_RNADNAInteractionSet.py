from __future__ import print_function
from __future__ import division
import sys
import unittest
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.RNADNAInteractionSet import *
import os
from rgt.Util import GenomeData
from rgt.Util import OverlapType

sample_txp = "/projects/lncRNA/data/dn_hotair_genes_test.txp"

"""Unit Test"""

class TestGenomicRegionSet(unittest.TestCase):
    
    def test_read_txp(self):
        txp = RNADNAInteractionSet(filename=sample_txp)
        print(txp)
    
    
   
if __name__ == "__main__":

    unittest.main()