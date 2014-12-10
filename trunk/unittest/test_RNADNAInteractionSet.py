from __future__ import print_function
from __future__ import division
import sys
#print(sys.path)
sys.path.append('/projects/reggen/reg-gen')
import unittest
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.RNADNAInteractionSet import *
import os
from rgt.Util import GenomeData
from rgt.Util import OverlapType

sample_txp = "/projects/lncRNA/data/sample.txp"

"""Unit Test"""

class TestGenomicRegionSet(unittest.TestCase):
    
    def test_read_txp(self):
        txp = RNADNAInteractionSet(organism="hg19", filename=sample_txp)
        for s in txp:
        	pass

        a = len(txp)
        b = txp[3]
        
    def test_sort_tfo(self):
    	txp = RNADNAInteractionSet(organism="hg19", filename=sample_txp)
    	txp.sort_tfo()
    	
    def test_sort_tts(self):
        txp = RNADNAInteractionSet(organism="hg19", filename=sample_txp)
    	txp.sort_tts()
    	

    def test_filter_tfo(self):
    	txp = RNADNAInteractionSet(organism="hg19", filename=sample_txp)
    	txp.filter_tfo(start=500,end=1000, output=False, adverse=True)

    def test_filter_tfo(self):
    	txp = RNADNAInteractionSet(organism="hg19", filename=sample_txp)
    	g = GenomicRegionSet("g")
    	s = GenomicRegion(chrom="chr2", initial=74000000, final=75000000)
    	g.add(s)
    	result = txp.filter_tts(g)
    	for i in result:
    	    print(i)
    
   # def test_
   
if __name__ == "__main__":

    unittest.main()