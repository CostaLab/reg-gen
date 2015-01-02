# Python Libraries
from __future__ import print_function
from collections import *
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from RNABindingSet import RNABinding, RNABindingSet
from SequenceSet import Sequence, SequenceSet

####################################################################################
####################################################################################

class RegionSequenceDict:
    """A structure to connect a region (GenomicRegion or RNABinding) to its sequence (Sequence)

    Author: Joseph Kuo 
    """
    def __init__(self, name):
        """Initiation"""
        self.name = name
        self.dict = OrderedDict()

    def add(self, a_region, its_sequence):
        """Add a region with its corresponding sequence"""
        self.dict[a_region] = its_sequence

    def __len__(self):
        """Return the number of region-sequece pair """
        return len(self.dict.keys())

    def __iter__(self):
        """Iteration for the dictionary"""
        return self.dict.iteritems()


if __name__ == '__main__':
    a = RNABinding("a",1,5)
    a_s = Sequence(name="a", seq="ATTGGC")
    b = RNABinding("b",4,8)
    b_s = Sequence(name="b", seq="ATTGGC")
    
    rs = RegionSequenceDict("test")
    rs.add(a, a_s)
    rs.add(b, b_s)

    print(rs)
    print(len(rs))
    