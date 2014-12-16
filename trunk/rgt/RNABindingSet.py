from __future__ import print_function
from rgt.RNABinding import RNABinding
from rgt.GenomicRegionSet import GenomicRegionSet

"""
Represent list of RNABinding.

Authors: Joseph Kuo

"""
class RNABindingSet(GenomicRegionSet):
    def __init__(self, name):
        """Initialize"""
        GenomicRegionSet.__init__(self, name = name)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = RNABindingSet.__cmp__)
        self.sorted = True
    
if __name__ == '__main__':
    a = RNABindingSet(name="a")
    a.add("a",1,5)
    a.add("a",10,15)
    
    b = RNABindingSet(name="b")
    b.add("a",4,8)
    
    print(len(a))
    print(len(b))
    
    a.subtract(b)
    
    print(a.sequences)
    print(b.sequences)
    