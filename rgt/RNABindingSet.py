from __future__ import print_function
from RNABinding import *
from rgt.GenomicRegionSet import GenomicRegionSet

"""
Represent list of RNABinding.

Authors: Joseph Kuo

"""
class RNABindingSet(GenomicRegionSet):
    """A collection of RNABinding with some functions"""
    
    def __init__(self, name):
        """Initialize"""
        GenomicRegionSet.__init__(self, name = name)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicRegion.__cmp__)
        self.sorted = True
    
    def write_rbs(self, filename):
        """Write the information into a file with .rbs """
        with open(filename+".rbs", "w") as f:
            print("# Sequence-ID\tStart\tEnd\tScore\tMotif\tError-rate\tSequence",file=f)
            for tfo in self.sequences:
                err_rate = 1 - tfo.score/(tfo.final - tfo.initial)
                print("\t".join([tfo.name, str(tfo.initial), str(tfo.final), str(tfo.score), 
                                 tfo.motif, "{0:.2f}".format(err_rate), tfo.seq]),file=f)
    def concatenate(self, another_RNABindingSet):
        """Concatenate another RNABindingSet without sorting"""
        self.sequences = self.sequences + another_RNABindingSet.sequences
        

if __name__ == '__main__':
    a = RNABindingSet(name="a")
    a.add(RNABinding("a",1,5))
    a.add(RNABinding("a",10,15))
    
    b = RNABindingSet(name="b")
    b.add(RNABinding("b",4,8))
    
    print(len(a))
    print(len(b))
    print(a.sequences)
    print(b.sequences)
    
    a.subtract(b)
    a.merge()
    print(a.sequences)
    print(b.sequences)
    