from __future__ import print_function

from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet

####################################################################################
####################################################################################

class RNABinding(GenomicRegion):
    """Describes a binding region on RNA including the information regarding to this region.

    Authors: Joseph Kuo
    """

    def __init__(self, name, initial, final, score=None, errors_bp=None, motif=None, 
                 orientation=None, guanine_rate=None):
        """Initialize"""
        GenomicRegion.__init__(self, chrom=name, initial=initial, final=final)
        
        self.name = name                      # RNA name
        self.score = score                    # Score for pattern matching
        self.errors_bp = errors_bp                  
        self.motif = motif
        self.orientation = orientation
        self.guanine_rate = str(guanine_rate)

    def __str__(self):
        """Give informal string representation"""
        infos = [ self.name, self.initial, self.final, self.score, self.errors_bp, 
                  self.motif, self.orientation, self.seq, self.guanine_rate ]
        return '-'.join( [str(x) for x in infos if x] )
        
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        infos = [ self.name, self.initial, self.final, self.score, self.errors_bp, 
                  self.motif, self.orientation, self.seq, self.guanine_rate ]
        return ','.join( [str(x) for x in infos if x] ) 
        
    def __len__(self):
        """Return the length of the binding site """
        return self.final - self.initial

    def __hash__(self):
        return hash(tuple([self.initial, self.final]))

####################################################################################
####################################################################################

class RNABindingSet(GenomicRegionSet):
    """Represent a collection of RNABinding with some functions

    Authors: Joseph Kuo
    """

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
        

####################################################################################
#   Test
####################################################################################

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
    