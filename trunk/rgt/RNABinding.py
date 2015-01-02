from rgt.GenomicRegion import GenomicRegion

"""
A RNABinding describes a binding region on RNA.

Authors: Joseph Kuo

"""


class RNABinding(GenomicRegion):

    def __init__(self, name, initial, final, score=None, errors=None, motif=None, 
                 orientation=None, seq=None, guanine_rate=None):
        """Initialize"""
        GenomicRegion.__init__(self, chrom=name, initial=initial, final=final)
        
        self.name = name                      # RNA name
        self.score = score                    # Score for pattern matching
        self.errors = errors                  
        self.motif = motif
        self.orientation = orientation
        self.seq = seq.upper()
        self.guanine_rate = str(guanine_rate)

    def __str__(self):
        """Give informal string representation"""
        infos = [ self.name, self.initial, self.final, self.score, self.errors, 
                  self.motif, self.orientation, self.seq, self.guanine_rate ]
        return '-'.join( [str(x) for x in infos if x] )
        
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        infos = [ self.name, self.initial, self.final, self.score, self.errors, 
                  self.motif, self.orientation, self.seq, self.guanine_rate ]
        return ','.join( [str(x) for x in infos if x] ) 
        
    def __len__(self):
        """Return the length of the binding site """
        return self.final - self.initial

    def __hash__(self):
        return hash(tuple([self.initial, self.final]))
"""
if __name__ == '__main__':
    a = RNABinding('rna_blablabla', 10, 35, 20, 'd6t8d12', 'M', 'P', 0.3)
    print(a)
    print(a.__repr__())

    print(len(a))
""" 