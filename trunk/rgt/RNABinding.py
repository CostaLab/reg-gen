"""
A RNABinding describes a binding region on RNA.

Authors: Joseph Kuo

"""

from rgt.GenomicRegion import GenomicRegion

class RNABinding(GenomicRegion):

    def __init__(self, name, initial, final, score, errors, motif, orientation, guanine_rate):
        """Initialize"""
        GenomicRegion.__init__(self, chrom=name, initial=initial, final=final)
        
        self.name = name
        self.score = score
        self.errors = errors
        self.motif = motif
        self.orientation = orientation
        self.guanine_rate = str(guanine_rate)

    def __str__(self):
        """Give informal string representation"""
        return '-'.join( [self.name, str(self.initial), str(self.final), str(self.score), self.errors, 
                          self.motif, self.orientation, self.guanine_rate] )
        
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        return ','.join( [self.name, str(self.initial), str(self.final), str(self.score), self.errors, 
                          self.motif, self.orientation, self.guanine_rate] )

    def __len__(self):
        """Return the length of the binding site """
        return self.final - self.initial

if __name__ == '__main__':
    a = RNABinding('rna_blablabla', 10, 35, 20, 'd6t8d12', 'M', 'P', 0.3)
    print(a)
    print(a.__repr__())

    print(len(a))
    