"""
Authors: Ivan G. Costa, Manuel Allhoff

A GenomicRegion describes a genomic region.

Methods:

extend(left, right):
Extend the region to the left and right side.

overlap(region):
Return true if region overlaps with argument-region, otherwise false.

"""

class GenomicRegion:

    def __init__(self, chrom, initial, final, name=None, orientation=None, data=None):
        """Initialize GenomicRegion"""
        self.chrom = chrom
        self.initial = initial
        self.final = final
        self.name = name
        self.orientation = orientation
        self.data = data

    def __len__(self):
        """Return length of GenomicRegion"""
        return self.final - self.initial + 1


    def __str__(self):
        """Give informal string representation"""
        s = '\t'.join( [self.chrom, str(self.initial), str(self.final)] )
        if self.name is not None:
            s += '\t' + self.name
        if self.orientation is not None:
            s += '\t' + self.orientation
        if self.data is not None:
            s += '\t' + self.data
        return s

    def extend(self, left, right):
        """Extend GenmocRegion both-sided"""
        #TODO: What is about orientation?
        self.initial -= left
        self.final += right    

    def overlap(self, region):
        """Return True, if GenomicRegion overlaps with region, else False."""
        if region.chrom == self.chrom:
            if self.initial < region.initial:
                if self.final >= region.initial:
                    return True
                else:
                    if self.initial <= region.final:
                        return True
        return False

                    
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        return ','.join( [self.chrom, str(self.initial), str(self.final)] )

    def __cmp__(self, x, y):
        """Return negative value if x < y, zero if x == y and strictly positive if x > y"""
        if x.chrom < y.chrom:
            return -1
        elif x.chrom > y.chrom:
            return 1
        else:
            if x.initial < y.initial:
                return -1
            elif x.initial > y.final:
                return 1
            else:
                if x.final < y.final:
                    return -1
                elif x.final > y.final:
                    return 1
                else:
                    return 0
