"""
A GenomicRegion describes a genomic region.

Authors: Ivan G. Costa, Manuel Allhoff

Define a GenomicRegion on [initial, final) on a particular chromosome.
The coordinates are 0-based.

Methods:

extend(left, right):
Extend the region to the left and right side.
Negative values are allowed.

overlap(region):
Return true if region overlaps with argument-region, otherwise false.

"""

class GenomicRegion:

    def __init__(self, chrom, initial, final, name=None, orientation=None, data=None):
        """Initialize GenomicRegion"""
        self.chrom = str(chrom) #chrom should be a string, not an integer
        self.initial = initial
        self.final = final
        self.name = name
        self.orientation = orientation
        self.data = data

    def __len__(self):
        """Return length of GenomicRegion"""
        return self.final - self.initial


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

    def toString(self):
        return self.chrom+":"+str(self.initial)+"-"+str(self.final)

    def extend(self, left, right):
        """Extend GenomicRegion both-sided"""
        self.initial -= left
        self.final += right 
        
        #if left, right are negative, switching the border may be necessary
        if self.initial > self.final:
            self.initial, self.final = self.final, self.initial
        
        self.initial = max(self.initial, 0)

    def overlap(self, region):
        """Return True, if GenomicRegion overlaps with region, else False."""
        if region.chrom == self.chrom:
            if self.initial < region.initial:
                if self.final > region.initial:
                    return True
            else:
                if self.initial < region.final:
                    return True
        return False
                    
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        return ','.join( [self.chrom, str(self.initial), str(self.final)] )

    def __cmp__(self, region):
        """Return negative value if x < y, zero if x == y and strictly positive if x > y"""
        if self.chrom < region.chrom:
            return -1
        elif self.chrom > region.chrom:
            return 1
        else:
            if self.initial < region.initial:
                return -1
            elif self.initial > region.final:
                return 1
            else:
                if self.final < region.final:
                    return -1
                elif self.final > region.final:
                    return 1
                else:
                    return 0
