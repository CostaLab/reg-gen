"""
GenomicRegion
===================
GenomicRegion describes a genomic region.

"""
class GenomicRegion:
    """*Keyword arguments:*

            - chrom -- Chromosome.
            - initial -- Start position
            - final -- End position
            - name -- Name of the region
            - orientation -- Orientation of the region, "+" or "-"
            - data -- Extra information
            - proximity -- Close genes
    """
    #__slots__ = ['chrom', 'initial', 'final', 'name', 'orientation', 'data', 'proximity']

    def __init__(self, chrom, initial, final, name=None, orientation=None, data=None, proximity = None):
        self.chrom = str(chrom) #chrom should be a string, not an integer
        if not isinstance( initial, int ) or not isinstance(final, int):
            raise ValueError('The initial and final input for GenomicRegion should be integer.')
        self.initial = initial
        self.final = final
        self.name = name
        self.orientation = orientation
        self.data = data # data can be integer, string, list
        self.proximity = proximity

    def get_data(self, as_list=False):
        """Return data as string (with special separating character (_$_)) or as list.

        *Keyword arguments:*
        
            - as_list -- Return a list instead of a string.
        """
        if not as_list:
            return self.data
        else:
            tmp = self.data.split("_$_")
            return tmp

    def __len__(self):
        """Return length of GenomicRegion."""
        return self.final - self.initial

    def __str__(self):
        """Give informal string representation."""
        s = '\t'.join( [self.chrom, str(self.initial), str(self.final)] )
        # Name
        if self.name is not None: s += '\t' + self.name
        else: s += '\t' + self.toString()
        # Score
        if not self.data: score = "."
        else:
            data = self.data.split("\t")
            try: 
                score = data[0]
            except: score = "."
        s += '\t' + score
        # Orientation
        if self.orientation: s += '\t' + self.orientation
        else: s += '\t' + "."
        # Else
        if self.data:
            s += '\t' + "\t".join(data[1:])
        return s

    def __hash__(self):
        return hash(tuple([self.chrom, self.initial, self.final, self.orientation]))

    def __eq__(self, other):
        return (self.chrom, self.initial, self.final, self.orientation) == (other.chrom, other.initial, other.final, other.orientation)

    def toString(self, space=False):
        """Return a string of GenomicRegion by its position.

        *Keyword arguments:*

            - space -- insert spaces between the values.
        """
        if space:
            return "chr "+self.chrom[3:]+": "+str(self.initial)+"-"+str(self.final)
        else:
            return self.chrom+":"+str(self.initial)+"-"+str(self.final)

    def extend(self, left, right):
        """Extend GenomicRegion on both sides.

        *Keyword arguments:*
        
            - left -- Define the length to extend on left.
            - right -- Define the length to extend on right.
        """
        self.initial -= left
        self.final += right 
        
        #if left, right are negative, switching the border may be necessary
        if self.initial > self.final:
            self.initial, self.final = self.final, self.initial
        
        self.initial = max(self.initial, 0)

    def overlap(self, region):
        """Return True, if GenomicRegion overlaps with region, else False.

        *Keyword arguments:*

            - region -- Given GenomicRegion to compare.
        """
        if region.chrom == self.chrom:
            if self.initial <= region.initial:
                if self.final > region.initial:
                    return True
                #elif self.initial == self.final:
                #    raise Exception(self.chrom+","+str(self.initial)+","+str(self.final)+"\tThe region length shouldn't be zero. Please extend the region.")
            else:
                if self.initial < region.final:
                    return True
        return False
                    
    def __repr__(self):
        """Return official representation of GenomicRegion."""
        return ','.join( [self.chrom, str(self.initial), str(self.final)] )

    def __cmp__(self, region):
        """Return negative value if x < y, zero if x == y and strictly positive if x > y.

        *Keyword arguments:*

            - region -- Given GenomicRegion to compare.
        """
        if self.chrom < region.chrom:
            return -1
        elif self.chrom > region.chrom:
            return 1
        else:
            if self.initial < region.initial:
                return -1
            elif self.initial > region.initial:
                return 1
            else:
                if self.final < region.final:
                    return -1
                elif self.final > region.final:
                    return 1
                else:
                    return 0

    def extract_blocks(self):
        """Extract the block information in self.data into a GenomicRegionSet."""
        z = []
        data = self.data.split("\t")
        nexon = int(data[4])
        width = data[5].split(",")
        posit = data[6].split(",")
        reverse = False
        for i in range(nexon):
            if self.orientation == "-":
                n = "_exon_"+str(nexon - i)
                reverse = True
            else:
                n = "_exon_"+str(i + 1)
                
            z.append(GenomicRegion(chrom=self.chrom, 
                                   initial=self.initial+int(posit[i]), 
                                   final=self.initial+int(posit[i])+int(width[i]), 
                                   name=self.name+n, orientation=self.orientation, data=data[0] ))
        if reverse: return z[::-1]
        else: return z

    def distance(self, y):
        """Return the distance between two GenomicRegions. If overlapping, return 0; if on different chromosomes, return None.

        *Keyword arguments:*

            - y -- Given GenomicRegion to compare.
        """
        
        if self.chrom == y.chrom:
            if self.overlap(y):
                return 0
            elif self < y:
                return y.initial - self.final
            else:
                return self.initial - y.final
        else:
            return None
