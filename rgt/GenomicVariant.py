"""
GenomicVariant
===================
GenomicVariant is a specialized GenomicRegion class and describes a SNP or InDel.

"""

from .GenomicRegion import GenomicRegion


class GenomicVariant(GenomicRegion):
    """*Keyword arguments:*

        - chrom -- chromosome
        - pos -- position
        - ref -- reference nucleotide
        - alt -- alternative nucleotide
        - qual -- quality
        - filter -- filter
        - id -- id
        - info -- further informaton
        - format -- SNP format
        - genotype -- genotype
        - samples -- sample
        
    .. note:: all necessary information are contained in a VCF file.
    """

    def __init__(self, chrom, pos, ref, alt, qual, filter=None, id=None, info=None, format=None, genotype=None,
                 samples=None):
        GenomicRegion.__init__(self, chrom, pos, pos + 1)

        self.chrom = str(chrom)
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = format
        self.genotype = genotype
        self.samples = samples

        self.name = self.__str__
        self.data = "_$_".join([str(x) for x in [self.id, self.ref, self.alt, self.qual, self.filter, self.info, self.format,
                                    self.genotype, self.samples]])

    def __str__(self):
        """Give informal string representation"""
        s = '-'.join([self.chrom, str(self.pos)])
        if self.id is not None:
            s += '-' + self.id

        return s

    def __repr__(self):
        """Return official representation of GenomicRegion"""
        return ','.join([self.chrom, str(self.pos)])


if __name__ == '__main__':
    a = GenomicVariant('chr1', 1, "T", "A", 20)
    print(a)
