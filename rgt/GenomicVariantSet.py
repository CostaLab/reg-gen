"""
GenomicVariantSet
===================
GenomicVariantSet represents list of GenomicVariant.

"""


from .GenomicVariant import GenomicVariant
from .GenomicRegionSet import GenomicRegionSet
import vcf


class GenomicVariantSet(GenomicRegionSet):
    """*Keyword arguments:*

        - vcf_path -- VCF file
        - name -- name
    """

    def __init__(self, vcf_path=None, name='GenomicVariantSet'):
        """Initialize"""
        GenomicRegionSet.__init__(self, name=name)
        if vcf_path:
            self.read_vcf(vcf_path)

    def sort(self):
        """Sort elements by criteria defined by GenomicVariant.
        
        .. note:: By default, the genomic position is used as sorting criteria.
        
        """
        self.sequences.sort()
        self.sorted = True

    def read_vcf(self, vcf_path):
        """
        Read SNPs and InDels from a VCF file.
        
        *Keyword arguments:*

        - vcf_path -- VCF file
        
        .. note:: vcf_path can also be defined in the initialization.
        
        *Example:*
            We load a VCF file::
            
            >>>from rgt.GenomicVariantSet import GenomicVariantSet
            >>>snps_sample1 = GenomicVariantSet('snps.vcf', name='sample1')
        """

        self.reader = vcf.Reader(open(vcf_path), 'r')

        self.metadata = self.reader.metadata
        self.infos = self.reader.infos
        self.filters = self.reader.filters
        self.formats = self.reader.formats
        self.samples = self.reader.samples

        for v in self.reader:
            variant = GenomicVariant(v.CHROM, v.POS, v.REF, v.ALT, v.QUAL, filter=v.FILTER, id=v.ID, \
                                     info=v.INFO, format=v.FORMAT, genotype=v.genotype, samples=v.samples)
            self.add(variant)

    def write_vcf(self, vcf_path):
        """
        Write VCF file.
        
        *Keyword arguments:*

        - vcf_path -- VCF file
        
        """

        if not self.reader:
            raise Exception("No data available")

        writer = vcf.Writer(open(vcf_path, 'w'), self.reader)
        for v in self.sequences:
            record = vcf.model._Record(v.chrom, v.pos, v.id, v.ref, v.alt, v.qual, [], v.info, v.format, [], v.samples)
            writer.write_record(record)

    def filter_dbSNP(self):
        """Filter for dbSNP.
        
        .. note:: the vcf file must already contain the dbSNP annotation.
        
        """

        self.sequences = [x for x in self.sequences if 'DB' not in list(x.info.keys())]

    def filter(self, at, op, t):
        """
        Filter for attributes. 
        
        *Keyword arguments:*

        - at -- VCF file
        - op -- operation to perform
        - t -- threshold
        
        *Example:*
            We load a VCF file::
            
            >>>from rgt.GenomicVariantSet import GenomicVariantSet
            >>>snps_sample1 = GenomicVariantSet('snps.vcf', name='sample1')
            
            And we filter by the mapping quality::
             
            >>>snps_sample1.filter(at='MQ', op'>', t=30)
            
            The mapping quality is tagged as MQ in the VCF file. We only want to keep SNPs that have a mapping quality higher than 30.
            
            .. note:: operation <op> and threhold <t> depend on the filtering tag <at>
        
        """
        self.sequences = [x for x in self.sequences if eval(str(x.info[at]) + op + str(t))]

    def _reconstruct_info(self, GenomicRegionSet):
        """Reconstruct all information for GenomicVariantSet that get lost when using a GenomicRegionSet method"""
        tmp_sequences = []
        for genomic_region in GenomicRegionSet:
            c, p = genomic_region.chrom, genomic_region.initial
            for var in self.sequences:
                if var.chrom == c and var.pos == p:  # 1:1 mapping
                    tmp_sequences.append(var)
                    break

        return tmp_sequences

    def subtract(self, x):
        """
        Subtract GenomicVariantSet.
        
        *Keyword arguments:*

        - x -- instance of GenomicVariantSet which is subtracted
        
        """
        tmp = GenomicRegionSet.subtract(self, x, whole_region=False)
        self.sequences = self._reconstruct_info(tmp)

    def intersect(self, x):
        """
        Intersect GenomicVariantSet.
        
        *Keyword arguments:*

        - x -- instance of GenomicVariantSet
        
        """

        tmp = self._intersect(x)
        self.sequences = self._reconstruct_info(tmp)

    def _intersect(self, y, rm_duplicates=False):
        """Return the overlapping regions with three different modes.
                   
        (mode = OverlapType.ORIGINAL)
        Return the regions of original GenomicRegionSet which have any intersections with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- the regions of original GenomicRegionSet which have any intersections with y
            
            Graphical explanation:
            self       ----------              ------
            y              ----------                    ----
            Result     ----------
            
        """
        a = self
        b = y

        z = GenomicRegionSet(a.name + ' + ' + b.name)
        # XXX - someone putted an special symbol and spaces in the name! this is used as file name, never use strange characters.
        if len(a) == 0 or len(b) == 0:
            return z

        else:
            # If there is overlap within self or y, they should be merged first. 
            if a.sorted == False: a.sort()
            if b.sorted == False: b.sort()

            iter_a = iter(a)
            s = next(iter_a)
            last_j = len(b) - 1
            j = 0
            cont_loop = True

            ########################### OverlapType.ORIGINAL ###################################
            while cont_loop:
                # print(str(s),"\t",str(b[j]))
                # When the regions overlap
                if s.overlap(b[j]):
                    z.add(s)

                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False

                elif s < b[j]:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False
                elif s > b[j]:
                    if j == last_j:
                        cont_loop = False
                    else:
                        j = j + 1
                else:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False

            if rm_duplicates:
                z.remove_duplicates()

            return z


if __name__ == '__main__':
    s = GenomicVariantSet('/home/manuel/data/humangenetics/01_S1_L001_R1_001.filtered.vcf')
    b = GenomicVariantSet('/home/manuel/data/humangenetics/K28_S8_L001_R1_001.filtered.vcf')
    print(len(s))
    s.subtract(b)

    print(len(s))
    # print(b.sequences[:10])
    # print(c.sequences[:10])
    print(len(s))
    s.filter('MQ', '>=', 40.0)
    print(len(s))
    # s.filter_DB()
    # print(len(s))
    # s.write_vcf('/home/manuel/1.vcf')
