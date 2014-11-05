from __future__ import print_function
from rgt.GenomicVariant import GenomicVariant
from rgt.GenomicRegionSet import GenomicRegionSet
import vcf

"""
Represent list of GenomicVariant.

Authors: Manuel Allhoff

"""

class GenomicVariantSet(GenomicRegionSet):

    def __init__(self, vcf_path = None, name='GenomicVariantSet'):
        """Initialize"""
        GenomicRegionSet.__init__(self, name = name)
        if vcf_path:
            self.read_vcf(vcf_path)
    
    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicVariant.__cmp__)
        self.sorted = True
         
    def read_vcf(self, vcf_path):
        """Read SNPs and InDels form VCF file"""
        self.reader = vcf.Reader(open(vcf_path), 'r')
         
        self.metadata = self.reader.metadata
        self.infos = self.reader.infos
        self.filters = self.reader.filters
        self.formats = self.reader.formats
        self.samples = self.reader.samples
         
        for v in self.reader:
            variant = GenomicVariant(v.CHROM, v.POS, v.REF, v.ALT, v.QUAL, filter = v.FILTER, id = v.ID, \
                                     info = v.INFO, format = v.FORMAT, genotype = v.genotype, samples = v.samples)
            self.add(variant)
    
    def write_vcf(self, vcf_path):
        """Write VCF file"""
        if not self.reader:
            raise Exception("No data available")
         
        writer = vcf.Writer(open(vcf_path, 'w'), self.reader)
        for v in self.sequences:
            record = vcf.model._Record(v.chrom, v.pos, v.id, v.ref, v.alt, v.qual, [], v.info, v.format, [], v.sample)
            writer.write_record(record)
 
    def filter_DB(self):
        """Filter for dbSNP"""
        self.sequences = filter(lambda x: 'DB' not in x.info.keys(), self.sequences)
    
    def filter(self, at, t, op = '>='):
        """Filter for Attributes evalute( at = <MQ, DP> <operator> <t> )"""
        self.sequences = filter(lambda x: eval(str(x.info[at]) + op + str(t)), self.sequences)
    
if __name__ == '__main__':
    s = GenomicVariantSet('/home/manuel/data/humangenetics/01_S1_L001_R1_001.filtered.vcf')
    #b = GenomicVariantSet('/home/manuel/data/humangenetics/K28_S8_L001_R1_001.filtered.vcf')
    #c = s.subtract(b)
    
    print(len(s))
    s.filter('MQ', 40.0)
    #s.filter_DB()
    print(len(s))
    s.write_vcf('/home/manuel/1.vcf')
    