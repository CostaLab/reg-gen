import unittest

from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegionSet import *

regions = GenomicRegionSet("test")
regions.add(GenomicRegion("chr1", 10000, 11000, "+"))
regions.add(GenomicRegion("chr1", 20000, 21000, "-"))

cov = CoverageSet("coverage", regions)

bamfile = "/projects/lncRNA/local/cardio/total_rna/bam/d4_1.bam"
bedfile = "~/rgtdata/hg38/genes_hg38.bed"


class CoverageSetTest(unittest.TestCase):
    def test_coverage_from_genomicset(self):
        cov.coverage_from_genomicset(bamfile)
        print(cov.coverage)
        self.assertEqual(cov.coverage, 4)
