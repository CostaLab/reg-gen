# Python 3 compatibility


# Python
import unittest
import os

# Internal
from rgt.motifanalysis.Util import Input, Result, is_bed, is_bb, bed_to_bb, bb_to_bed
from rgt.GeneSet import GeneSet


class UtilTest(unittest.TestCase):

    def test_is_bed(self):
        self.assertTrue(is_bed(os.path.join(os.path.dirname(__file__), "target_regions_mpbs.bed")))
        self.assertFalse(is_bed(os.path.join(os.path.dirname(__file__), "test.bb")))

    def test_is_bb(self):
        self.assertTrue(is_bb("test.bb"))
        self.assertFalse(is_bb(os.path.join(os.path.dirname(__file__), "target_regions_mpbs.bed")))

    def test_bed_to_bb(self):
        self.assertEqual("test.bb", bed_to_bb("test.bb", "chrom_sizes.csv"))
        with self.assertRaises(ValueError):
            bed_to_bb("test.txt", "chrom_sizes.csv")

    def test_bb_to_bed(self):
        self.assertEqual("target_regions_mpbs.bed", bb_to_bed("target_regions_mpbs.bed"))
        with self.assertRaises(ValueError):
            bb_to_bed("test.txt")

    def test_result(self):
        res = Result()
        self.assertEqual(res.name, "")
        self.assertEqual(res.p_value, 0.0)
        self.assertEqual(res.corr_p_value, 0.0)
        self.assertEqual(res.a, 0)
        self.assertEqual(res.b, 0)
        self.assertEqual(res.c, 0)
        self.assertEqual(res.d, 0)
        self.assertEqual(res.percent, 0.0)
        self.assertEqual(res.back_percent, 0.0)
        self.assertEqual(res.genes, [])
        self.assertEqual(str(res), "\t".join(["", "0.0", "0.0", "0", "0", "0", "0", "0.0", "0.0", ""]))

    def test_input(self):
        gene_set = GeneSet("gene_set")
        region_list = []
        test_input = Input(gene_set, region_list)

        self.assertEqual(test_input.gene_set, gene_set)
        self.assertEqual(test_input.region_list, region_list)