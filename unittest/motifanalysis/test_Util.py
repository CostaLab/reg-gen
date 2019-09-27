# Python 3 compatibility
from __future__ import print_function

# Python
import unittest
import os

# Internal
from rgt.motifanalysis.Util import Input, Result, is_bed, is_bb, bed_to_bb, bb_to_bed, parse_filter
from rgt.GeneSet import GeneSet


class UtilTest(unittest.TestCase):

    def test_is_bed(self):
        self.assertTrue(is_bed(os.path.join(os.path.dirname(__file__), "target_regions_mpbs.bed")))
        self.assertFalse(is_bed(os.path.join(os.path.dirname(__file__), "test.bb")))

    def test_is_bb(self):
        self.assertTrue(is_bb("test.bb"))
        self.assertFalse(is_bb(os.path.join(os.path.dirname(__file__), "target_regions_mpbs.bed")))

    def test_bed_to_bb(self):
        # 1. input is already bigBed, 2. input is bed, 3. input is neither bed nor bigBed (Error)
        self.assertEqual("test.bb", bed_to_bb("test.bb", "chrom.sizes.hg19"))
        self.assertEqual("test.bb", bed_to_bb("test.bed", "chrom.sizes.hg19"))
        with self.assertRaises(ValueError):
            bed_to_bb("test.txt", "chrom_sizes.csv")

    def test_bb_to_bed(self):
        # 1. input is already bed, 2. input is bigBed, 3. input is neither bed nor bigBed (Error)
        self.assertEqual("target_regions_mpbs.bed", bb_to_bed("target_regions_mpbs.bed"))
        self.assertEqual("test.bed", bb_to_bed("test.bb"))
        with self.assertRaises(ValueError):
            bb_to_bed("test.txt")

    def test_parse_filter(self):
        # 1. empty filter, 2. normal case, 3. invalid key, 4. invalid format, 5. gene names, 6. names
        self.assertEqual({}, parse_filter(""))
        self.assertEqual({'data_source': ['selex'], 'species': ['mus', 'human'], 'name': ['ERF1', 'CEBPB', 'REST'],
                          'database': ['meme']},
                         parse_filter("species:mus,human;database:meme;name:ERF1,CEBPB,REST;data_source:selex"))
        with self.assertRaises(ValueError):
            parse_filter("wrong_key:test")
        with self.assertRaises(ValueError):
            parse_filter("species:mus,human;database:meme;name:ERF1;CEBPB,REST;data_source:selex")
        self.assertEqual({'gene_names': ['FOXC1', 'HINFP']}, parse_filter("gene_names:FOXC1,HINFP,ZNF282;"
                                                                          "gene_names_file=genes.txt"))
        self.assertEqual({'name': ['MA0014.3.PAX5', 'MA0075.2.Prrx2']}, parse_filter("name:MA0014.3.PAX5,"
                                                                                     "MA0075.2.Prrx2,MA1154.1.ZNF282;"
                                                                                     "name_file=names.txt"))

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