
# Python 3 compatibility
from __future__ import print_function

# Python
import unittest
import os

# Internal
from rgt.MotifSet import MotifSet


# NB: test based on hocomoco
class MotifSetTest(unittest.TestCase):
    def setUp(self):
        dirname = os.path.dirname(__file__)
        mtf_file = os.path.join(dirname, "../data/motifs/hocomoco.mtf")

        # we must enforce the use hocomoco as database
        self.motif_set = MotifSet(preload_motifs=False)
        self.motif_set.read_mtf(mtf_file)

    def test_create_default(self):

        ms = MotifSet()
        self.assertEqual(len(ms.motifs_map), 0, msg="motif dictionary must be empty by default (no preload)")

    def test_create_empty(self):
        ms = MotifSet(preload_motifs=False)
        self.assertEqual(len(ms.motifs_map), 0, msg="motif dictionary must be empty")

    def test_create_non_empty(self):
        ms = MotifSet(preload_motifs=True)
        self.assertGreater(len(ms.motifs_map), 0, msg="motif dictionary must be non empty")

    def test_filter_keys_not_list(self):
        with self.assertRaises(ValueError):
            self.motif_set.filter("test")

    def test_filter_wrong_key_type(self):
        with self.assertRaises(ValueError):
            self.motif_set.filter([], key_type="test")

    def test_filter_names(self):
        ms2 = self.motif_set.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="name", search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter(["ALX1"], key_type="name", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2 = self.motif_set.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter(["ALX1"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter(["ALX"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)

        ms2 = self.motif_set.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="name", search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter(["ALX1.*"], key_type="name", search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter(["ALX[134]_.*"], key_type="name", search="regex")
        self.assertEqual(len(ms2.motifs_map), 3)

    def test_filter_genes(self):
        ms2 = self.motif_set.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter(["ALX1"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["ALX"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter(["ALX1"], key_type="gene_names", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["ALX"], key_type="gene_names", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

        ms2 = self.motif_set.filter(["ALX1.*"], key_type="gene_names", search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["ALX[134]"], key_type="gene_names", search="regex")
        self.assertEqual(len(ms2.motifs_map), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

    def test_filter_family(self):
        ms2 = self.motif_set.filter(["Paired-related HD factors"], key_type="family", search="exact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["factors"], key_type="family", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter(["Paired-related HD factors"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["Paired-related HD"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["factors"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

        ms2 = self.motif_set.filter(["Paired.*factors"], key_type="family", search="regex")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["Paired-related.*"], key_type="family", search="regex")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter([".*factors"], key_type="family", search="regex")
        self.assertEqual(len(ms2.motifs_map), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

    def test_filter_uniprot(self):
        ms2 = self.motif_set.filter(["Q9H3D4"], key_type="uniprot_ids", search="exact")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["Q9H"], key_type="uniprot_ids", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter(["Q9H3D4"], key_type="uniprot_ids", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["Q9H"], key_type="uniprot_ids", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

        ms2 = self.motif_set.filter(["Q9H3D4"], key_type="uniprot_ids", search="regex")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["Q9H.*"], key_type="uniprot_ids", search="regex")
        self.assertEqual(len(ms2.motifs_map), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

    def test_filter_data_source(self):
        # implicitly, we are also testing the case insensitiveness of the string matching of all three types

        ms2 = self.motif_set.filter(["chip-seq"], key_type="data_source", search="exact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["chip"], key_type="data_source", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter(["chip-seq"], key_type="data_source", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["chip"], key_type="data_source", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["chip-seq"], key_type="data_source", search="regex")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter(["(chip|selex)"], key_type="data_source", search="regex")
        self.assertEqual(len(ms2.motifs_map), 588)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 588)
        self.assertEqual(len(k2m), 2)
