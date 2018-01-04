
# Python 3 compatibility
from __future__ import print_function

# Python
import unittest

# Internal
from rgt.MotifSet import MotifSet


# NB: test based on hocomoco
# TODO: must enforce it by only loading hocomoco mtf
class MotifSetTest(unittest.TestCase):
    def test_create_default(self):
        ms = MotifSet()
        self.assertEqual(len(ms.motifs_map), 0, msg="motif dictionary must be empty")

    def test_create_empty(self):
        ms = MotifSet(preload_motifs=False)
        self.assertEqual(len(ms.motifs_map), 0, msg="motif dictionary must be empty")

    def test_create_non_empty(self):
        ms = MotifSet(preload_motifs=True)
        self.assertGreater(len(ms.motifs_map), 0, msg="motif dictionary must be non empty")

    def test_filter_keys_not_list(self):
        ms = MotifSet()

        with self.assertRaises(ValueError):
            ms.filter("test")

    def test_filter_wrong_key_type(self):
        ms = MotifSet()

        with self.assertRaises(ValueError):
            ms.filter([], key_type="test")

    def test_filter_names(self):
        ms = MotifSet(preload_motifs=True)

        ms2, _, _ = ms.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="name", search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2, _, _ = ms.filter(["ALX1"], key_type="name", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2, _, _ = ms.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2, _, _ = ms.filter(["ALX1"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2, _, _ = ms.filter(["ALX"], key_type="name", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)

    def test_filter_genes(self):
        ms = MotifSet(preload_motifs=True)

        ms2, _, _ = ms.filter(["ALX1_HUMAN.H11MO.0.B"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2, _, _ = ms.filter(["ALX1"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2, _, _ = ms.filter(["ALX"], key_type="gene_names", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2, _, _ = ms.filter(["ALX1"], key_type="gene_names", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2, _, _ = ms.filter(["ALX"], key_type="gene_names", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)

    def test_filter_family(self):
        ms = MotifSet(preload_motifs=True)

        ms2, _, _ = ms.filter(["Paired-related HD factors"], key_type="family", search="exact")
        self.assertEqual(len(ms2.motifs_map), 35)

        ms2, _, _ = ms.filter(["factors"], key_type="family", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2, _, _ = ms.filter(["Paired-related HD factors"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)

        ms2, _, _ = ms.filter(["Paired-related HD"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)

        ms2, _, _ = ms.filter(["factors"], key_type="family", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 637)

    def test_filter_uniprot(self):
        ms = MotifSet(preload_motifs=True)

        ms2, _, _ = ms.filter(["Q9H3D4"], key_type="uniprot_ids", search="exact")
        self.assertEqual(len(ms2.motifs_map), 2)

        ms2, _, _ = ms.filter(["Q9H"], key_type="uniprot_ids", search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2, _, _ = ms.filter(["Q9H3D4"], key_type="uniprot_ids", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 2)

        ms2, _, _ = ms.filter(["Q9H"], key_type="uniprot_ids", search="inexact")
        self.assertEqual(len(ms2.motifs_map), 20)
