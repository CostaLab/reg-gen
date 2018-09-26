
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

    def test_filter_values_not_dict(self):
        with self.assertRaises(ValueError):
            self.motif_set.filter("test")

    def test_filter_wrong_key_type(self):
        with self.assertRaises(ValueError):
            self.motif_set.filter({'test': []})

    def test_get_mappings_wrong_key_type(self):
        ms = MotifSet()
        with self.assertRaises(ValueError):
            ms.get_mappings("test")

    def test_filter_names(self):
        ms2 = self.motif_set.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)

        ms2 = self.motif_set.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter({'name': ["ALX"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)

        ms2 = self.motif_set.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)

        ms2 = self.motif_set.filter({'name': ["ALX[134]_.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 3)

    def test_filter_genes(self):
        ms2 = self.motif_set.filter({'gene_names': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX[134]"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

    def test_filter_family(self):
        ms2 = self.motif_set.filter({'family': ["Paired-related HD factors"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["factors"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'family': ["Paired-related HD factors"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["Paired-related HD"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["factors"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

        ms2 = self.motif_set.filter({'family': ["Paired.*factors"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["Paired-related.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': [".*factors"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

    def test_filter_uniprot(self):
        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

    def test_filter_data_source(self):
        # implicitly, we are also testing the case insensitiveness of the string matching of all three types

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip"], 'species': ["Homo sapiens", "Mus musculus"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 0)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["(chip|selex)"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2.motifs_map), 588)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 588)
        self.assertEqual(len(k2m), 2)

    def test_filter(self):
        ms2 = self.motif_set.filter({'data_source': ["chip-seq", "integrative"]}, search="exact")
        self.assertEqual(len(ms2.motifs_map), 1138)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 1138)
        self.assertEqual(len(k2m), 2)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq", "integrative"], 'family': ["Steroid hormone receptors (NR3)"], 'species': ["Mus musculus"]}, search="exact")
        self.assertEqual(len(ms2.motif_map), 1140)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 1140)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'family': ["Steroid hormone receptors (NR3)"], 'tax_group': ["vertebrates", "plants"]}, search="exact")
        self.assertEqual(len(ms2.motif_map), 1296)
        m2k, k2m = ms2.get_mappings(key_type="tax_group")
        self.assertEqual(len(m2k), 1296)
        self.assertEqual(len(k2m), 1)