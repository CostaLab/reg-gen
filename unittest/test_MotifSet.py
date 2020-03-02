
# Python
import unittest
import os

# Internal
from rgt.MotifSet import MotifSet


# NB: test based on hocomoco
class MotifSetTest(unittest.TestCase):
    def setUp(self):
        # we must enforce the use hocomoco as database
        self.motif_set = MotifSet(preload_motifs="hocomoco")

    def test_built_in_functions(self):
        ms = MotifSet(preload_motifs="hocomoco")
        self.assertTrue(str(ms).startswith("MotifSet:{"), msg="str(ms): wrong format")
        self.assertTrue(repr(ms) == str(ms), msg="MotifSet: repr does not equal str")
        ms2 = ms.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertTrue("'name': 'ALX1_HUMAN.H11MO.0.B'" in str(ms2), msg="str(ms2): wrong MotifMap")
        self.assertTrue(str(ms2).startswith("MotifSet:{"), msg="str(ms2): wrong format")
        ma = ms2.__getitem__("ALX1_HUMAN.H11MO.0.B")
        self.assertTrue("'thresholds': {0.005: 3.1595, 0.001: 6.52, 0.0005: 7.778, 0.0001: 10.3565, "
                        "5e-05: 11.318, 1e-05: 13.4015}" in str(ma), msg="str(ma): threshold missing")
        self.assertTrue("'name': 'ALX1_HUMAN.H11MO.0.B'" in str(ma), msg="str(ma): wrong Motif")
        self.assertTrue(repr(ma) == str(ma), msg="MotifAnnotation: repr does not equal str")

    def test_create_default(self):
        ms = MotifSet()
        self.assertEqual(len(ms), 0, msg="motif dictionary must be empty by default (no preload)")
        motif_list = ms.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 0)

    def test_create_multiple(self):
        ms = MotifSet(preload_motifs=["hocomoco", "jaspar_vertebrates"])
        self.assertEqual(len(ms), 2042, msg="motif dictionary must contain sum of motifs in files from jaspar"
                                            "_vertebrates and hocomoco")

    def test_create_empty(self):
        ms = MotifSet(preload_motifs=None)
        self.assertEqual(len(ms), 0, msg="motif dictionary must be empty")
        motif_list = ms.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 0)

    def test_create_non_empty(self):
        ms = MotifSet(preload_motifs="hocomoco")
        self.assertGreater(len(ms), 0, msg="motif dictionary must be non empty")

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
        self.assertEqual(len(ms2), 1)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 0)

        ms2 = self.motif_set.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 1)

        ms2 = self.motif_set.filter({'name': ["ALX"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 3)

        ms2 = self.motif_set.filter({'name': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 1)

        ms2 = self.motif_set.filter({'name': ["ALX1.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 1)

        ms2 = self.motif_set.filter({'name': ["ALX[134]_.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 3)

    def test_filter_genes(self):
        ms2 = self.motif_set.filter({'gene_names': ["ALX1_HUMAN.H11MO.0.B"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

        ms2 = self.motif_set.filter({'gene_names': ["ALX1.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 1)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 1)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'gene_names': ["ALX[134]"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 3)
        m2k, k2m = ms2.get_mappings(key_type="gene_names")
        self.assertEqual(len(m2k), 3)
        self.assertEqual(len(k2m), 3)

    def test_filter_family(self):
        ms2 = self.motif_set.filter({'family': ["Paired-related HD factors"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["factors"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'family': ["Paired-related HD factors"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["Paired-related HD"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["factors"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

        ms2 = self.motif_set.filter({'family': ["Paired.*factors"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': ["Paired-related.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 35)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 35)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': [".*factors"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 673)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 673)
        self.assertEqual(len(k2m), 59)

    def test_filter_uniprot(self):
        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H3D4"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 2)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 2)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'uniprot_ids': ["Q9H.*"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 20)
        m2k, k2m = ms2.get_mappings(key_type="uniprot_ids")
        self.assertEqual(len(m2k), 20)
        self.assertEqual(len(k2m), 16)

    def test_filter_data_source(self):
        # implicitly, we are also testing the case insensitiveness of the string matching of all three types

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip"], 'species': ["Homo sapiens", "Mus musculus"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip"], 'species': ["Homo sapiens"]}, search="inexact")
        self.assertEqual(len(ms2), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 431)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 431)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["(chip|selex)"], 'species': ["Homo sapiens"]}, search="regex")
        self.assertEqual(len(ms2), 588)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 588)
        self.assertEqual(len(k2m), 2)

    def test_filter_database(self):
        ms2 = self.motif_set.filter({'database': ["hocomoco"]}, search="exact")
        self.assertEqual(len(ms2), 1296)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 1296)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'database': ["jaspar_vertebrates"]}, search="exact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'database': ["jaspar_vertebrates"]}, search="inexact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'database': ["uniprobe"]}, search="inexact")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'database': ["jaspar"]}, search="regex")
        self.assertEqual(len(ms2), 0)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 0)
        self.assertEqual(len(k2m), 0)

        ms2 = self.motif_set.filter({'database': ["(hocomoco|jaspar)"]}, search="regex")
        self.assertEqual(len(ms2), 1296)
        m2k, k2m = ms2.get_mappings(key_type="database")
        self.assertEqual(len(m2k), 1296)
        self.assertEqual(len(k2m), 1)

    def test_filter(self):
        #test different combinations of key_types and keys

        ms2 = self.motif_set.filter({'data_source': ["chip-seq", "integrative"]}, search="exact")
        self.assertEqual(len(ms2), 1138)
        m2k, k2m = ms2.get_mappings(key_type="data_source")
        self.assertEqual(len(m2k), 1138)
        self.assertEqual(len(k2m), 2)

        ms2 = self.motif_set.filter(
            {'data_source': ["chip-seq", "integrative"], 'family': ["Steroid hormone receptors (NR3)"],
             'species': ["Mus musculus"]}, search="exact")
        self.assertEqual(len(ms2), 14)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 14)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'family': ["Steroid hormone receptors (NR3)"],
                                     'tax_group': ["vertebrates", "plants"]}, search="exact")
        self.assertEqual(len(ms2), 25)
        m2k, k2m = ms2.get_mappings(key_type="tax_group")
        self.assertEqual(len(m2k), 25)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'data_source': ["chip-seq"], 'family': ["Steroid hormone receptors (NR3)"],
                                     'species': ["Mus musculus", "Homo sapiens"]}, search="exact")
        self.assertEqual(len(ms2), 25)
        m2k, k2m = ms2.get_mappings(key_type="species")
        self.assertEqual(len(m2k), 25)
        self.assertEqual(len(k2m), 2)

        ms2 = self.motif_set.filter({'data_source': ["chip"], 'family': ["NR3"], 'tax_group': ["brates"]},
                                    search="inexact")
        self.assertEqual(len(ms2), 25)
        m2k, k2m = ms2.get_mappings(key_type="tax_group")
        self.assertEqual(len(m2k), 25)
        self.assertEqual(len(k2m), 1)

        ms2 = self.motif_set.filter({'family': [".*related factors"]}, search="regex")
        self.assertEqual(len(ms2), 587)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 587)
        self.assertEqual(len(k2m), 36)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 587)

        ms2 = self.motif_set.filter({'data_source': ["(chip|integr)"], 'family': ["multiple"],
                                     'species': ["(musculus|sapiens)"]}, search="regex")
        self.assertEqual(len(ms2), 57)
        m2k, k2m = ms2.get_mappings(key_type="family")
        self.assertEqual(len(m2k), 57)
        self.assertEqual(len(k2m), 1)

    def test_filter_with_empty_dict(self):
        ms2 = self.motif_set.filter({}, search="exact")
        self.assertEqual(len(ms2), 1296)

    def test_create_motif_list(self):
        ms2 = self.motif_set.filter({'name': ["PITX"]}, search="inexact")  # 5 Motifs
        threshold = ms2["PITX2_HUMAN.H11MO.0.D"].thresholds[0.0001]
        # we remove the pre-calculated thresholds so we can test whether the calculation works
        for ma in iter(ms2):
            for fpr in [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]:
                ma.thresholds[fpr] = []
        # is the new threshold equal to the mtf one?
        ml = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(ml), len(ms2))
        self.assertEqual(ml[2].threshold, threshold, msg="create_motif_list calculates threshold incorrectly")
        # is the threshold calculated for non-standard fpr?
        for ma in iter(ms2):
            ma.thresholds = {}
        ml = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(ml[2].threshold, threshold, msg="create_motif_list doesn't work for empty thresholds")
        self.assertEqual(len(ml), len(ms2))


# Test is based on Custom DB
class CustomDBTest(unittest.TestCase):
    def setUp(self):
        # use CustomDB
        self.motif_set = MotifSet(preload_motifs=[os.path.join(os.path.dirname(__file__), "TestCustomDB")],
                                  motif_dbs=True)

    def test_loading(self):
        self.assertEqual(len(self.motif_set.motifs_map), 3, msg="loaded wrong number of motifs")
        self.assertIsNone(self.motif_set.motifs_map["firstMotif_5.0.B"].gene_names, msg="gene_names not None")
        self.assertIsNone(self.motif_set.motifs_map["secondMotif_5.0.B"].data_source, msg="data_source not None")
        self.assertEqual(len(self.motif_set.motifs_map["thirdMotif_5.0.B"].thresholds), 0, msg="thresholds is not an empty dict")

    def test_built_in_functions(self):
        self.assertTrue(str(self.motif_set).startswith("MotifSet:{"), msg="str(ms): wrong format")
        self.assertTrue(repr(self.motif_set) == str(self.motif_set), msg="MotifSet: repr does not equal str")
        ms2 = self.motif_set.filter({'name': ['firstMotif_5.0.B']}, search="exact")
        self.assertTrue("'name': 'firstMotif_5.0.B'" in str(ms2), msg="str(ms2): wrong MotifMap")
        self.assertTrue(str(ms2).startswith("MotifSet:{"), msg="str(ms2): wrong format")
        ma = ms2.__getitem__("firstMotif_5.0.B")
        self.assertTrue("'name': 'firstMotif_5.0.B'" in str(ma), msg="str(ma): wrong Motif")
        self.assertTrue(repr(ma) == str(ma), msg="MotifAnnotation: repr does not equal str")

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
        ms2 = self.motif_set.filter({'name': ["firstMotif_5.0.B"]}, search="exact")
        self.assertEqual(len(ms2), 1)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 1)

        ms2 = self.motif_set.filter({'name': ["secondMotif_5.0.B"]}, search="inexact")
        self.assertEqual(len(ms2), 1)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 1)

        ms2 = self.motif_set.filter({'name': ["thirdMotif_5.0.B"]}, search="regex")
        self.assertEqual(len(ms2), 1)
        motif_list = ms2.get_motif_list(1.0, 0.0001)
        self.assertEqual(len(motif_list), 1)

    def test_filter_database(self):
        ms2 = self.motif_set.filter({'database': ["hocomoco"]}, search="exact")
        self.assertEqual(len(ms2), 0)

        ms2 = self.motif_set.filter({'database': ["TestCustomDB"]}, search="exact")
        self.assertEqual(len(ms2), 3)

        ms2 = self.motif_set.filter({'database': ["TestCustomDB"]}, search="inexact")
        self.assertEqual(len(ms2), 3)

        ms2 = self.motif_set.filter({'database': ["TestCustomDB"]}, search="regex")
        self.assertEqual(len(ms2), 3)

    def test_filter_with_empty_dict(self):
        ms2 = self.motif_set.filter({}, search="exact")
        self.assertEqual(len(ms2), 3)
