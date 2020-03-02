
# Python
import unittest
import numpy as np
import os

# Internal
from rgt.motifanalysis.Statistics import ecdf, multiple_test_correction, fisher_table, get_fisher_dict
from rgt.GenomicRegionSet import GenomicRegionSet


class StatisticsTest(unittest.TestCase):

    def test_ecdf(self):
        pvals = np.arange(0, 200)
        test = ecdf(pvals)
        compare = np.arange(0.005, 1.005, 0.005)
        self.assertTrue(np.allclose(test, compare), msg="ecdf produces wrong output")

    def test_multiple_test_correction(self):
        pvalues = [0.001, 0.005, 0.1, 0.5, 0.01]
        rejected_hypotheses = np.array([True, True, False, False, True])
        corrected_pvalues = np.array([0.005, 0.0125, 0.125, 0.5, 0.01666667])

        rej, cor = multiple_test_correction(pvalues)

        self.assertTrue(np.array_equal(rej, rejected_hypotheses), msg="multiple_test_correction returns wrong list of "
                                                                      "rejected hypotheses")
        self.assertTrue(np.allclose(cor, corrected_pvalues), msg="multiple_test_correction returns wrong list of"
                                                                 " corrected pvalues")

        pvalues = [0.01, 0.05, 0.1, 0.5]
        rejected_hypotheses = np.array([False, False, False, False])
        corrected_pvalues = np.array([0.08333333, 0.20833333, 0.27777777, 1])

        rej, cor = multiple_test_correction(pvalues, method="negcorr")

        self.assertTrue(np.array_equal(rej, rejected_hypotheses), msg="multiple_test_correction(negcorr) returns wrong"
                                                                      " list of rejected hypotheses")
        self.assertTrue(np.allclose(cor, corrected_pvalues), msg="multiple_test_correction(negcorr) returns wrong"
                                                                 " list of corrected pvalues")

        with self.assertRaises(ValueError):
            multiple_test_correction(pvalues, method="some_method")

    def test_multiple_test_correction_using_R(self):
        pvalues_list = [[0.001, 0.005, 0.1, 0.5, 0.01], [0.03, 0.8, 0.47,0.1], [0.0003, 0.4, 0.002]]

        # corrected pvalues calculated by using R function p.adjust
        corrected_pvalues = [np.array([0.005, 0.0125, 0.125, 0.5, 0.01666667]), np.array([0.12, 0.8, 0.62666667, 0.2]),
                             np.array([0.0009, 0.4, 0.003])]

        # get corrected pvalues from own function
        for i in range(0, len(pvalues_list)):
            res = multiple_test_correction(pvalues_list[i])
            self.assertTrue(np.allclose(res[1], corrected_pvalues[i]), msg="multiple_test_correction returns wrong list"
                                                                           "of corrected pvalues")

    def test_fisher_table(self):
        regions = GenomicRegionSet("regions")
        regions.read(os.path.join(os.path.dirname(__file__), "test.bed"))
        mpbs = GenomicRegionSet("mpbs")
        mpbs.read(os.path.join(os.path.dirname(__file__), "test_mpbs.bed"))

        i, ni, gs, ms = fisher_table("GGT1", regions, mpbs)
        self.assertEqual(i, 0)
        self.assertEqual(ni, 36)

        i, ni, gs, ms = fisher_table("HIC2", regions, mpbs)
        self.assertEqual(i, 8)
        self.assertEqual(ni, 28)

        i, ni, gs, ms = fisher_table("RAC2", regions, mpbs, gene_set=True, mpbs_set=True)
        self.assertEqual(len(gs), 0)
        self.assertEqual(len(ms), 0)

    def test_get_fisher_dict(self):
        motif_names = ["TFIP11", "ACO2", "HIC2", "HAT5"]

        regions = GenomicRegionSet("regions")
        regions.read(os.path.join(os.path.dirname(__file__), "test.bed"))
        mpbs = GenomicRegionSet("mpbs")
        mpbs.read(os.path.join(os.path.dirname(__file__), "test_mpbs.bed"))

        regions2 = GenomicRegionSet("regions2")
        regions2.read(os.path.join(os.path.dirname(__file__), "test2.bed"))
        mpbs2 = GenomicRegionSet("mpbs2")
        mpbs2.read(os.path.join(os.path.dirname(__file__), "test2_mpbs.bed"))

        regions3 = GenomicRegionSet("regions3")
        regions3.read(os.path.join(os.path.dirname(__file__), "test3.bed"))
        mpbs3 = GenomicRegionSet("mpbs3")
        mpbs3.read(os.path.join(os.path.dirname(__file__), "test3_mpbs.bed"))

        result = get_fisher_dict(motif_names, regions, mpbs)
        intersecting = result[0]
        not_intersecting = result[1]

        for mn in ["TFIP11", "ACO2", "HAT5"]:
            self.assertEqual(intersecting[mn], 0)

        self.assertEqual(intersecting["HIC2"], 8)

        for mn in motif_names:
            self.assertEqual(intersecting[mn]+not_intersecting[mn], 36)

        result = get_fisher_dict(motif_names, regions, mpbs, gene_set=True, mpbs_set=True)
        gs = result[2]
        ms = result[3]

        for mn in ["TFIP11", "ACO2", "HAT5"]:
            self.assertEqual(len(gs[mn]), 0)

        self.assertEqual(len(gs["HIC2"]), 8)
        self.assertEqual(len(ms["HIC2"]), 8)

        # test whether both approaches lead to the same result:
        # safe: for each input file subtract target regions(manually or with GenomicRegionSet.subtract()) from
        # background regions and compute fisher dict for smaller background region (background_tmp)
        # fast: compute fisher dict for background region and subtract respective counts for each input file
        # from counts inside background fisher dict

        # TODO does not work because there are target regions which are not part of the background regions

        # motif_names = ["LARGE", "MCAT", "PDXP", "HIC2"]
        #
        # for target_regions in [regions2]:
        #
        #     background_tmp = regions.subtract(target_regions, whole_region=False, merge=False)
        #
        #     self.assertEqual(len(regions.sequences)-len(target_regions.sequences), len(background_tmp.sequences),
        #                      msg=str(len(regions.sequences)-len(target_regions.sequences)) + " != " +
        #                      str(len(background_tmp.sequences)) + ", " + str(target_regions.name))
        #
        #     bg_c_dict, bg_d_dict, _, _ = get_fisher_dict(motif_names, regions, mpbs)  # regions = background
        #     bg_a_dict, bg_b_dict, _, _ = get_fisher_dict(motif_names, target_regions, mpbs)  # target regions
        #     bg_c_tmp_dict, bg_d_tmp_dict, _, _ = get_fisher_dict(motif_names, background_tmp, mpbs)
        #
        #     for gene in bg_c_tmp_dict:
        #         self.assertEqual(bg_c_tmp_dict[gene], bg_c_dict[gene]-bg_a_dict[gene])
        #
        #     for gene in bg_d_tmp_dict:
        #         self.assertEqual(bg_d_tmp_dict[gene], bg_d_dict[gene]-bg_b_dict[gene])
