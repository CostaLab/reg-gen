# Python 3 compatibility
from __future__ import print_function

# Python
import unittest
import os

# Internal
from rgt.motifanalysis.Statistics import ecdf, multiple_test_correction, get_fisher_dict
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet

# External
from numpy import any, all


class StatisticsTest(unittest.TestCase):

    def test_ecdf(self):
        p_vals = [0.01, 0.05, 0.15, 0.005]
        self.assertTrue(all(ecdf(p_vals)) in [0.25, 0.5, 0.75, 1.0], msg="ecdf does not work correctly")
        self.assertTrue(any(ecdf(p_vals) == 0.25), msg="ecdf result is not complete")
        self.assertTrue(any(ecdf(p_vals) == 0.5), msg="ecdf result is not complete")
        self.assertTrue(any(ecdf(p_vals) == 0.75), msg="ecdf result is not complete")
        self.assertTrue(any(ecdf(p_vals) == 1.0), msg="ecdf result is not complete")

    def test_multiple_test_correction(self):
        # TODO finish test
        p_vals = [0.01, 0.05]
        a, b = multiple_test_correction(p_vals)
        self.assertTrue(any(a), msg="multiple test correlation returns wrong list of rejected hypotheses")
        self.assertFalse(all(a), msg="multiple test correlation returns wrong list of rejected hypotheses")
