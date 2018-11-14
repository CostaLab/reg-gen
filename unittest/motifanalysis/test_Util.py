# Python 3 compatibility
from __future__ import print_function

# Python
import unittest
import os

# Internal
from rgt.motifanalysis.Util import is_bed, is_bb


class StatisticsTest(unittest.TestCase):

    def test_is_bed(self):
        self.assertTrue(is_bed(os.path.join(os.path.dirname(__file__), "target_regions_mpbs.bed")))

    def test_is_bb(self):
        self.assertTrue(is_bb("test.bb"))