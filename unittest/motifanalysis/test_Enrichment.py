
# Python
import unittest
import os

# Internal
from rgt.GenomicRegionSet import GenomicRegionSet


class EnrichmentTest(unittest.TestCase):

    def test_subtract_exact(self):
        reference = GenomicRegionSet("reference")
        reference.read(os.path.join(os.path.dirname(__file__), "test_result.bed"))
        background = GenomicRegionSet("background")
        background.read(os.path.join(os.path.dirname(__file__), "test_background.bed"))
        target = GenomicRegionSet("target")
        target.read(os.path.join(os.path.dirname(__file__), "test_target.bed"))

        background_tmp = background.subtract(target, exact=True)

        reference.sort()

        self.assertEqual(len(background_tmp.sequences), len(reference.sequences))

        for region, region_ref in zip(background_tmp.sequences, reference.sequences):
            self.assertEqual(region, region_ref)
