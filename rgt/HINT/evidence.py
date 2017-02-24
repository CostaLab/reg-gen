# Import
import os
import sys

# Internal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import OverlapType

"""
Creates a bed file with MPBSs with (in green) and without (in red) evidence.
Also, the name of the instances will be Y for evidence, N for non evidence.

Authors: Eduardo G. Gusmao, Zhijian Li
"""
class Evidence:

    def __init__(self, peak_ext, mpbs_name, tfbs_summit_fname, mpbs_fname, output_location):
        self.peak_ext = peak_ext
        self.mpbs_name = mpbs_name
        self.tfbs_summit_fname = tfbs_summit_fname
        self.mpbs_fname = mpbs_fname
        self.output_location = output_location

    def create_file(self):
        # Expanding summits
        tfbs_summit_regions = GenomicRegionSet("TFBS Summit Regions")
        tfbs_summit_regions.read_bed(self.tfbs_summit_fname)

        for region in iter(tfbs_summit_regions):
            summit = int(region.data.split()[-1]) + region.initial
            region.initial = max(summit - (self.peak_ext / 2), 0)
            region.final = summit + (self.peak_ext / 2)

        # Calculating intersections
        mpbs_regions = GenomicRegionSet("MPBS Regions")
        mpbs_regions.read_bed(self.mpbs_fname)

        tfbs_summit_regions.sort()
        mpbs_regions.sort()

        with_overlap_regions = mpbs_regions.intersect(tfbs_summit_regions, mode=OverlapType.ORIGINAL)
        without_overlap_regions = mpbs_regions.subtract(tfbs_summit_regions, whole_region=True)
        tfbs_regions = GenomicRegionSet("TFBS Regions")

        for region in iter(with_overlap_regions):
            region.name = region.name.split(":")[0] + ":Y"
            tfbs_regions.add(region)

        for region in iter(without_overlap_regions):
            region.name = region.name.split(":")[0] + ":N"
            tfbs_regions.add(region)

        tfbs_regions.sort()

        tfbs_fname = os.path.join(self.output_location, "{}.bed".format(self.mpbs_name))
        tfbs_regions.write_bed(tfbs_fname)