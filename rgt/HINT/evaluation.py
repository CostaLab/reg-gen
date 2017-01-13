###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import os
import sys
import math
import operator
import numpy as np
from pickle import load
from sklearn.metrics import auc
from scipy.integrate import simps, trapz
from optparse import OptionParser, BadOptionError, AmbiguousOptionError


# Internal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.SequenceSet import *
from rgt.Util import GenomeData, OverlapType, AuxiliaryFunctions, Library_path



"""
Evaluate the footprints prediction using TF ChIP-seq or expression data.

Authors: Eduardo G. Gusmao, Zhijian Li
"""


class Evaluation:
    """
    Contains two different methodologies: TF ChIP-seq Based Evaluation
    Usage:
    1. Initialize class.
    2. Call load_sg_coefs once.
    3. Call get_signal as many times as needed.

    """

    def __init__(self, pred_footprints_fname, mpbs_fname):
        self.pred_footprints_fname = pred_footprints_fname
        self.mpbs_fname = mpbs_fname

    # Function standardize
    def standardize(vec):
        maxN = max(vec)
        minN = min(vec)
        return [(e - minN) / (maxN - minN) for e in vec]

    def chip_evaluate(self, output_location):
        """
        This evaluation methodology uses motif-predicted binding sites (MPBSs) together with TF ChIP-seq data
        to evaluate the footprint predictions.

        return:
        """

        pred_footprints_gen_regions = GenomicRegionSet("Footprints Prediction")
        pred_footprints_gen_regions.read_bed(self.pred_footprints_fname)
        mpbs_gen_regions = GenomicRegionSet("MPBS")
        mpbs_gen_regions.read_bed(self.mpbs_fname)


        # Verifying the maximum score of the MPBS file
        max_score = -99999999
        for region in iter(mpbs_gen_regions):
            score = int(str(region).split("\t")[4])
            if score > max_score:
                max_score = score
        max_score = max_score + 1

        # Sort footprint prediction and mpbs bed files
        pred_footprints_gen_regions.sort()
        mpbs_gen_regions.sort()

        # Increasing the score of MPBS entry once if any overlaps found in the predicted footprints.
        increased_score_mpbs_regions = GenomicRegionSet("Increased score mpbs")
        intersect_mpbs_regions = mpbs_gen_regions.intersect(pred_footprints_gen_regions, mode = OverlapType.ORIGINAL)
        for region in iter(mpbs_gen_regions):
            if intersect_mpbs_regions.include(region):
                region.data = str(int(region.data) + max_score)
            increased_score_mpbs_regions.add(region)
        increased_score_mpbs_regions.write_bed("score.bed")


