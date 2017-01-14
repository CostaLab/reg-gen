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
from sklearn import metrics
import matplotlib.pyplot as plt
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

        maxPoints = 10000
        fpr_auc_1 = 0.1
        fpr_auc_2 = 0.01

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
        increased_score_mpbs_regions = GenomicRegionSet("Increased Regions")
        intersect_mpbs_regions = mpbs_gen_regions.intersect(pred_footprints_gen_regions, mode=OverlapType.ORIGINAL)
        for region in iter(intersect_mpbs_regions):
            region.data = str(int(region.data) + max_score)
            increased_score_mpbs_regions.add(region)

        without_intersect_regions = GenomicRegionSet("Without Increased Regions")
        without_intersect_regions = mpbs_gen_regions.subtract(pred_footprints_gen_regions, whole_region=True)
        for region in iter(without_intersect_regions):
            increased_score_mpbs_regions.add(region)

        # Evaluate Statistics
        mpbs_name = self.mpbs_fname.split("/")[-1].split(".")[-2]

        increased_score_mpbs_regions.sort_score()
        increased_score_mpbs_regions.write_bed("score.sorted.bed")
        gs = GenomicRegionSet("gs")
        gs.read_bed("gs_o.bed")
        # Calculating receiver operating characteristic curve (ROC)
        # Reading data points
        y_true = list()
        y_score = list()
        for region in iter(gs):
            if str(region.name).split(":")[-1] == "N":
                y_true.append(0)
            else:
                y_true.append(1)
            y_score.append(int(region.data))
        fpr, tpr, thresholds = metrics.roc_curve(np.array(y_true), np.array(y_score))
        roc_auc = metrics.auc(fpr, tpr)
        pr_auc = metrics.average_precision_score(np.array(y_true), np.array(y_score))
        print(pr_auc)