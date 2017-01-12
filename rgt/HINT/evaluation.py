
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import math
import operator
import numpy as np
from pickle import load
from sklearn.metrics import auc
from scipy.integrate import simps, trapz
from optparse import OptionParser,BadOptionError,AmbiguousOptionError

"""
Evaluate the footprints prediction using TF ChIP-seq or expression data.

Authors: Eduardo G. Gusmao.
"""

class Evaluation:
    """
    Contains two different methodologies: TF ChIP-seq Based Evaluation
    Usage:
    1. Initialize class.
    2. Call load_sg_coefs once.
    3. Call get_signal as many times as needed.

    """

    def __init__(self, mpbs_fname, footprints_fname):
        self.mpbs_fname = mpbs_fname
        self.footprints_fname = footprints_fname

    # Function standardize
    def standardize(vec):
        maxN = max(vec)
        minN = min(vec)
        return [(e-minN)/(maxN-minN) for e in vec]

    def chip_evaluate(self):
        """
        This evaluation methodology uses motif-predicted binding sites (MPBSs) together with TF ChIP-seq data
        to evaluate the footprint predictions.

        return:
        """

        # Verifying the maximum score of the MPBS file
        max_score = -99999999
        mpbs_file = open(self.mpbs_fname, "r")
        for line in mpbs_file:
            ll = line.strip().split("\t")
            score = int(ll[4])
            if score > max_score:
                max_score = score
        max_score = max_score + 1

        # Sort prediction and mpbs bed files


