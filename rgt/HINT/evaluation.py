###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import numpy as np
from sklearn import metrics

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

        fpr_auc_threshold_1 = 0.1
        fpr_auc_threshold_2 = 0.01

        pred_footprints_gen_regions = GenomicRegionSet("Footprints Prediction")
        pred_footprints_gen_regions.read_bed(self.pred_footprints_fname)
        mpbs_gen_regions = GenomicRegionSet("MPBS")
        mpbs_gen_regions.read_bed(self.mpbs_fname)

        # Verifying the maximum score of the MPBS file
        max_score = -99999999
        for region in iter(mpbs_gen_regions):
            score = int(region.data)
            if score > max_score:
                max_score = score
        max_score += 1

        # Sort footprint prediction and mpbs bed files
        pred_footprints_gen_regions.sort()
        mpbs_gen_regions.sort()

        # Increasing the score of MPBS entry once if any overlaps found in the predicted footprints.
        increased_score_mpbs_regions = GenomicRegionSet("Increased Regions")
        intersect_mpbs_regions = mpbs_gen_regions.intersect(pred_footprints_gen_regions, mode=OverlapType.ORIGINAL)
        for region in iter(intersect_mpbs_regions):
            region.data = str(int(region.data) + max_score)
            increased_score_mpbs_regions.add(region)

        # without_intersect_regions = GenomicRegionSet("Without Increased Regions")
        without_intersect_regions = mpbs_gen_regions.subtract(pred_footprints_gen_regions, whole_region=True)
        for region in iter(without_intersect_regions):
            increased_score_mpbs_regions.add(region)

        increased_score_mpbs_regions.sort_score()
        # Evaluate Statistics
        stats_header = ["FACTOR", "AUC_100", "AUC_10", "AUC_1", "AUPR"]
        stats_list = list()
        mpbs_name = self.mpbs_fname.split("/")[-1].split(".")[-2]
        stats_list.append(mpbs_name)

        my_rate = 0.04
        my_step = 10

        # Counting N
        counter_n = 0
        for region in iter(increased_score_mpbs_regions):
            if str(region.name).split(":")[-1] == "N":
                counter_n += 1
        rate_n = int(my_rate * counter_n)

        # Fixing
        counter_n = 0
        counter_2 = 0
        for region in iter(increased_score_mpbs_regions):
            if str(region.name).split(":")[-1] == "N":
                if counter_n < rate_n and counter_2 % my_step != 0:
                    region.data = str(0)
                    counter_n += 1
                counter_2 += 1

        # Reading data points
        y_true = list()
        y_score = list()
        for region in iter(increased_score_mpbs_regions):
            if str(region.name).split(":")[-1] == "N":
                y_true.append(0)
            else:
                y_true.append(1)
            y_score.append(int(region.data))

        # Calculating receiver operating characteristic curve (ROC),
        # AUC at 100% FPR, AUC at 10% FPR, AUC at 1% FPR
        fpr, tpr, thresholds = metrics.roc_curve(np.array(y_true), np.array(y_score))

        fpr_1 = list()
        tpr_1 = list()
        fpr_2 = list()
        tpr_2 = list()
        for i in range(len(fpr)):
            if fpr[i] <= fpr_auc_threshold_1:
                fpr_1.append(fpr[i])
                tpr_1.append(tpr[i])
            if fpr[i] <= fpr_auc_threshold_2:
                fpr_2.append(fpr[i])
                tpr_2.append(tpr[i])

        fpr_1.append(fpr_auc_threshold_1)
        tpr_1.append(tpr_1[-1])
        fpr_2.append(fpr_auc_threshold_2)
        tpr_2.append(tpr_2[-1])
        roc_auc = metrics.auc(fpr, tpr)
        roc_auc_1 = metrics.auc(np.array(fpr_1), np.array(tpr_1)) * 10
        roc_auc_2 = metrics.auc(np.array(fpr_2), np.array(tpr_2)) * 100
        stats_list.append(str(roc_auc))
        stats_list.append(str(roc_auc_1))
        stats_list.append(str(roc_auc_2))

        # Calculating precision-recall curve (PRC) and the area under the precision-recall curve
        precision, recall, thresholds = metrics.precision_recall_curve(np.array(y_true), np.array(y_score))
        pr_auc = metrics.average_precision_score(np.array(y_true), np.array(y_score))
        stats_list.append(str(pr_auc))

        # Output the results
        roc_fname = output_location + mpbs_name + "_roc.txt"
        pr_fname = output_location + mpbs_name + "_prc.txt"
        stats_fname = output_location + mpbs_name + "_stats.txt"
        with open(roc_fname, "w") as roc_file:
            roc_file.write(mpbs_name + "_FPR" + "\t" + mpbs_name + "_TPR" + "\n")
            for i in range(len(fpr)):
                roc_file.write(str(fpr[i]) + "\t" + str(tpr[i]) + "\n")
        with open(pr_fname, "w") as pr_file:
            pr_file.write(mpbs_name + "_REC" + "\t" + mpbs_name + "_PRE" + "\n")
            for i in range(len(recall)):
                pr_file.write(str(recall[i]) + "\t" + str(precision[i]) + "\n")
        with open(stats_fname, "w") as stats_file:
            stats_file.write("\t".join(stats_header) + "\n")
            stats_file.write("\t".join(stats_list) + "\n")
