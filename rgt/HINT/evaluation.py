###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import numpy as np
from sklearn import metrics
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pylab
# Internal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import OverlapType

"""
Evaluate the footprints prediction using TF ChIP-seq or expression data.

Authors: Eduardo G. Gusmao, Zhijian Li
"""


class Evaluation:
    """
    Contains two different methodologies: TF ChIP-seq Based Evaluation

    """

    def __init__(self, mpbs_file, footprint_file, footprint_name, print_roc_curve, print_pr_curve, output_location):
        self.mpbs_fname = mpbs_file
        self.footprint_fname = footprint_file.split(",")
        self.footprint_name = footprint_name.split(",")
        self.print_roc_curve = print_roc_curve
        self.print_pr_curve = print_pr_curve
        self.output_location = output_location
        if self.output_location[-1] != "/":
            self.output_location += "/"

    def chip_evaluate(self, alignment_file):
        """
        This evaluation methodology uses motif-predicted binding sites (MPBSs) together with TF ChIP-seq data
        to evaluate the footprint predictions.

        return:
        """

        fpr_auc_threshold_1 = 0.1
        fpr_auc_threshold_2 = 0.01

        mpbs_gen_regions = GenomicRegionSet("MPBS")
        mpbs_gen_regions.read_bed(self.mpbs_fname)
        mpbs_gen_regions.sort()

        mpbs_name = self.mpbs_fname.split("/")[-1].split(".")[-2]

        # Verifying the maximum score of the MPBS file
        max_score = -99999999
        for region in iter(mpbs_gen_regions):
            score = int(region.data)
            if score > max_score:
                max_score = score
        max_score += 1

        # Evaluate Statistics
        stats_header = ["METHOD", "AUC_100", "AUC_10", "AUC_1", "AUPR"]
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        roc_auc_1 = dict()
        roc_auc_2 = dict()
        recall = dict()
        precision = dict()
        prc_auc = dict()
        for i in range(len(self.footprint_fname)):
            footprints_gen_regions = GenomicRegionSet("Footprints Prediction")
            footprints_gen_regions.read_bed(self.footprint_fname[i])

            # Sort footprint prediction bed files
            footprints_gen_regions.sort()

            ################################################
            # Extend 10 bp for all methods
            for region in iter(footprints_gen_regions):
                mid = (region.initial + region.final) / 2
                region.initial = max(0, mid - 5)
                region.final = mid + 5
            #################################################

            # Increasing the score of MPBS entry once if any overlaps found in the predicted footprints.
            increased_score_mpbs_regions = GenomicRegionSet("Increased Regions")
            intersect_regions = mpbs_gen_regions.intersect(footprints_gen_regions, mode=OverlapType.ORIGINAL)
            for region in iter(intersect_regions):
                region.data = str(int(region.data) + max_score)
                increased_score_mpbs_regions.add(region)

            # Keep the score of remained MPBS entry unchanged
            without_intersect_regions = mpbs_gen_regions.subtract(footprints_gen_regions, whole_region=True)
            for region in iter(without_intersect_regions):
                increased_score_mpbs_regions.add(region)

            increased_score_mpbs_regions.sort_score()

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
            fpr[i], tpr[i], roc_auc[i] = self.compute_fpr_tpr(y_true, y_score, 1.0)
            _, _, roc_auc_1[i] = self.compute_fpr_tpr(y_true, y_score, fpr_auc_threshold_1)
            _, _, roc_auc_2[i] = self.compute_fpr_tpr(y_true, y_score, fpr_auc_threshold_2)

            # Calculating precision-recall curve (PRC) and the area under the precision-recall curve
            recall[i], precision[i], prc_auc[i] = self.compute_precision_recall(y_true, y_score, 1.0)

        # Output the statistics results into text
        stats_fname = self.output_location + mpbs_name + "_stats.txt"
        with open(stats_fname, "w") as stats_file:
            stats_file.write("\t".join(stats_header) + "\n")
            for i in range(len(self.footprint_name)):
                stats_file.write(self.footprint_name[i] + "\t" + str(roc_auc[i]) + "\t" + str(roc_auc_1[i]) + "\t"
                                 + str(roc_auc_2[i]) + "\t" + str(prc_auc[i]) + "\n")

        # Output the curves
        if self.print_roc_curve:
            label_x = "False Positive Rate"
            label_y = "True Positive Rate"
            curve_name = "ROC"
            self.plot_curve(fpr, tpr, roc_auc, label_x, label_y, mpbs_name, curve_name)
        if self.print_pr_curve:
            label_x = "Recall"
            label_y = "Precision"
            curve_name = "PRC"
            self.plot_curve(recall, precision, prc_auc, label_x, label_y, mpbs_name, curve_name)

    def plot_curve(self, data_x, data_y, stats, label_x, label_y, mpbs_name, curve_name):
        color_list = ["#000000", "#000099", "#006600", "#990000", "#660099", "#CC00CC", "#222222", "#CC9900",
                      "#FF6600", "#0000CC", "#336633", "#CC0000", "#6600CC", "#FF00FF", "#555555", "#CCCC00",
                      "#FF9900", "#0000FF", "#33CC33", "#FF0000", "#663399", "#FF33FF", "#888888", "#FFCC00",
                      "#663300", "#009999", "#66CC66", "#FF3333", "#9933FF", "#FF66FF", "#AAAAAA", "#FFCC33",
                      "#993300", "#00FFFF", "#99FF33", "#FF6666", "#CC99FF", "#FF99FF", "#CCCCCC", "#FFFF00"]

        matplotlib.use('Agg')
        # Creating figure
        fig = plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
        ax = fig.add_subplot(111)
        for i in range(len(self.footprint_name)):
            ax.plot(data_x[i], data_y[i], color=color_list[i],
                    label=self.footprint_name[i] + ": " + str(round(stats[i], 4)))

        # Plot red diagonal line
        ax.plot([0, 1.0], [0, 1.0], color="#CCCCCC", linestyle="--", alpha=1.0)

        # Line legend
        leg = ax.legend(bbox_to_anchor=(1.0, 0.0, 1.0, 1.0), loc=2, ncol=2, borderaxespad=0., title="AUC",
                        prop={'size': 4})
        for e in leg.legendHandles:
            e.set_linewidth(2.0)

        # Titles and Axis Labels
        ax.set_title(mpbs_name)
        ax.set_xlabel(label_x)
        ax.set_ylabel(label_y)

        # Ticks
        ax.grid(True, which='both')
        ax.set_xticks(np.arange(0.0, 1.01, 0.1))
        ax.set_yticks(np.arange(0.0, 1.01, 0.1))
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(10)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(10)

        # Axis limits
        pylab.xlim([0, 1.0])
        pylab.ylim([0, 1.0])

        # Saving figure
        figure_name = self.output_location + mpbs_name + "_" + curve_name + ".png"
        fig.savefig(figure_name, format="png", dpi=300, bbox_inches='tight', bbox_extra_artists=[leg])

    def compute_fpr_tpr(self, y_true, y_score, fpr_cutoff):
        fpr, tpr, _ = metrics.roc_curve(np.array(y_true), np.array(y_score))
        fpr_at_cutoff = list()
        tpr_at_cutoff = list()
        for idx in range(len(fpr)):
            if fpr[idx] <= fpr_cutoff:
                fpr_at_cutoff.append(fpr[idx])
                tpr_at_cutoff.append(tpr[idx])
        fpr_at_cutoff.append(fpr_cutoff)
        tpr_at_cutoff.append(tpr_at_cutoff[-1])
        scale = 1 / fpr_cutoff
        auc_at_cutoff = metrics.auc(np.array(fpr_at_cutoff), np.array(tpr_at_cutoff)) * scale
        return fpr, tpr, auc_at_cutoff

    def compute_precision_recall(self, y_true, y_score, fdr_cutoff):
        precision, recall, _ = metrics.precision_recall_curve(np.array(y_true), np.array(y_score))
        precision_at_cutoff = list()
        recall_at_cutoff = list()
        for idx in range(len(precision)):
            fdr = 1 - precision[idx]
            if fdr <= fdr_cutoff:
                precision_at_cutoff.append(precision[idx])
                recall_at_cutoff.append(recall[idx])
        precision_at_cutoff.append(fdr_cutoff)
        recall_at_cutoff.append(recall[-1])
        auc_at_cutoff = metrics.auc(np.array(precision_at_cutoff), np.array(recall_at_cutoff), reorder=True)
        return recall, precision, auc_at_cutoff
