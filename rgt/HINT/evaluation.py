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
from pysam import Samfile

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

    def chip_evaluate(self, bam_fname):
        """
        This evaluation methodology uses motif-predicted binding sites (MPBSs) together with TF ChIP-seq data
        to evaluate the footprint predictions.

        return:
        """
        bam = Samfile(bam_fname, "rb")

        fpr_auc_threshold_1 = 0.1
        fpr_auc_threshold_2 = 0.01

        mpbs_regions = GenomicRegionSet("MPBS")
        mpbs_regions.read_bed(self.mpbs_fname)

        mpbs_name = self.mpbs_fname.split("/")[-1].split(".")[-2]

        # Evaluate Statistics
        stats_header = ["METHOD", "AUC_100", "AUC_10", "AUC_1", "AUPR"]
        fpr = dict()
        tpr = dict()
        roc_auc = dict()
        roc_auc_1 = dict()
        roc_auc_2 = dict()
        precision = dict()
        recall = dict()
        average_precision = dict()
        for i in range(len(self.footprint_fname)):
            footprints_regions = GenomicRegionSet("Footprints Prediction")
            footprints_regions.read_bed(self.footprint_fname[i])

            # Sort footprint prediction bed files
            footprints_regions.sort()

            ################################################
            ## Extend 10 bp for HINT, HINTBC, HINTBCN
            extend_list = ["HINT", "HINTBC", "HINTBCN"]
            if self.footprint_fname[i] in extend_list:
                for region in iter(footprints_regions):
                    if len(region) < 10:
                        region.initial = max(0, region.initial - 10)
                        region.final = region.final + 10
            #################################################

            # if MPBS entry overlaps with the predicted footprints, their score will be TC of the footprints
            score_mpbs_regions = GenomicRegionSet("Increased Regions")
            intersect_mpbs_regions = mpbs_regions.intersect(footprints_regions, mode=OverlapType.ORIGINAL)
            intersect_footprints_regions = footprints_regions.intersect(mpbs_regions, mode=OverlapType.ORIGINAL)
            for region in iter(intersect_footprints_regions):
                if not region.data:
                    region.data = str(region.name)
            intersect_footprints_regions.sort_score()

            for mpbs_region in iter(intersect_mpbs_regions):
                for footprints_region in iter(intersect_footprints_regions):
                    if mpbs_region.overlap(footprints_region):
                        mpbs_region.data = footprints_region.data
                        score_mpbs_regions.add(mpbs_region)
                        break

            # if without overlap, score equals TC of themselves
            # without_intersect_regions = mpbs_regions.subtract(footprints_regions, whole_region=True)
            # for region in iter(without_intersect_regions):
            #    region.data = str(bam.count(reference=region.chrom, start=region.initial, end=region.final))
            #    score_mpbs_regions.add(region)

            score_mpbs_regions.sort_score()

            # Reading data points
            y_true = list()
            y_score = list()
            for region in iter(score_mpbs_regions):
                if str(region.name).split(":")[-1] == "N":
                    y_true.append(0)
                else:
                    y_true.append(1)
                y_score.append(int(region.data))

            # Calculating receiver operating characteristic curve (ROC),
            # AUC at 100% FPR, AUC at 10% FPR, AUC at 1% FPR
            fpr[i], tpr[i], _ = metrics.roc_curve(np.array(y_true), np.array(y_score))
            fpr_1 = list()
            tpr_1 = list()
            fpr_2 = list()
            tpr_2 = list()
            for index in range(len(fpr[i])):
                if fpr[i][index] <= fpr_auc_threshold_1:
                    fpr_1.append(fpr[i][index])
                    tpr_1.append(tpr[i][index])
                if fpr[i][index] <= fpr_auc_threshold_2:
                    fpr_2.append(fpr[i][index])
                    tpr_2.append(tpr[i][index])

            fpr_1.append(fpr_auc_threshold_1)
            tpr_1.append(tpr_1[-1])
            fpr_2.append(fpr_auc_threshold_2)
            tpr_2.append(tpr_2[-1])
            roc_auc[i] = metrics.auc(fpr[i], tpr[i])
            roc_auc_1[i] = metrics.auc(np.array(fpr_1), np.array(tpr_1)) * 10
            roc_auc_2[i] = metrics.auc(np.array(fpr_2), np.array(tpr_2)) * 100

            # Calculating precision-recall curve (PRC) and the area under the precision-recall curve
            precision[i], recall[i], _ = metrics.precision_recall_curve(np.array(y_true), np.array(y_score))
            average_precision[i] = metrics.average_precision_score(np.array(y_true), np.array(y_score))

        # Output the statistics results into text
        stats_fname = self.output_location + mpbs_name + "_stats.txt"
        with open(stats_fname, "w") as stats_file:
            stats_file.write("\t".join(stats_header) + "\n")
            for i in range(len(self.footprint_name)):
                stats_file.write(self.footprint_name[i] + "\t" + str(roc_auc[i]) + "\t" + str(roc_auc_1[i]) + "\t"
                                 + str(roc_auc_2[i]) + "\t" + str(average_precision[i]) + "\n")

        # Output the curves
        if self.print_roc_curve:
            label_x = "False Positive Rate"
            label_y = "True Positive Rate"
            curve_type = "ROC"
            self.plot_curve(fpr, tpr, roc_auc, label_x, label_y, mpbs_name, curve_type)
        if self.print_pr_curve:
            label_x = "Recall"
            label_y = "Precision"
            curve_type = "PRC"
            self.plot_curve(recall, precision, average_precision,
                            label_x, label_y, mpbs_name, curve_type)

    def plot_curve(self, data_x, data_y, stats, label_x, label_y, mpbs_name, curve_type):
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
        figure_name = self.output_location + mpbs_name + "_" + curve_type + ".png"
        fig.savefig(figure_name, format="png", dpi=300, bbox_inches='tight', bbox_extra_artists=[leg])
