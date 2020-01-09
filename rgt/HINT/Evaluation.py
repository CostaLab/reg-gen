
import os
import numpy as np
import math
from sklearn import metrics
from scipy.integrate import trapz
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


def evaluation_args(parser):
    # Input Options
    parser.add_argument("--tfbs-file", type=str,
                        metavar="FILE", default=None,
                        help="A bed file containing all motif-predicted binding sites (MPBSs)."
                             "The values in the bed SCORE field will be used to rank the MPBSs."
                             "The extension must be (.bed).")
    parser.add_argument("--footprint-file", type=str,
                        metavar="FILE1,FILE2,FILE3,FILE4...", default=None,
                        help="The bed files containing the footprint predictions."
                             "The extension must be (.bed).")
    parser.add_argument("--footprint-name", type=str,
                        metavar="NAME1,NAME2,NAME3,NAME4...", default=None,
                        help="The name used to distinguish different footprint files."
                             "The number of methods name must be consistent with that of footprint files")
    parser.add_argument("--footprint-type", type=str,
                        metavar="TYPE1,TYPE2,TYPE3,TYPE4...", default=None,
                        help="The methods type used to predicted the footprint."
                             "Must be SC (site centric) or SEG (segmentation approach)")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd())
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default=None)

    parser.add_argument("--print-roc-curve", action="store_true", default=False,
                        help="If set, HINT will print the receiver operating characteristic curve.")
    parser.add_argument("--print-pr-curve", action="store_true", default=False,
                        help="If set, HINT will print the precision recall curve.")


def evaluation_run(args):
    chip_evaluate(args)


def chip_evaluate(args):
    # Evaluate Statistics
    fpr = dict()
    tpr = dict()
    roc_auc_1 = dict()
    roc_auc_10 = dict()
    roc_auc_50 = dict()
    roc_auc_100 = dict()
    recall = dict()
    precision = dict()
    prc_auc_1 = dict()
    prc_auc_10 = dict()
    prc_auc_50 = dict()
    prc_auc_100 = dict()

    footprint_file = args.footprint_file.split(",")
    footprint_name = args.footprint_name.split(",")
    footprint_type = args.footprint_type.split(",")

    max_score = 0
    if "SEG" in footprint_type:
        mpbs_regions = GenomicRegionSet("TFBS")
        mpbs_regions.read(args.tfbs_file)

        # Verifying the maximum score of the MPBS file
        for region in iter(mpbs_regions):
            score = int(region.data.split("\t")[0])
            if score > max_score:
                max_score = score
        max_score += 1

    max_points = []
    for i in range(len(footprint_file)):
        footprints_regions = GenomicRegionSet("Footprints Prediction")
        footprints_regions.read(footprint_file[i])
        footprints_regions.sort()

        if footprint_type[i] == "SEG":
            # Increasing the score of MPBS entry once if any overlaps found in the predicted footprints.
            increased_score_mpbs_regions = GenomicRegionSet("Increased Regions")
            intersect_regions = mpbs_regions.intersect(footprints_regions, mode=OverlapType.ORIGINAL)
            for region in iter(intersect_regions):
                region.data = str(int(region.data.split("\t")[0]) + max_score)
                increased_score_mpbs_regions.add(region)

            # Keep the score of remained MPBS entry unchanged
            without_intersect_regions = mpbs_regions.subtract(footprints_regions, whole_region=True)
            for region in iter(without_intersect_regions):
                increased_score_mpbs_regions.add(region)

            increased_score_mpbs_regions.sort_score()

            fpr[i], tpr[i], roc_auc_1[i], roc_auc_10[i], roc_auc_50[i], roc_auc_100[i] = \
                roc_curve(increased_score_mpbs_regions)
            recall[i], precision[i], prc_auc_1[i], prc_auc_10[i], prc_auc_50[i], prc_auc_100[i] = \
                precision_recall_curve(increased_score_mpbs_regions)

            max_points.append(len(intersect_regions))

        elif footprint_type[i] == "SC":
            footprints_regions.sort_score()
            fpr[i], tpr[i], roc_auc_1[i], roc_auc_10[i], roc_auc_50[i], roc_auc_100[i] = \
                roc_curve(footprints_regions)
            recall[i], precision[i], prc_auc_1[i], prc_auc_10[i], prc_auc_50[i], prc_auc_100[i] = \
                precision_recall_curve(footprints_regions)

            max_points.append(len(footprints_regions))

    # Output the statistics results into text
    stats_fname = os.path.join(args.output_location, "{}_stats.txt".format(args.output_prefix))
    stats_header = ["METHOD", "AUC_100", "AUC_50", "AUC_10", "AUC_1", "AUPR_100", "AUPR_50", "AUPR_10", "AUPR_1"]
    with open(stats_fname, "w") as stats_file:
        stats_file.write("\t".join(stats_header) + "\n")
        for i in range(len(footprint_name)):
            stats_file.write(footprint_name[i] + "\t" +
                             str(roc_auc_100[i]) + "\t" + str(roc_auc_50[i]) + "\t" + str(roc_auc_10[i]) + "\t" +
                             str(roc_auc_1[i]) + "\t" + str(prc_auc_100[i]) + "\t" + str(prc_auc_50[i]) + "\t" +
                             str(prc_auc_10[i]) + "\t" + str(prc_auc_1[i]) + "\n")

    # Output the curves
    if args.print_roc_curve:
        label_x = "False Positive Rate"
        label_y = "True Positive Rate"
        curve_name = "ROC"
        plot_curve(footprint_name, args.output_location, fpr, tpr, roc_auc_100, label_x, label_y, args.output_prefix,
                   curve_name, max_points=max_points)
    if args.print_pr_curve:
        label_x = "Recall"
        label_y = "Precision"
        curve_name = "PRC"
        plot_curve(footprint_name, args.output_location, recall, precision, prc_auc_100, label_x, label_y,
                   args.output_prefix, curve_name, max_points=max_points)

    output_points(footprint_name, args.output_location, args.output_prefix, fpr, tpr, recall, precision)


def plot_curve(footprint_name, output_location, data_x, data_y, stats, label_x, label_y,
               tf_name, curve_name, max_points=None):
    color_list = ["#000000", "#000099", "#006600", "#990000", "#660099", "#CC00CC", "#222222", "#CC9900",
                  "#FF6600", "#0000CC", "#336633", "#CC0000", "#6600CC", "#FF00FF", "#555555", "#CCCC00",
                  "#FF9900", "#0000FF", "#33CC33", "#FF0000", "#663399", "#FF33FF", "#888888", "#FFCC00",
                  "#663300", "#009999", "#66CC66", "#FF3333", "#9933FF", "#FF66FF", "#AAAAAA", "#FFCC33",
                  "#993300", "#00FFFF", "#99FF33", "#FF6666", "#CC99FF", "#FF99FF", "#CCCCCC", "#FFFF00"]

    matplotlib.use('Agg')
    # Creating figure
    fig = plt.figure(figsize=(8, 5), facecolor='w', edgecolor='k')
    ax = fig.add_subplot(111)
    for idx, name in enumerate(footprint_name):
        max_point = max_points[idx]
        ax.plot(data_x[idx][:max_point], data_y[idx][:max_point], color=color_list[idx],
                label=name + ": " + str(round(stats[idx], 4)))


    # Plot red diagonal line
    ax.plot([0, 1.0], [0, 1.0], color="#CCCCCC", linestyle="--", alpha=1.0)

    # Line legend
    leg = ax.legend(bbox_to_anchor=(1.0, 0.0, 1.0, 1.0), loc=2, ncol=2, borderaxespad=0., title="AUC",
                    prop={'size': 4})
    for e in leg.legendHandles:
        e.set_linewidth(2.0)

    # Titles and Axis Labels
    ax.set_title(tf_name)
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
    figure_name = os.path.join(output_location, "{}_{}.png".format(tf_name, curve_name))
    fig.savefig(figure_name, format="png", dpi=300, bbox_inches='tight', bbox_extra_artists=[leg])


def roc_curve(sort_score_regions):
    count_x = 0
    count_y = 0
    fpr = [count_x]
    tpr = [count_y]

    for region in iter(sort_score_regions):
        if str(region.name).split(":")[-1] == "Y":
            count_y += 1
            fpr.append(count_x)
            tpr.append(count_y)
        else:
            count_x += 1
            fpr.append(count_x)
            tpr.append(count_y)
    fpr = [e * (1.0 / count_x) for e in fpr]
    tpr = [e * (1.0 / count_y) for e in tpr]

    roc_auc_100 = metrics.auc(fpr, tpr)

    # Evaluating 1% AUC
    fpr_auc = 0.01
    fpr_1 = list()
    tpr_1 = list()
    for idx in range(0, len(fpr)):
        if fpr[idx] > fpr_auc:
            break
        fpr_1.append(fpr[idx])
        tpr_1.append(tpr[idx])
    if len(fpr_1) < 2:
        roc_auc_1 = 0  # At least 2 points are needed to compute area under curve
    else:
        roc_auc_1 = metrics.auc(standardize(fpr_1), tpr_1)

    # Evaluating 10% AUC
    fpr_auc = 0.1
    fpr_10 = list()
    tpr_10 = list()
    for idx in range(0, len(fpr)):
        if fpr[idx] > fpr_auc:
            break
        fpr_10.append(fpr[idx])
        tpr_10.append(tpr[idx])
    if len(fpr_10) < 2:
        roc_auc_10 = 0  # At least 2 points are needed to compute area under curve
    else:
        roc_auc_10 = metrics.auc(standardize(fpr_10), tpr_10)

    # Evaluating 10% AUC
    fpr_auc = 0.5
    fpr_50 = list()
    tpr_50 = list()
    for idx in range(0, len(fpr)):
        if fpr[idx] > fpr_auc:
            break
        fpr_50.append(fpr[idx])
        tpr_50.append(tpr[idx])
    if len(fpr_50) < 2:
        roc_auc_50 = 0  # At least 2 points are needed to compute area under curve
    else:
        roc_auc_50 = metrics.auc(standardize(fpr_50), tpr_50)

    return fpr, tpr, roc_auc_1, roc_auc_10, roc_auc_50, roc_auc_100


def precision_recall_curve(sort_score_regions):
    count_x = 0
    count_y = 0
    precision = [0.0]
    recall = [0.0]
    for region in iter(sort_score_regions):
        if str(region.name).split(":")[-1] == "Y":
            count_y += 1
            precision.append(float(count_y) / (count_x + count_y))
            recall.append(count_y)
        else:
            count_x += 1
            precision.append(float(count_y) / (count_x + count_y))
            recall.append(count_y)

    recall = [e * (1.0 / count_y) for e in recall]

    precision.append(0.0)
    recall.append(1.0)

    # Evaluating 100% AUPR
    pr_auc_100 = (abs(trapz(recall, precision)))

    # Evaluating 1% AUPR
    recall_auc = 0.01
    recall_1 = list()
    precision_1 = list()
    for idx in range(0, len(recall)):
        if recall[idx] > recall_auc:
            break
        recall_1.append(recall[idx])
        precision_1.append(precision[idx])
    if len(recall_1) < 2:
        pr_auc_1 = 0  # At least 2 points are needed to compute area under curve
    else:
        pr_auc_1 = metrics.auc(standardize(recall_1), precision_1)

    # Evaluating 10% AUPR
    recall_auc = 0.1
    recall_10 = list()
    precision_10 = list()
    for idx in range(0, len(recall)):
        if recall[idx] > recall_auc:
            break
        recall_10.append(recall[idx])
        precision_10.append(precision[idx])
    if len(recall_10) < 2:
        pr_auc_10 = 0  # At least 2 points are needed to compute area under curve
    else:
        pr_auc_10 = metrics.auc(standardize(recall_10), precision_10)

    # Evaluating 50% AUPR
    recall_auc = 0.5
    recall_50 = list()
    precision_50 = list()
    for idx in range(0, len(recall)):
        if recall[idx] > recall_auc:
            break
        recall_50.append(recall[idx])
        precision_50.append(precision[idx])
    if len(recall_50) < 2:
        pr_auc_50 = 0  # At least 2 points are needed to compute area under curve
    else:
        pr_auc_50 = metrics.auc(standardize(recall_50), precision_50)

    return recall, precision, pr_auc_1, pr_auc_10, pr_auc_50, pr_auc_100


def standardize(vector):
    max_num = max(vector)
    min_num = min(vector)
    if max_num == min_num:
        return vector
    else:
        return [(e - min_num) / (max_num - min_num) for e in vector]


def optimize_roc_points(footprint_name, fpr, tpr, max_points=1000):
    new_fpr = dict()
    new_tpr = dict()
    for i in range(len(footprint_name)):
        if len(fpr[i]) > max_points:
            new_idx_list = [int(math.ceil(e)) for e in np.linspace(0, len(fpr[i]) - 1, max_points)]
            new_idx_fpr = []
            new_idx_tpr = []
            for j in new_idx_list:
                new_idx_fpr.append(fpr[i][j])
                new_idx_tpr.append(tpr[i][j])
            new_fpr[i] = new_idx_fpr
            new_tpr[i] = new_idx_tpr
        else:
            new_fpr[i] = fpr[i]
            new_tpr[i] = tpr[i]

    return new_fpr, new_tpr


def optimize_pr_points(footprint_name, recall, precision, max_points=1000):
    new_recall = dict()
    new_precision = dict()

    for i in range(len(footprint_name)):
        data_recall = []
        data_precision = []
        for j in range(min(max_points, len(recall[i]))):
            data_recall.append(recall[i].pop(0))
            data_precision.append(precision[i].pop(0))
        if len(recall[i]) > max_points:
            new_idx_list = [int(math.ceil(e)) for e in np.linspace(0, len(recall[i]) - 1, max_points)]
            for j in new_idx_list:
                data_recall.append(recall[i][j])
                data_precision.append(precision[i][j])
            new_recall[i] = data_recall
            new_precision[i] = data_precision
        else:
            new_recall[i] = data_recall + recall[i]
            new_precision[i] = data_precision + precision[i]

    return new_recall, new_precision


def output_points(footprint_name, output_location, output_prefix, fpr, tpr, recall, precision):
    roc_fname = os.path.join(output_location, "{}_roc.txt".format(output_prefix))
    #new_fpr, new_tpr = optimize_roc_points(footprint_name, fpr, tpr)
    new_fpr, new_tpr = fpr, tpr
    header = list()
    len_vec = list()
    for i in range(len(footprint_name)):
        header.append(footprint_name[i] + "_FPR")
        header.append(footprint_name[i] + "_TPR")
        len_vec.append(len(new_fpr[i]))

    with open(roc_fname, "w") as roc_file:
        roc_file.write("\t".join(header) + "\n")
        max_idx = len_vec.index(max(len_vec))
        for j in range(len(new_fpr[max_idx])):
            to_write = list()
            for i in range(len(footprint_name)):
                if j >= len(new_fpr[i]):
                    to_write.append("NA")
                else:
                    to_write.append(str(new_fpr[i][j]))
                if j >= len(new_tpr[i]):
                    to_write.append("NA")
                else:
                    to_write.append(str(new_tpr[i][j]))
            roc_file.write("\t".join(to_write) + "\n")

    prc_fname = os.path.join(output_location, "{}_prc.txt".format(output_prefix))
    #new_recall, new_precision = optimize_pr_points(footprint_name, recall, precision)
    new_recall, new_precision = recall, precision
    header = list()
    len_vec = list()
    for i in range(len(footprint_name)):
        header.append(footprint_name[i] + "_REC")
        header.append(footprint_name[i] + "_PRE")
        len_vec.append(len(new_recall[i]))

    with open(prc_fname, "w") as prc_file:
        prc_file.write("\t".join(header) + "\n")
        max_idx = len_vec.index(max(len_vec))
        for j in range(len(new_recall[max_idx])):
            to_write = list()
            for i in range(len(footprint_name)):
                if j >= len(new_recall[i]):
                    to_write.append("NA")
                else:
                    to_write.append(str(new_recall[i][j]))
                if j >= len(new_precision[i]):
                    to_write.append("NA")
                else:
                    to_write.append(str(new_precision[i][j]))
            prc_file.write("\t".join(to_write) + "\n")
