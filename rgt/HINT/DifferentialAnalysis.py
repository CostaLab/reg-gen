import os
import numpy as np
from pysam import Samfile, Fastafile
from math import ceil, floor
from Bio import motifs
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pyx
from scipy.stats.mvn import mvnun
from argparse import SUPPRESS

from multiprocessing import Pool, cpu_count

# Internal
from rgt.Util import ErrorHandler, AuxiliaryFunctions, GenomeData, HmmData
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.HINT.biasTable import BiasTable

"""
Perform differential footprints analysis based on the prediction of transcription factor binding sites.

Authors: Eduardo G. Gusmao, Zhijian Li
"""

dic = {"A": 0, "C": 1, "G": 2, "T": 3}


def diff_analysis_args(parser):
    # Input Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19, hg38. mm9, and mm10. DEFAULT: hg19")
    parser.add_argument("--mpbs-file1", type=str, metavar="FILE", default=None,
                        help="motif predicted binding sites file for condition 1, must be .bed file. DEFAULT: None")
    parser.add_argument("--mpbs-file2", type=str, metavar="FILE", default=None,
                        help="motif predicted binding sites file for condition 2, must be .bed file. DEFAULT: None")
    parser.add_argument("--reads-file1", type=str, metavar="FILE", default=None,
                        help="The BAM file containing the DNase-seq or ATAC-seq reads for condition 1. DEFAULT: None")
    parser.add_argument("--reads-file2", type=str, metavar="FILE", default=None,
                        help="The BAM file containing the DNase-seq or ATAC-seq reads for condition 2. DEFAULT: None")

    parser.add_argument("--window-size", type=int, metavar="INT", default=200,
                        help="The window size for differential analysis. DEFAULT: 200")
    parser.add_argument("--factor1", type=float, metavar="FLOAT", default=None,
                        help="The normalization factor for condition 1. DEFAULT: None")
    parser.add_argument("--factor2", type=float, metavar="FLOAT", default=None,
                        help="The normalization factor for condition 1. DEFAULT: None")

    parser.add_argument("--forward-shift", type=int, metavar="INT", default=5, help=SUPPRESS)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=-4, help=SUPPRESS)
    parser.add_argument("--bias-table1", type=str, metavar="FILE1_F,FILE1_R", default=None, help=SUPPRESS)
    parser.add_argument("--bias-table2", type=str, metavar="FILE2_F,FILE2_R", default=None, help=SUPPRESS)

    parser.add_argument("--condition1", type=str, metavar="STRING", default="condition1",
                        help="The name of condition1. DEFAULT: condition1")
    parser.add_argument("--condition2", type=str, metavar="STRING", default="condition1",
                        help="The name of condition2. DEFAULT: condition2")
    parser.add_argument("--fdr", type=float, metavar="FLOAT", default=0.05,
                        help="The false discovery rate. DEFAULT: 0.05")
    parser.add_argument("--bc", action="store_true", default=False,
                        help="If set, all analysis will be based on bias corrected signal. DEFAULT: False")
    parser.add_argument("--nc", type=int, metavar="INT", default=cpu_count(),
                        help="The number of cores. DEFAULT: 1")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="differential",
                        help="The prefix for results files. DEFAULT: differential")
    parser.add_argument("--standardize", action="store_true", default=False,
                        help="If set, the signal will be rescaled to (0, 1) for plotting.")
    parser.add_argument("--output-profiles", default=False, action='store_true',
                        help="If set, the footprint profiles will be writen into a text, in which each row is a "
                             "specific instance of the given motif. DEFAULT: False")


def get_raw_signal(arguments):
    (mpbs_name, mpbs_file1, mpbs_file2, reads_file1, reads_file2, organism,
     window_size, forward_shift, reverse_shift) = arguments

    mpbs1 = GenomicRegionSet("Motif Predicted Binding Sites of Condition1")
    mpbs1.read(mpbs_file1)

    mpbs2 = GenomicRegionSet("Motif Predicted Binding Sites of Condition2")
    mpbs2.read(mpbs_file2)

    mpbs = mpbs1.combine(mpbs2, output=True)
    mpbs.sort()

    bam1 = Samfile(reads_file1, "rb")
    bam2 = Samfile(reads_file2, "rb")

    genome_data = GenomeData(organism)
    fasta = Fastafile(genome_data.get_genome())

    signal_1 = np.zeros(window_size)
    signal_2 = np.zeros(window_size)
    motif_len = None
    pwm = dict([("A", [0.0] * window_size), ("C", [0.0] * window_size),
                ("G", [0.0] * window_size), ("T", [0.0] * window_size),
                ("N", [0.0] * window_size)])

    mpbs_regions = mpbs.by_names([mpbs_name])
    num_motif = len(mpbs_regions)

    for region in mpbs_regions:
        if motif_len is None:
            motif_len = region.final - region.initial

        mid = (region.final + region.initial) / 2
        p1 = mid - window_size / 2
        p2 = mid + window_size / 2

        if p1 <= 0:
            continue

        # Fetch raw signal
        for read in bam1.fetch(region.chrom, p1, p2):
            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1 <= cut_site < p2:
                    signal_1[cut_site - p1] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1 <= cut_site < p2:
                    signal_1[cut_site - p1] += 1.0

        for read in bam2.fetch(region.chrom, p1, p2):
            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1 <= cut_site < p2:
                    signal_2[cut_site - p1] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1 <= cut_site < p2:
                    signal_2[cut_site - p1] += 1.0
        update_pwm(pwm, fasta, region, p1, p2)

    return signal_1, signal_2, motif_len, pwm, num_motif


def get_bc_signal(arguments):
    (mpbs_name, mpbs_file1, mpbs_file2, reads_file1, reads_file2, organism,
     window_size, forward_shift, reverse_shift, bias_table1, bias_table2) = arguments

    mpbs1 = GenomicRegionSet("Motif Predicted Binding Sites of Condition1")
    mpbs1.read(mpbs_file1)

    mpbs2 = GenomicRegionSet("Motif Predicted Binding Sites of Condition2")
    mpbs2.read(mpbs_file2)

    mpbs = mpbs1.combine(mpbs2, output=True)
    mpbs.sort()

    bam1 = Samfile(reads_file1, "rb")
    bam2 = Samfile(reads_file2, "rb")

    genome_data = GenomeData(organism)
    fasta = Fastafile(genome_data.get_genome())

    signal_1 = np.zeros(window_size)
    signal_2 = np.zeros(window_size)
    motif_len = None
    pwm = dict([("A", [0.0] * window_size), ("C", [0.0] * window_size),
                ("G", [0.0] * window_size), ("T", [0.0] * window_size),
                ("N", [0.0] * window_size)])

    mpbs_regions = mpbs.by_names([mpbs_name])
    num_motif = len(mpbs_regions)

    # Fetch bias corrected signal
    for region in mpbs_regions:
        if motif_len is None:
            motif_len = region.final - region.initial

        mid = (region.final + region.initial) / 2
        p1 = mid - window_size / 2
        p2 = mid + window_size / 2

        if p1 <= 0:
            continue
        # Fetch raw signal
        signal1 = bias_correction(chrom=region.chrom, start=p1, end=p2, bam=bam1,
                                  bias_table=bias_table1, genome_file_name=genome_data.get_genome(),
                                  forward_shift=forward_shift, reverse_shift=reverse_shift)

        signal2 = bias_correction(chrom=region.chrom, start=p1, end=p2, bam=bam2,
                                  bias_table=bias_table2, genome_file_name=genome_data.get_genome(),
                                  forward_shift=forward_shift, reverse_shift=reverse_shift)

        if len(signal1) != len(signal_1) or len(signal2) != len(signal_2):
            continue

        signal_1 = np.add(signal_1, np.array(signal1))
        signal_2 = np.add(signal_2, np.array(signal2))

        update_pwm(pwm, fasta, region, p1, p2)

    return signal_1, signal_2, motif_len, pwm, num_motif


def diff_analysis_run(args):
    # Initializing Error Handler
    err = ErrorHandler()

    output_location = os.path.join(args.output_location, "Lineplots")
    try:
        if not os.path.isdir(output_location):
            os.makedirs(output_location)
    except Exception:
        err.throw_error("MM_OUT_FOLDER_CREATION")

    mpbs1 = GenomicRegionSet("Motif Predicted Binding Sites of Condition1")
    mpbs1.read(args.mpbs_file1)

    mpbs2 = GenomicRegionSet("Motif Predicted Binding Sites of Condition2")
    mpbs2.read(args.mpbs_file2)

    mpbs = mpbs1.combine(mpbs2, output=True)
    mpbs.sort()
    mpbs_name_list = list(set(mpbs.get_names()))

    signal_dict_by_tf_1 = dict()
    signal_dict_by_tf_2 = dict()
    motif_len_dict = dict()
    motif_num_dict = dict()
    pwm_dict_by_tf = dict()

    pool = Pool(processes=args.nc)
    # differential analysis using bias corrected signal
    if args.bc:
        hmm_data = HmmData()
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table1 = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)
        bias_table2 = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)

        mpbs_list = list()
        for mpbs_name in mpbs_name_list:
            mpbs_list.append((mpbs_name, args.mpbs_file1, args.mpbs_file2, args.reads_file1, args.reads_file2,
                              args.organism, args.window_size, args.forward_shift, args.reverse_shift,
                              bias_table1, bias_table2))
        res = pool.map(get_bc_signal, mpbs_list)

    # differential analysis using raw signal
    else:
        mpbs_list = list()
        for mpbs_name in mpbs_name_list:
            mpbs_list.append((mpbs_name, args.mpbs_file1, args.mpbs_file2, args.reads_file1, args.reads_file2,
                              args.organism, args.window_size, args.forward_shift, args.reverse_shift))
        res = pool.map(get_raw_signal, mpbs_list)

    for idx, mpbs_name in enumerate(mpbs_name_list):
        signal_dict_by_tf_1[mpbs_name] = res[idx][0]
        signal_dict_by_tf_2[mpbs_name] = res[idx][1]
        motif_len_dict[mpbs_name] = res[idx][2]
        pwm_dict_by_tf[mpbs_name] = res[idx][3]
        motif_num_dict[mpbs_name] = res[idx][4]

    if args.factor1 is None or args.factor2 is None:
        args.factor1, args.factor2 = compute_factors(signal_dict_by_tf_1, signal_dict_by_tf_2)
        output_factor(args, args.factor1, args.factor2)

    if args.output_profiles:
        output_profiles(mpbs_name_list, signal_dict_by_tf_1, output_location, args.condition1)
        output_profiles(mpbs_name_list, signal_dict_by_tf_2, output_location, args.condition2)

    ps_tc_results_by_tf = dict()

    plots_list = list()
    for mpbs_name in mpbs_name_list:
        plots_list.append((mpbs_name, motif_num_dict[mpbs_name], signal_dict_by_tf_1[mpbs_name],
                           signal_dict_by_tf_2[mpbs_name], args.factor1, args.factor2, args.condition1,
                           args.condition2, pwm_dict_by_tf[mpbs_name], output_location, args.window_size,
                           args.standardize))

    pool.map(line_plot, plots_list)

    for mpbs_name in mpbs_name_list:
        res = get_ps_tc_results(signal_dict_by_tf_1[mpbs_name], signal_dict_by_tf_2[mpbs_name],
                                args.factor1, args.factor2, motif_num_dict[mpbs_name], motif_len_dict[mpbs_name])
        #
        #     # only use the factors whose protection scores are greater than 0
        #     if res[0] > 0 and res[1] < 0:
        ps_tc_results_by_tf[mpbs_name] = res
    #
    stat_results_by_tf = get_stat_results(ps_tc_results_by_tf)
    scatter_plot(args, stat_results_by_tf)
    output_stat_results(args, stat_results_by_tf)


def bias_correction(chrom, start, end, bam, bias_table, genome_file_name, forward_shift, reverse_shift):
    # Parameters
    window = 50
    defaultKmerValue = 1.0

    # Initialization
    fastaFile = Fastafile(genome_file_name)
    fBiasDict = bias_table[0]
    rBiasDict = bias_table[1]
    k_nb = len(fBiasDict.keys()[0])
    p1 = start
    p2 = end
    p1_w = p1 - (window / 2)
    p2_w = p2 + (window / 2)
    p1_wk = p1_w - int(floor(k_nb / 2.))
    p2_wk = p2_w + int(ceil(k_nb / 2.))
    if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
        # Return raw counts
        bc_signal = [0.0] * (p2 - p1)
        for read in bam.fetch(chrom, p1, p2):
            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1 <= cut_site < p2:
                    bc_signal[cut_site - p1] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1 <= cut_site < p2:
                    bc_signal[cut_site - p1] += 1.0

        return bc_signal

    # Raw counts
    nf = [0.0] * (p2_w - p1_w)
    nr = [0.0] * (p2_w - p1_w)
    for read in bam.fetch(chrom, p1_w, p2_w):
        if not read.is_reverse:
            cut_site = read.pos + forward_shift
            if p1_w <= cut_site < p2_w:
                nf[cut_site - p1_w] += 1.0
        else:
            cut_site = read.aend + reverse_shift - 1
            if p1_w <= cut_site < p2_w:
                nr[cut_site - p1_w] += 1.0

    # Smoothed counts
    Nf = []
    Nr = []
    f_sum = sum(nf[:window])
    r_sum = sum(nr[:window])
    f_last = nf[0]
    r_last = nr[0]
    for i in range((window / 2), len(nf) - (window / 2)):
        Nf.append(f_sum)
        Nr.append(r_sum)
        f_sum -= f_last
        f_sum += nf[i + (window / 2)]
        f_last = nf[i - (window / 2) + 1]
        r_sum -= r_last
        r_sum += nr[i + (window / 2)]
        r_last = nr[i - (window / 2) + 1]

    # Fetching sequence
    currStr = str(fastaFile.fetch(chrom, p1_wk, p2_wk - 1)).upper()
    currRevComp = AuxiliaryFunctions.revcomp(str(fastaFile.fetch(chrom, p1_wk + 1, p2_wk)).upper())

    # Iterating on sequence to create signal
    af = []
    ar = []
    for i in range(int(ceil(k_nb / 2.)), len(currStr) - int(floor(k_nb / 2)) + 1):
        fseq = currStr[i - int(floor(k_nb / 2.)):i + int(ceil(k_nb / 2.))]
        rseq = currRevComp[len(currStr) - int(ceil(k_nb / 2.)) - i:len(currStr) + int(floor(k_nb / 2.)) - i]
        try:
            af.append(fBiasDict[fseq])
        except Exception:
            af.append(defaultKmerValue)
        try:
            ar.append(rBiasDict[rseq])
        except Exception:
            ar.append(defaultKmerValue)

    # Calculating bias and writing to wig file
    f_sum = sum(af[:window])
    r_sum = sum(ar[:window])
    f_last = af[0]
    r_last = ar[0]
    bc_signal = []
    for i in range((window / 2), len(af) - (window / 2)):
        nhatf = Nf[i - (window / 2)] * (af[i] / f_sum)
        nhatr = Nr[i - (window / 2)] * (ar[i] / r_sum)
        bc_signal.append(nhatf + nhatr)
        f_sum -= f_last
        f_sum += af[i + (window / 2)]
        f_last = af[i - (window / 2) + 1]
        r_sum -= r_last
        r_sum += ar[i + (window / 2)]
        r_last = ar[i - (window / 2) + 1]

    # Termination
    fastaFile.close()
    return bc_signal


def get_ps_tc_results(signal_1, signal_2, factor1, factor2, num_motif, motif_len):
    signal_1 = (signal_1 / factor1) / num_motif
    signal_2 = (signal_2 / factor2) / num_motif

    # signal_1, signal_2 = standard(signal_1, signal_2)

    signal_half_len = len(signal_1) / 2

    nc = sum(signal_1[signal_half_len - motif_len / 2:signal_half_len + motif_len / 2])
    nr = sum(signal_1[signal_half_len + motif_len / 2:signal_half_len + motif_len / 2 + motif_len])
    nl = sum(signal_1[signal_half_len - motif_len / 2 - motif_len:signal_half_len - motif_len / 2])

    protect_score1 = (nr - nc) / motif_len + (nl - nc) / motif_len
    tc1 = (sum(signal_1) - nc) / (len(signal_1) - motif_len)

    nc = sum(signal_2[signal_half_len - motif_len / 2:signal_half_len + motif_len / 2])
    nr = sum(signal_2[signal_half_len + motif_len / 2:signal_half_len + motif_len / 2 + motif_len])
    nl = sum(signal_2[signal_half_len - motif_len / 2 - motif_len:signal_half_len - motif_len / 2])

    protect_score2 = (nr - nc) / motif_len + (nl - nc) / motif_len
    tc2 = (sum(signal_2) - nc) / (len(signal_2) - motif_len)

    protect_diff = protect_score2 - protect_score1
    tc_diff = tc2 - tc1
    return [protect_score1, protect_score2, protect_diff, tc1, tc2, tc_diff]


def update_pwm(pwm, fasta, region, p1, p2):
    # Update pwm
    aux_plus = 1
    dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
    if (region.final - region.initial) % 2 == 0:
        aux_plus = 0
    dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                             p1 + aux_plus, p2 + aux_plus)).upper())
    if region.orientation == "+":
        for i in range(0, len(dna_seq)):
            pwm[dna_seq[i]][i] += 1
    elif region.orientation == "-":
        for i in range(0, len(dna_seq_rev)):
            pwm[dna_seq_rev[i]][i] += 1


def compute_factors(signal_dict_by_tf_1, signal_dict_by_tf_2):
    keys = signal_dict_by_tf_1.keys()

    signal_1 = np.zeros(len(keys))
    signal_2 = np.zeros(len(keys))

    for idx, key in enumerate(keys):
        signal_1[idx] = sum(signal_dict_by_tf_1[key])
        signal_2[idx] = sum(signal_dict_by_tf_2[key])

    # Take log
    log_tc1 = np.log(signal_1)
    log_tc2 = np.log(signal_2)

    # Average
    average_log_tc = np.add(log_tc1, log_tc2) / 2

    # Filter
    filter_log_tc1 = log_tc1[~np.isnan(log_tc1)]
    filter_log_tc2 = log_tc2[~np.isnan(log_tc2)]
    filter_log_tc = average_log_tc[~np.isnan(average_log_tc)]

    # Subtract
    sub_tc1 = np.subtract(filter_log_tc1, filter_log_tc)
    sub_tc2 = np.subtract(filter_log_tc2, filter_log_tc)

    median_tc1 = np.median(sub_tc1)
    median_tc2 = np.median(sub_tc2)

    factor1 = np.exp(median_tc1)
    factor2 = np.exp(median_tc2)

    return factor1, factor2


def line_plot(arguments):
    (mpbs_name, num_fp, signal_1, signal_2, factor1, factor2, condition1, condition2,
     pwm_dict, output_location, window_size, standardize) = arguments

    mpbs_name = mpbs_name.replace("(", "_")
    mpbs_name = mpbs_name.replace(")", "")
    mean_signal_1 = (signal_1 / num_fp) / factor1
    mean_signal_2 = (signal_2 / num_fp) / factor2

    if standardize:
        mean_signal_1, mean_signal_2 = standard(mean_signal_1, mean_signal_2)

    # Output PWM and create logo
    pwm_fname = os.path.join(output_location, "{}.pwm".format(mpbs_name))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(output_location, "{}.logo.eps".format(mpbs_name))
    pwm = motifs.read(open(pwm_fname), "pfm")
    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    start = -(window_size / 2)
    end = (window_size / 2) - 1
    x = np.linspace(start, end, num=window_size)

    plt.close('all')
    fig, ax = plt.subplots()
    ax.plot(x, mean_signal_1, color='red', label=condition1)
    ax.plot(x, mean_signal_2, color='blue', label=condition2)
    ax.text(0.15, 0.9, 'n = {}'.format(num_fp), verticalalignment='bottom', horizontalalignment='right',
            transform=ax.transAxes, fontweight='bold')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    min_signal = min(min(mean_signal_1), min(mean_signal_2))
    max_signal = max(max(mean_signal_1), max(mean_signal_2))
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)

    ax.set_title(mpbs_name, fontweight='bold')
    ax.set_xlim(start, end)
    ax.set_ylim([min_signal, max_signal])
    ax.legend(loc="upper right", frameon=False)
    ax.spines['bottom'].set_position(('outward', 70))

    figure_name = os.path.join(output_location, "{}.line.eps".format(mpbs_name))
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(output_location, "{}.eps".format(mpbs_name))

    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(0.45, 0.8, logo_fname, width=16.5, height=3))
    c.writeEPSfile(output_fname)
    os.system(" ".join(["epstopdf", output_fname]))

    os.remove(figure_name)
    os.remove(logo_fname)
    os.remove(output_fname)
    os.remove(pwm_fname)


def scatter_plot(args, stat_results_by_tf):
    tc_diff = list()
    ps_diff = list()
    mpbs_name_list = stat_results_by_tf.keys()
    P_values = list()
    for mpbs_name in mpbs_name_list:
        ps_diff.append(float(stat_results_by_tf[mpbs_name][2]))
        tc_diff.append(float(stat_results_by_tf[mpbs_name][-3]))
        P_values.append(np.log10(float(stat_results_by_tf[mpbs_name][-1])))

    fig, ax = plt.subplots(figsize=(12, 12))
    for i, mpbs_name in enumerate(mpbs_name_list):
        if stat_results_by_tf[mpbs_name][-1] < args.fdr:
            ax.scatter(tc_diff[i], ps_diff[i], c="red")
            ax.annotate(mpbs_name, (tc_diff[i], ps_diff[i]), alpha=0.6)
        else:
            ax.scatter(tc_diff[i], ps_diff[i], c="black", alpha=0.6)
    ax.margins(0.05)

    tc_diff_mean = np.mean(tc_diff)
    ps_diff_mean = np.mean(ps_diff)
    ax.axvline(x=tc_diff_mean, linewidth=2, linestyle='dashed')
    ax.axhline(y=ps_diff_mean, linewidth=2, linestyle='dashed')

    ax.set_xlabel("{} $\longrightarrow$ {} \n $\Delta$ Open Chromatin Score".format(args.condition1, args.condition2),
                  fontweight='bold', fontsize=20)
    ax.set_ylabel("$\Delta$ Protection Score \n {} $\longrightarrow$ {}".format(args.condition1, args.condition2),
                  fontweight='bold', rotation=90, fontsize=20)

    figure_name = os.path.join(args.output_location, "{}_{}_statistics.pdf".format(args.condition1, args.condition2))
    fig.savefig(figure_name, format="pdf", dpi=300)


def output_results(args, ps_tc_results_by_tf):
    mpbs_name_list = ps_tc_results_by_tf.keys()
    header = ["Motif",
              "Protection_Score_{}".format(args.condition1), "Protection_Score_{}".format(args.condition2),
              "Protection_Diff_{}_{}".format(args.condition1, args.condition2),
              "TC_{}".format(args.condition1), "TC_{}".format(args.condition2),
              "TC_Diff_{}_{}".format(args.condition1, args.condition2)]
    output_fname = os.path.join(args.output_location, "{}_{}_results.txt".format(args.condition1, args.condition2))
    with open(output_fname, "w") as f:
        f.write("\t".join(header) + "\n")
        for mpbs_name in mpbs_name_list:
            f.write(mpbs_name + "\t" + "\t".join(map(str, ps_tc_results_by_tf[mpbs_name])) + "\n")


def output_stat_results(args, stat_results_by_tf):
    output_fname = os.path.join(args.output_location, "{}_{}_statistics.txt".format(args.condition1, args.condition2))
    header = ["Motif",
              "Protection_Score_{}".format(args.condition1), "Protection_Score_{}".format(args.condition2),
              "Protection_Diff_{}_{}".format(args.condition1, args.condition2),
              "TC_{}".format(args.condition1), "TC_{}".format(args.condition2),
              "TC_Diff_{}_{}".format(args.condition1, args.condition2), "P_values", "Adjust_p_values"]
    with open(output_fname, "w") as f:
        f.write("\t".join(header) + "\n")
        for mpbs_name in stat_results_by_tf.keys():
            f.write(mpbs_name + "\t" + "\t".join(map(str, stat_results_by_tf[mpbs_name])) + "\n")


def output_factor(args, factor1, factor2):
    output_file = os.path.join(args.output_location, "{}_{}_factor.txt".format(args.condition1, args.condition2))
    f = open(output_file, "w")
    f.write("Factor1: " + str(factor1) + "\n")
    f.write("Factor2: " + str(factor2) + "\n")
    f.close()


def output_mu(args, median_diff_prot, median_diff_tc):
    output_file = os.path.join(args.output_location, "{}_{}_mu.txt".format(args.condition1, args.condition2))
    f = open(output_file, "w")
    f.write("median_diff_prot: " + str(median_diff_prot) + "\n")
    f.write("median_diff_tc: " + str(median_diff_tc) + "\n")
    f.close()


def get_stat_results(ps_tc_results_by_tf):
    ps_diff = list()
    tc_diff = list()
    mpbs_name_list = ps_tc_results_by_tf.keys()
    for mpbs_name in mpbs_name_list:
        ps_diff.append(ps_tc_results_by_tf[mpbs_name][2])
        tc_diff.append(ps_tc_results_by_tf[mpbs_name][-1])

    ps_tc_diff = np.array([ps_diff, tc_diff]).T
    mu = np.mean(ps_tc_diff, axis=0)
    cov_ps_tc_diff = np.cov(ps_tc_diff.T)

    low = np.zeros(2)
    upp = np.zeros(2)
    p_values = list()
    for idx, mpbs_name in enumerate(mpbs_name_list):
        if ps_diff[idx] >= mu[0]:
            low[0] = ps_diff[idx]
            upp[0] = float('inf')
        else:
            low[0] = -float('inf')
            upp[0] = ps_diff[idx]

        if tc_diff[idx] >= mu[1]:
            low[1] = tc_diff[idx]
            upp[1] = float('inf')
        else:
            low[1] = -float('inf')
            upp[1] = tc_diff[idx]

        p_value, i = mvnun(low, upp, mu, cov_ps_tc_diff)
        ps_tc_results_by_tf[mpbs_name].append(p_value)
        p_values.append(p_value)

    adjusted_p_values = adjust_p_values(p_values)
    for idx, mpbs_name in enumerate(mpbs_name_list):
        ps_tc_results_by_tf[mpbs_name].append(adjusted_p_values[idx])

    return ps_tc_results_by_tf


def standard(vector1, vector2):
    max_ = max(max(vector1), max(vector2))
    min_ = min(min(vector1), min(vector2))
    if max_ > min_:
        return [(e - min_) / (max_ - min_) for e in vector1], [(e - min_) / (max_ - min_) for e in vector2]
    else:
        return vector1, vector2


def adjust_p_values(p_values):
    p = np.asfarray(p_values)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]


def output_profiles(mpbs_name_list, signal_dict_by_tf, output_location, condition):
    for mpbs_name in mpbs_name_list:
        output_fname = os.path.join(output_location, "{}_{}.txt".format(mpbs_name, condition))
        with open(output_fname, "w") as f:
            f.write("\t".join(map(str, signal_dict_by_tf[mpbs_name])) + "\n")
