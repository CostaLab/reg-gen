# Import
import os
import numpy as np
from pysam import Samfile, Fastafile
from math import log, ceil, floor
from Bio import motifs
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

import pyx
from scipy import stats

# Internal
from ..Util import AuxiliaryFunctions, GenomeData
from rgt.GenomicRegionSet import GenomicRegionSet
from biasTable import BiasTable
from signalProcessing import GenomicSignal

"""
Perform differential footprints analysis based on the prediction.

Authors: Eduardo G. Gusmao, Zhijian Li
"""

dic = {"A": 0, "C": 1, "G": 2, "T": 3}


###################################################################################################
# Classes
###################################################################################################

class DiffFootprints:
    def __init__(self, organism, mpbs_file, reads_file1, reads_file2, bias_table1, bias_table2,
                 factor1, factor2, window_size, motif_ext, min_value, initial_clip,
                 downstream_ext, upstream_ext, forward_shift, reverse_shift, k_nb, output_location,
                 output_prefix):
        self.organism = organism
        self.mpbs_file = mpbs_file
        self.reads_file1 = reads_file1
        self.reads_file2 = reads_file2
        self.bias_table1 = bias_table1
        self.bias_table2 = bias_table2
        self.factor1 = factor1
        self.factor2 = factor2
        self.window_size = window_size
        self.motif_ext = motif_ext
        self.min_value = min_value
        self.initial_clip = initial_clip
        self.downstream_ext = downstream_ext
        self.upstream_ext = upstream_ext
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift
        self.k_nb = k_nb
        self.output_location = output_location
        self.output_prefix = output_prefix

    def get_bc_signal(self, chrom, start, end, bam, bias_table, genome_file_name,
                      forward_shift, reverse_shift):
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
        fSum = sum(nf[:window])
        rSum = sum(nr[:window])
        fLast = nf[0]
        rLast = nr[0]
        for i in range((window / 2), len(nf) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += nf[i + (window / 2)]
            fLast = nf[i - (window / 2) + 1]
            rSum -= rLast
            rSum += nr[i + (window / 2)]
            rLast = nr[i - (window / 2) + 1]

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
        fSum = sum(af[:window])
        rSum = sum(ar[:window])
        fLast = af[0]
        rLast = ar[0]
        bc_signal = []
        for i in range((window / 2), len(af) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (ar[i] / rSum)
            bc_signal.append(nhatf + nhatr)
            fSum -= fLast
            fSum += af[i + (window / 2)]
            fLast = af[i - (window / 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window / 2)]
            rLast = ar[i - (window / 2) + 1]

        # Termination
        fastaFile.close()
        return bc_signal

    def get_stat_results(self, mean_signal_1, mean_signal_2, motif_len):
        signal_half_len = len(mean_signal_1) / 2

        nc = sum(mean_signal_1[signal_half_len - motif_len / 2:signal_half_len + motif_len / 2])
        nr = sum(mean_signal_1[signal_half_len + motif_len / 2:signal_half_len + motif_len / 2 + motif_len])
        nl = sum(mean_signal_1[signal_half_len - motif_len / 2 - motif_len:signal_half_len - motif_len / 2])

        protect_score1 = (nr - nc) / motif_len + (nl - nc) / motif_len
        fos1 = (nr + 1) / (nc + 1) + (nl + 1) / (nc + 1)
        tc1 = sum(mean_signal_1) / len(mean_signal_1)

        nc = sum(mean_signal_2[signal_half_len - motif_len / 2:signal_half_len + motif_len / 2])
        nr = sum(mean_signal_2[signal_half_len + motif_len / 2:signal_half_len + motif_len / 2 + motif_len])
        nl = sum(mean_signal_2[signal_half_len - motif_len / 2 - motif_len:signal_half_len - motif_len / 2])

        protect_score2 = (nr - nc) / motif_len + (nl - nc) / motif_len
        fos2 = (nr + 1) / (nc + 1) + (nl + 1) / (nc + 1)
        tc2 = sum(mean_signal_2) / len(mean_signal_2)

        protect_diff = protect_score1 - protect_score2
        tc_diff = tc1 - tc2
        fos_diff = fos1 - fos2
        return [protect_score1, protect_score2, protect_diff, tc1, tc2, tc_diff, fos1, fos2, fos_diff]

    def update_pwm(self, pwm, fasta, region, p1, p2):
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

    def run(self):
        mpbs = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs.read(self.mpbs_file)

        mpbs_name_list = list(set(mpbs.get_names()))

        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        bam1 = Samfile(self.reads_file1, "rb")
        bam2 = Samfile(self.reads_file2, "rb")

        signal_dict_by_tf_1 = dict()
        signal_dict_by_tf_2 = dict()
        motif_len_dict = dict()
        pwm_dict_by_tf = dict()

        if self.bias_table1 is None or self.bias_table2 is None:
            # differential analysis using raw reads number
            for mpbs_name in mpbs_name_list:
                signal_dict_by_tf_1[mpbs_name] = list()
                signal_dict_by_tf_2[mpbs_name] = list()
                pwm_dict_by_tf[mpbs_name] = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                                                  ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                                                  ("N", [0.0] * self.window_size)])
                motif_len_dict[mpbs_name] = 0

                mpbs_regions = mpbs.by_names([mpbs_name])
                for region in mpbs_regions:
                    if motif_len_dict[mpbs_name] == 0:
                        motif_len_dict[mpbs_name] = region.final - region.initial

                    mid = (region.final + region.initial) / 2
                    p1 = max(mid - self.window_size / 2, 0)
                    p2 = mid + self.window_size / 2

                    # Fetch raw signal
                    tc1 = np.zeros(self.window_size)
                    for read in bam1.fetch(region.chrom, p1, p2):
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                tc1[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                tc1[cut_site - p1] += 1.0
                    signal_dict_by_tf_1[mpbs_name].append(tc1.tolist())

                    tc2 = np.zeros(self.window_size)
                    for read in bam2.fetch(region.chrom, p1, p2):
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                tc2[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                tc2[cut_site - p1] += 1.0
                    signal_dict_by_tf_2[mpbs_name].append(tc2.tolist())
                    self.update_pwm(pwm_dict_by_tf[mpbs_name], fasta, region, p1, p2)
        else:
            # using bias corrected signal
            bias_table1 = None
            bias_table2 = None
            if self.bias_table1:
                table_list = self.bias_table1.split(",")
                bias_table1 = BiasTable().load_table(table_file_name_F=table_list[0], table_file_name_R=table_list[1])
            if self.bias_table2:
                table_list = self.bias_table2.split(",")
                bias_table2 = BiasTable().load_table(table_file_name_F=table_list[0], table_file_name_R=table_list[1])

            for mpbs_name in mpbs_name_list:
                signal_dict_by_tf_1[mpbs_name] = list()
                signal_dict_by_tf_2[mpbs_name] = list()
                pwm_dict_by_tf[mpbs_name] = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                                                  ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                                                  ("N", [0.0] * self.window_size)])
                motif_len_dict[mpbs_name] = 0

                mpbs_regions = mpbs.by_names([mpbs_name])
                for region in mpbs_regions:
                    if motif_len_dict[mpbs_name] == 0:
                        motif_len_dict[mpbs_name] = region.final - region.initial

                    mid = (region.final + region.initial) / 2
                    p1 = max(mid - self.window_size / 2, 0)
                    p2 = mid + self.window_size / 2

                    # Fetch bias corrected signal
                    signal_1 = self.get_bc_signal(chrom=region.chrom, start=p1, end=p2, bam=bam1,
                                                  bias_table=bias_table1, genome_file_name=genome_data.get_genome(),
                                                  forward_shift=self.forward_shift, reverse_shift=self.reverse_shift)
                    signal_dict_by_tf_1[mpbs_name].append(signal_1)

                    signal_2 = self.get_bc_signal(chrom=region.chrom, start=p1, end=p2, bam=bam2,
                                                  bias_table=bias_table2, genome_file_name=genome_data.get_genome(),
                                                  forward_shift=self.forward_shift, reverse_shift=self.reverse_shift)
                    signal_dict_by_tf_2[mpbs_name].append(signal_2)

                    self.update_pwm(pwm_dict_by_tf[mpbs_name], fasta, region, p1, p2)

        if self.factor1 is None or self.factor2 is None:
            self.factor1, self.factor2 = self.compute_factors(signal_dict_by_tf_1, signal_dict_by_tf_2)
            self.output_factor(self.factor1, self.factor2)

        stat_results_by_tf = dict()
        for mpbs_name in mpbs_name_list:
            num_fp = len(signal_dict_by_tf_1[mpbs_name])
            mean_signal_1 = np.zeros(self.window_size)
            mean_signal_2 = np.zeros(self.window_size)

            for i in range(num_fp):
                mean_signal_1 = np.add(mean_signal_1, np.array(signal_dict_by_tf_1[mpbs_name][i]))
                mean_signal_2 = np.add(mean_signal_2, np.array(signal_dict_by_tf_2[mpbs_name][i]))

            mean_signal_1 = mean_signal_1 / (num_fp * self.factor1)
            mean_signal_2 = mean_signal_2 / (num_fp * self.factor2)

            self.plot(mpbs_name, num_fp, mean_signal_1, mean_signal_2, pwm_dict_by_tf[mpbs_name], self.output_location,
                      self.output_prefix)

            stat_results_by_tf[mpbs_name] = self.get_stat_results(mean_signal_1, mean_signal_2,
                                                                  motif_len_dict[mpbs_name])

        self.output_stat_results(stat_results_by_tf)

    def compute_factors(self, signal_dict_by_tf_1, signal_dict_by_tf_2):
        keys = signal_dict_by_tf_1.keys()

        signal_1 = np.zeros(len(keys))
        signal_2 = np.zeros(len(keys))

        for key in keys:
            for i in range(len(signal_dict_by_tf_1[key])):
                idx = keys.index(key)
                signal_1[idx] += sum(signal_dict_by_tf_1[key][i])
                signal_2[idx] += sum(signal_dict_by_tf_2[key][i])

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

    def plot(self, mpbs_name, num_fp, mean_signal_1, mean_signal_2, pwm_dict, output_location, output_prefix):
        # Output PWM and create logo
        pwm_fname = os.path.join(output_location, "{}.pwm".format(mpbs_name))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(output_location, "{}.logo.eps".format(mpbs_name))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig, ax = plt.subplots()

        cell_names = output_prefix.split("_")
        ax.plot(x, mean_signal_1, color='red', label=cell_names[0])
        ax.plot(x, mean_signal_2, color='blue', label=cell_names[1])
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
        ax.spines['bottom'].set_position(('outward', 40))

        figure_name = os.path.join(output_location, "{}.line.eps".format(mpbs_name))
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(output_location, "{}.eps".format(mpbs_name))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(0.45, 0.8, logo_fname, width=16.5, height=3))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + output_fname)
        os.remove(figure_name)
        os.remove(logo_fname)
        os.remove(output_fname)

    def output_stat_results(self, stat_results_by_tf):
        output_fname = os.path.join(self.output_location, "footprint_statistics.txt")
        conditions = self.output_prefix.split("_")
        header = ["Motif",
                  "Protection_Score_{}".format(conditions[0]), "Protection_Score_{}".format(conditions[1]),
                  "Protection_Diff_{}_{}".format(conditions[0], conditions[1]),
                  "TC_{}".format(conditions[0]), "TC_{}".format(conditions[1]),
                  "TC_Diff_{}_{}".format(conditions[0], conditions[1]),
                  "FOS_{}".format(conditions[0]), "FOS_{}".format(conditions[1]),
                  "FOS_Diff_{}_{}".format(conditions[0], conditions[1])]
        with open(output_fname, "w") as f:
            f.write("\t".join(header) + "\n")
            for mpbs_name in stat_results_by_tf.keys():
                f.write(mpbs_name + "\t" + "\t".join(map(str, stat_results_by_tf[mpbs_name])) + "\n")

    def output_factor(self, factor1, factor2):
        output_file = os.path.join(self.output_location, "{}_factor.txt".format(self.output_prefix))
        f = open(output_file, "w")
        f.write("Factor1: " + str(factor1) + "\n")
        f.write("Factor2: " + str(factor2) + "\n")
        f.close()

    def output_mu(self, median_diff_prot, median_diff_tc):
        output_file = os.path.join(self.output_location, "{}_mu.txt".format(self.output_prefix))
        f = open(output_file, "w")
        f.write("median_diff_prot: " + str(median_diff_prot) + "\n")
        f.write("median_diff_tc: " + str(median_diff_tc) + "\n")
        f.close()

    def hotellings(self, X, mu):
        '''
        One-sample Hotelling's T2 test.
        :param X:
        :param mu:
        :return:
        '''
        x = X - mu
        n = x.shape[0]
        k = x.shape[1]
        m = x.mean(axis=0)  # mean vector
        S = np.cov(x.T)  # covariance
        t2 = n * m * np.linalg.inv(S) * m.T
        pvalue = stats.chi2.sf(np.sqrt(t2), k)

        return t2.tolist()[0][0], pvalue[0][0]
