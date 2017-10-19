###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import numpy as np
from pysam import Fastafile, Samfile
from Bio import motifs
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyx
from math import log, ceil, floor, isnan
# Internal
from ..Util import GenomeData
from signalProcessing import GenomicSignal
from rgt.GenomicRegionSet import GenomicRegionSet
from biasTable import BiasTable
from ..Util import AuxiliaryFunctions
import collections
from scipy.signal import savgol_filter
from scipy.stats import scoreatpercentile

class Plot:
    """

    """

    def __init__(self, organism, reads_file, motif_file, window_size,
                 downstream_ext, upstream_ext, forward_shift, reverse_shift,
                 initial_clip, bias_table, k_nb, output_loc, output_prefix):
        self.organism = organism
        self.reads_file = reads_file
        self.motif_file = motif_file
        self.window_size = window_size
        self.downstream_ext = downstream_ext
        self.upstream_ext = upstream_ext
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift
        self.initial_clip = initial_clip
        self.bias_table = bias_table
        self.k_nb = k_nb
        self.output_loc = output_loc
        self.output_prefix = output_prefix

    def line1(self):
        signal = GenomicSignal(self.reads_file)
        signal.load_sg_coefs(slope_window_size=9)
        bias_table = BiasTable()
        bias_table_list = self.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        num_sites = 0
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        mean_norm_signal = None
        mean_norm_signal_f = None
        mean_norm_signal_r = None

        pwm_dict = None

        size = 0
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                # mid = (region.initial + region.final) / 2
                # p1 = mid - (self.window_size / 2)
                # p2 = mid + (self.window_size / 2)

                p1 = region.initial - (self.window_size / 2)
                p2 = region.final + (self.window_size / 2)

                size = p2 - p1
                # Fetch raw signal
                norm_signal, norm_signal_f, norm_signal_r = \
                    self.get_signal1(ref=region.chrom, start=p1, end=p2, bam=bam, fasta=fasta, bias_table=table,
                                     signal=signal)

                num_sites += 1

                if mean_norm_signal is None:
                    mean_norm_signal = np.zeros(p2 - p1)

                if mean_norm_signal_f is None:
                    mean_norm_signal_f = np.zeros(p2 - p1)

                if mean_norm_signal_r is None:
                    mean_norm_signal_r = np.zeros(p2 - p1)

                mean_norm_signal = np.add(mean_norm_signal, norm_signal)
                mean_norm_signal_f = np.add(mean_norm_signal_f, norm_signal_f)
                mean_norm_signal_r = np.add(mean_norm_signal_r, norm_signal_r)

                # Update pwm

                if pwm_dict is None:
                    pwm_dict = dict([("A", [0.0] * (p2 - p1)), ("C", [0.0] * (p2 - p1)),
                                     ("G", [0.0] * (p2 - p1)), ("T", [0.0] * (p2 - p1)),
                                     ("N", [0.0] * (p2 - p1))])

                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        # mean_norm_signal = mean_norm_signal / num_sites
        # mean_norm_signal_f = mean_norm_signal_f / num_sites
        # mean_norm_signal_r = mean_norm_signal_r / num_sites

        mean_norm_signal_f = self.rescaling(mean_norm_signal_f)
        mean_norm_signal_r = self.rescaling(mean_norm_signal_r)

        # mean_norm_signal = signal.boyle_norm(mean_norm_signal)
        # perc = scoreatpercentile(mean_norm_signal, 98)
        # std = np.std(mean_norm_signal)
        # mean_norm_signal = signal.hon_norm_atac(mean_norm_signal, perc, std)

        # mean_norm_signal_f = signal.boyle_norm(mean_norm_signal_f)
        # perc = scoreatpercentile(mean_norm_signal_f, 98)
        # std = np.std(mean_norm_signal_f)
        # mean_norm_signal_f = signal.hon_norm_atac(mean_norm_signal_f, perc, std)

        # mean_norm_signal_r = signal.boyle_norm(mean_norm_signal_r)
        # perc = scoreatpercentile(mean_norm_signal_r, 98)
        # std = np.std(mean_norm_signal_r)
        # mean_norm_signal_r = signal.hon_norm_atac(mean_norm_signal_r, perc, std)

        # mean_slope_signal = signal.slope(mean_norm_signal, signal.sg_coefs)
        # mean_slope_signal_f = signal.slope(mean_norm_signal_f, signal.sg_coefs)
        # mean_slope_signal_r = signal.slope(mean_norm_signal_r, signal.sg_coefs)

        # mean_slope_signal = signal.boyle_norm(mean_slope_signal)
        # perc = scoreatpercentile(mean_slope_signal, 98)
        # std = np.std(mean_slope_signal)
        # mean_slope_signal = signal.hon_norm_atac(mean_slope_signal, perc, std)

        # mean_slope_signal_f = signal.boyle_norm(mean_slope_signal_f)
        # perc = scoreatpercentile(mean_slope_signal_f, 98)
        # std = np.std(mean_slope_signal_f)
        # mean_slope_signal_f = signal.hon_norm_atac(mean_slope_signal_f, perc, std)

        # mean_slope_signal_r = signal.boyle_norm(mean_slope_signal_r)
        # perc = scoreatpercentile(mean_slope_signal_r, 98)
        # std = np.std(mean_slope_signal_r)
        # mean_slope_signal_r = signal.hon_norm_atac(mean_slope_signal_r, perc, std)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        # f.write("\t".join((map(str, mean_norm_signal))) + "\n")
        # f.write("\t".join((map(str, mean_slope_signal))) + "\n")
        f.write("\t".join((map(str, mean_norm_signal_f))) + "\n")
        # f.write("\t".join((map(str, mean_slope_signal_f))) + "\n")
        f.write("\t".join((map(str, mean_norm_signal_r))) + "\n")
        # f.write("\t".join((map(str, mean_slope_signal_r))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=True, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(size / 2)
        end = (size / 2) - 1
        x = np.linspace(start, end, num=size)

        fig = plt.figure(figsize=(8, 4))
        ax2 = fig.add_subplot(111)

        min_signal = min(min(mean_norm_signal_f), min(mean_norm_signal_r))
        max_signal = max(max(mean_norm_signal_f), max(mean_norm_signal_r))
        ax2.plot(x, mean_norm_signal_f, color='red', label='Forward')
        ax2.plot(x, mean_norm_signal_r, color='green', label='Reverse')
        # ax2.plot(x, mean_norm_signal_f, color='red')
        # ax2.plot(x, mean_norm_signal_r, color='green')
        ax2.set_title(self.output_prefix, fontweight='bold')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])
        ax2.legend(loc="upper right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        # ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        # ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(1.48, 0.92, logo_fname, width=18.5, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        # os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        # os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        # os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def get_signal1(self, ref, start, end, bam, fasta, bias_table, signal):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window / 2)
        p2_w = p2 + (window / 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        signal_raw_f = [0.0] * (p2_w - p1_w)
        signal_raw_r = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + self.forward_shift
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(signal_raw_f[:window])
        rSum = sum(signal_raw_r[:window])
        fLast = signal_raw_f[0]
        rLast = signal_raw_r[0]
        for i in range((window / 2), len(signal_raw_f) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += signal_raw_f[i + (window / 2)]
            fLast = signal_raw_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_raw_r[i + (window / 2)]
            rLast = signal_raw_r[i - (window / 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        signal_bc = []
        signal_bc_f = []
        signal_bc_r = []
        for i in range((window / 2), len(signal_bias_f) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (signal_bias_r[i] / rSum)
            signal_bc.append(nhatf + nhatr)
            signal_bc_f.append(nhatf)
            signal_bc_r.append(nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window / 2)]
            fLast = signal_bias_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window / 2)]
            rLast = signal_bias_r[i - (window / 2) + 1]

        return signal_bc, signal_bc_f, signal_bc_r

    def line2(self):
        signal = GenomicSignal(self.reads_file)
        signal.load_sg_coefs(slope_window_size=9)
        bias_table = BiasTable()
        bias_table_list = self.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        num_sites = 0

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)

        bam = Samfile(self.reads_file, "rb")

        mean_signal_bias_f = np.zeros(self.window_size)
        mean_signal_bias_r = np.zeros(self.window_size)
        mean_signal_raw = np.zeros(self.window_size)
        mean_signal_bc = np.zeros(self.window_size)
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by window_size
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                signal_bias_f, signal_bias_r, signal_raw, signal_bc = \
                    self.get_signal2(ref=region.chrom, start=p1, end=p2, bam=bam, fasta=fasta, bias_table=table)

                num_sites += 1
                mean_signal_bias_f = np.add(mean_signal_bias_f, np.array(signal_bias_f))
                mean_signal_bias_r = np.add(mean_signal_bias_r, np.array(signal_bias_r))
                mean_signal_raw = np.add(mean_signal_raw, np.array(signal_raw))
                mean_signal_bc = np.add(mean_signal_bc, np.array(signal_bc))

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        mean_signal_bias_f = mean_signal_bias_f / num_sites
        mean_signal_bias_r = mean_signal_bias_r / num_sites
        mean_signal_raw = mean_signal_raw / num_sites
        mean_signal_bc = mean_signal_bc / num_sites

        mean_signal_raw = self.rescaling(mean_signal_raw)
        mean_signal_bc = self.rescaling(mean_signal_bc)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_signal_bias_f))) + "\n")
        f.write("\t".join((map(str, mean_signal_bias_r))) + "\n")
        f.write("\t".join((map(str, mean_signal_raw))) + "\n")
        f.write("\t".join((map(str, mean_signal_bc))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")

        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5))

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        ############################################################
        # bias signal per strand
        min_ = min(min(mean_signal_bias_f), min(mean_signal_bias_r))
        max_ = max(max(mean_signal_bias_f), max(mean_signal_bias_r))
        ax1.plot(x, mean_signal_bias_f, color='purple', label='Forward')
        ax1.plot(x, mean_signal_bias_r, color='green', label='Reverse')
        # ax1.plot(x, mean_signal_bias_f, color='purple')
        # ax1.plot(x, mean_signal_bias_r, color='green')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.spines['bottom'].set_position(('outward', 5))
        ax1.tick_params(direction='out')
        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        ax1.set_yticks([min_, max_])
        ax1.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
        ax1.set_title(self.output_prefix, fontweight='bold')
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_, max_])
        ax1.legend(loc="upper right", frameon=False)
        # ax1.set_ylabel("Bias Signal", rotation=90, fontweight='bold')
        ####################################################################

        #####################################################################
        # Bias corrected, non-bias corrected (not strand specific)
        min_ = min(min(mean_signal_raw), min(mean_signal_bc))
        max_ = max(max(mean_signal_raw), max(mean_signal_bc))
        ax2.plot(x, mean_signal_raw, color='blue', label='Uncorrected')
        ax2.plot(x, mean_signal_bc, color='red', label='Corrected')
        # ax2.plot(x, mean_signal_raw, color='blue')
        # ax2.plot(x, mean_signal_bc, color='red')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_, max_])
        ax2.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_, max_])
        ax2.legend(loc="upper right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        # ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        # ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')
        ###################################################################################

        ###############################################################################
        # merge the above figures
        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        # c.insert(pyx.epsfile.epsfile(1.10, 0.92, logo_fname, width=23.35, height=1.75))
        c.insert(pyx.epsfile.epsfile(1.60, 0.92, logo_fname, width=23.35, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def get_signal2(self, ref, start, end, bam, fasta, bias_table):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window / 2)
        p2_w = p2 + (window / 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        signal_raw_f = [0.0] * (p2_w - p1_w)
        signal_raw_r = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + self.forward_shift
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(signal_raw_f[:window])
        rSum = sum(signal_raw_r[:window])
        fLast = signal_raw_f[0]
        rLast = signal_raw_r[0]
        for i in range((window / 2), len(signal_raw_f) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += signal_raw_f[i + (window / 2)]
            fLast = signal_raw_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_raw_r[i + (window / 2)]
            rLast = signal_raw_r[i - (window / 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        signal_bias_forward = []
        signal_bias_reverse = []
        signal_raw = []
        signal_bc = []
        for i in range((window / 2), len(signal_bias_f) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (signal_bias_r[i] / rSum)
            signal_bias_forward.append(signal_bias_f[i])
            signal_bias_reverse.append(signal_bias_r[i])
            signal_raw.append(signal_raw_f[i] + signal_raw_r[i])
            signal_bc.append(nhatf + nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window / 2)]
            fLast = signal_bias_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window / 2)]
            rLast = signal_bias_r[i - (window / 2) + 1]

        return signal_bias_forward, signal_bias_reverse, signal_raw, signal_bc

    def line3(self, bias_table1, bias_table2):
        signal = GenomicSignal(self.reads_file)
        signal.load_sg_coefs(slope_window_size=9)
        bias_table = BiasTable()
        bias_table_list = bias_table1.split(",")
        table1 = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                       table_file_name_R=bias_table_list[1])

        bias_table = BiasTable()
        bias_table_list = bias_table2.split(",")
        table2 = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                       table_file_name_R=bias_table_list[1])

        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        num_sites = 0
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        mean_signal_raw = np.zeros(self.window_size)
        mean_signal_bc1 = np.zeros(self.window_size)
        mean_signal_bc2 = np.zeros(self.window_size)
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by window_size
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                signal_raw, signal_bc1 = \
                    self.get_signal3(ref=region.chrom, start=p1, end=p2, bam=bam, fasta=fasta, bias_table=table1)

                signal_raw, signal_bc2 = \
                    self.get_signal3(ref=region.chrom, start=p1, end=p2, bam=bam, fasta=fasta, bias_table=table2)

                num_sites += 1
                mean_signal_raw = np.add(mean_signal_raw, np.array(signal_raw))
                mean_signal_bc1 = np.add(mean_signal_bc1, np.array(signal_bc1))
                mean_signal_bc2 = np.add(mean_signal_bc2, np.array(signal_bc2))

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(
                    str(fasta.fetch(region.chrom, p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        mean_signal_raw = mean_signal_raw / num_sites
        mean_signal_bc1 = mean_signal_bc1 / num_sites
        mean_signal_bc2 = mean_signal_bc2 / num_sites

        mean_signal_raw = self.rescaling(mean_signal_raw)
        mean_signal_bc1 = self.rescaling(mean_signal_bc1)
        mean_signal_bc2 = self.rescaling(mean_signal_bc2)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_signal_raw))) + "\n")
        f.write("\t".join((map(str, mean_signal_bc1))) + "\n")
        f.write("\t".join((map(str, mean_signal_bc2))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")

        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        fig, (ax1, ax2) = plt.subplots(2)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        ############################################################
        # bias signal per strand
        min_ = min(mean_signal_raw)
        max_ = max(mean_signal_raw)
        ax1.plot(x, mean_signal_raw, color='red')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.spines['bottom'].set_position(('outward', 5))
        ax1.tick_params(direction='out')
        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        ax1.set_yticks([min_, max_])
        ax1.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
        ax1.set_title(self.output_prefix, fontweight='bold')
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_, max_])
        ax1.legend(loc="upper right", frameon=False)
        ax1.set_ylabel("Raw Signal", rotation=90, fontweight='bold')
        ####################################################################

        #####################################################################
        # Bias corrected, non-bias corrected (not strand specific)
        min_ = min(min(mean_signal_bc1), min(mean_signal_bc2))
        max_ = max(max(mean_signal_bc1), max(mean_signal_bc2))
        ax2.plot(x, mean_signal_bc1, color='red', label='VOM_8_NonNaked')
        ax2.plot(x, mean_signal_bc2, color='green', label='KMER_6_NonNaked')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_, max_])
        ax2.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_, max_])
        ax2.legend(loc="lower right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        ax2.set_ylabel("Bias Corrected Signal", rotation=90, fontweight='bold')
        ###################################################################################

        ###############################################################################
        # merge the above figures
        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(2.0, 1.39, logo_fname, width=13.8, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def get_signal3(self, ref, start, end, bam, fasta, bias_table):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window / 2)
        p2_w = p2 + (window / 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        signal_raw_f = [0.0] * (p2_w - p1_w)
        signal_raw_r = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + self.forward_shift
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(signal_raw_f[:window])
        rSum = sum(signal_raw_r[:window])
        fLast = signal_raw_f[0]
        rLast = signal_raw_r[0]
        for i in range((window / 2), len(signal_raw_f) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += signal_raw_f[i + (window / 2)]
            fLast = signal_raw_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_raw_r[i + (window / 2)]
            rLast = signal_raw_r[i - (window / 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        signal_raw = []
        signal_bc = []
        for i in range((window / 2), len(signal_bias_f) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (signal_bias_r[i] / rSum)
            signal_raw.append(signal_raw_f[i] + signal_raw_r[i])
            signal_bc.append(nhatf + nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window / 2)]
            fLast = signal_bias_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window / 2)]
            rLast = signal_bias_r[i - (window / 2) + 1]

        return signal_raw, signal_bc

    def line4(self):
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        pwm_dict = None
        signal_raw_f = None
        signal_raw_r = None

        num_sites = 0
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # p1 = region.initial - (self.window_size / 2)
                # p2 = region.final + (self.window_size / 2)

                size = p2 - p1

                if pwm_dict is None:
                    pwm_dict = dict([("A", [0.0] * size), ("C", [0.0] * size),
                                     ("G", [0.0] * size), ("T", [0.0] * size),
                                     ("N", [0.0] * size)])
                if signal_raw_f is None:
                    signal_raw_f = np.zeros(size)

                if signal_raw_r is None:
                    signal_raw_r = np.zeros(size)

                # Fetch raw signal
                for read in bam.fetch(region.chrom, p1, p2):
                    if (not read.is_reverse):
                        cut_site = read.pos + self.forward_shift
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_r[cut_site - p1] += 1.0

                num_sites += 1

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        # mean_signal_raw_f = self.rescaling(signal_raw_f)
        # mean_signal_raw_r = self.rescaling(signal_raw_r)
        mean_signal_raw_f = signal_raw_f
        mean_signal_raw_r = signal_raw_r
        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_signal_raw_f))) + "\n")
        f.write("\t".join((map(str, mean_signal_raw_r))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(size / 2)
        end = (size / 2) - 1
        x = np.linspace(start, end, num=size)

        fig = plt.figure(figsize=(8, 4))
        ax2 = fig.add_subplot(111)

        min_signal = min(min(mean_signal_raw_f), min(mean_signal_raw_r))
        max_signal = max(max(mean_signal_raw_f), max(mean_signal_raw_r))
        ax2.plot(x, mean_signal_raw_f, color='red', label='Forward')
        ax2.plot(x, mean_signal_raw_r, color='green', label='Reverse')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_title(self.output_prefix, fontweight='bold')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])
        ax2.legend(loc="upper right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        # ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        # ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(1.48, 0.92, logo_fname, width=18.5, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def line5(self, reads_file1, reads_file2, bias_table1, bias_table2):
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        bias_table = BiasTable()
        bias_table_list = bias_table1.split(",")
        table1 = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                       table_file_name_R=bias_table_list[1])

        bias_table_list = bias_table2.split(",")
        table2 = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                       table_file_name_R=bias_table_list[1])

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam_atac = Samfile(reads_file1, "rb")
        bam_dnase = Samfile(reads_file2, "rb")

        mean_signal_atac = np.zeros(self.window_size)
        mean_signal_dnase = np.zeros(self.window_size)

        num_sites = 0
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                signal_atac = self.get_signal5(ref=region.chrom, start=p1, end=p2, bam=bam_atac,
                                               fasta=fasta, bias_table=table1,
                                               forward_shift=5, reverse_shift=-4)

                signal_dnase = self.get_signal5(ref=region.chrom, start=p1, end=p2, bam=bam_dnase,
                                                fasta=fasta, bias_table=table2,
                                                forward_shift=0, reverse_shift=0)

                num_sites += 1

                mean_signal_atac = np.add(mean_signal_atac, np.array(signal_atac))
                mean_signal_dnase = np.add(mean_signal_dnase, np.array(signal_dnase))

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        mean_signal_atac = self.rescaling(mean_signal_atac)
        mean_signal_dnase = self.rescaling(mean_signal_dnase)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_signal_atac))) + "\n")
        f.write("\t".join((map(str, mean_signal_dnase))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig = plt.figure(figsize=(8, 4))
        ax2 = fig.add_subplot(111)

        min_signal = min(min(mean_signal_atac), min(mean_signal_dnase))
        max_signal = max(max(mean_signal_atac), max(mean_signal_dnase))
        ax2.plot(x, mean_signal_atac, color='red', label='ATAC-seq')
        ax2.plot(x, mean_signal_dnase, color='green', label='DNase-seq')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_title(self.output_prefix, fontweight='bold')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])
        ax2.legend(loc="upper right", frameon=False)
        # ax2.legend(loc="lower right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        # ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        # ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(1.51, 0.89, logo_fname, width=18.3, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def get_signal5(self, ref, start, end, bam, fasta, bias_table, forward_shift, reverse_shift):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window / 2)
        p2_w = p2 + (window / 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        signal_raw_f = [0.0] * (p2_w - p1_w)
        signal_raw_r = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + self.forward_shift
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if cut_site >= p1_w and cut_site < p2_w:
                    signal_raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(signal_raw_f[:window])
        rSum = sum(signal_raw_r[:window])
        fLast = signal_raw_f[0]
        rLast = signal_raw_r[0]
        for i in range((window / 2), len(signal_raw_f) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += signal_raw_f[i + (window / 2)]
            fLast = signal_raw_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_raw_r[i + (window / 2)]
            rLast = signal_raw_r[i - (window / 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        signal_bc = []
        for i in range((window / 2), len(signal_bias_f) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (signal_bias_r[i] / rSum)
            signal_bc.append(nhatf + nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window / 2)]
            fLast = signal_bias_f[i - (window / 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window / 2)]
            rLast = signal_bias_r[i - (window / 2) + 1]

        return signal_bc

    def line6(self, reads_file1, reads_file2):
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam_atac = Samfile(reads_file1, "rb")
        bam_dnase = Samfile(reads_file2, "rb")

        mean_signal_atac = np.zeros(self.window_size)
        mean_signal_dnase = np.zeros(self.window_size)

        num_sites = 0
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                for read in bam_atac.fetch(region.chrom, p1, p2):
                    if (not read.is_reverse):
                        cut_site = read.pos + self.forward_shift
                        if cut_site >= p1 and cut_site < p2:
                            mean_signal_atac[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if cut_site >= p1 and cut_site < p2:
                            mean_signal_atac[cut_site - p1] += 1.0

                # Fetch raw signal
                for read in bam_dnase.fetch(region.chrom, p1, p2):
                    if (not read.is_reverse):
                        cut_site = read.pos
                        if cut_site >= p1 and cut_site < p2:
                            mean_signal_dnase[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend - 1
                        if cut_site >= p1 and cut_site < p2:
                            mean_signal_dnase[cut_site - p1] += 1.0

                num_sites += 1

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        mean_signal_atac = self.rescaling(mean_signal_atac)
        mean_signal_dnase = self.rescaling(mean_signal_dnase)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_signal_atac))) + "\n")
        f.write("\t".join((map(str, mean_signal_dnase))) + "\n")
        f.close()

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
        pwm_file = open(pwm_fname, "w")
        for e in ["A", "C", "G", "T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig = plt.figure(figsize=(8, 4))
        ax2 = fig.add_subplot(111)

        min_signal = min(min(mean_signal_atac), min(mean_signal_dnase))
        max_signal = max(max(mean_signal_atac), max(mean_signal_dnase))
        ax2.plot(x, mean_signal_atac, color='red', label='ATAC-seq')
        ax2.plot(x, mean_signal_dnase, color='green', label='DNase-seq')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_title(self.output_prefix, fontweight='bold')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])
        ax2.legend(loc="upper right", frameon=False)
        # ax2.legend(loc="lower right", frameon=False)

        ax2.spines['bottom'].set_position(('outward', 40))
        # ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        # ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(1.51, 0.89, logo_fname, width=18.3, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    def line7(self):
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        signal_raw_f_145 = np.zeros(self.window_size)
        signal_raw_r_145 = np.zeros(self.window_size)

        signal_raw_f_146_307 = np.zeros(self.window_size)
        signal_raw_r_146_307 = np.zeros(self.window_size)

        signal_raw_f_308_500 = np.zeros(self.window_size)
        signal_raw_r_308_500 = np.zeros(self.window_size)

        signal_raw_f_501 = np.zeros(self.window_size)
        signal_raw_r_501 = np.zeros(self.window_size)

        signal_raw_f_read_1 = np.zeros(self.window_size)
        signal_raw_r_read_1 = np.zeros(self.window_size)

        signal_raw_f_read_2 = np.zeros(self.window_size)
        signal_raw_r_read_2 = np.zeros(self.window_size)

        signal_raw_f = np.zeros(self.window_size)
        signal_raw_r = np.zeros(self.window_size)

        signal_raw = np.zeros(self.window_size)

        num_sites = 0
        for region in mpbs_regions:
            if region.orientation == "+": continue
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                for read in bam.fetch(region.chrom, p1, p2):
                    if (not read.is_reverse):
                        cut_site = read.pos + self.forward_shift
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                            signal_raw[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_r[cut_site - p1] += 1.0
                            signal_raw[cut_site - p1] += 1.0

                    if abs(read.template_length) <= 145:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_145[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_145[cut_site - p1] += 1.0

                    if abs(read.template_length) > 145 and abs(read.template_length) <= 307:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_146_307[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_146_307[cut_site - p1] += 1.0

                    if abs(read.template_length) > 307 and abs(read.template_length) <= 500:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_308_500[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_308_500[cut_site - p1] += 1.0

                    if abs(read.template_length) > 500:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_501[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_501[cut_site - p1] += 1.0

                    if read.is_read1:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_read_1[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_read_1[cut_site - p1] += 1.0

                    if read.is_read2:
                        if (not read.is_reverse):
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_f_read_2[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_r_read_2[cut_site - p1] += 1.0

                num_sites += 1

                # Update pwm
                aux_plus = 1
                dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
                if (region.final - region.initial) % 2 == 0:
                    aux_plus = 0
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
                                                                         p1 + aux_plus, p2 + aux_plus)).upper())
                if region.orientation == "+":
                    for i in range(0, len(dna_seq)):
                        pwm_dict[dna_seq[i]][i] += 1
                elif region.orientation == "-":
                    for i in range(0, len(dna_seq_rev)):
                        pwm_dict[dna_seq_rev[i]][i] += 1

        # Apply a Savitzky-Golay filter to an array.
        # signal_raw = savgol_filter(signal_raw, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f = savgol_filter(signal_raw_f, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r = savgol_filter(signal_raw_r, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_145 = savgol_filter(signal_raw_f_145, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_145 = savgol_filter(signal_raw_r_145, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_146_307 = savgol_filter(signal_raw_f_146_307, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_146_307 = savgol_filter(signal_raw_r_146_307, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_308_500 = savgol_filter(signal_raw_f_308_500, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_308_500 = savgol_filter(signal_raw_r_308_500, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_501 = savgol_filter(signal_raw_f_501, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_501 = savgol_filter(signal_raw_r_501, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_read_1 = savgol_filter(signal_raw_f_read_1, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_read_1 = savgol_filter(signal_raw_r_read_1, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_f_read_2 = savgol_filter(signal_raw_f_read_2, window_length=11, polyorder=2, mode='nearest')
        # signal_raw_r_read_2 = savgol_filter(signal_raw_r_read_2, window_length=11, polyorder=2, mode='nearest')

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, signal_raw))) + "\n")
        f.write("\t".join((map(str, signal_raw_f))) + "\n")
        f.write("\t".join((map(str, signal_raw_r))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_145))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_145))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_146_307))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_146_307))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_308_500))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_308_500))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_501))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_501))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_read_1))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_read_1))) + "\n")
        f.write("\t".join((map(str, signal_raw_f_read_2))) + "\n")
        f.write("\t".join((map(str, signal_raw_r_read_2))) + "\n")
        f.close()

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6), (ax7, ax8)) = plt.subplots(4, 2, figsize=(12, 12))

        min_signal = min(signal_raw)
        max_signal = max(signal_raw)
        ax1.plot(x, signal_raw, color='red')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.tick_params(direction='out')
        ax1.set_title(self.output_prefix + " of all reads", fontweight='bold')
        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        ax1.set_yticks([min_signal, max_signal])
        ax1.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_signal, max_signal])

        min_signal = min(min(signal_raw_f), min(signal_raw_r))
        max_signal = max(max(signal_raw_f), max(signal_raw_r))
        ax2.plot(x, signal_raw_f, color='red', label='Forward')
        ax2.plot(x, signal_raw_r, color='green', label='Reverse')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_title(self.output_prefix + " of all reads", fontweight='bold')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])
        ax2.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_145), min(signal_raw_r_145))
        max_signal = max(max(signal_raw_f_145), max(signal_raw_r_145))
        max_f_index = np.argmax(signal_raw_f_145) - 500
        max_r_index = np.argmax(signal_raw_r_145) - 500
        ax3.vlines(max_r_index, min_signal, max_signal)
        ax3.vlines(max_f_index, min_signal, max_signal)
        ax3.plot(x, signal_raw_f_145, color='red', label='Forward')
        ax3.plot(x, signal_raw_r_145, color='green', label='Reverse')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_position(('outward', 15))
        ax3.tick_params(direction='out')
        ax3.set_title(self.output_prefix + " of length <= 145", fontweight='bold')
        ax3.set_xticks([start, max_f_index, 0, max_r_index, end])
        ax3.set_xticklabels([str(start), str(max_f_index), 0, str(max_r_index), str(end)])
        ax3.set_yticks([min_signal, max_signal])
        ax3.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax3.set_xlim(start, end)
        ax3.set_ylim([min_signal, max_signal])

        ax3.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_146_307), min(signal_raw_r_146_307))
        max_signal = max(max(signal_raw_f_146_307), max(signal_raw_r_146_307))
        max_f_index_1 = np.argmax(signal_raw_f_146_307[:self.window_size / 2]) - 500
        max_f_index_2 = np.argmax(signal_raw_f_146_307[self.window_size / 2:])
        max_r_index_1 = np.argmax(signal_raw_r_146_307[:self.window_size / 2]) - 500
        max_r_index_2 = np.argmax(signal_raw_r_146_307[self.window_size / 2:])
        ax4.vlines(max_f_index_1, min_signal, max_signal)
        ax4.vlines(max_f_index_2, min_signal, max_signal)
        ax4.vlines(max_r_index_1, min_signal, max_signal)
        ax4.vlines(max_r_index_2, min_signal, max_signal)

        ax4.plot(x, signal_raw_f_146_307, color='red', label='Forward')
        ax4.plot(x, signal_raw_r_146_307, color='green', label='Reverse')
        ax4.xaxis.set_ticks_position('bottom')
        ax4.yaxis.set_ticks_position('left')
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_position(('outward', 15))
        ax4.tick_params(direction='out')
        ax4.set_title(self.output_prefix + " of length > 145 and <=307", fontweight='bold')
        ax4.set_xticks([start, max_f_index_1, max_r_index_1, 0, max_f_index_2, max_r_index_2, end])
        ax4.set_xticklabels(
            [str(start), str(max_f_index_1), str(max_r_index_1), 0, str(max_f_index_2), str(max_r_index_2), str(end)])
        ax4.set_yticks([min_signal, max_signal])
        ax4.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax4.set_xlim(start, end)
        ax4.set_ylim([min_signal, max_signal])
        ax4.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_308_500), min(signal_raw_r_308_500))
        max_signal = max(max(signal_raw_f_308_500), max(signal_raw_r_308_500))
        max_f_index_1 = np.argmax(signal_raw_f_308_500[:self.window_size / 2]) - 500
        max_f_index_2 = np.argmax(signal_raw_f_308_500[self.window_size / 2:])
        max_r_index_1 = np.argmax(signal_raw_r_308_500[:self.window_size / 2]) - 500
        max_r_index_2 = np.argmax(signal_raw_r_308_500[self.window_size / 2:])
        ax5.vlines(max_f_index_1, min_signal, max_signal)
        ax5.vlines(max_f_index_2, min_signal, max_signal)
        ax5.vlines(max_r_index_1, min_signal, max_signal)
        ax5.vlines(max_r_index_2, min_signal, max_signal)

        ax5.plot(x, signal_raw_f_308_500, color='red', label='Forward')
        ax5.plot(x, signal_raw_r_308_500, color='green', label='Reverse')
        ax5.xaxis.set_ticks_position('bottom')
        ax5.yaxis.set_ticks_position('left')
        ax5.spines['top'].set_visible(False)
        ax5.spines['right'].set_visible(False)
        ax5.spines['left'].set_position(('outward', 15))
        ax5.tick_params(direction='out')
        ax5.set_title(self.output_prefix + " of length > 307 and <= 500", fontweight='bold')
        ax5.set_xticks([start, max_f_index_1, max_r_index_1, 0, max_f_index_2, max_r_index_2, end])
        ax5.set_xticklabels(
            [str(start), str(max_f_index_1), str(max_r_index_1), 0, str(max_f_index_2), str(max_r_index_2), str(end)])
        ax5.set_yticks([min_signal, max_signal])
        ax5.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax5.set_xlim(start, end)
        ax5.set_ylim([min_signal, max_signal])
        ax5.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_501), min(signal_raw_r_501))
        max_signal = max(max(signal_raw_f_501), max(signal_raw_r_501))
        ax6.plot(x, signal_raw_f_501, color='red', label='Forward')
        ax6.plot(x, signal_raw_r_501, color='green', label='Reverse')
        ax6.xaxis.set_ticks_position('bottom')
        ax6.yaxis.set_ticks_position('left')
        ax6.spines['top'].set_visible(False)
        ax6.spines['right'].set_visible(False)
        ax6.spines['left'].set_position(('outward', 15))
        ax6.tick_params(direction='out')
        ax6.set_title(self.output_prefix + " of length > 500", fontweight='bold')
        ax6.set_xticks([start, 0, end])
        ax6.set_xticklabels([str(start), 0, str(end)])
        ax6.set_yticks([min_signal, max_signal])
        ax6.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax6.set_xlim(start, end)
        ax6.set_ylim([min_signal, max_signal])
        ax6.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_read_1), min(signal_raw_r_read_1))
        max_signal = max(max(signal_raw_f_read_1), max(signal_raw_r_read_1))
        ax7.plot(x, signal_raw_f_read_1, color='red', label='Forward')
        ax7.plot(x, signal_raw_r_read_1, color='green', label='Reverse')
        ax7.xaxis.set_ticks_position('bottom')
        ax7.yaxis.set_ticks_position('left')
        ax7.spines['top'].set_visible(False)
        ax7.spines['right'].set_visible(False)
        ax7.spines['left'].set_position(('outward', 15))
        ax7.tick_params(direction='out')
        ax7.set_title(self.output_prefix + " of read 1", fontweight='bold')
        ax7.set_xticks([start, 0, end])
        ax7.set_xticklabels([str(start), 0, str(end)])
        ax7.set_yticks([min_signal, max_signal])
        ax7.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax7.set_xlim(start, end)
        ax7.set_ylim([min_signal, max_signal])
        ax7.legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_raw_f_read_2), min(signal_raw_r_read_2))
        max_signal = max(max(signal_raw_f_read_2), max(signal_raw_r_read_2))
        ax8.plot(x, signal_raw_f_read_2, color='red', label='Forward')
        ax8.plot(x, signal_raw_r_read_2, color='green', label='Reverse')
        ax8.xaxis.set_ticks_position('bottom')
        ax8.yaxis.set_ticks_position('left')
        ax8.spines['top'].set_visible(False)
        ax8.spines['right'].set_visible(False)
        ax8.spines['left'].set_position(('outward', 15))
        ax8.tick_params(direction='out')
        ax8.set_title(self.output_prefix + " of read 2", fontweight='bold')
        ax8.set_xticks([start, 0, end])
        ax8.set_xticklabels([str(start), 0, str(end)])
        ax8.set_yticks([min_signal, max_signal])
        ax8.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax8.set_xlim(start, end)
        ax8.set_ylim([min_signal, max_signal])
        ax8.legend(loc="upper right", frameon=False)

        figure_name = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="pdf", dpi=300)

    def line8(self, bias_tables):
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        table = BiasTable()
        bias_table_list = bias_tables.split(",")
        bias_table = table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

        signal = GenomicSignal()
        signal.load_sg_coefs(9)

        signal_bc_f_max_145 = np.zeros(self.window_size)
        signal_bc_r_max_145 = np.zeros(self.window_size)

        signal_bc_f_min_145 = np.zeros(self.window_size)
        signal_bc_r_min_145 = np.zeros(self.window_size)

        signal_bc_f_145_307 = np.zeros(self.window_size)
        signal_bc_r_145_307 = np.zeros(self.window_size)

        signal_bc_f_min_307 = np.zeros(self.window_size)
        signal_bc_r_min_307 = np.zeros(self.window_size)

        signal_bc_f_307_500 = np.zeros(self.window_size)
        signal_bc_r_307_500 = np.zeros(self.window_size)

        signal_bc_f_min_500 = np.zeros(self.window_size)
        signal_bc_r_min_500 = np.zeros(self.window_size)

        signal_bc_f = np.zeros(self.window_size)
        signal_bc_r = np.zeros(self.window_size)

        signal_bc = np.zeros(self.window_size)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                bc_f_max_145, bc_r_max_145 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=None, max_length=145, strand=True)
                bc_f_min_145, bc_r_min_145 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=145, max_length=None, strand=True)
                bc_f_145_307, bc_r_145_307 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=145, max_length=307, strand=True)
                bc_f_min_307, bc_r_min_307 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=307, max_length=None, strand=True)
                bc_f_307_500, bc_r_307_500 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=307, max_length=500, strand=True)
                bc_f_min_500, bc_r_min_500 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=500, max_length=None, strand=True)
                bc_f, bc_r = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=None, max_length=None, strand=True)
                bc = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=None, max_length=None, strand=False)

                signal_bc_f_max_145 = np.add(signal_bc_f_max_145, bc_f_max_145)
                signal_bc_r_max_145 = np.add(signal_bc_r_max_145, bc_r_max_145)
                signal_bc_f_min_145 = np.add(signal_bc_f_min_145, bc_f_min_145)
                signal_bc_r_min_145 = np.add(signal_bc_r_min_145, bc_r_min_145)
                signal_bc_f_145_307 = np.add(signal_bc_f_145_307, bc_f_145_307)
                signal_bc_r_145_307 = np.add(signal_bc_r_145_307, bc_r_145_307)
                signal_bc_f_min_307 = np.add(signal_bc_f_min_307, bc_f_min_307)
                signal_bc_r_min_307 = np.add(signal_bc_r_min_307, bc_r_min_307)
                signal_bc_f_307_500 = np.add(signal_bc_f_307_500, bc_f_307_500)
                signal_bc_r_307_500 = np.add(signal_bc_r_307_500, bc_r_307_500)
                signal_bc_f_min_500 = np.add(signal_bc_f_min_500, bc_f_min_500)
                signal_bc_r_min_500 = np.add(signal_bc_r_min_500, bc_r_min_500)
                signal_bc_f = np.add(signal_bc_f, bc_f)
                signal_bc_r = np.add(signal_bc_r, bc_r)
                signal_bc = np.add(signal_bc, bc)

        # Pre-process signal
        signal_bc = signal.boyle_norm(signal_bc)
        perc = scoreatpercentile(signal_bc, 98)
        std = np.array(signal_bc).std()
        signal_bc = signal.hon_norm_atac(signal_bc, perc, std)
        signal_bc_slope = signal.slope(signal_bc, signal.sg_coefs)

        signal_bc_f = signal.boyle_norm(signal_bc_f)
        perc = scoreatpercentile(signal_bc_f, 98)
        std = np.array(signal_bc_f).std()
        signal_bc_f = signal.hon_norm_atac(signal_bc_f, perc, std)
        signal_bc_f_slope = signal.slope(signal_bc_f, signal.sg_coefs)

        signal_bc_r = signal.boyle_norm(signal_bc_r)
        perc = scoreatpercentile(signal_bc_r, 98)
        std = np.array(signal_bc_r).std()
        signal_bc_r = signal.hon_norm_atac(signal_bc_r, perc, std)
        signal_bc_r_slope = signal.slope(signal_bc_r, signal.sg_coefs)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc))) + "\n")
            f.write("\t".join((map(str, signal_bc_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_slope))) + "\n")

        # Output signal of class a
        signal_bc_f_max_145 = signal.boyle_norm(signal_bc_f_max_145)
        perc = scoreatpercentile(signal_bc_f_max_145, 98)
        std = np.array(signal_bc_f_max_145).std()
        signal_bc_f_max_145 = signal.hon_norm_atac(signal_bc_f_max_145, perc, std)
        signal_bc_f_max_145_slope = signal.slope(signal_bc_f_max_145, signal.sg_coefs)

        signal_bc_r_max_145 = signal.boyle_norm(signal_bc_r_max_145)
        perc = scoreatpercentile(signal_bc_r_max_145, 98)
        std = np.array(signal_bc_r_max_145).std()
        signal_bc_r_max_145 = signal.hon_norm_atac(signal_bc_r_max_145, perc, std)
        signal_bc_r_max_145_slope = signal.slope(signal_bc_r_max_145, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_a.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")

        # Output signal of class b
        signal_bc_f_min_145 = signal.boyle_norm(signal_bc_f_min_145)
        perc = scoreatpercentile(signal_bc_f_min_145, 98)
        std = np.array(signal_bc_f_min_145).std()
        signal_bc_f_min_145 = signal.hon_norm_atac(signal_bc_f_min_145, perc, std)
        signal_bc_f_min_145_slope = signal.slope(signal_bc_f_min_145, signal.sg_coefs)

        signal_bc_r_min_145 = signal.boyle_norm(signal_bc_r_min_145)
        perc = scoreatpercentile(signal_bc_r_min_145, 98)
        std = np.array(signal_bc_r_min_145).std()
        signal_bc_r_min_145 = signal.hon_norm_atac(signal_bc_r_min_145, perc, std)
        signal_bc_r_min_145_slope = signal.slope(signal_bc_r_min_145, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_b.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_145_slope))) + "\n")

        # Output signal of class c
        signal_bc_f_145_307 = signal.boyle_norm(signal_bc_f_145_307)
        perc = scoreatpercentile(signal_bc_f_145_307, 98)
        std = np.array(signal_bc_f_145_307).std()
        signal_bc_f_145_307 = signal.hon_norm_atac(signal_bc_f_145_307, perc, std)
        signal_bc_f_145_307_slope = signal.slope(signal_bc_f_145_307, signal.sg_coefs)

        signal_bc_r_145_307 = signal.boyle_norm(signal_bc_r_145_307)
        perc = scoreatpercentile(signal_bc_r_145_307, 98)
        std = np.array(signal_bc_r_145_307).std()
        signal_bc_r_145_307 = signal.hon_norm_atac(signal_bc_r_145_307, perc, std)
        signal_bc_r_145_307_slope = signal.slope(signal_bc_r_145_307, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_c.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307_slope))) + "\n")

        # Output signal of class d
        signal_bc_f_min_307 = signal.boyle_norm(signal_bc_f_min_307)
        perc = scoreatpercentile(signal_bc_f_min_307, 98)
        std = np.array(signal_bc_f_min_307).std()
        signal_bc_f_min_307 = signal.hon_norm_atac(signal_bc_f_min_307, perc, std)
        signal_bc_f_min_307_slope = signal.slope(signal_bc_f_min_307, signal.sg_coefs)

        signal_bc_r_min_307 = signal.boyle_norm(signal_bc_r_min_307)
        perc = scoreatpercentile(signal_bc_r_min_307, 98)
        std = np.array(signal_bc_r_min_307).std()
        signal_bc_r_min_307 = signal.hon_norm_atac(signal_bc_r_min_307, perc, std)
        signal_bc_r_min_307_slope = signal.slope(signal_bc_r_min_307, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_d.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_307_slope))) + "\n")

        # Output signal of class e
        signal_bc_f_307_500 = signal.boyle_norm(signal_bc_f_307_500)
        perc = scoreatpercentile(signal_bc_f_307_500, 98)
        std = np.array(signal_bc_f_307_500).std()
        signal_bc_f_307_500 = signal.hon_norm_atac(signal_bc_f_307_500, perc, std)
        signal_bc_f_307_500_slope = signal.slope(signal_bc_f_307_500, signal.sg_coefs)

        signal_bc_r_307_500 = signal.boyle_norm(signal_bc_r_307_500)
        perc = scoreatpercentile(signal_bc_r_307_500, 98)
        std = np.array(signal_bc_r_307_500).std()
        signal_bc_r_307_500 = signal.hon_norm_atac(signal_bc_r_307_500, perc, std)
        signal_bc_r_307_500_slope = signal.slope(signal_bc_r_307_500, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_e.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_307_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_307_500_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_307_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_307_500_slope))) + "\n")

        # Output signal of class f
        signal_bc_f_min_500 = signal.boyle_norm(signal_bc_f_min_500)
        perc = scoreatpercentile(signal_bc_f_min_500, 98)
        std = np.array(signal_bc_f_min_500).std()
        signal_bc_f_min_500 = signal.hon_norm_atac(signal_bc_f_min_500, perc, std)
        signal_bc_f_min_500_slope = signal.slope(signal_bc_f_min_500, signal.sg_coefs)

        signal_bc_r_min_500 = signal.boyle_norm(signal_bc_r_min_500)
        perc = scoreatpercentile(signal_bc_r_min_500, 98)
        std = np.array(signal_bc_r_min_500).std()
        signal_bc_r_min_500 = signal.hon_norm_atac(signal_bc_r_min_500, perc, std)
        signal_bc_r_min_500_slope = signal.slope(signal_bc_r_min_500, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_f.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((map(str, signal_bc_f_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_max_145_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_145_307_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_307_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_307_500_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_307_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_307_500_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_f_min_500_slope))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_500))) + "\n")
            f.write("\t".join((map(str, signal_bc_r_min_500_slope))) + "\n")


        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig, ax = plt.subplots(4, 2, figsize=(12, 12))

        min_signal = min(signal_bc)
        max_signal = max(signal_bc)
        ax[0, 0].plot(x, signal_bc, color='red')
        ax[0, 0].xaxis.set_ticks_position('bottom')
        ax[0, 0].yaxis.set_ticks_position('left')
        ax[0, 0].spines['top'].set_visible(False)
        ax[0, 0].spines['right'].set_visible(False)
        ax[0, 0].spines['left'].set_position(('outward', 15))
        ax[0, 0].tick_params(direction='out')
        ax[0, 0].set_title(self.output_prefix + " of all reads", fontweight='bold')
        ax[0, 0].set_xticks([start, 0, end])
        ax[0, 0].set_xticklabels([str(start), 0, str(end)])
        ax[0, 0].set_yticks([min_signal, max_signal])
        ax[0, 0].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[0, 0].set_xlim(start, end)
        ax[0, 0].set_ylim([min_signal, max_signal])

        min_signal = min(min(signal_bc_f), min(signal_bc_r))
        max_signal = max(max(signal_bc_f), max(signal_bc_r))
        ax[0, 1].plot(x, signal_bc_f, color='red', label='Forward')
        ax[0, 1].plot(x, signal_bc_r, color='green', label='Reverse')
        ax[0, 1].xaxis.set_ticks_position('bottom')
        ax[0, 1].yaxis.set_ticks_position('left')
        ax[0, 1].spines['top'].set_visible(False)
        ax[0, 1].spines['right'].set_visible(False)
        ax[0, 1].spines['left'].set_position(('outward', 15))
        ax[0, 1].tick_params(direction='out')
        ax[0, 1].set_title(self.output_prefix + " of all reads", fontweight='bold')
        ax[0, 1].set_xticks([start, 0, end])
        ax[0, 1].set_xticklabels([str(start), 0, str(end)])
        ax[0, 1].set_yticks([min_signal, max_signal])
        ax[0, 1].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[0, 1].set_xlim(start, end)
        ax[0, 1].set_ylim([min_signal, max_signal])
        ax[0, 1].legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_bc_f_max_145), min(signal_bc_r_max_145))
        max_signal = max(max(signal_bc_f_max_145), max(signal_bc_r_max_145))
        ax[1, 0].plot(x, signal_bc_f_max_145, color='red', label='Forward')
        ax[1, 0].plot(x, signal_bc_r_max_145, color='green', label='Reverse')
        ax[1, 0].xaxis.set_ticks_position('bottom')
        ax[1, 0].yaxis.set_ticks_position('left')
        ax[1, 0].spines['top'].set_visible(False)
        ax[1, 0].spines['right'].set_visible(False)
        ax[1, 0].spines['left'].set_position(('outward', 15))
        ax[1, 0].tick_params(direction='out')
        ax[1, 0].set_title(self.output_prefix + " from fragment length <= 145 ", fontweight='bold')
        ax[1, 0].set_xticks([start, 0, end])
        ax[1, 0].set_xticklabels([str(start), 0, str(end)])
        ax[1, 0].set_yticks([min_signal, max_signal])
        ax[1, 0].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[1, 0].set_xlim(start, end)
        ax[1, 0].set_ylim([min_signal, max_signal])
        ax[1, 0].legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_bc_f_min_145), min(signal_bc_r_min_145))
        max_signal = max(max(signal_bc_f_min_145), max(signal_bc_r_min_145))
        ax[1, 1].plot(x, signal_bc_f_min_145, color='red', label='Forward')
        ax[1, 1].plot(x, signal_bc_r_min_145, color='green', label='Reverse')
        ax[1, 1].xaxis.set_ticks_position('bottom')
        ax[1, 1].yaxis.set_ticks_position('left')
        ax[1, 1].spines['top'].set_visible(False)
        ax[1, 1].spines['right'].set_visible(False)
        ax[1, 1].spines['left'].set_position(('outward', 15))
        ax[1, 1].tick_params(direction='out')
        ax[1, 1].set_title(self.output_prefix + " from fragment length > 145 ", fontweight='bold')
        ax[1, 1].set_xticks([start, 0, end])
        ax[1, 1].set_xticklabels([str(start), 0, str(end)])
        ax[1, 1].set_yticks([min_signal, max_signal])
        ax[1, 1].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[1, 1].set_xlim(start, end)
        ax[1, 1].set_ylim([min_signal, max_signal])
        ax[1, 1].legend(loc="upper right", frameon=False)


        min_signal = min(min(signal_bc_f_145_307), min(signal_bc_r_145_307))
        max_signal = max(max(signal_bc_f_145_307), max(signal_bc_r_145_307))
        ax[2, 0].plot(x, signal_bc_f_145_307, color='red', label='Forward')
        ax[2, 0].plot(x, signal_bc_r_145_307, color='green', label='Reverse')
        ax[2, 0].xaxis.set_ticks_position('bottom')
        ax[2, 0].yaxis.set_ticks_position('left')
        ax[2, 0].spines['top'].set_visible(False)
        ax[2, 0].spines['right'].set_visible(False)
        ax[2, 0].spines['left'].set_position(('outward', 15))
        ax[2, 0].tick_params(direction='out')
        ax[2, 0].set_title(self.output_prefix + " from fragment length > 145 and <= 307 ", fontweight='bold')
        ax[2, 0].set_xticks([start, 0, end])
        ax[2, 0].set_xticklabels([str(start), 0, str(end)])
        ax[2, 0].set_yticks([min_signal, max_signal])
        ax[2, 0].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[2, 0].set_xlim(start, end)
        ax[2, 0].set_ylim([min_signal, max_signal])
        ax[2, 0].legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_bc_f_min_307), min(signal_bc_r_min_307))
        max_signal = max(max(signal_bc_f_min_307), max(signal_bc_r_min_307))
        ax[2, 1].plot(x, signal_bc_f_min_307, color='red', label='Forward')
        ax[2, 1].plot(x, signal_bc_r_min_307, color='green', label='Reverse')
        ax[2, 1].xaxis.set_ticks_position('bottom')
        ax[2, 1].yaxis.set_ticks_position('left')
        ax[2, 1].spines['top'].set_visible(False)
        ax[2, 1].spines['right'].set_visible(False)
        ax[2, 1].spines['left'].set_position(('outward', 15))
        ax[2, 1].tick_params(direction='out')
        ax[2, 1].set_title(self.output_prefix + " from fragment length > 307 ", fontweight='bold')
        ax[2, 1].set_xticks([start, 0, end])
        ax[2, 1].set_xticklabels([str(start), 0, str(end)])
        ax[2, 1].set_yticks([min_signal, max_signal])
        ax[2, 1].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[2, 1].set_xlim(start, end)
        ax[2, 1].set_ylim([min_signal, max_signal])
        ax[2, 1].legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_bc_f_307_500), min(signal_bc_r_307_500))
        max_signal = max(max(signal_bc_f_307_500), max(signal_bc_r_307_500))
        ax[3, 0].plot(x, signal_bc_f_307_500, color='red', label='Forward')
        ax[3, 0].plot(x, signal_bc_r_307_500, color='green', label='Reverse')
        ax[3, 0].xaxis.set_ticks_position('bottom')
        ax[3, 0].yaxis.set_ticks_position('left')
        ax[3, 0].spines['top'].set_visible(False)
        ax[3, 0].spines['right'].set_visible(False)
        ax[3, 0].spines['left'].set_position(('outward', 15))
        ax[3, 0].tick_params(direction='out')
        ax[3, 0].set_title(self.output_prefix + " from fragment length > 307 and <= 500 ", fontweight='bold')
        ax[3, 0].set_xticks([start, 0, end])
        ax[3, 0].set_xticklabels([str(start), 0, str(end)])
        ax[3, 0].set_yticks([min_signal, max_signal])
        ax[3, 0].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[3, 0].set_xlim(start, end)
        ax[3, 0].set_ylim([min_signal, max_signal])
        ax[3, 0].legend(loc="upper right", frameon=False)

        min_signal = min(min(signal_bc_f_min_500), min(signal_bc_r_min_500))
        max_signal = max(max(signal_bc_f_min_500), max(signal_bc_r_min_500))
        ax[3, 1].plot(x, signal_bc_f_min_500, color='red', label='Forward')
        ax[3, 1].plot(x, signal_bc_r_min_500, color='green', label='Reverse')
        ax[3, 1].xaxis.set_ticks_position('bottom')
        ax[3, 1].yaxis.set_ticks_position('left')
        ax[3, 1].spines['top'].set_visible(False)
        ax[3, 1].spines['right'].set_visible(False)
        ax[3, 1].spines['left'].set_position(('outward', 15))
        ax[3, 1].tick_params(direction='out')
        ax[3, 1].set_title(self.output_prefix + " from fragment length > 500 ", fontweight='bold')
        ax[3, 1].set_xticks([start, 0, end])
        ax[3, 1].set_xticklabels([str(start), 0, str(end)])
        ax[3, 1].set_yticks([min_signal, max_signal])
        ax[3, 1].set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax[3, 1].set_xlim(start, end)
        ax[3, 1].set_ylim([min_signal, max_signal])
        ax[3, 1].legend(loc="upper right", frameon=False)

        figure_name = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="pdf", dpi=300)

    def distribution_of_frag_length(self, reads_file=None, motif_file=None):
        bam = Samfile(self.reads_file, "rb")
        lengths = list()
        lengths_shorter_than_145 = list()
        lengths_between_146_307 = list()
        lengths_between_307_500 = list()
        lengths_longer_than_500 = list()
        if motif_file is None:
            # Use all fragments to plot the distribution
            for read in bam.fetch():
                if read.is_read1:
                    length = abs(read.template_length)
                    lengths.append(length)
                    if length <= 145: lengths_shorter_than_145.append(length)
                    if length > 145 and length <= 307: lengths_between_146_307.append(length)
                    if length > 307 and length <= 500: lengths_between_307_500.append(length)
                    if length > 500: lengths_longer_than_500.append(length)

        elif motif_file is not None:
            mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
            mpbs_regions.read(self.motif_file)
            for region in mpbs_regions:
                if str(region.name).split(":")[-1] == "Y":
                    # Extend by 50 bp
                    mid = (region.initial + region.final) / 2
                    p1 = mid - (self.window_size / 2)
                    p2 = mid + (self.window_size / 2)

                    for read in bam.fetch(region.chrom, p1, p2):
                        if read.is_read1:
                            length = abs(read.template_length)
                            lengths.append(length)
                            if length <= 145: lengths_shorter_than_145.append(length)
                            if length > 145 and length <= 307: lengths_between_146_307.append(length)
                            if length > 307 and length <= 500: lengths_between_307_500.append(length)
                            if length > 500: lengths_longer_than_500.append(length)

        # output the plot
        unique, counts = np.unique(lengths, return_counts=True)
        count_list = list()
        length_dict = dict(zip(unique, counts))
        sort_keys = sorted(length_dict.keys())
        txt_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(txt_fname, "w") as f:
            for key in sort_keys:
                f.write(str(key) + "\t" + str(length_dict[key]) + "\n")
                count_list.append(length_dict[key])

        fig, ax = plt.subplots(5, 1, figsize=(6, 15))
        mean = np.mean(lengths)
        std = np.std(lengths)
        ax[0].plot(unique, counts)
        ax[0].set_title("Histogram of all fragment length", fontweight='bold')
        ax[0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0].transAxes, fontsize=12)
        ax[0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0].transAxes, fontsize=12)

        mean = np.mean(lengths_shorter_than_145)
        std = np.std(lengths_shorter_than_145)
        unique, counts = np.unique(lengths_shorter_than_145, return_counts=True)
        ax[1].plot(unique, counts)
        ax[1].set_title("Histogram of fragment length <= 145", fontweight='bold')
        ax[1].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_shorter_than_145))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[1].transAxes, fontsize=12)
        ax[1].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[1].transAxes, fontsize=12)
        ax[1].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[1].transAxes, fontsize=12)

        mean = np.mean(lengths_between_146_307)
        std = np.std(lengths_between_146_307)
        unique, counts = np.unique(lengths_between_146_307, return_counts=True)
        ax[2].plot(unique, counts)
        ax[2].set_title("Histogram of fragment length > 145 and <= 307", fontweight='bold')
        ax[2].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_between_146_307))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[2].transAxes, fontsize=12)
        ax[2].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[2].transAxes, fontsize=12)
        ax[2].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[2].transAxes, fontsize=12)

        mean = np.mean(lengths_between_307_500)
        std = np.std(lengths_between_307_500)
        unique, counts = np.unique(lengths_between_307_500, return_counts=True)
        ax[3].plot(unique, counts)
        ax[3].set_title("Histogram of fragment length > 307 and <= 500", fontweight='bold')
        ax[3].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_between_307_500))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[3].transAxes, fontsize=12)
        ax[3].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[3].transAxes, fontsize=12)
        ax[3].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[3].transAxes, fontsize=12)

        mean = np.mean(lengths_longer_than_500)
        std = np.std(lengths_longer_than_500)
        unique, counts = np.unique(lengths_longer_than_500, return_counts=True)
        ax[4].plot(unique, counts)
        ax[4].set_title("Histogram of fragment length > 500", fontweight='bold')
        ax[4].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_longer_than_500))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[4].transAxes, fontsize=12)
        ax[4].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[4].transAxes, fontsize=12)
        ax[4].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[4].transAxes, fontsize=12)

        figure_fname = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_fname, format="pdf", dpi=300)

    def distribution_of_frag_length_by_strand(self, reads_file=None, motif_file=None):
        bam = Samfile(self.reads_file, "rb")
        lengths_forward = list()
        lengths_reverse = list()
        if motif_file is None:
            # Use all fragments to plot the distribution
            for read in bam.fetch():
                if read.is_read1 and not read.is_reverse:
                    lengths_forward.append(abs(read.template_length))
                elif read.is_read1 and read.is_reverse:
                    lengths_reverse.append(abs(read.template_length))

        elif motif_file is not None:
            mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
            mpbs_regions.read(self.motif_file)
            for region in mpbs_regions:
                if str(region.name).split(":")[-1] == "Y":
                    # Extend by 50 bp
                    mid = (region.initial + region.final) / 2
                    p1 = mid - (self.window_size / 2)
                    p2 = mid + (self.window_size / 2)

                    for read in bam.fetch(region.chrom, p1, p2):
                        if read.is_read1:
                            if read.is_read1 and not read.is_reverse:
                                lengths_forward.append(abs(read.template_length))
                            elif read.is_read1 and read.is_reverse:
                                lengths_reverse.append(abs(read.template_length))

        # output the plot
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6))

        ax1.hist(np.asarray(lengths_forward), bins='auto', histtype='stepfilled', color='red')
        ax1.set_title("Histogram for Fragment length of Forward", fontweight='bold')

        ax2.hist(np.asarray(lengths_reverse), bins='auto', histtype='stepfilled', color='blue')
        ax2.set_title("Histogram for Fragment length of Reverse", fontweight='bold')

        figure_fname = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_fname, format="pdf", dpi=300)

        unique, counts = np.unique(lengths_forward, return_counts=True)
        length_f_dict = dict(zip(unique, counts))
        unique, counts = np.unique(lengths_reverse, return_counts=True)
        length_r_dict = dict(zip(unique, counts))
        sort_keys = sorted(length_f_dict.keys())
        txt_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(txt_fname, "w") as f:
            f.write("Length\tCount_Forward\tCount_Reverse\n")
            for key in sort_keys:
                f.write(str(key) + "\t" + str(length_f_dict[key]) + "\t" + str(length_r_dict[key]) + "\n")

    def distribution_of_frag_length_by_position(self, reads_file=None, motif_file=None):
        bam = Samfile(self.reads_file, "rb")
        lengths = list()
        lengths_class_1 = list()
        lengths_class_1_f = list()
        lengths_class_1_r = list()
        lengths_class_2 = list()
        lengths_class_2_f = list()
        lengths_class_2_r = list()
        lengths_class_3 = list()
        lengths_class_3_f = list()
        lengths_class_3_r = list()
        lengths_class_4 = list()
        lengths_class_4_f = list()
        lengths_class_4_r = list()

        signal_raw_f = np.zeros(self.window_size)
        signal_raw_r = np.zeros(self.window_size)
        signal_raw_class_1_f = np.zeros(self.window_size)
        signal_raw_class_1_r = np.zeros(self.window_size)
        signal_raw_class_2_f = np.zeros(self.window_size)
        signal_raw_class_2_r = np.zeros(self.window_size)
        signal_raw_class_3_f = np.zeros(self.window_size)
        signal_raw_class_3_r = np.zeros(self.window_size)
        signal_raw_class_4_f = np.zeros(self.window_size)
        signal_raw_class_4_r = np.zeros(self.window_size)

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                for read in bam.fetch(region.chrom, p1, p2):
                    length = abs(read.template_length)
                    lengths.append(length)

                    if (not read.is_reverse):
                        cut_site = read.pos + self.forward_shift
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_r[cut_site - p1] += 1.0

                    if length <= 145:
                        # class 1
                        if not read.is_reverse:
                            if mid - 73 < read.pos < mid + 73 and mid - 73 < read.pos + length < mid + 73:
                                lengths_class_1.append(length)
                                lengths_class_1_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_1_f[cut_site - p1] += 1.0
                        else:
                            if mid - 73 < read.aend < mid + 73 and mid - 73 < read.aend - length < mid + 73:
                                lengths_class_1.append(length)
                                lengths_class_1_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_1_r[cut_site - p1] += 1.0

                    elif 145 < length <= 307:
                        # class 2
                        if not read.is_reverse:
                            if read.pos < mid - 73 and mid - 73 <= read.pos + length <= mid:
                                lengths_class_2.append(length)
                                lengths_class_2_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_2_f[cut_site - p1] += 1.0
                            elif mid <= read.pos <= mid + 73 and read.pos + length > mid + 73:
                                lengths_class_2.append(length)
                                lengths_class_2_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_2_f[cut_site - p1] += 1.0
                        else:
                            if mid - 73 <= read.aend <= mid and read.aend - length < mid - 73:
                                lengths_class_2.append(length)
                                lengths_class_2_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_2_r[cut_site - p1] += 1.0
                            elif mid <= read.aend - length <= mid + 73 and read.aend > mid + 73:
                                lengths_class_2.append(length)
                                lengths_class_2_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_2_r[cut_site - p1] += 1.0

                        # class 3
                        if not read.is_reverse:
                            if read.pos <= mid - 73 and mid <= read.pos + length <= mid + 73:
                                lengths_class_3.append(length)
                                lengths_class_3_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_3_f[cut_site - p1] += 1.0
                            elif mid - 73 <= read.pos <= mid and read.pos + length >= mid + 73:
                                lengths_class_3.append(length)
                                lengths_class_3_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_3_f[cut_site - p1] += 1.0
                        else:
                            if mid <= read.aend <= mid + 73 and read.aend - length <= mid - 73:
                                lengths_class_3.append(length)
                                lengths_class_3_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_3_r[cut_site - p1] += 1.0
                            elif read.aend >= mid + 73 and mid - 73 <= read.aend - length <= mid:
                                lengths_class_3.append(length)
                                lengths_class_3_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if cut_site >= p1 and cut_site < p2:
                                    signal_raw_class_3_r[cut_site - p1] += 1.0

                    elif length > 307:
                        lengths_class_4.append(length)
                        if (not read.is_reverse):
                            lengths_class_4_f.append(length)
                            cut_site = read.pos + self.forward_shift
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_class_4_f[cut_site - p1] += 1.0
                        else:
                            lengths_class_4_r.append(length)
                            cut_site = read.aend + self.reverse_shift - 1
                            if cut_site >= p1 and cut_site < p2:
                                signal_raw_class_4_r[cut_site - p1] += 1.0



        # output the plot
        unique, counts = np.unique(lengths, return_counts=True)
        count_list = list()
        length_dict = dict(zip(unique, counts))
        sort_keys = sorted(length_dict.keys())
        txt_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(txt_fname, "w") as f:
            for key in sort_keys:
                f.write(str(key) + "\t" + str(length_dict[key]) + "\n")
                count_list.append(length_dict[key])

        fig, ax = plt.subplots(5, 2, figsize=(10, 15))
        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)


        mean = np.mean(lengths)
        std = np.std(lengths)
        ax[0, 0].plot(unique, counts)
        ax[0, 0].set_title("Histogram of all fragment length", fontweight='bold')
        ax[0, 0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0, 0].transAxes, fontsize=12)
        ax[0, 0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0, 0].transAxes, fontsize=12)
        ax[0, 0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[0, 0].transAxes, fontsize=12)
        ax[0, 1].plot(x, signal_raw_f, color='red', label='Forward')
        ax[0, 1].plot(x, signal_raw_r, color='green', label='Reverse')
        ax[0, 1].text(0.2, 0.7, "#forward:{}".format(str(sum(signal_raw_f[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[0, 1].transAxes, fontsize=12)
        ax[0, 1].text(0.2, 0.6, "#reverse:{}".format(str(sum(signal_raw_r[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[0, 1].transAxes, fontsize=12)
        ax[0, 1].text(0.8, 0.7, "#forward:{}".format(str(sum(signal_raw_f[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[0, 1].transAxes, fontsize=12)
        ax[0, 1].text(0.8, 0.6, "#reverse:{}".format(str(sum(signal_raw_r[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[0, 1].transAxes, fontsize=12)
        ax[0, 1].legend(loc="upper right", frameon=False)

        mean = np.mean(lengths_class_1)
        std = np.std(lengths_class_1)

        unique, counts = np.unique(lengths_class_1, return_counts=True)
        ax[1, 0].plot(unique, counts, color='red', label='All')
        unique, counts = np.unique(lengths_class_1_f, return_counts=True)
        ax[1, 0].plot(unique, counts, color='blue', label='Forward')
        unique, counts = np.unique(lengths_class_1_r, return_counts=True)
        ax[1, 0].plot(unique, counts, color='green', label='Reverse')
        ax[1, 0].legend(loc="lower left", frameon=False)

        ax[1, 0].set_title("Histogram of Class 1", fontweight='bold')
        ax[1, 0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_class_1))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[1, 0].transAxes, fontsize=12)
        ax[1, 0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[1, 0].transAxes, fontsize=12)
        ax[1, 0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[1, 0].transAxes, fontsize=12)
        ax[1, 1].plot(x, signal_raw_class_1_f, color='red', label='Forward')
        ax[1, 1].plot(x, signal_raw_class_1_r, color='green', label='Reverse')
        ax[1, 1].text(0.2, 0.7, "#forward:{}".format(str(sum(signal_raw_class_1_f[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[1, 1].transAxes, fontsize=12)
        ax[1, 1].text(0.2, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_1_r[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[1, 1].transAxes, fontsize=12)
        ax[1, 1].text(0.8, 0.7, "#forward:{}".format(str(sum(signal_raw_class_1_f[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[1, 1].transAxes, fontsize=12)
        ax[1, 1].text(0.8, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_1_r[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[1, 1].transAxes, fontsize=12)
        ax[1, 1].legend(loc="upper right", frameon=False)

        mean = np.mean(lengths_class_2)
        std = np.std(lengths_class_2)

        unique, counts = np.unique(lengths_class_2, return_counts=True)
        ax[2, 0].plot(unique, counts, color='red', label='All')
        unique, counts = np.unique(lengths_class_2_f, return_counts=True)
        ax[2, 0].plot(unique, counts, color='blue', label='Forward')
        unique, counts = np.unique(lengths_class_2_r, return_counts=True)
        ax[2, 0].plot(unique, counts, color='green', label='Reverse')
        ax[2, 0].legend(loc="lower left", frameon=False)

        ax[2, 0].set_title("Histogram of Class 2", fontweight='bold')
        ax[2, 0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_class_2))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[2, 0].transAxes, fontsize=12)
        ax[2, 0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[2, 0].transAxes, fontsize=12)
        ax[2, 0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[2, 0].transAxes, fontsize=12)
        ax[2, 1].plot(x, signal_raw_class_2_f, color='red', label='Forward')
        ax[2, 1].plot(x, signal_raw_class_2_r, color='green', label='Reverse')
        ax[2, 1].text(0.2, 0.7, "#forward:{}".format(str(sum(signal_raw_class_2_f[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[2, 1].transAxes, fontsize=12)
        ax[2, 1].text(0.2, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_2_r[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[2, 1].transAxes, fontsize=12)
        ax[2, 1].text(0.8, 0.7, "#forward:{}".format(str(sum(signal_raw_class_2_f[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[2, 1].transAxes, fontsize=12)
        ax[2, 1].text(0.8, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_2_r[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[2, 1].transAxes, fontsize=12)
        ax[1, 1].legend(loc="upper right", frameon=False)
        ax[2, 1].legend(loc="upper right", frameon=False)

        mean = np.mean(lengths_class_3)
        std = np.std(lengths_class_3)

        unique, counts = np.unique(lengths_class_3, return_counts=True)
        ax[3, 0].plot(unique, counts, color='red', label='All')
        unique, counts = np.unique(lengths_class_3_f, return_counts=True)
        ax[3, 0].plot(unique, counts, color='blue', label='Forward')
        unique, counts = np.unique(lengths_class_3_r, return_counts=True)
        ax[3, 0].plot(unique, counts, color='green', label='Reverse')
        ax[3, 0].legend(loc="lower left", frameon=False)

        ax[3, 0].set_title("Histogram of Class 3", fontweight='bold')
        ax[3, 0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_class_3))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[3, 0].transAxes, fontsize=12)
        ax[3, 0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[3, 0].transAxes, fontsize=12)
        ax[3, 0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[3, 0].transAxes, fontsize=12)
        ax[3, 1].plot(x, signal_raw_class_3_f, color='red', label='Forward')
        ax[3, 1].plot(x, signal_raw_class_3_r, color='green', label='Reverse')
        ax[3, 1].text(0.2, 0.7, "#forward:{}".format(str(sum(signal_raw_class_3_f[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[3, 1].transAxes, fontsize=12)
        ax[3, 1].text(0.2, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_3_r[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[3, 1].transAxes, fontsize=12)
        ax[3, 1].text(0.8, 0.7, "#forward:{}".format(str(sum(signal_raw_class_3_f[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[3, 1].transAxes, fontsize=12)
        ax[3, 1].text(0.8, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_3_r[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[3, 1].transAxes, fontsize=12)
        ax[3, 1].legend(loc="upper right", frameon=False)


        mean = np.mean(lengths_class_4)
        std = np.std(lengths_class_4)

        unique, counts = np.unique(lengths_class_4, return_counts=True)
        ax[4, 0].plot(unique, counts, color='red', label='All')
        unique, counts = np.unique(lengths_class_4_f, return_counts=True)
        ax[4, 0].plot(unique, counts, color='blue', label='Forward')
        unique, counts = np.unique(lengths_class_4_r, return_counts=True)
        ax[4, 0].plot(unique, counts, color='green', label='Reverse')
        ax[4, 0].legend(loc="lower left", frameon=False)

        ax[4, 0].set_title("Histogram of Class 4", fontweight='bold')
        ax[4, 0].text(0.7, 0.9, "number of fragments:{}".format(str(len(lengths_class_4))),
                   verticalalignment='center',
                   horizontalalignment='center', transform=ax[4, 0].transAxes, fontsize=12)
        ax[4, 0].text(0.8, 0.8, "mean:{}".format(str(round(mean, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[4, 0].transAxes, fontsize=12)
        ax[4, 0].text(0.8, 0.7, "std:{}".format(str(round(std, 2))), verticalalignment='center',
                   horizontalalignment='center', transform=ax[4, 0].transAxes, fontsize=12)
        ax[4, 1].plot(x, signal_raw_class_4_f, color='red', label='Forward')
        ax[4, 1].plot(x, signal_raw_class_4_r, color='green', label='Reverse')
        ax[4, 1].text(0.2, 0.7, "#forward:{}".format(str(sum(signal_raw_class_4_f[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[4, 1].transAxes, fontsize=12)
        ax[4, 1].text(0.2, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_4_r[:self.window_size / 2]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[4, 1].transAxes, fontsize=12)
        ax[4, 1].text(0.8, 0.7, "#forward:{}".format(str(sum(signal_raw_class_4_f[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[4, 1].transAxes, fontsize=12)
        ax[4, 1].text(0.8, 0.6, "#reverse:{}".format(str(sum(signal_raw_class_4_r[self.window_size / 2:]))),
                      verticalalignment='center',
                      horizontalalignment='center', transform=ax[4, 1].transAxes, fontsize=12)
        ax[4, 1].legend(loc="upper right", frameon=False)

        figure_fname = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_fname, format="pdf", dpi=300)

    def rescaling(self, vector):
        maxN = max(vector)
        minN = min(vector)
        return [(e - minN) / (maxN - minN) for e in vector]
