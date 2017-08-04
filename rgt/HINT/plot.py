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
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        mean_norm_signal = np.zeros(self.window_size)
        mean_norm_signal_f = np.zeros(self.window_size)
        mean_norm_signal_r = np.zeros(self.window_size)

        num_sites = 0

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)

        bam = Samfile(self.reads_file, "rb")

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                norm_signal, norm_signal_f, norm_signal_r = \
                    self.get_signal1(ref=region.chrom, start=p1, end=p2, bam=bam, fasta=fasta, bias_table=table,
                                     signal=signal)

                num_sites += 1
                mean_norm_signal = np.add(mean_norm_signal, norm_signal)
                mean_norm_signal_f = np.add(mean_norm_signal_f, norm_signal_f)
                mean_norm_signal_r = np.add(mean_norm_signal_r, norm_signal_r)

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

        mean_norm_signal = mean_norm_signal / num_sites
        mean_norm_signal_f = mean_norm_signal_f / num_sites
        mean_norm_signal_r = mean_norm_signal_r / num_sites

        mean_norm_signal = signal.boyle_norm(mean_norm_signal)
        perc = scoreatpercentile(mean_norm_signal, 98)
        std = np.std(mean_norm_signal)
        mean_norm_signal = signal.hon_norm_atac(mean_norm_signal, perc, std)

        mean_norm_signal_f = signal.boyle_norm(mean_norm_signal_f)
        perc = scoreatpercentile(mean_norm_signal_f, 98)
        std = np.std(mean_norm_signal_f)
        mean_norm_signal_f = signal.hon_norm_atac(mean_norm_signal_f, perc, std)

        mean_norm_signal_r = signal.boyle_norm(mean_norm_signal_r)
        perc = scoreatpercentile(mean_norm_signal_r, 98)
        std = np.std(mean_norm_signal_r)
        mean_norm_signal_r = signal.hon_norm_atac(mean_norm_signal_r, perc, std)

        mean_slope_signal = signal.slope(mean_norm_signal, signal.sg_coefs)
        mean_slope_signal_f = signal.slope(mean_norm_signal_f, signal.sg_coefs)
        mean_slope_signal_r = signal.slope(mean_norm_signal_r, signal.sg_coefs)

        mean_slope_signal = signal.boyle_norm(mean_slope_signal)
        perc = scoreatpercentile(mean_slope_signal, 98)
        std = np.std(mean_slope_signal)
        mean_slope_signal = signal.hon_norm_atac(mean_slope_signal, perc, std)

        mean_slope_signal_f = signal.boyle_norm(mean_slope_signal_f)
        perc = scoreatpercentile(mean_slope_signal_f, 98)
        std = np.std(mean_slope_signal_f)
        mean_slope_signal_f = signal.hon_norm_atac(mean_slope_signal_f, perc, std)

        mean_slope_signal_r = signal.boyle_norm(mean_slope_signal_r)
        perc = scoreatpercentile(mean_slope_signal_r, 98)
        std = np.std(mean_slope_signal_r)
        mean_slope_signal_r = signal.hon_norm_atac(mean_slope_signal_r, perc, std)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((map(str, mean_norm_signal))) + "\n")
        f.write("\t".join((map(str, mean_slope_signal))) + "\n")
        f.write("\t".join((map(str, mean_norm_signal_f))) + "\n")
        f.write("\t".join((map(str, mean_slope_signal_f))) + "\n")
        f.write("\t".join((map(str, mean_norm_signal_r))) + "\n")
        f.write("\t".join((map(str, mean_slope_signal_r))) + "\n")
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

        min_signal = min(min(mean_norm_signal_f), min(mean_norm_signal_r))
        max_signal = max(max(mean_norm_signal_f), max(mean_norm_signal_r))
        ax2.plot(x, mean_norm_signal_f, color='red', label='Forward')
        ax2.plot(x, mean_norm_signal_r, color='green', label='Reverse')

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
        ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(2, 1.4, logo_fname, width=17.8, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

        os.remove(pwm_fname)
        os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
        #os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
        os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

    # def line1(self):
    #     signal = GenomicSignal(self.reads_file)
    #     signal.load_sg_coefs(slope_window_size=9)
    #     bias_table = BiasTable()
    #     bias_table_list = self.bias_table.split(",")
    #     table = bias_table.load_table(table_file_name_F=bias_table_list[0],
    #                                   table_file_name_R=bias_table_list[1])
    #
    #     genome_data = GenomeData(self.organism)
    #     fasta = Fastafile(genome_data.get_genome())
    #     pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
    #                      ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
    #                      ("N", [0.0] * self.window_size)])
    #
    #     mean_norm_signal = np.zeros(self.window_size)
    #     mean_norm_signal_f = np.zeros(self.window_size)
    #     mean_norm_signal_r = np.zeros(self.window_size)
    #
    #     num_sites = 0
    #
    #     mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    #     mpbs_regions.read(self.motif_file)
    #
    #     bam = Samfile(self.reads_file, "rb")
    #
    #     for region in mpbs_regions:
    #         if str(region.name).split(":")[-1] == "Y":
    #             # Extend by 50 bp
    #             mid = (region.initial + region.final) / 2
    #             p1 = mid - (self.window_size / 2)
    #             p2 = mid + (self.window_size / 2)
    #
    #             # Fetch raw signal
    #             norm_signal, norm_signal_f, norm_signal_r = \
    #                 self.get_signal1(ref=region.chrom, start=p1,end=p2, bam=bam, fasta=fasta, bias_table=table,
    #                                  signal=signal)
    #
    #             num_sites += 1
    #             mean_norm_signal = np.add(mean_norm_signal, norm_signal)
    #             mean_norm_signal_f = np.add(mean_norm_signal_f, norm_signal_f)
    #             mean_norm_signal_r = np.add(mean_norm_signal_r, norm_signal_r)
    #
    #             # Update pwm
    #             aux_plus = 1
    #             dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
    #             if (region.final - region.initial) % 2 == 0:
    #                 aux_plus = 0
    #             dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom,
    #                                                                      p1 + aux_plus, p2 + aux_plus)).upper())
    #             if region.orientation == "+":
    #                 for i in range(0, len(dna_seq)):
    #                     pwm_dict[dna_seq[i]][i] += 1
    #             elif region.orientation == "-":
    #                 for i in range(0, len(dna_seq_rev)):
    #                     pwm_dict[dna_seq_rev[i]][i] += 1
    #
    #     mean_norm_signal = mean_norm_signal / num_sites
    #     mean_norm_signal_f = mean_norm_signal_f / num_sites
    #     mean_norm_signal_r = mean_norm_signal_r / num_sites
    #
    #     mean_norm_signal = signal.boyle_norm(mean_norm_signal)
    #     perc = scoreatpercentile(mean_norm_signal, 98)
    #     std = np.std(mean_norm_signal)
    #     mean_norm_signal = signal.hon_norm_atac(mean_norm_signal, perc, std)
    #
    #     mean_norm_signal_f = signal.boyle_norm(mean_norm_signal_f)
    #     perc = scoreatpercentile(mean_norm_signal_f, 98)
    #     std = np.std(mean_norm_signal_f)
    #     mean_norm_signal_f = signal.hon_norm_atac(mean_norm_signal_f, perc, std)
    #
    #     mean_norm_signal_r = signal.boyle_norm(mean_norm_signal_r)
    #     perc = scoreatpercentile(mean_norm_signal_r, 98)
    #     std = np.std(mean_norm_signal_r)
    #     mean_norm_signal_r = signal.hon_norm_atac(mean_norm_signal_r, perc, std)
    #
    #     mean_slope_signal = signal.slope(mean_norm_signal, signal.sg_coefs)
    #     mean_slope_signal_f = signal.slope(mean_norm_signal_f, signal.sg_coefs)
    #     mean_slope_signal_r = signal.slope(mean_norm_signal_r, signal.sg_coefs)
    #
    #     mean_slope_signal = signal.boyle_norm(mean_slope_signal)
    #     perc = scoreatpercentile(mean_slope_signal, 98)
    #     std = np.std(mean_slope_signal)
    #     mean_slope_signal = signal.hon_norm_atac(mean_slope_signal, perc, std)
    #
    #     mean_slope_signal_f = signal.boyle_norm(mean_slope_signal_f)
    #     perc = scoreatpercentile(mean_slope_signal_f, 98)
    #     std = np.std(mean_slope_signal_f)
    #     mean_slope_signal_f = signal.hon_norm_atac(mean_slope_signal_f, perc, std)
    #
    #     mean_slope_signal_r = signal.boyle_norm(mean_slope_signal_r)
    #     perc = scoreatpercentile(mean_slope_signal_r, 98)
    #     std = np.std(mean_slope_signal_r)
    #     mean_slope_signal_r = signal.hon_norm_atac(mean_slope_signal_r, perc, std)
    #
    #     # Output the norm and slope signal
    #     output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
    #     f = open(output_fname, "w")
    #     f.write("\t".join((map(str, mean_norm_signal))) + "\n")
    #     f.write("\t".join((map(str, mean_slope_signal))) + "\n")
    #     f.write("\t".join((map(str, mean_norm_signal_f))) + "\n")
    #     f.write("\t".join((map(str, mean_slope_signal_f))) + "\n")
    #     f.write("\t".join((map(str, mean_norm_signal_r))) + "\n")
    #     f.write("\t".join((map(str, mean_slope_signal_r))) + "\n")
    #     f.close()
    #
    #     # Output PWM and create logo
    #     pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.output_prefix))
    #     pwm_file = open(pwm_fname, "w")
    #     for e in ["A", "C", "G", "T"]:
    #         pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    #     pwm_file.close()
    #
    #     logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix))
    #     pwm = motifs.read(open(pwm_fname), "pfm")
    #     pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
    #                 color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
    #                 show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
    #                 show_fineprint=False, show_ends=False)
    #
    #     fig, (ax1, ax2) = plt.subplots(2)
    #
    #     start = -(self.window_size / 2)
    #     end = (self.window_size / 2) - 1
    #     x = np.linspace(start, end, num=self.window_size)
    #
    #     ax1.plot(x, mean_norm_signal, color='red', label='ATAC-seq')
    #
    #     ax1.xaxis.set_ticks_position('bottom')
    #     ax1.yaxis.set_ticks_position('left')
    #     ax1.spines['top'].set_visible(False)
    #     ax1.spines['right'].set_visible(False)
    #     ax1.spines['left'].set_position(('outward', 15))
    #     ax1.spines['bottom'].set_position(('outward', 5))
    #     ax1.tick_params(direction='out')
    #
    #     ax1.set_xticks([start, 0, end])
    #     ax1.set_xticklabels([str(start), 0, str(end)])
    #     min_signal = min(mean_norm_signal)
    #     max_signal = max(mean_norm_signal)
    #     ax1.set_yticks([min_signal, max_signal])
    #     ax1.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
    #
    #     ax1.set_title(self.output_prefix, fontweight='bold')
    #     ax1.set_xlim(start, end)
    #     ax1.set_ylim([min_signal, max_signal])
    #     ax1.legend(loc="upper right", frameon=False)
    #     ax1.set_ylabel("Average Signal", rotation=90, fontweight='bold')
    #
    #     min_signal = min(min(mean_norm_signal_f), min(mean_norm_signal_r))
    #     max_signal = max(max(mean_norm_signal_f), max(mean_norm_signal_r))
    #     ax2.plot(x, mean_norm_signal_f, color='red', label='Forward')
    #     ax2.plot(x, mean_norm_signal_r, color='green', label='Reverse')
    #
    #     ax2.xaxis.set_ticks_position('bottom')
    #     ax2.yaxis.set_ticks_position('left')
    #     ax2.spines['top'].set_visible(False)
    #     ax2.spines['right'].set_visible(False)
    #     ax2.spines['left'].set_position(('outward', 15))
    #     ax2.tick_params(direction='out')
    #     ax2.set_xticks([start, 0, end])
    #     ax2.set_xticklabels([str(start), 0, str(end)])
    #     ax2.set_yticks([min_signal, max_signal])
    #     ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
    #     ax2.set_xlim(start, end)
    #     ax2.set_ylim([min_signal, max_signal])
    #     ax2.legend(loc="upper right", frameon=False)
    #
    #     ax2.spines['bottom'].set_position(('outward', 40))
    #     ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
    #     ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')
    #
    #     figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
    #     fig.subplots_adjust(bottom=.2, hspace=.5)
    #     fig.tight_layout()
    #     fig.savefig(figure_name, format="eps", dpi=300)
    #
    #     # Creating canvas and printing eps / pdf with merged results
    #     output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
    #     c = pyx.canvas.canvas()
    #     c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    #     c.insert(pyx.epsfile.epsfile(0, 1.55, logo_fname, width=17.5, height=3))
    #     c.writeEPSfile(output_fname)
    #     os.system("epstopdf " + figure_name)
    #     os.system("epstopdf " + logo_fname)
    #     os.system("epstopdf " + output_fname)
    #
    #     os.remove(pwm_fname)
    #     os.remove(os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix)))
    #     os.remove(os.path.join(self.output_loc, "{}.logo.eps".format(self.output_prefix)))
    #     os.remove(os.path.join(self.output_loc, "{}.line.pdf".format(self.output_prefix)))
    #     os.remove(os.path.join(self.output_loc, "{}.logo.pdf".format(self.output_prefix)))
    #     os.remove(os.path.join(self.output_loc, "{}.eps".format(self.output_prefix)))

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
        # for read in bam.fetch(ref, p1_w, p2_w):
        #     if (not read.is_reverse):
        #         for i in range(max(read.pos + self.forward_shift - self.upstream_ext, start),
        #                       min(read.pos + self.forward_shift + self.downstream_ext, end - 1)):
        #             signal_raw_f[i - start] += 1.0
        #
        #     else:
        #         for i in range(max(read.aend + self.reverse_shift - self.downstream_ext, start),
        #                        min(read.aend + self.reverse_shift + self.upstream_ext, end - 1)):
        #             signal_raw_r[i - start] += 1.0

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

        fig, (ax1, ax2) = plt.subplots(2)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        ############################################################
        # bias signal per strand
        min_ = min(min(mean_signal_bias_f), min(mean_signal_bias_r))
        max_ = max(max(mean_signal_bias_f), max(mean_signal_bias_r))
        ax1.plot(x, mean_signal_bias_f, color='red', label='Forward')
        ax1.plot(x, mean_signal_bias_r, color='blue', label='Reverse')
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
        ax1.set_ylabel("Bias Signal", rotation=90, fontweight='bold')
        ####################################################################

        #####################################################################
        # Bias corrected, non-bias corrected (not strand specific)
        min_ = min(min(mean_signal_raw), min(mean_signal_bc))
        max_ = max(max(mean_signal_raw), max(mean_signal_bc))
        ax2.plot(x, mean_signal_raw, color='red', label='Uncorrected')
        ax2.plot(x, mean_signal_bc, color='green', label='Corrected')
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
        ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')
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
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom, p1 + aux_plus, p2 + aux_plus)).upper())
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
        pwm_dict = dict([("A", [0.0] * self.window_size), ("C", [0.0] * self.window_size),
                         ("G", [0.0] * self.window_size), ("T", [0.0] * self.window_size),
                         ("N", [0.0] * self.window_size)])

        signal_raw_f = np.zeros(self.window_size)
        signal_raw_r = np.zeros(self.window_size)
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        num_sites = 0
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                for read in bam.fetch(region.chrom, p1, p2):
                    # if (not read.is_reverse):
                    #     cut_site = read.pos + self.forward_shift
                    #     for i in range(max(cut_site - self.upstream_ext, p1),
                    #                    min(cut_site + self.downstream_ext, p2)):
                    #         signal_raw_f[i - p1] += 1.0
                    #
                    # else:
                    #     cut_site = read.aend + self.reverse_shift - 1
                    #     for i in range(max(cut_site - self.downstream_ext, p1),
                    #                    min(cut_site + self.upstream_ext, p2)):
                    #         signal_raw_r[i - p1] += 1.0
                    if (not read.is_reverse):
                        cut_site = read.pos + self.forward_shift
                        if cut_site >= p1 and cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift
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

        mean_signal_raw_f = signal_raw_f / num_sites
        mean_signal_raw_r = signal_raw_r / num_sites

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
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(self.window_size),
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

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
        ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        ax2.set_ylabel("Average Signal", rotation=90, fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(2, 1.4, logo_fname, width=17.8, height=1.75))
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

    def rescaling(self, vector):
        maxN = max(vector)
        minN = min(vector)
        return [(e - minN) / (maxN - minN) for e in vector]