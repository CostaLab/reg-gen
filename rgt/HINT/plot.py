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

# Internal
from ..Util import GenomeData
from signalProcessing import GenomicSignal
from rgt.GenomicRegionSet import GenomicRegionSet
from biasTable import BiasTable
from ..Util import AuxiliaryFunctions


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

    def line(self):
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

        mean_raw_signal = np.zeros(self.window_size)
        mean_bc_signal = np.zeros(self.window_size)

        mean_bias_signal_f = np.zeros(self.window_size)
        mean_bias_signal_r = np.zeros(self.window_size)
        num_sites = 0

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read_bed(self.motif_file)

        bam = Samfile(self.reads_file, "rb")
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                num_sites += 1
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                raw_signal, bc_signal = self.get_signal(ref=region.chrom, start=p1, end=p2,
                                                        bam=bam, fasta=fasta, bias_table=table)

                mean_raw_signal = np.add(mean_raw_signal, raw_signal)
                mean_bc_signal = np.add(mean_bc_signal, bc_signal)

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

                # Create bias signal
                bias_table_f = table[0]
                bias_table_r = table[1]
                self.k_nb = len(bias_table_f.keys()[0])
                bias_signal_f = []
                bias_signal_r = []
                p1_wk = p1 - int(self.k_nb / 2)
                p2_wk = p2 + int(self.k_nb / 2)
                dna_seq = str(fasta.fetch(region.chrom, p1_wk, p2_wk - 1)).upper()
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom, p1_wk, p2_wk + 1)).upper())
                for i in range(int(self.k_nb / 2), len(dna_seq) - int(self.k_nb / 2) + 1):
                    fseq = dna_seq[i - int(self.k_nb / 2):i + int(self.k_nb / 2)]
                    rseq = dna_seq_rev[len(dna_seq) - int(self.k_nb / 2) - i:len(dna_seq) + int(self.k_nb / 2) - i]
                    try:
                        bias_signal_f.append(bias_table_f[fseq])
                    except Exception:
                        bias_signal_f.append(1)
                    try:
                        bias_signal_r.append(bias_table_r[rseq])
                    except Exception:
                        bias_signal_r.append(1)

                mean_bias_signal_f = np.add(mean_bias_signal_f, np.array(bias_signal_f))
                mean_bias_signal_r = np.add(mean_bias_signal_r, np.array(bias_signal_r))

        mean_raw_signal = mean_raw_signal / num_sites
        mean_bc_signal = mean_bc_signal / num_sites

        mean_raw_signal = self.rescaling(mean_raw_signal)
        mean_bc_signal = self.rescaling(mean_bc_signal)

        mean_bias_signal_f = mean_bias_signal_f / num_sites
        mean_bias_signal_r = mean_bias_signal_r / num_sites

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

        # Output the raw, bias corrected signal and protection score
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        output_file = open(output_fname, "w")

        output_file.write("raw signal: \n" + np.array_str(np.array(mean_raw_signal)) + "\n")
        output_file.write("bias corrected signal: \n" + np.array_str(np.array(mean_bc_signal)) + "\n")

        output_file.write("forward bias signal: \n" + np.array_str(mean_bias_signal_f) + "\n")
        output_file.write("reverse bias signal: \n" + np.array_str(mean_bias_signal_r) + "\n")
        output_file.close()

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1

        fig, (ax1, ax2) = plt.subplots(2)
        x = np.linspace(start, end, num=self.window_size)

        ax1.plot(x, mean_bias_signal_f, color='red', label='Forward')
        ax1.plot(x, mean_bias_signal_r, color='blue', label='Reverse')

        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.spines['bottom'].set_position(('outward', 5))
        ax1.tick_params(direction='out')

        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        min_bias_signal = min(min(mean_bias_signal_f), min(mean_bias_signal_r))
        max_bias_signal = max(max(mean_bias_signal_f), max(mean_bias_signal_r))
        ax1.set_yticks([min_bias_signal, max_bias_signal])
        ax1.set_yticklabels([str(round(min_bias_signal, 2)), str(round(max_bias_signal, 2))], rotation=90)

        ax1.text(start + 2, max_bias_signal, '# Sites = {}'.format(str(num_sites)), fontweight='bold')
        ax1.set_title(self.output_prefix, fontweight='bold')
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_bias_signal, max_bias_signal])
        ax1.legend(loc="upper right", frameon=False)
        ax1.set_ylabel("Bias Signal", rotation=90, fontweight='bold')

        ax2.plot(x, mean_raw_signal, color='red', label='Uncorrected')
        ax2.plot(x, mean_bc_signal, color='green', label='Corrected')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_xticks([start, 0, end])
        ax2.set_xticklabels([str(start), 0, str(end)])
        min_signal = min(min(mean_raw_signal), min(mean_bc_signal))
        max_signal = max(max(mean_raw_signal), max(mean_bc_signal))
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])

        ax2.spines['bottom'].set_position(('outward', 40))
        ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
        ax2.set_ylabel("ATAC-seq Signal", rotation=90, fontweight='bold')
        ax2.legend(loc="center", frameon=False, bbox_to_anchor=(0.85, 0.06))

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.output_prefix))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        c.insert(pyx.epsfile.epsfile(2.68, 1.35, logo_fname, width=17, height=2.45))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

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
        mpbs_regions.read_bed(self.motif_file)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                norm_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                   bias_table=table,
                                                   downstream_ext=self.downstream_ext,
                                                   upstream_ext=self.upstream_ext,
                                                   forward_shift=self.forward_shift,
                                                   reverse_shift=self.reverse_shift,
                                                   genome_file_name=genome_data.get_genome())
                norm_signal_f, _, norm_signal_r, _ = signal.get_signal1(ref=region.chrom, start=p1, end=p2,
                                                                        bias_table=table,
                                                                        downstream_ext=self.downstream_ext,
                                                                        upstream_ext=self.upstream_ext,
                                                                        forward_shift=self.forward_shift,
                                                                        reverse_shift=self.reverse_shift,
                                                                        genome_file_name=genome_data.get_genome())

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

        mean_norm_signal = self.rescaling(mean_norm_signal)
        mean_norm_signal_f = self.rescaling(mean_norm_signal_f)
        mean_norm_signal_r = self.rescaling(mean_norm_signal_r)

        mean_slope_signal = signal.slope(mean_norm_signal, signal.sg_coefs)
        mean_slope_signal_f = signal.slope(mean_norm_signal_f, signal.sg_coefs)
        mean_slope_signal_r = signal.slope(mean_norm_signal_r, signal.sg_coefs)

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

        fig, (ax1, ax2) = plt.subplots(2)

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        ax1.plot(x, mean_norm_signal, color='red', label='ATAC-seq')

        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.spines['bottom'].set_position(('outward', 5))
        ax1.tick_params(direction='out')

        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        min_signal = min(mean_norm_signal)
        max_signal = max(mean_norm_signal)
        ax1.set_yticks([min_signal, max_signal])
        ax1.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)

        ax1.set_title(self.output_prefix, fontweight='bold')
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_signal, max_signal])
        ax1.legend(loc="upper right", frameon=False)
        ax1.set_ylabel("Average Signal", rotation=90, fontweight='bold')

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
        c.insert(pyx.epsfile.epsfile(2.10, 1.55, logo_fname, width=17.5, height=3))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

    def get_signal(self, ref, start, end, bam, fasta, bias_table):
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

        # Raw counts
        nf = [0.0] * (p2_w - p1_w)
        nr = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + self.forward_shift
                if cut_site >= p1_w and cut_site < p2_w:
                    nf[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if cut_site >= p1_w and cut_site < p2_w:
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

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create signal
        af = []
        ar = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
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
        raw_signal = []
        for i in range((window / 2), len(af) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (ar[i] / rSum)
            raw_signal.append(nf[i] + nr[i])
            bc_signal.append(nhatf + nhatr)
            fSum -= fLast
            fSum += af[i + (window / 2)]
            fLast = af[i - (window / 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window / 2)]
            rLast = ar[i - (window / 2) + 1]

        return np.array(raw_signal), np.array(bc_signal)

    def rescaling(self, vector):
        maxN = max(vector)
        minN = min(vector)
        return [(e - minN) / (maxN - minN) for e in vector]
