###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import numpy as np
import math
from pysam import Samfile
from pysam import Fastafile
from Bio import motifs
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Internal
from ..Util import GenomeData
from signalProcessing import GenomicSignal
from rgt.GenomicRegionSet import GenomicRegionSet
from biasTable import BiasTable
from ..Util import AuxiliaryFunctions



class Plot:
    """

    """
    def __init__(self, bam_file, motif_file, motif_name, window_size,
                atac_downstream_ext, atac_upstream_ext, atac_forward_shift, atac_reverse_shift,
                initial_clip, organism, bias_table, k_nb, output_loc):
        self.bam_file = bam_file
        self.motif_file = motif_file
        self.motif_name = motif_name
        self.window_size = window_size
        self.atac_downstream_ext = atac_downstream_ext
        self.atac_upstream_ext = atac_upstream_ext
        self.atac_forward_shift = atac_forward_shift
        self.atac_reverse_shift = atac_reverse_shift
        self.initial_clip = initial_clip
        self.organism = organism
        self.bias_table = bias_table
        self.k_nb = k_nb
        self.output_loc = output_loc

    def line(self):
        signal = GenomicSignal(self.bam_file)
        signal.load_sg_coefs(slope_window_size=9)
        bias_table = BiasTable()
        bias_table_list = self.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])
        genome_data = GenomeData(self.organism)
        bam = Samfile(self.bam_file, "rb")
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
        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                num_sites += 1
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                raw_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                  downstream_ext=self.atac_downstream_ext,
                                                  upstream_ext=self.atac_upstream_ext,
                                                  forward_shift=self.atac_forward_shift,
                                                  reverse_shift=self.atac_reverse_shift,
                                                  genome_file_name=genome_data.get_genome(),)

                mean_raw_signal = np.add(mean_raw_signal, raw_signal)

                # Fetch bias correction signal
                bc_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                 bias_table=table,
                                                 downstream_ext=self.atac_downstream_ext,
                                                 upstream_ext=self.atac_upstream_ext,
                                                 forward_shift=self.atac_forward_shift,
                                                 reverse_shift=self.atac_reverse_shift,
                                                 genome_file_name=genome_data.get_genome(),)

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
                bias_signal_f = []
                bias_signal_r = []
                p1_wk = p1 - int(math.floor(self.k_nb / 2.))
                p2_wk = p2 + int(math.ceil(self.k_nb / 2.))
                dna_seq = str(fasta.fetch(region.chrom, p1_wk, p2_wk)).upper()
                dna_seq_rev = AuxiliaryFunctions.revcomp(str(fasta.fetch(region.chrom, p1_wk + 1, p2_wk)).upper())

                for i in range(int(math.ceil(self.k_nb / 2.)), len(dna_seq) - int(math.floor(self.k_nb / 2))):
                    fseq = dna_seq[i - int(math.floor(self.k_nb / 2.)):i + int(math.ceil(self.k_nb / 2.))]
                    rseq = dna_seq_rev[len(dna_seq) - int(math.ceil(self.k_nb / 2.))
                                       - i:len(dna_seq) + int(math.floor(self.k_nb / 2.)) - i]
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
        mean_bias_signal_f = mean_bias_signal_f / num_sites
        mean_bias_signal_r = mean_bias_signal_r / num_sites

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.motif_name))
        pwm_file = open(pwm_fname,"w")
        for e in ["A","C","G","T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]])+"\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.eps".format(self.motif_name))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line="50",
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        # Output the raw and bias corrected signal
        signal_fname = os.path.join(self.output_loc, "{}.signal.txt".format(self.motif_name))
        signal_file = open(signal_fname, "w")
        signal_file.write("raw signal: " + np.array_str(mean_raw_signal) + "\n")
        signal_file.write("bias corrected signal: " + np.array_str(mean_bc_signal) + "\n")
        signal_file.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.grid(True, which='both')
        x = np.linspace(0, self.window_size, num=self.window_size)
        ax.plot(x, mean_raw_signal, color='red')
        ax.plot(x, mean_bc_signal, color='green')
        figure_name = os.path.join(self.output_loc, "{}.signal.pdf".format(self.motif_name))
        fig.savefig(figure_name, format="pdf", dpi=300)

        # Output the forward and reverse bias signal
        bias_fname = os.path.join(self.output_loc, "{}.bias.txt".format(self.motif_name))
        bias_file = open(bias_fname, "w")
        bias_file.write("forward bias signal: " + np.array_str(mean_bias_signal_f) + "\n")
        bias_file.write("reverse bias signal: " + np.array_str(mean_bias_signal_r) + "\n")
        bias_file.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.xaxis.set_ticklabels([])
        ax.yaxis.set_ticklabels([])
        ax.grid(True, which='both')
        x = np.linspace(0, self.window_size, num=self.window_size)
        ax.plot(x, mean_bias_signal_f, color='red')
        ax.plot(x, mean_bias_signal_r, color='blue')
        figure_name = os.path.join(self.output_loc, "{}.bias.pdf".format(self.motif_name))
        fig.savefig(figure_name, format="pdf", dpi=300)






