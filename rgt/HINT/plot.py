###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import numpy as np
from pysam import Fastafile
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
    def __init__(self, bam_file, motif_file, motif_name, window_size,
                 atac_downstream_ext, atac_upstream_ext, atac_forward_shift, atac_reverse_shift,
                 initial_clip, organism, bias_table, k_nb, protection_score,
                 strands_specific, output_loc):
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
        self.strands_specific = strands_specific
        self.protection_score = protection_score
        self.output_loc = output_loc

    def line(self):
        signal = GenomicSignal(self.bam_file)
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
        mean_raw_signal_f = np.zeros(self.window_size)
        mean_bc_signal_f = np.zeros(self.window_size)
        mean_raw_signal_r = np.zeros(self.window_size)
        mean_bc_signal_r = np.zeros(self.window_size)

        mean_bias_signal_f = np.zeros(self.window_size)
        mean_bias_signal_r = np.zeros(self.window_size)
        num_sites = 0

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read_bed(self.motif_file)

        total_nc_signal = 0
        total_nl_signal = 0
        total_nr_signal = 0

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                num_sites += 1
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                if not self.strands_specific:
                    # Fetch raw signal
                    raw_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                  downstream_ext=self.atac_downstream_ext,
                                                  upstream_ext=self.atac_upstream_ext,
                                                  forward_shift=self.atac_forward_shift,
                                                  reverse_shift=self.atac_reverse_shift,
                                                  genome_file_name=genome_data.get_genome())

                    mean_raw_signal = np.add(mean_raw_signal, raw_signal)

                    # Fetch bias correction signal
                    bc_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                     bias_table=table,
                                                     downstream_ext=self.atac_downstream_ext,
                                                     upstream_ext=self.atac_upstream_ext,
                                                     forward_shift=self.atac_forward_shift,
                                                     reverse_shift=self.atac_reverse_shift,
                                                     genome_file_name=genome_data.get_genome())

                    mean_bc_signal = np.add(mean_bc_signal, bc_signal)
                else:
                    raw_signal_f, _, raw_signal_r, _ =  signal.get_signal_per_strand(ref=region.chrom, start=p1, end=p2,
                                                                        downstream_ext=self.atac_downstream_ext,
                                                                        upstream_ext=self.atac_upstream_ext,
                                                                        forward_shift=self.atac_forward_shift,
                                                                        reverse_shift=self.atac_reverse_shift,
                                                                        genome_file_name=genome_data.get_genome())
                    mean_raw_signal_f = np.add(mean_raw_signal_f, raw_signal_f)
                    mean_raw_signal_r = np.add(mean_raw_signal_r, raw_signal_r)

                    bc_signal_f, _, bc_signal_r, _ = signal.get_signal_per_strand(ref=region.chrom, start=p1, end=p2,
                                                                                  bias_table=table,
                                                                                  downstream_ext=self.atac_downstream_ext,
                                                                                  upstream_ext=self.atac_upstream_ext,
                                                                                  forward_shift=self.atac_forward_shift,
                                                                                  reverse_shift=self.atac_reverse_shift,
                                                                                  genome_file_name=genome_data.get_genome())
                    mean_bc_signal_f = np.add(mean_bc_signal_f, bc_signal_f)
                    mean_bc_signal_r = np.add(mean_bc_signal_r, bc_signal_r)

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

                if self.protection_score:
                    # signal in the center of the MPBS
                    p1 = region.initial
                    p2 = region.final
                    nc_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                     bias_table=table,
                                                     downstream_ext=self.atac_downstream_ext,
                                                     upstream_ext=self.atac_upstream_ext,
                                                     forward_shift=self.atac_forward_shift,
                                                     reverse_shift=self.atac_reverse_shift,
                                                     genome_file_name=genome_data.get_genome())
                    total_nc_signal += sum(nc_signal)
                    p1 = region.final
                    p2 = 2 * region.final - region.initial
                    nr_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                     bias_table=table,
                                                     downstream_ext=self.atac_downstream_ext,
                                                     upstream_ext=self.atac_upstream_ext,
                                                     forward_shift=self.atac_forward_shift,
                                                     reverse_shift=self.atac_reverse_shift,
                                                     genome_file_name=genome_data.get_genome())
                    total_nr_signal += sum(nr_signal)
                    p1 = 2 * region.initial - region.final
                    p2 = region.final
                    nl_signal, _ = signal.get_signal(ref=region.chrom, start=p1, end=p2,
                                                     bias_table=table,
                                                     downstream_ext=self.atac_downstream_ext,
                                                     upstream_ext=self.atac_upstream_ext,
                                                     forward_shift=self.atac_forward_shift,
                                                     reverse_shift=self.atac_reverse_shift,
                                                     genome_file_name=genome_data.get_genome())
                    total_nl_signal += sum(nl_signal)


        mean_raw_signal = mean_raw_signal / num_sites
        mean_bc_signal = mean_bc_signal / num_sites

        mean_raw_signal_f = mean_raw_signal_f / num_sites
        mean_raw_signal_r = mean_raw_signal_r / num_sites
        mean_bc_signal_f = mean_bc_signal_f / num_sites
        mean_bc_signal_r = mean_bc_signal_r / num_sites

        mean_bias_signal_f = mean_bias_signal_f / num_sites
        mean_bias_signal_r = mean_bias_signal_r / num_sites

        protection_score = (total_nl_signal + total_nr_signal - 2 * total_nc_signal) / (2 * num_sites)

        # Output PWM and create logo
        pwm_fname = os.path.join(self.output_loc, "{}.pwm".format(self.motif_name))
        pwm_file = open(pwm_fname,"w")
        for e in ["A","C","G","T"]:
            pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]])+"\n")
        pwm_file.close()

        logo_fname = os.path.join(self.output_loc, "{}.logo.eps".format(self.motif_name))
        pwm = motifs.read(open(pwm_fname), "pfm")
        pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line="100",
                    color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                    show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                    show_fineprint=False, show_ends=False)

        # Output the raw, bias corrected signal and protection score
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.motif_name))
        output_file = open(output_fname, "w")
        if not self.strands_specific:
            output_file.write("raw signal: \n" + np.array_str(mean_raw_signal) + "\n")
            output_file.write("bias corrected signal: \n" + np.array_str(mean_bc_signal) + "\n")
        else:
            output_file.write("raw forward signal: \n" + np.array_str(mean_raw_signal_f) + "\n")
            output_file.write("bias corrected forward signal: \n" + np.array_str(mean_bc_signal_f) + "\n")
            output_file.write("raw reverse signal: \n" + np.array_str(mean_raw_signal_r) + "\n")
            output_file.write("bias reverse corrected signal: \n" + np.array_str(mean_bc_signal_r) + "\n")
        output_file.write("forward bias signal: \n" + np.array_str(mean_bias_signal_f) + "\n")
        output_file.write("reverse bias signal: \n" + np.array_str(mean_bias_signal_r) + "\n")
        if self.protection_score:
            output_file.write("protection score: \n" + str(protection_score) + "\n")
        output_file.close()

        if self.strands_specific:
            fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(12.0, 10.0))
        else:
            fig, (ax1, ax2) = plt.subplots(2)
        x = np.linspace(-50, 49, num=self.window_size)

        ax1.plot(x, mean_bias_signal_f, color='red', label='Forward')
        ax1.plot(x, mean_bias_signal_r, color='blue', label='Reverse')

        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.spines['bottom'].set_position(('outward', 5))
        ax1.tick_params(direction='out')

        ax1.set_xticks([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 49])
        ax1.set_xticklabels(['-50', '-40', '-30', '-20', '-10', '0', '10', '20', '30', '40', '49'])
        min_bias_signal = min(min(mean_bias_signal_f), min(mean_bias_signal_r))
        max_bias_signal = max(max(mean_bias_signal_f), max(mean_bias_signal_r))
        ax1.set_yticks([min_bias_signal, max_bias_signal])
        ax1.set_yticklabels([str(round(min_bias_signal,2)), str(round(max_bias_signal,2))], rotation=90)

        ax1.text(-48, max_bias_signal, '# Sites = {}'.format(str(num_sites)), fontweight='bold')
        ax1.set_title(self.motif_name, fontweight='bold')
        ax1.set_xlim(-50, 49)
        ax1.set_ylim([min_bias_signal, max_bias_signal])
        ax1.legend(loc="upper right", frameon=False)
        ax1.set_ylabel("Average Bias \nSignal", rotation=90, fontweight='bold')

        if not self.strands_specific:
            mean_raw_signal = self.standardize(mean_raw_signal)
            mean_bc_signal = self.standardize(mean_bc_signal)
            ax2.plot(x, mean_raw_signal, color='red', label='Uncorrected')
            ax2.plot(x, mean_bc_signal, color='green', label='Corrected')
        else:
            mean_raw_signal_f = self.standardize(mean_raw_signal_f)
            mean_raw_signal_r = self.standardize(mean_raw_signal_r)
            mean_bc_signal_f = self.standardize(mean_bc_signal_f)
            mean_bc_signal_r = self.standardize(mean_bc_signal_r)
            ax2.plot(x, mean_raw_signal_f, color='red', label='Forward')
            ax2.plot(x, mean_raw_signal_r, color='green', label='Reverse')
            ax3.plot(x, mean_bc_signal_f, color='red', label='Forward')
            ax3.plot(x, mean_bc_signal_r, color='green', label='Reverse')

        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_xticks([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 49])
        ax2.set_xticklabels(['-50', '-40', '-30', '-20', '-10', '0', '10', '20', '30', '40', '49'])
        ax2.set_yticks([0, 1])
        ax2.set_yticklabels([str(0), str(1)], rotation=90)
        ax2.set_xlim(-50, 49)
        ax2.set_ylim([0, 1])

        if not self.strands_specific:
            ax2.spines['bottom'].set_position(('outward', 40))
            ax2.set_xlabel("Coordinates from Motif Center", fontweight='bold')
            ax2.set_ylabel("Average ATAC-seq \nSignal", rotation=90, fontweight='bold')
            ax2.legend(loc="center", frameon=False, bbox_to_anchor=(0.85, 0.06))
        else:
            ax2.spines['bottom'].set_position(('outward', 5))
            ax2.set_ylabel("Average ATAC-seq \n Uncorrected Signal", rotation=90, fontweight='bold')
            ax2.legend(loc="lower right", frameon=False)

            ax3.xaxis.set_ticks_position('bottom')
            ax3.yaxis.set_ticks_position('left')
            ax3.spines['top'].set_visible(False)
            ax3.spines['right'].set_visible(False)
            ax3.spines['left'].set_position(('outward', 15))
            ax3.tick_params(direction='out')
            ax3.set_xticks([-50, -40, -30, -20, -10, 0, 10, 20, 30, 40, 49])
            ax3.set_xticklabels(['-50', '-40', '-30', '-20', '-10', '0', '10', '20', '30', '40', '49'])
            ax3.set_yticks([0, 1])
            ax3.set_yticklabels([str(0), str(1)], rotation=90)
            ax3.set_xlim(-50, 49)
            ax3.set_ylim([0, 1])
            ax3.legend(loc="lower right", frameon=False)
            ax3.spines['bottom'].set_position(('outward', 40))
            ax3.set_xlabel("Coordinates from Motif Center", fontweight='bold')
            ax3.set_ylabel("Average ATAC-seq \n Corrected Signal", rotation=90, fontweight='bold')
            ax3.text(-48, 0.05, '# K-mer = {}\n# Forward Shift = {}'.format(str(self.k_nb), str(self.atac_forward_shift)),
                     fontweight='bold')

        figure_name = os.path.join(self.output_loc, "{}.line.eps".format(self.motif_name))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="eps", dpi=300)

        # Creating canvas and printing eps / pdf with merged results
        output_fname = os.path.join(self.output_loc, "{}.eps".format(self.motif_name))
        c = pyx.canvas.canvas()
        c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
        if self.strands_specific:
            c.insert(pyx.epsfile.epsfile(2.76, 1.58, logo_fname, width=27.2, height=2.45))
        else:
            c.insert(pyx.epsfile.epsfile(2.5, 1.54, logo_fname, width=16, height=1.75))
        c.writeEPSfile(output_fname)
        os.system("epstopdf " + figure_name)
        os.system("epstopdf " + logo_fname)
        os.system("epstopdf " + output_fname)

    def standardize(self, vector):
        maxN = max(vector)
        minN = min(vector)
        return [(e - minN) / (maxN - minN) for e in vector]