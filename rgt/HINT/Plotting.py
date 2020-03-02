

###################################################################################################
# Libraries
###################################################################################################
import os
import numpy as np
from pysam import Samfile, Fastafile
from Bio import motifs
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.signal import savgol_filter
from scipy.stats import scoreatpercentile
from argparse import SUPPRESS
import pyx

# Internal
from rgt.Util import GenomeData, AuxiliaryFunctions
from rgt.HINT.signalProcessing import GenomicSignal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.HINT.biasTable import BiasTable


def plotting_args(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help=("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file."))
    parser.add_argument("--reads-file", type=str, metavar="FILE", default=None)
    parser.add_argument("--region-file", type=str, metavar="FILE", default=None)
    parser.add_argument("--reads-file1", type=str, metavar="FILE", default=None)
    parser.add_argument("--reads-file2", type=str, metavar="FILE", default=None)
    parser.add_argument("--motif-file", type=str, metavar="FILE", default=None)
    parser.add_argument("--bias-table", type=str, metavar="FILE1_F,FILE1_R", default=None)
    parser.add_argument("--bias-table1", type=str, metavar="FILE1_F,FILE1_R", default=None)
    parser.add_argument("--bias-table2", type=str, metavar="FILE1_F,FILE1_R", default=None)
    parser.add_argument("--window-size", type=int, metavar="INT", default=400)

    # Hidden Options
    parser.add_argument("--initial-clip", type=int, metavar="INT", default=50, help=SUPPRESS)
    parser.add_argument("--downstream-ext", type=int, metavar="INT", default=1, help=SUPPRESS)
    parser.add_argument("--upstream-ext", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--forward-shift", type=int, metavar="INT", default=5, help=SUPPRESS)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=-5, help=SUPPRESS)
    parser.add_argument("--k-nb", type=int, metavar="INT", default=6, help=SUPPRESS)
    parser.add_argument("--y-lim", type=float, metavar="FLOAT", default=0.3, help=SUPPRESS)

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written.")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default=None,
                        help="The prefix for results files.")

    # plot type
    parser.add_argument("--seq-logo", default=False, action='store_true')
    parser.add_argument("--bias-raw-bc-line", default=False, action='store_true')
    parser.add_argument("--raw-bc-line", default=False, action='store_true')
    parser.add_argument("--strand-line", default=False, action='store_true')
    parser.add_argument("--unstrand-line", default=False, action='store_true')
    parser.add_argument("--bias-line", default=False, action='store_true')
    parser.add_argument("--atac-dnase-line", default=False, action='store_true')
    parser.add_argument("--bias-raw-bc-strand-line2", default=False, action='store_true')
    parser.add_argument("--fragment-raw-size-line", default=False, action='store_true')
    parser.add_argument("--fragment-bc-size-line", default=False, action='store_true')


def plotting_run(args):
    if args.seq_logo:
        seq_logo(args)

    if args.bias_raw_bc_line:
        bias_raw_bc_strand_line(args)

    if args.strand_line:
        strand_line(args)

    if args.unstrand_line:
        unstrand_line(args)

    if args.raw_bc_line:
        raw_bc_line(args)

    if args.bias_raw_bc_strand_line2:
        bias_raw_bc_strand_line2(args)

    if args.fragment_raw_size_line:
        fragment_size_raw_line(args)

    if args.fragment_bc_size_line:
        fragment_size_bc_line(args)


def seq_logo(args):
    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm_file = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_dict = dict(
        [("A", [0.0] * args.window_size), ("C", [0.0] * args.window_size), ("G", [0.0] * args.window_size),
         ("T", [0.0] * args.window_size), ("N", [0.0] * args.window_size)])

    genome_data = GenomeData(args.organism)
    fasta_file = Fastafile(genome_data.get_genome())
    bam = Samfile(args.reads_file, "rb")
    regions = GenomicRegionSet("Peaks")
    regions.read(args.region_file)

    for region in regions:
        for r in bam.fetch(region.chrom, region.initial, region.final):
            if not r.is_reverse:
                cut_site = r.pos - 1
                p1 = cut_site - int(args.window_size / 2)
            else:
                cut_site = r.aend + 1
                p1 = cut_site - int(args.window_size / 2)
            p2 = p1 + args.window_size

            # Fetching k-mer
            currStr = str(fasta_file.fetch(region.chrom, p1, p2)).upper()
            if r.is_reverse: continue
            for i in range(0, len(currStr)):
                pwm_dict[currStr[i]][i] += 1

    with open(pwm_file, "w") as f:
        for e in ["A", "C", "G", "T"]:
            f.write(" ".join([str(int(c)) for c in pwm_dict[e]]) + "\n")

    pwm = motifs.read(open(pwm_file), "pfm")
    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False, yaxis_scale=args.y_lim)

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size).tolist()

    fig = plt.figure(figsize=(8, 2))
    ax = fig.add_subplot(111)

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')

    ax.xaxis.set_ticks(list(map(int, x)))
    x1 = list(map(int, x))
    ax.set_xticklabels(list(map(str, x1)), rotation=90)
    ax.set_xlabel("Coordinates from Read Start", fontweight='bold')

    ax.set_ylim([0, args.y_lim])
    ax.yaxis.set_ticks([0, args.y_lim])
    ax.set_yticklabels([str(0), str(args.y_lim)], rotation=90)
    ax.set_ylabel("bits", rotation=90)

    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.5, 1.5, logo_fname, width=18.8, height=3.5))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + output_fname)

    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def bias_raw_bc_line(args):
    signal = GenomicSignal(args.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    bias_table = BiasTable()
    bias_table_list = args.bias_table.split(",")
    table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                  table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())
    pwm_dict = dict([("A", [0.0] * args.window_size), ("C", [0.0] * args.window_size),
                     ("G", [0.0] * args.window_size), ("T", [0.0] * args.window_size),
                     ("N", [0.0] * args.window_size)])

    num_sites = 0

    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)

    bam = Samfile(args.reads_file, "rb")

    mean_signal_bias_f = np.zeros(args.window_size)
    mean_signal_bias_r = np.zeros(args.window_size)
    mean_signal_raw = np.zeros(args.window_size)
    mean_signal_raw_f = np.zeros(args.window_size)
    mean_signal_raw_r = np.zeros(args.window_size)
    mean_signal_bc = np.zeros(args.window_size)
    mean_signal_bc_f = np.zeros(args.window_size)
    mean_signal_bc_r = np.zeros(args.window_size)

    motif_len = 0

    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            # Extend by window_size
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)
            motif_len = region.final - region.initial

            signal_bias_f, signal_bias_r, raw, raw_f, raw_r, bc, bc_f, bc_r = \
                signal.get_bias_raw_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam,
                                              fasta=fasta, bias_table=table,
                                              forward_shift=args.forward_shift,
                                              reverse_shift=args.reverse_shift,
                                              strand=True)

            num_sites += 1
            mean_signal_bias_f = np.add(mean_signal_bias_f, np.array(signal_bias_f))
            mean_signal_bias_r = np.add(mean_signal_bias_r, np.array(signal_bias_r))
            mean_signal_raw = np.add(mean_signal_raw, np.array(raw))
            mean_signal_raw_f = np.add(mean_signal_raw_f, np.array(raw_f))
            mean_signal_raw_r = np.add(mean_signal_raw_r, np.array(raw_r))
            mean_signal_bc = np.add(mean_signal_bc, np.array(bc))
            mean_signal_bc_f = np.add(mean_signal_bc_f, np.array(bc_f))
            mean_signal_bc_r = np.add(mean_signal_bc_r, np.array(bc_r))

            # Update pwm
            aux_plus = 1
            dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
            if (region.final - region.initial) % 2 == 0:
                aux_plus = 0

            if region.orientation == "+":
                for i in range(len(dna_seq)):
                    pwm_dict[dna_seq[i]][i] += 1

    mean_signal_bias_f = mean_signal_bias_f / num_sites
    mean_signal_bias_r = mean_signal_bias_r / num_sites
    mean_signal_raw = mean_signal_raw / num_sites
    mean_signal_bc = mean_signal_bc / num_sites
    mean_signal_bc_f = mean_signal_bc_f / num_sites
    mean_signal_bc_r = mean_signal_bc_r / num_sites

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, mean_signal_bias_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bias_r)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_raw)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc)))) + "\n")
    f.close()

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")

    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 6))

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    if motif_len % 2 == 0:
        x1 = int(- (motif_len / 2))
        x2 = int(motif_len / 2)
    else:
        x1 = int(-(motif_len / 2) - 1)
        x2 = int((motif_len / 2) + 1)

    ############################################################
    # bias signal per strand
    fp_score = sum(mean_signal_raw[args.window_size / 2 + x1: args.window_size / 2 + x2])
    shoulder_l = sum(mean_signal_raw[args.window_size / 2 + x1 - motif_len:args.window_size / 2 + x1])
    shoulder_r = sum(mean_signal_raw[args.window_size / 2 + x2:args.window_size / 2 + x2 + motif_len])
    sfr = (shoulder_l + shoulder_r) / (2 * fp_score)
    min_ax1 = min(mean_signal_raw)
    max_ax1 = max(mean_signal_raw)
    ax1.plot(x, mean_signal_raw, color='blue', label='Uncorrected')
    ax1.text(0.15, 0.9, 'n = {}'.format(num_sites), verticalalignment='bottom',
             horizontalalignment='right', transform=ax1.transAxes, fontweight='bold')
    ax1.text(0.35, 0.15, 'SFR = {}'.format(round(sfr, 2)), verticalalignment='bottom',
             horizontalalignment='right', transform=ax1.transAxes, fontweight='bold')
    ax1.xaxis.set_ticks_position('bottom')
    ax1.yaxis.set_ticks_position('left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.spines['left'].set_position(('outward', 15))
    ax1.spines['bottom'].set_position(('outward', 5))
    ax1.tick_params(direction='out')
    ax1.set_xticks([start, 0, end])
    ax1.set_xticklabels([str(start), 0, str(end)])
    ax1.set_yticks([min_ax1, max_ax1])
    ax1.set_yticklabels([str(round(min_ax1, 2)), str(round(max_ax1, 2))], rotation=90)
    ax1.set_title(args.output_prefix, fontweight='bold')
    ax1.set_xlim(start, end)
    ax1.set_ylim([min_ax1, max_ax1])
    ax1.legend(loc="lower right", frameon=False)
    ####################################################################

    #####################################################################
    # Bias corrected, non-bias corrected (not strand specific)
    fp_score = sum(mean_signal_bc[args.window_size / 2 + x1: args.window_size / 2 + x2])
    shoulder_l = sum(mean_signal_bc[args.window_size / 2 + x1 - motif_len:args.window_size / 2 + x1])
    shoulder_r = sum(mean_signal_bc[args.window_size / 2 + x2:args.window_size / 2 + x2 + motif_len])
    sfr = (shoulder_l + shoulder_r) / (2 * fp_score)
    min_ax2 = min(mean_signal_bc)
    max_ax2 = max(mean_signal_bc)
    ax2.plot(x, mean_signal_bc, color='red', label='Corrected')
    ax2.text(0.35, 0.15, 'SFR = {}'.format(round(sfr, 2)), verticalalignment='bottom',
             horizontalalignment='right', transform=ax2.transAxes, fontweight='bold')
    ax2.xaxis.set_ticks_position('bottom')
    ax2.yaxis.set_ticks_position('left')
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.spines['left'].set_position(('outward', 15))
    ax2.tick_params(direction='out')
    ax2.set_xticks([start, 0, end])
    ax2.set_xticklabels([str(start), 0, str(end)])
    ax2.set_yticks([min_ax2, max_ax2])
    ax2.set_yticklabels([str(round(min_ax2, 2)), str(round(max_ax2, 2))], rotation=90)
    ax2.set_xlim(start, end)
    ax2.set_ylim([min_ax2, max_ax2])
    ax2.legend(loc="lower right", frameon=False)

    fp_score_f = sum(mean_signal_bc_f[args.window_size / 2 + x1: args.window_size / 2 + x2])
    shoulder_l_f = sum(mean_signal_bc_f[args.window_size / 2 + x1 - motif_len:args.window_size / 2 + x1])
    shoulder_r_f = sum(mean_signal_bc_f[args.window_size / 2 + x2:args.window_size / 2 + x2 + motif_len])
    sfr_f = (shoulder_l_f + shoulder_r_f) / (2 * fp_score_f)
    fp_score_r = sum(mean_signal_bc_r[args.window_size / 2 + x1: args.window_size / 2 + x2])
    shoulder_l_r = sum(mean_signal_bc_r[args.window_size / 2 + x1 - motif_len:args.window_size / 2 + x1])
    shoulder_r_r = sum(mean_signal_bc_r[args.window_size / 2 + x2:args.window_size / 2 + x2 + motif_len])
    sfr_r = (shoulder_l_r + shoulder_r_r) / (2 * fp_score_r)
    min_ax3 = min(min(mean_signal_bc_f), min(mean_signal_bc_r))
    max_ax3 = max(max(mean_signal_bc_f), max(mean_signal_bc_r))
    ax3.plot(x, mean_signal_bc_f, color='purple', label='Forward')
    ax3.plot(x, mean_signal_bc_r, color='green', label='Reverse')
    ax3.text(0.35, 0.15, 'SFR_f = {}'.format(round(sfr_f, 2)), verticalalignment='bottom',
             horizontalalignment='right', transform=ax3.transAxes, fontweight='bold')
    ax3.text(0.35, 0.05, 'SFR_r = {}'.format(round(sfr_r, 2)), verticalalignment='bottom',
             horizontalalignment='right', transform=ax3.transAxes, fontweight='bold')
    ax3.xaxis.set_ticks_position('bottom')
    ax3.yaxis.set_ticks_position('left')
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.spines['left'].set_position(('outward', 15))
    ax3.tick_params(direction='out')
    ax3.set_xticks([start, 0, end])
    ax3.set_xticklabels([str(start), 0, str(end)])
    ax3.set_yticks([min_ax3, max_ax3])
    ax3.set_yticklabels([str(round(min_ax3, 2)), str(round(max_ax3, 2))], rotation=90)
    ax3.set_xlim(start, end)
    ax3.set_ylim([min_ax3, max_ax3])
    ax3.legend(loc="lower right", frameon=False)

    ax3.spines['bottom'].set_position(('outward', 40))

    ax1.axvline(x=x1, ymin=-0.3, ymax=1, c="black", lw=0.5, ls='dashed', zorder=0, clip_on=False)
    ax1.axvline(x=x2, ymin=-0.3, ymax=1, c="black", lw=0.5, ls='dashed', zorder=0, clip_on=False)
    ax2.axvline(x=x1, ymin=-0.5, ymax=1.2, c="black", lw=0.5, ls='dashed', zorder=0, clip_on=False)
    ax2.axvline(x=x2, ymin=-0.5, ymax=1.2, c="black", lw=0.5, ls='dashed', zorder=0, clip_on=False)
    ###############################################################################
    # merge the above figures
    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.45, 0.89, logo_fname, width=18.3, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def raw_bc_line(args):
    signal = GenomicSignal(args.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    bias_table = BiasTable()
    bias_table_list = args.bias_table.split(",")
    table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                  table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())
    pwm_dict = dict([("A", [0.0] * args.window_size), ("C", [0.0] * args.window_size),
                     ("G", [0.0] * args.window_size), ("T", [0.0] * args.window_size),
                     ("N", [0.0] * args.window_size)])

    num_sites = 0

    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)

    bam = Samfile(args.reads_file, "rb")

    mean_signal_raw = np.zeros(args.window_size)
    mean_signal_bc = np.zeros(args.window_size)

    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            # Extend by window_size
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)
            signal_bias_f, signal_bias_r, raw, raw_f, raw_r, bc, bc_f, bc_r = \
                signal.get_bias_raw_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam,
                                              fasta=fasta, bias_table=table,
                                              forward_shift=args.forward_shift,
                                              reverse_shift=args.reverse_shift,
                                              strand=True)

            num_sites += 1
            mean_signal_raw = np.add(mean_signal_raw, np.array(raw))
            mean_signal_bc = np.add(mean_signal_bc, np.array(bc))

            # Update pwm
            aux_plus = 1
            dna_seq = str(fasta.fetch(region.chrom, p1, p2)).upper()
            if (region.final - region.initial) % 2 == 0:
                aux_plus = 0

            if region.orientation == "+":
                for i in range(len(dna_seq)):
                    pwm_dict[dna_seq[i]][i] += 1

    mean_signal_raw = mean_signal_raw / num_sites
    mean_signal_bc = mean_signal_bc / num_sites

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, mean_signal_raw)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc)))) + "\n")
    f.close()

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")

    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    ############################################################
    min_ = min(min(mean_signal_raw), min(mean_signal_bc))
    max_ = max(max(mean_signal_raw), max(mean_signal_bc))
    ax.plot(x, mean_signal_raw, color='red', label='Uncorrected')
    ax.plot(x, mean_signal_bc, color='blue', label='Corrected')
    ax.text(0.15, 0.9, 'n = {}'.format(num_sites), verticalalignment='bottom',
            horizontalalignment='right', transform=ax.transAxes, fontweight='bold')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.tick_params(direction='out')
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    ax.set_yticks([min_, max_])
    ax.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
    ax.set_title(args.output_prefix, fontweight='bold')
    ax.set_xlim(start, end)
    ax.set_ylim(min_, max_)
    ax.legend(loc="lower right", frameon=False)

    ax.spines['bottom'].set_position(('outward', 40))

    ###############################################################################
    # merge the above figures
    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.45, 0.89, logo_fname, width=18.3, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def bias_raw_bc_strand_line(args):
    signal = GenomicSignal(args.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    bias_table = BiasTable()
    bias_table_list = args.bias_table.split(",")
    table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                  table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())
    pwm_dict = dict([("A", [0.0] * args.window_size), ("C", [0.0] * args.window_size),
                     ("G", [0.0] * args.window_size), ("T", [0.0] * args.window_size),
                     ("N", [0.0] * args.window_size)])

    num_sites = 0

    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)

    bam = Samfile(args.reads_file, "rb")

    mean_signal_bias_f = np.zeros(args.window_size)
    mean_signal_bias_r = np.zeros(args.window_size)
    mean_signal_raw = np.zeros(args.window_size)
    mean_signal_bc = np.zeros(args.window_size)
    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            signal_bias_f, signal_bias_r, signal_raw, signal_bc = \
                signal.get_bias_raw_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam,
                                              fasta=fasta, bias_table=table,
                                              forward_shift=args.forward_shift, reverse_shift=args.reverse_shift)

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

    mean_signal_bias_f = mean_signal_bias_f / num_sites
    mean_signal_bias_r = mean_signal_bias_r / num_sites
    mean_signal_raw = mean_signal_raw / num_sites
    mean_signal_bc = mean_signal_bc / num_sites

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, mean_signal_bias_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bias_r)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_raw)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc)))) + "\n")
    f.close()

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")

    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4))

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    ############################################################
    # bias signal per strand
    min_ = min(min(mean_signal_bias_f), min(mean_signal_bias_r))
    max_ = max(max(mean_signal_bias_f), max(mean_signal_bias_r))
    ax1.plot(x, mean_signal_bias_f, color='purple', label='Forward')
    ax1.plot(x, mean_signal_bias_r, color='green', label='Reverse')
    ax1.text(0.15, 0.9, 'n = {}'.format(num_sites), verticalalignment='bottom',
             horizontalalignment='right', transform=ax1.transAxes, fontweight='bold')
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
    ax1.set_title(args.output_prefix, fontweight='bold')
    ax1.set_xlim(start, end)
    ax1.set_ylim([min_, max_])
    ax1.legend(loc="upper right", frameon=False)
    ####################################################################

    #####################################################################
    # Bias corrected, non-bias corrected (not strand specific)
    min_ = min(min(mean_signal_raw), min(mean_signal_bc))
    max_ = max(max(mean_signal_raw), max(mean_signal_bc))
    ax2.plot(x, mean_signal_raw, color='blue', label='Uncorrected')
    ax2.plot(x, mean_signal_bc, color='red', label='Corrected')
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

    ###############################################################################
    # merge the above figures
    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.51, 0.89, logo_fname, width=18.3, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def bias_raw_bc_strand_line2(args):
    signal = GenomicSignal(args.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    bias_table = BiasTable()
    bias_table_list = args.bias_table.split(",")
    table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                  table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())
    pwm_dict = dict([("A", [0.0] * args.window_size), ("C", [0.0] * args.window_size),
                     ("G", [0.0] * args.window_size), ("T", [0.0] * args.window_size),
                     ("N", [0.0] * args.window_size)])

    num_sites = 0

    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)

    bam = Samfile(args.reads_file, "rb")

    mean_signal_bias_f = np.zeros(args.window_size)
    mean_signal_bias_r = np.zeros(args.window_size)
    mean_signal_raw = np.zeros(args.window_size)
    mean_signal_raw_f = np.zeros(args.window_size)
    mean_signal_raw_r = np.zeros(args.window_size)
    mean_signal_bc = np.zeros(args.window_size)
    mean_signal_bc_f = np.zeros(args.window_size)
    mean_signal_bc_r = np.zeros(args.window_size)
    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            signal_bias_f, signal_bias_r, signal_raw, signal_raw_f, signal_raw_r, signal_bc, signal_bc_f, signal_bc_r = \
                signal.get_bias_raw_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam,
                                              fasta=fasta, bias_table=table,
                                              forward_shift=args.forward_shift,
                                              reverse_shift=args.reverse_shift,
                                              strand=True)

            num_sites += 1
            mean_signal_bias_f = np.add(mean_signal_bias_f, np.array(signal_bias_f))
            mean_signal_bias_r = np.add(mean_signal_bias_r, np.array(signal_bias_r))
            mean_signal_raw = np.add(mean_signal_raw, np.array(signal_raw))
            mean_signal_raw_f = np.add(mean_signal_raw_f, np.array(signal_raw_f))
            mean_signal_raw_r = np.add(mean_signal_raw_r, np.array(signal_raw_r))
            mean_signal_bc = np.add(mean_signal_bc, np.array(signal_bc))
            mean_signal_bc_f = np.add(mean_signal_bc_f, np.array(signal_bc_f))
            mean_signal_bc_r = np.add(mean_signal_bc_r, np.array(signal_bc_r))

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

    mean_signal_bias_f = mean_signal_bias_f / num_sites
    mean_signal_bias_r = mean_signal_bias_r / num_sites
    mean_signal_raw = mean_signal_raw / num_sites
    mean_signal_raw_f = mean_signal_raw_f / num_sites
    mean_signal_raw_r = mean_signal_raw_r / num_sites
    mean_signal_bc = mean_signal_bc / num_sites
    mean_signal_bc_f = mean_signal_bc_f / num_sites
    mean_signal_bc_r = mean_signal_bc_r / num_sites

    # mean_signal_raw = rescaling(mean_signal_raw)
    # mean_signal_bc = rescaling(mean_signal_bc)

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, mean_signal_bias_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bias_r)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_raw)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_raw_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_raw_r)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_signal_bc_r)))) + "\n")
    f.close()

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")

    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    # fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, figsize=(8, 8))
    fig, (ax1, ax4) = plt.subplots(2, figsize=(8, 4))
    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    ############################################################
    # bias signal per strand
    min_ = min(min(mean_signal_bias_f), min(mean_signal_bias_r))
    max_ = max(max(mean_signal_bias_f), max(mean_signal_bias_r))
    ax1.plot(x, mean_signal_bias_f, color='purple', label='Forward')
    ax1.plot(x, mean_signal_bias_r, color='green', label='Reverse')
    ax1.text(0.15, 0.9, 'n = {}'.format(num_sites), verticalalignment='bottom',
             horizontalalignment='right', transform=ax1.transAxes, fontweight='bold')
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
    ax1.set_title(args.output_prefix, fontweight='bold')
    ax1.set_xlim(start, end)
    ax1.set_ylim([min_, max_])
    ax1.legend(loc="upper right", frameon=False)
    ####################################################################

    #####################################################################
    # Bias corrected, non-bias corrected (not strand specific)
    # min_ = min(min(mean_signal_raw_f), min(mean_signal_raw_r))
    # max_ = max(max(mean_signal_raw_f), max(mean_signal_raw_r))
    # ax2.plot(x, mean_signal_raw_f, color='red', label='Forward')
    # ax2.plot(x, mean_signal_raw_r, color='green', label='Reverse')
    # ax2.xaxis.set_ticks_position('bottom')
    # ax2.yaxis.set_ticks_position('left')
    # ax2.spines['top'].set_visible(False)
    # ax2.spines['right'].set_visible(False)
    # ax2.spines['left'].set_position(('outward', 15))
    # ax2.tick_params(direction='out')
    # ax2.set_xticks([start, -1, 0, 1, end])
    # ax2.set_xticklabels([str(start), -1, 0,1, str(end)])
    # ax2.set_yticks([min_, max_])
    # ax2.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
    # ax2.set_xlim(start, end)
    # ax2.set_ylim([min_, max_])
    # ax2.legend(loc="upper right", frameon=False)

    #####################################################################
    # Bias corrected and strand specific
    # min_ = min(min(mean_signal_bc_f), min(mean_signal_bc_r))
    # max_ = max(max(mean_signal_bc_f), max(mean_signal_bc_r))
    # ax3.plot(x, mean_signal_bc_f, color='red', label='Forward')
    # ax3.plot(x, mean_signal_bc_r, color='green', label='Reverse')
    # ax3.xaxis.set_ticks_position('bottom')
    # ax3.yaxis.set_ticks_position('left')
    # ax3.spines['top'].set_visible(False)
    # ax3.spines['right'].set_visible(False)
    # ax3.spines['left'].set_position(('outward', 15))
    # ax3.tick_params(direction='out')
    # ax3.set_xticks([start, 0, end])
    # ax3.set_xticklabels([str(start), 0, str(end)])
    # ax3.set_yticks([min_, max_])
    # ax3.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
    # ax3.set_xlim(start, end)
    # ax3.set_ylim([min_, max_])
    # ax3.legend(loc="upper right", frameon=False)

    #####################################################################
    # Bias corrected, non-bias corrected (not strand specific)
    min_ = min(min(mean_signal_raw), min(mean_signal_bc))
    max_ = max(max(mean_signal_raw), max(mean_signal_bc))
    ax4.plot(x, mean_signal_raw, color='blue', label='Uncorrected')
    ax4.plot(x, mean_signal_bc, color='red', label='Corrected')
    ax4.xaxis.set_ticks_position('bottom')
    ax4.yaxis.set_ticks_position('left')
    ax4.spines['top'].set_visible(False)
    ax4.spines['right'].set_visible(False)
    ax4.spines['left'].set_position(('outward', 15))
    ax4.tick_params(direction='out')
    ax4.set_xticks([start, 0, end])
    ax4.set_xticklabels([str(start), 0, str(end)])
    ax4.set_yticks([min_, max_])
    ax4.set_yticklabels([str(round(min_, 2)), str(round(max_, 2))], rotation=90)
    ax4.set_xlim(start, end)
    ax4.set_ylim([min_, max_])
    ax4.legend(loc="upper right", frameon=False)

    ax4.spines['bottom'].set_position(('outward', 40))

    ###############################################################################
    # merge the above figures
    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.45, 0.89, logo_fname, width=18.3, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def strand_line(args):
    genomic_signal = GenomicSignal(args.reads_file)
    genomic_signal.load_sg_coefs(slope_window_size=9)

    table = None
    if args.bias_table is not None:
        bias_table = BiasTable()
        bias_table_list = args.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())

    num_sites = 0
    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)
    bam = Samfile(args.reads_file, "rb")

    mean_signal_f = np.zeros(args.window_size)
    mean_signal_r = np.zeros(args.window_size)

    pwm_dict = None
    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            # Extend by 50 bp
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            if args.bias_table is not None:
                signal_f, signal_r = genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1,
                                                                                     end=p2, bam=bam, fasta=fasta,
                                                                                     bias_table=table,
                                                                                     forward_shift=args.forward_shift,
                                                                                     reverse_shift=args.reverse_shift)
            else:
                signal_f, signal_r = genomic_signal.get_raw_signal_by_fragment_length(ref=region.chrom, start=p1,
                                                                                      end=p2,
                                                                                      bam=bam,
                                                                                      forward_shift=args.forward_shift,
                                                                                      reverse_shift=args.reverse_shift)

            num_sites += 1

            mean_signal_f = np.add(mean_signal_f, signal_f)
            mean_signal_r = np.add(mean_signal_r, signal_r)

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

    mean_norm_signal_f = genomic_signal.boyle_norm(mean_signal_f)
    perc = scoreatpercentile(mean_norm_signal_f, 98)
    std = np.std(mean_norm_signal_f)
    mean_norm_signal_f = genomic_signal.hon_norm_atac(mean_norm_signal_f, perc, std)

    mean_norm_signal_r = genomic_signal.boyle_norm(mean_signal_r)
    perc = scoreatpercentile(mean_norm_signal_r, 98)
    std = np.std(mean_norm_signal_r)
    mean_norm_signal_r = genomic_signal.hon_norm_atac(mean_norm_signal_r, perc, std)

    mean_slope_signal_f = genomic_signal.slope(mean_norm_signal_f, genomic_signal.sg_coefs)
    mean_slope_signal_r = genomic_signal.slope(mean_norm_signal_r, genomic_signal.sg_coefs)

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, mean_norm_signal_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_slope_signal_f)))) + "\n")
    f.write("\t".join((list(map(str, mean_norm_signal_r)))) + "\n")
    f.write("\t".join((list(map(str, mean_slope_signal_r)))) + "\n")
    f.close()

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")
    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)

    min_signal = min(min(mean_signal_f), min(mean_signal_r))
    max_signal = max(max(mean_signal_f), max(mean_signal_r))
    ax.plot(x, mean_signal_f, color='red', label='Forward')
    ax.plot(x, mean_signal_r, color='green', label='Reverse')
    ax.set_title(args.output_prefix, fontweight='bold')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
    ax.set_xlim(start, end)
    ax.set_ylim([min_signal, max_signal])
    ax.legend(loc="upper right", frameon=False)
    ax.spines['bottom'].set_position(('outward', 40))

    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.37, 0.89, logo_fname, width=18.5, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    # os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def unstrand_line(args):
    genomic_signal = GenomicSignal(args.reads_file)
    genomic_signal.load_sg_coefs(slope_window_size=9)

    table = None
    if args.bias_table is not None:
        bias_table = BiasTable()
        bias_table_list = args.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0], table_file_name_R=bias_table_list[1])

    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())

    num_sites = 0
    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)
    bam = Samfile(args.reads_file, "rb")

    mean_signal = np.zeros(args.window_size)

    pwm_dict = None
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))

    with open(output_fname, "w") as output_f:
        for region in mpbs_regions:
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            if args.bias_table is not None:
                signal = genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                         bam=bam, fasta=fasta,
                                                                         bias_table=table,
                                                                         forward_shift=args.forward_shift,
                                                                         reverse_shift=args.reverse_shift,
                                                                         strand=False)
            else:
                signal = genomic_signal.get_raw_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                          bam=bam,
                                                                          forward_shift=args.forward_shift,
                                                                          reverse_shift=args.reverse_shift,
                                                                          strand=False)

            if region.orientation == "-":
                signal = np.flip(signal)

            name = "{}_{}_{}".format(region.chrom, str(region.initial), str(region.final))
            output_f.write(name + "\t" + "\t".join(map(str, list(map(int, signal)))) + "\n")
            num_sites += 1

            mean_signal = np.add(mean_signal, signal)

            # Update pwm
            if pwm_dict is None:
                pwm_dict = dict([("A", [0.0] * (p2 - p1)), ("C", [0.0] * (p2 - p1)),
                                 ("G", [0.0] * (p2 - p1)), ("T", [0.0] * (p2 - p1)),
                                 ("N", [0.0] * (p2 - p1))])

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

    mean_signal = mean_signal / num_sites

    # Output PWM and create logo
    pwm_fname = os.path.join(args.output_location, "{}.pwm".format(args.output_prefix))
    pwm_file = open(pwm_fname, "w")
    for e in ["A", "C", "G", "T"]:
        pwm_file.write(" ".join([str(int(f)) for f in pwm_dict[e]]) + "\n")
    pwm_file.close()

    logo_fname = os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix))
    pwm = motifs.read(open(pwm_fname), "pfm")
    pwm.weblogo(logo_fname, format="eps", stack_width="large", stacks_per_line=str(args.window_size),
                color_scheme="color_classic", unit_name="", show_errorbars=False, logo_title="",
                show_xaxis=False, xaxis_label="", show_yaxis=False, yaxis_label="",
                show_fineprint=False, show_ends=False)

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)

    fig = plt.figure(figsize=(8, 4))
    ax = fig.add_subplot(111)

    min_signal = min(mean_signal)
    max_signal = max(mean_signal)
    ax.plot(x, mean_signal, color='red')
    ax.set_title(args.output_prefix, fontweight='bold')

    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')
    ax.set_xticks([start, 0, end])
    ax.set_xticklabels([str(start), 0, str(end)])
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels([str(round(min_signal, 2)), str(round(max_signal, 2))], rotation=90)
    ax.set_xlim(start, end)
    ax.set_ylim([min_signal, max_signal])
    ax.legend(loc="upper right", frameon=False)

    ax.spines['bottom'].set_position(('outward', 40))

    figure_name = os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="eps", dpi=300)

    # Creating canvas and printing eps / pdf with merged results
    output_fname = os.path.join(args.output_location, "{}.eps".format(args.output_prefix))
    c = pyx.canvas.canvas()
    c.insert(pyx.epsfile.epsfile(0, 0, figure_name, scale=1.0))
    c.insert(pyx.epsfile.epsfile(1.31, 0.89, logo_fname, width=18.5, height=1.75))
    c.writeEPSfile(output_fname)
    os.system("epstopdf " + figure_name)
    os.system("epstopdf " + logo_fname)
    os.system("epstopdf " + output_fname)

    os.remove(pwm_fname)
    os.remove(os.path.join(args.output_location, "{}.line.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.logo.eps".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.line.pdf".format(args.output_prefix)))
    # os.remove(os.path.join(args.output_location, "{}.logo.pdf".format(args.output_prefix)))
    os.remove(os.path.join(args.output_location, "{}.eps".format(args.output_prefix)))


def fragment_size_raw_line(args):
    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)
    bam = Samfile(args.reads_file, "rb")

    signal_f_max_145 = np.zeros(args.window_size)
    signal_r_max_145 = np.zeros(args.window_size)

    signal_f_146_307 = np.zeros(args.window_size)
    signal_r_146_307 = np.zeros(args.window_size)

    signal_f_min_307 = np.zeros(args.window_size)
    signal_r_min_307 = np.zeros(args.window_size)

    signal_f = np.zeros(args.window_size)
    signal_r = np.zeros(args.window_size)

    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            # Extend by 50 bp
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            # Fetch raw signal
            for read in bam.fetch(region.chrom, p1, p2):
                # All reads
                if not read.is_reverse:
                    cut_site = read.pos + args.forward_shift
                    if p1 <= cut_site < p2:
                        signal_f[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + args.reverse_shift - 1
                    if p1 <= cut_site < p2:
                        signal_r[cut_site - p1] += 1.0

                # length <= 145
                if abs(read.template_length) <= 145:
                    if not read.is_reverse:
                        cut_site = read.pos + args.forward_shift
                        if p1 <= cut_site < p2:
                            signal_f_max_145[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + args.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            signal_r_max_145[cut_site - p1] += 1.0

                # length > 145 and <= 307
                if 145 < abs(read.template_length) <= 307:
                    if not read.is_reverse:
                        cut_site = read.pos + args.forward_shift
                        if p1 <= cut_site < p2:
                            signal_f_146_307[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + args.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            signal_r_146_307[cut_site - p1] += 1.0

                # length > 307
                if abs(read.template_length) > 307:
                    if not read.is_reverse:
                        cut_site = read.pos + args.forward_shift
                        if p1 <= cut_site < p2:
                            signal_f_min_307[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + args.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            signal_r_min_307[cut_site - p1] += 1.0

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, signal_f)))) + "\n")
    f.write("\t".join((list(map(str, signal_r)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_max_145)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_max_145)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_146_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_146_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_min_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_min_307)))) + "\n")
    f.close()

    # find out the linker position
    pos_f_1, pos_r_1, pos_f_2, pos_r_2 = get_linkers_position(signal_f_146_307,
                                                              signal_r_146_307,
                                                              signal_f_min_307,
                                                              signal_r_min_307)
    p1 = (pos_f_1 - pos_f_2) / 2 + pos_f_2
    p2 = p1 + 180
    p3 = args.window_size - p2
    p4 = args.window_size - p1

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)
    x_ticks = [start, p1 - 500, p2 - 500, 0, p3 - 500, p4 - 500, end]

    update_axes_for_fragment_size_line(ax1, x, x_ticks, start, end, signal_f, signal_r, p1, p2,
                                       p3, p4)
    update_axes_for_fragment_size_line(ax2, x, x_ticks, start, end, signal_f_max_145, signal_r_max_145, p1, p2,
                                       p3, p4)
    update_axes_for_fragment_size_line(ax3, x, x_ticks, start, end, signal_f_146_307, signal_r_146_307, p1, p2,
                                       p3, p4)
    update_axes_for_fragment_size_line(ax4, x, x_ticks, start, end, signal_f_min_307, signal_r_min_307, p1, p2,
                                       p3, p4)

    figure_name = os.path.join(args.output_location, "{}.pdf".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="pdf", dpi=300)


def fragment_size_bc_line(args):
    mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
    mpbs_regions.read(args.motif_file)

    genomic_signal = GenomicSignal(args.reads_file)
    genomic_signal.load_sg_coefs(11)
    bam = Samfile(args.reads_file, "rb")
    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())

    bias_table = BiasTable()
    bias_table_list = args.bias_table.split(",")
    table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                  table_file_name_R=bias_table_list[1])

    signal_f_max_145 = np.zeros(args.window_size)
    signal_r_max_145 = np.zeros(args.window_size)

    signal_f_146_307 = np.zeros(args.window_size)
    signal_r_146_307 = np.zeros(args.window_size)

    signal_f_min_307 = np.zeros(args.window_size)
    signal_r_min_307 = np.zeros(args.window_size)

    signal_f = np.zeros(args.window_size)
    signal_r = np.zeros(args.window_size)

    for region in mpbs_regions:
        if str(region.name).split(":")[-1] == "Y":
            mid = (region.initial + region.final) / 2
            p1 = mid - (args.window_size / 2)
            p2 = mid + (args.window_size / 2)

            # All reads
            signal_bc_f, signal_bc_r = \
                genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta,
                                                                bias_table=table,
                                                                forward_shift=args.forward_shift,
                                                                reverse_shift=args.reverse_shift,
                                                                min_length=None, max_length=None,
                                                                strand=True)
            # length <= 145
            signal_bc_max_145_f, signal_bc_max_145_r = \
                genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta,
                                                                bias_table=table,
                                                                forward_shift=args.forward_shift,
                                                                reverse_shift=args.reverse_shift,
                                                                min_length=None, max_length=145,
                                                                strand=True)
            # length > 145 and <= 307
            signal_bc_146_307_f, signal_bc_146_307_r = \
                genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta,
                                                                bias_table=table,
                                                                forward_shift=args.forward_shift,
                                                                reverse_shift=args.reverse_shift,
                                                                min_length=145, max_length=307,
                                                                strand=True)
            # length > 307
            signal_bc_min_307_f, signal_bc_min_307_r = \
                genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta,
                                                                bias_table=table,
                                                                forward_shift=args.forward_shift,
                                                                reverse_shift=args.reverse_shift,
                                                                min_length=307, max_length=None,
                                                                strand=True)

            signal_f = np.add(signal_f, np.array(signal_bc_f))
            signal_r = np.add(signal_r, np.array(signal_bc_r))
            signal_f_max_145 = np.add(signal_f_max_145, np.array(signal_bc_max_145_f))
            signal_r_max_145 = np.add(signal_r_max_145, np.array(signal_bc_max_145_r))
            signal_f_146_307 = np.add(signal_f_146_307, np.array(signal_bc_146_307_f))
            signal_r_146_307 = np.add(signal_r_146_307, np.array(signal_bc_146_307_r))
            signal_f_min_307 = np.add(signal_f_min_307, np.array(signal_bc_min_307_f))
            signal_r_min_307 = np.add(signal_r_min_307, np.array(signal_bc_min_307_r))

    # Output the norm and slope signal
    output_fname = os.path.join(args.output_location, "{}.txt".format(args.output_prefix))
    f = open(output_fname, "w")
    f.write("\t".join((list(map(str, signal_f)))) + "\n")
    f.write("\t".join((list(map(str, signal_r)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_max_145)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_max_145)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_146_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_146_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_f_min_307)))) + "\n")
    f.write("\t".join((list(map(str, signal_r_min_307)))) + "\n")
    f.close()

    # find out the linker position
    pos_f_1, pos_r_1, pos_f_2, pos_r_2 = get_linkers_position(signal_f_146_307,
                                                              signal_r_146_307,
                                                              signal_f_min_307,
                                                              signal_r_min_307)
    p1 = (pos_f_1 - pos_f_2) / 2 + pos_f_2
    p2 = p1 + 180
    p3 = args.window_size - p2
    p4 = args.window_size - p1

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

    start = -(args.window_size / 2)
    end = (args.window_size / 2) - 1
    x = np.linspace(start, end, num=args.window_size)
    x_ticks = [start, p1 - 500, p2 - 500, 0, p3 - 500, p4 - 500, end]

    update_axes_for_fragment_size_line(ax1, x, x_ticks, start, end, signal_f, signal_r, p1, p2, p3, p4)
    update_axes_for_fragment_size_line(ax2, x, x_ticks, start, end, signal_f_max_145, signal_r_max_145, p1, p2,
                                       p3, p4)
    update_axes_for_fragment_size_line(ax3, x, x_ticks, start, end, signal_f_146_307, signal_r_146_307, p1, p2,
                                       p3, p4)
    update_axes_for_fragment_size_line(ax4, x, x_ticks, start, end, signal_f_min_307, signal_r_min_307, p1, p2,
                                       p3, p4)

    figure_name = os.path.join(args.output_location, "{}.pdf".format(args.output_prefix))
    fig.subplots_adjust(bottom=.2, hspace=.5)
    fig.tight_layout()
    fig.savefig(figure_name, format="pdf", dpi=300)


def update_axes_for_fragment_size_line(ax, x, x_ticks, start, end, signal_f, signal_r, p1, p2, p3, p4):
    max_signal = max(max(signal_f), max(signal_r))
    min_signal = min(min(signal_f), min(signal_r))
    ax.plot(x, signal_f, color='red', label='Forward')
    ax.plot(x, signal_r, color='green', label='Reverse')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 15))
    ax.tick_params(direction='out')
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(list(map(str, x_ticks)))
    ax.set_xlim(start, end)
    ax.set_yticks([min_signal, max_signal])
    ax.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
    ax.set_ylim([min_signal, max_signal])
    ax.legend().set_visible(False)
    f_1, r_1 = sum(signal_f[:p1]) / sum(signal_r), sum(signal_r[:p1]) / sum(signal_r)
    f_2, r_2 = sum(signal_f[p1:p2]) / sum(signal_r), sum(signal_r[p1:p2]) / sum(signal_r)
    f_3, r_3 = sum(signal_f[p2:500]) / sum(signal_r), sum(signal_r[p2:500]) / sum(signal_r)
    f_4, r_4 = sum(signal_f[500:p3]) / sum(signal_r), sum(signal_r[500:p3]) / sum(signal_r)
    f_5, r_5 = sum(signal_f[p3:p4]) / sum(signal_r), sum(signal_r[p3:p4]) / sum(signal_r)
    f_6, r_6 = sum(signal_f[p4:]) / sum(signal_r), sum(signal_r[p4:]) / sum(signal_r)
    text_x_1 = ((p1 - 0) / 2.0 + 0) / 1000
    text_x_2 = ((p2 - p1) / 2.0 + p1) / 1000
    text_x_3 = ((500 - p2) / 2.0 + p2) / 1000
    text_x_4 = ((p3 - 500) / 2.0 + 500) / 1000
    text_x_5 = ((p4 - p3) / 2.0 + p3) / 1000
    text_x_6 = ((1000 - p4) / 2.0 + p4) / 1000
    ax.text(text_x_1, 1.0, str(round(f_1, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_1, 0.9, str(round(r_1, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_2, 1.0, str(round(f_2, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_2, 0.9, str(round(r_2, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_3, 1.0, str(round(f_3, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_3, 0.9, str(round(r_3, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_4, 1.0, str(round(f_4, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_4, 0.9, str(round(r_4, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_5, 1.0, str(round(f_5, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_5, 0.9, str(round(r_5, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_6, 1.0, str(round(f_6, 2)), verticalalignment='center', color='red',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)
    ax.text(text_x_6, 0.9, str(round(r_6, 2)), verticalalignment='center', color='green',
            horizontalalignment='center', transform=ax.transAxes, fontsize=12)


def get_linkers_position(signal_f_146_307, signal_r_146_307, signal_f_min_307, signal_r_min_307):
    smooth_signal_f_146_307 = savgol_filter(signal_f_146_307, window_length=51, polyorder=2)
    smooth_signal_r_146_307 = savgol_filter(signal_r_146_307, window_length=51, polyorder=2)
    smooth_signal_f_min_307 = savgol_filter(signal_f_min_307, window_length=51, polyorder=2)
    smooth_signal_r_min_307 = savgol_filter(signal_r_min_307, window_length=51, polyorder=2)

    position_f_1 = np.argmax(smooth_signal_f_146_307[:400])
    position_f_2 = np.argmax(smooth_signal_f_min_307[:position_f_1])

    position_r_1 = np.argmax(smooth_signal_r_146_307[600:]) + 600
    position_r_2 = np.argmax(smooth_signal_r_min_307[position_r_1:]) + position_r_1

    return position_f_1, position_r_1, position_f_2, position_r_2


def rescaling(vector):
    maxN = max(vector)
    minN = min(vector)
    return [(e - minN) / (maxN - minN) for e in vector]


class Plot:

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
        f.write("\t".join((list(map(str, mean_signal_raw)))) + "\n")
        f.write("\t".join((list(map(str, mean_signal_bc1)))) + "\n")
        f.write("\t".join((list(map(str, mean_signal_bc2)))) + "\n")
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
                    if not read.is_reverse:
                        cut_site = read.pos + self.forward_shift
                        if p1 <= cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if p1 <= cut_site < p2:
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
        f.write("\t".join((list(map(str, mean_signal_raw_f)))) + "\n")
        f.write("\t".join((list(map(str, mean_signal_raw_r)))) + "\n")
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

    def atac_dnase_bc_line(self, reads_file1, reads_file2, bias_table1, bias_table2):
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
                signal_atac = self.get_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam_atac,
                                                 fasta=fasta, bias_table=table1,
                                                 forward_shift=5, reverse_shift=-4)

                signal_dnase = self.get_bc_signal(ref=region.chrom, start=p1, end=p2, bam=bam_dnase,
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
        f.write("\t".join((list(map(str, mean_signal_atac)))) + "\n")
        f.write("\t".join((list(map(str, mean_signal_dnase)))) + "\n")
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

    def get_bc_signal(self, ref, start, end, bam, fasta, bias_table, forward_shift, reverse_shift):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
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
            if not read.is_reverse:
                cut_site = read.pos + self.forward_shift
                if p1_w <= cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + self.reverse_shift - 1
                if p1_w <= cut_site < p2_w:
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

    def atac_dnase_raw_line(self, reads_file1, reads_file2):
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
                    if not read.is_reverse:
                        cut_site = read.pos + self.forward_shift
                        if p1 <= cut_site < p2:
                            mean_signal_atac[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            mean_signal_atac[cut_site - p1] += 1.0

                # Fetch raw signal
                for read in bam_dnase.fetch(region.chrom, p1, p2):
                    if not read.is_reverse:
                        cut_site = read.pos
                        if p1 <= cut_site < p2:
                            mean_signal_dnase[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend - 1
                        if p1 <= cut_site < p2:
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
        f.write("\t".join((list(map(str, mean_signal_atac)))) + "\n")
        f.write("\t".join((list(map(str, mean_signal_dnase)))) + "\n")
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

        ax2.spines['bottom'].set_position(('outward', 40))

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

    def fragment_size_raw_line(self):
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        signal_f_max_145 = np.zeros(self.window_size)
        signal_r_max_145 = np.zeros(self.window_size)

        signal_f_146_307 = np.zeros(self.window_size)
        signal_r_146_307 = np.zeros(self.window_size)

        signal_f_min_307 = np.zeros(self.window_size)
        signal_r_min_307 = np.zeros(self.window_size)

        signal_f = np.zeros(self.window_size)
        signal_r = np.zeros(self.window_size)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                for read in bam.fetch(region.chrom, p1, p2):
                    # All reads
                    if not read.is_reverse:
                        cut_site = read.pos + self.forward_shift
                        if p1 <= cut_site < p2:
                            signal_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            signal_r[cut_site - p1] += 1.0

                    # length <= 145
                    if abs(read.template_length) <= 145:
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_f_max_145[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_r_max_145[cut_site - p1] += 1.0

                    # length > 145 and <= 307
                    if 145 < abs(read.template_length) <= 307:
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_f_146_307[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_r_146_307[cut_site - p1] += 1.0

                    # length > 307
                    if abs(read.template_length) > 307:
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_f_min_307[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_r_min_307[cut_site - p1] += 1.0

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((list(map(str, signal_f)))) + "\n")
        f.write("\t".join((list(map(str, signal_r)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_max_145)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_max_145)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_min_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_min_307)))) + "\n")
        f.close()

        # find out the linker position
        pos_f_1, pos_r_1, pos_f_2, pos_r_2 = self.get_linkers_position(signal_f_146_307,
                                                                       signal_r_146_307,
                                                                       signal_f_min_307,
                                                                       signal_r_min_307)
        p1 = (pos_f_1 - pos_f_2) / 2 + pos_f_2
        p2 = p1 + 180
        p3 = self.window_size - p2
        p4 = self.window_size - p1

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)
        x_ticks = [start, p1 - 500, p2 - 500, 0, p3 - 500, p4 - 500, end]

        self.update_axes_for_fragment_size_line(ax1, x, x_ticks, start, end, signal_f, signal_r, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax2, x, x_ticks, start, end, signal_f_max_145, signal_r_max_145, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax3, x, x_ticks, start, end, signal_f_146_307, signal_r_146_307, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax4, x, x_ticks, start, end, signal_f_min_307, signal_r_min_307, p1, p2,
                                                p3, p4)

        figure_name = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="pdf", dpi=300)

    def fragment_size_bc_line(self, bias_table_files):
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)

        genomic_signal = GenomicSignal(self.reads_file)
        genomic_signal.load_sg_coefs(11)
        bam = Samfile(self.reads_file, "rb")
        genome_data = GenomeData(self.organism)
        fasta = Fastafile(genome_data.get_genome())

        bias_table = BiasTable()
        bias_table_list = bias_table_files.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

        signal_f_max_145 = np.zeros(self.window_size)
        signal_r_max_145 = np.zeros(self.window_size)

        signal_f_146_307 = np.zeros(self.window_size)
        signal_r_146_307 = np.zeros(self.window_size)

        signal_f_min_307 = np.zeros(self.window_size)
        signal_r_min_307 = np.zeros(self.window_size)

        signal_f = np.zeros(self.window_size)
        signal_r = np.zeros(self.window_size)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # All reads
                signal_bc_f, signal_bc_r = \
                    genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                    bam=bam, fasta=fasta,
                                                                    bias_table=table,
                                                                    forward_shift=self.forward_shift,
                                                                    reverse_shift=self.reverse_shift,
                                                                    min_length=None, max_length=None,
                                                                    strand=True)
                # length <= 145
                signal_bc_max_145_f, signal_bc_max_145_r = \
                    genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                    bam=bam, fasta=fasta,
                                                                    bias_table=table,
                                                                    forward_shift=self.forward_shift,
                                                                    reverse_shift=self.reverse_shift,
                                                                    min_length=None, max_length=145,
                                                                    strand=True)
                # length > 145 and <= 307
                signal_bc_146_307_f, signal_bc_146_307_r = \
                    genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                    bam=bam, fasta=fasta,
                                                                    bias_table=table,
                                                                    forward_shift=self.forward_shift,
                                                                    reverse_shift=self.reverse_shift,
                                                                    min_length=145, max_length=307,
                                                                    strand=True)
                # length > 307
                signal_bc_min_307_f, signal_bc_min_307_r = \
                    genomic_signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                    bam=bam, fasta=fasta,
                                                                    bias_table=table,
                                                                    forward_shift=self.forward_shift,
                                                                    reverse_shift=self.reverse_shift,
                                                                    min_length=307, max_length=None,
                                                                    strand=True)

                signal_f = np.add(signal_f, np.array(signal_bc_f))
                signal_r = np.add(signal_r, np.array(signal_bc_r))
                signal_f_max_145 = np.add(signal_f_max_145, np.array(signal_bc_max_145_f))
                signal_r_max_145 = np.add(signal_r_max_145, np.array(signal_bc_max_145_r))
                signal_f_146_307 = np.add(signal_f_146_307, np.array(signal_bc_146_307_f))
                signal_r_146_307 = np.add(signal_r_146_307, np.array(signal_bc_146_307_r))
                signal_f_min_307 = np.add(signal_f_min_307, np.array(signal_bc_min_307_f))
                signal_r_min_307 = np.add(signal_r_min_307, np.array(signal_bc_min_307_r))

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((list(map(str, signal_f)))) + "\n")
        f.write("\t".join((list(map(str, signal_r)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_max_145)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_max_145)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_min_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_min_307)))) + "\n")
        f.close()

        # find out the linker position
        pos_f_1, pos_r_1, pos_f_2, pos_r_2 = self.get_linkers_position(signal_f_146_307,
                                                                       signal_r_146_307,
                                                                       signal_f_min_307,
                                                                       signal_r_min_307)
        p1 = (pos_f_1 - pos_f_2) / 2 + pos_f_2
        p2 = p1 + 180
        p3 = self.window_size - p2
        p4 = self.window_size - p1

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)
        x_ticks = [start, p1 - 500, p2 - 500, 0, p3 - 500, p4 - 500, end]

        self.update_axes_for_fragment_size_line(ax1, x, x_ticks, start, end, signal_f, signal_r, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax2, x, x_ticks, start, end, signal_f_max_145, signal_r_max_145, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax3, x, x_ticks, start, end, signal_f_146_307, signal_r_146_307, p1, p2,
                                                p3, p4)
        self.update_axes_for_fragment_size_line(ax4, x, x_ticks, start, end, signal_f_min_307, signal_r_min_307, p1, p2,
                                                p3, p4)

        figure_name = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="pdf", dpi=300)

    def slope_of_reads_number(self):
        genome_data = GenomeData(self.organism)
        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)
        bam = Samfile(self.reads_file, "rb")

        signal_f_146_307 = np.zeros(self.window_size)
        signal_r_146_307 = np.zeros(self.window_size)

        signal_f_min_307 = np.zeros(self.window_size)
        signal_r_min_307 = np.zeros(self.window_size)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                # Extend by 50 bp
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                # Fetch raw signal
                for read in bam.fetch(region.chrom, p1, p2):
                    # length > 145 and <= 307
                    if 145 < abs(read.template_length) <= 307:
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_f_146_307[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_r_146_307[cut_site - p1] += 1.0

                    # length > 307
                    if abs(read.template_length) > 307:
                        if not read.is_reverse:
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_f_min_307[cut_site - p1] += 1.0
                        else:
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_r_min_307[cut_site - p1] += 1.0

        # Apply a Savitzky-Golay filter to an array.
        smooth_signal_f_146_307 = savgol_filter(signal_f_146_307, window_length=51, polyorder=2)
        smooth_signal_r_146_307 = savgol_filter(signal_r_146_307, window_length=51, polyorder=2)
        smooth_signal_f_min_307 = savgol_filter(signal_f_min_307, window_length=51, polyorder=2)
        smooth_signal_r_min_307 = savgol_filter(signal_r_min_307, window_length=51, polyorder=2)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        f = open(output_fname, "w")
        f.write("\t".join((list(map(str, signal_f_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_146_307)))) + "\n")
        f.write("\t".join((list(map(str, smooth_signal_f_146_307)))) + "\n")
        f.write("\t".join((list(map(str, smooth_signal_r_146_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_f_min_307)))) + "\n")
        f.write("\t".join((list(map(str, signal_r_min_307)))) + "\n")
        f.write("\t".join((list(map(str, smooth_signal_f_min_307)))) + "\n")
        f.write("\t".join((list(map(str, smooth_signal_r_min_307)))) + "\n")
        f.close()

        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)

        fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(8, 8))

        min_signal = min(min(signal_f_146_307), min(signal_r_146_307))
        max_signal = max(max(signal_f_146_307), max(signal_r_146_307))
        ax1.plot(x, signal_f_146_307, color='red', label='Forward')
        ax1.plot(x, signal_r_146_307, color='green', label='Reverse')
        ax1.xaxis.set_ticks_position('bottom')
        ax1.yaxis.set_ticks_position('left')
        ax1.spines['top'].set_visible(False)
        ax1.spines['right'].set_visible(False)
        ax1.spines['left'].set_position(('outward', 15))
        ax1.tick_params(direction='out')
        ax1.set_title("Raw reads number of 1N", fontweight='bold')
        ax1.set_xticks([start, 0, end])
        ax1.set_xticklabels([str(start), 0, str(end)])
        ax1.set_yticks([min_signal, max_signal])
        ax1.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax1.set_xlim(start, end)
        ax1.set_ylim([min_signal, max_signal])

        min_signal = min(min(smooth_signal_f_146_307), min(smooth_signal_r_146_307))
        max_signal = max(max(smooth_signal_f_146_307), max(smooth_signal_r_146_307))
        ax2.plot(x, smooth_signal_f_146_307, color='red', label='Forward')
        ax2.plot(x, smooth_signal_r_146_307, color='green', label='Reverse')
        ax2.xaxis.set_ticks_position('bottom')
        ax2.yaxis.set_ticks_position('left')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)
        ax2.spines['left'].set_position(('outward', 15))
        ax2.tick_params(direction='out')
        ax2.set_title("Smooth reads number of 1N", fontweight='bold')
        ax2.set_yticks([min_signal, max_signal])
        ax2.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax2.set_xlim(start, end)
        ax2.set_ylim([min_signal, max_signal])

        id_max_f_1 = np.argmax(smooth_signal_f_146_307[:500]) - 500
        id_max_r_1 = np.argmax(smooth_signal_r_146_307[500:])
        ax2.axvline(x=id_max_f_1)
        ax2.axvline(x=id_max_r_1)
        x_tick_labels = [str(start), str(id_max_f_1), 0, str(id_max_r_1), str(end)]
        ax2.set_xticks(list(map(int, x_tick_labels)))
        ax2.set_xticklabels(x_tick_labels, rotation=90)

        min_signal = min(min(signal_f_min_307), min(signal_r_min_307))
        max_signal = max(max(signal_f_min_307), max(signal_r_min_307))
        ax3.plot(x, signal_f_min_307, color='red', label='Forward')
        ax3.plot(x, signal_r_min_307, color='green', label='Reverse')
        ax3.xaxis.set_ticks_position('bottom')
        ax3.yaxis.set_ticks_position('left')
        ax3.spines['top'].set_visible(False)
        ax3.spines['right'].set_visible(False)
        ax3.spines['left'].set_position(('outward', 15))
        ax3.tick_params(direction='out')
        ax3.set_title("Raw reads number of +2N", fontweight='bold')
        ax3.set_xticks([start, 0, end])
        ax3.set_xticklabels([str(start), 0, str(end)])
        ax3.set_yticks([min_signal, max_signal])
        ax3.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax3.set_xlim(start, end)
        ax3.set_ylim([min_signal, max_signal])

        min_signal = min(min(smooth_signal_f_min_307), min(smooth_signal_r_min_307))
        max_signal = max(max(smooth_signal_f_min_307), max(smooth_signal_r_min_307))
        ax4.plot(x, smooth_signal_f_min_307, color='red', label='Forward')
        ax4.plot(x, smooth_signal_r_min_307, color='green', label='Reverse')
        ax4.xaxis.set_ticks_position('bottom')
        ax4.yaxis.set_ticks_position('left')
        ax4.spines['top'].set_visible(False)
        ax4.spines['right'].set_visible(False)
        ax4.spines['left'].set_position(('outward', 15))
        ax4.tick_params(direction='out')
        ax4.set_title("Smooth reads number of +2N", fontweight='bold')
        ax4.set_yticks([min_signal, max_signal])
        ax4.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax4.set_xlim(start, end)
        ax4.set_ylim([min_signal, max_signal])

        id_max_f = np.argmax(smooth_signal_f_146_307[:500])
        id_max_r = np.argmax(smooth_signal_r_146_307[500:]) + 500
        id_max_f_2 = np.argmax(smooth_signal_f_min_307[:id_max_f]) - 500
        id_max_r_2 = np.argmax(smooth_signal_r_min_307[id_max_r:]) + id_max_r - 500
        ax4.axvline(x=id_max_f_2)
        ax4.axvline(x=id_max_r_2)
        x_tick_labels = [str(start), str(id_max_f_2), 0, str(id_max_r_2), str(end)]
        ax4.set_xticks(list(map(int, x_tick_labels)))
        ax4.set_xticklabels(x_tick_labels, rotation=90)

        figure_name = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_name, format="pdf", dpi=300)

    def strand_line_by_size(self, bias_tables):
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
            f.write("\t".join((list(map(str, signal_bc_f)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_145_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_307_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_307_500_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_307_500_slope)))) + "\n")

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
            f.write("\t".join((list(map(str, signal_bc_f_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_307_500_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_307_500_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_f_min_500_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_r_min_500_slope)))) + "\n")

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

    def unstrand_line_by_size(self, bias_tables):
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

        signal_bc_max_145 = np.zeros(self.window_size)

        signal_bc_min_145 = np.zeros(self.window_size)

        signal_bc_145_307 = np.zeros(self.window_size)

        signal_bc_min_307 = np.zeros(self.window_size)

        signal_bc_307_500 = np.zeros(self.window_size)

        signal_bc_min_500 = np.zeros(self.window_size)

        signal_bc = np.zeros(self.window_size)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                bc_max_145 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=None, max_length=145, strand=False)
                bc_min_145 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=145, max_length=None, strand=False)
                bc_145_307 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=145, max_length=307, strand=False)
                bc_min_307 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=307, max_length=None, strand=False)
                bc_307_500 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=307, max_length=500, strand=False)
                bc_min_500 = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=500, max_length=None, strand=False)

                bc = \
                    signal.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                            bam=bam, fasta=fasta, bias_table=bias_table,
                                                            forward_shift=self.forward_shift,
                                                            reverse_shift=self.reverse_shift,
                                                            min_length=None, max_length=None, strand=False)

                signal_bc_max_145 = np.add(signal_bc_max_145, bc_max_145)
                signal_bc_min_145 = np.add(signal_bc_min_145, bc_min_145)
                signal_bc_145_307 = np.add(signal_bc_145_307, bc_145_307)
                signal_bc_min_307 = np.add(signal_bc_min_307, bc_min_307)
                signal_bc_307_500 = np.add(signal_bc_307_500, bc_307_500)
                signal_bc_min_500 = np.add(signal_bc_min_500, bc_min_500)
                signal_bc = np.add(signal_bc, bc)

        # Pre-process signal
        signal_bc = signal.boyle_norm(signal_bc)
        perc = scoreatpercentile(signal_bc, 98)
        std = np.array(signal_bc).std()
        signal_bc = signal.hon_norm_atac(signal_bc, perc, std)
        signal_bc_slope = signal.slope(signal_bc, signal.sg_coefs)

        # Output the norm and slope signal
        output_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_slope)))) + "\n")

        # Output signal of class a
        signal_bc_max_145 = signal.boyle_norm(signal_bc_max_145)
        perc = scoreatpercentile(signal_bc_max_145, 98)
        std = np.array(signal_bc_max_145).std()
        signal_bc_max_145 = signal.hon_norm_atac(signal_bc_max_145, perc, std)
        signal_bc_max_145_slope = signal.slope(signal_bc_max_145, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_a.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")

        # Output signal of class b
        signal_bc_min_145 = signal.boyle_norm(signal_bc_min_145)
        perc = scoreatpercentile(signal_bc_min_145, 98)
        std = np.array(signal_bc_min_145).std()
        signal_bc_min_145 = signal.hon_norm_atac(signal_bc_min_145, perc, std)
        signal_bc_min_145_slope = signal.slope(signal_bc_min_145, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_b.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_145_slope)))) + "\n")

        # Output signal of class c
        signal_bc_145_307 = signal.boyle_norm(signal_bc_145_307)
        perc = scoreatpercentile(signal_bc_145_307, 98)
        std = np.array(signal_bc_145_307).std()
        signal_bc_145_307 = signal.hon_norm_atac(signal_bc_145_307, perc, std)
        signal_bc_145_307_slope = signal.slope(signal_bc_145_307, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_c.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307_slope)))) + "\n")

        # Output signal of class d
        signal_bc_min_307 = signal.boyle_norm(signal_bc_min_307)
        perc = scoreatpercentile(signal_bc_min_307, 98)
        std = np.array(signal_bc_min_307).std()
        signal_bc_min_307 = signal.hon_norm_atac(signal_bc_min_307, perc, std)
        signal_bc_min_307_slope = signal.slope(signal_bc_min_307, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_d.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_307_slope)))) + "\n")

        # Output signal of class e
        signal_bc_307_500 = signal.boyle_norm(signal_bc_307_500)
        perc = scoreatpercentile(signal_bc_307_500, 98)
        std = np.array(signal_bc_307_500).std()
        signal_bc_307_500 = signal.hon_norm_atac(signal_bc_307_500, perc, std)
        signal_bc_307_500_slope = signal.slope(signal_bc_307_500, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_e.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_307_500_slope)))) + "\n")

        # Output signal of class f
        signal_bc_min_500 = signal.boyle_norm(signal_bc_min_500)
        perc = scoreatpercentile(signal_bc_min_500, 98)
        std = np.array(signal_bc_min_500).std()
        signal_bc_min_500 = signal.hon_norm_atac(signal_bc_min_500, perc, std)
        signal_bc_min_500_slope = signal.slope(signal_bc_min_500, signal.sg_coefs)

        output_fname = os.path.join(self.output_loc, "{}_f.txt".format(self.output_prefix))
        with open(output_fname, "w") as f:
            f.write("\t".join((list(map(str, signal_bc_max_145)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_max_145_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_145_307_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_307_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_307_500_slope)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_500)))) + "\n")
            f.write("\t".join((list(map(str, signal_bc_min_500_slope)))) + "\n")

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

        min_signal = min(signal_bc)
        max_signal = max(signal_bc)
        ax[0, 1].plot(x, signal_bc, color='red')
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

        min_signal = min(signal_bc_max_145)
        max_signal = max(signal_bc_max_145)
        ax[1, 0].plot(x, signal_bc_max_145, color='red')
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

        min_signal = min(signal_bc_min_145)
        max_signal = max(signal_bc_min_145)
        ax[1, 1].plot(x, signal_bc_min_145, color='red')
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

        min_signal = min(signal_bc_145_307)
        max_signal = max(signal_bc_145_307)
        ax[2, 0].plot(x, signal_bc_145_307, color='red')
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

        min_signal = min(signal_bc_min_307)
        max_signal = max(signal_bc_min_307)
        ax[2, 1].plot(x, signal_bc_min_307, color='red')
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

        min_signal = min(signal_bc_307_500)
        max_signal = max(signal_bc_307_500)
        ax[3, 0].plot(x, signal_bc_307_500, color='red')
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

        min_signal = min(signal_bc_min_500)
        max_signal = max(signal_bc_min_500)
        ax[3, 1].plot(x, signal_bc_min_500, color='red')
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
                    if 145 < length <= 307: lengths_between_146_307.append(length)
                    if 307 < length <= 500: lengths_between_307_500.append(length)
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
                            if 145 < length <= 307: lengths_between_146_307.append(length)
                            if 307 < length <= 500: lengths_between_307_500.append(length)
                            if length > 500: lengths_longer_than_500.append(length)

        # output the plot
        unique, counts = np.unique(lengths, return_counts=True)
        count_list = list()
        length_dict = dict(list(zip(unique, counts)))
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
        length_f_dict = dict(list(zip(unique, counts)))
        unique, counts = np.unique(lengths_reverse, return_counts=True)
        length_r_dict = dict(list(zip(unique, counts)))
        sort_keys = sorted(length_f_dict.keys())
        txt_fname = os.path.join(self.output_loc, "{}.txt".format(self.output_prefix))
        with open(txt_fname, "w") as f:
            f.write("Length\tCount_Forward\tCount_Reverse\n")
            for key in sort_keys:
                f.write(str(key) + "\t" + str(length_f_dict[key]) + "\t" + str(length_r_dict[key]) + "\n")

    def signal_by_position(self, reads_file=None, motif_file=None):
        bam = Samfile(self.reads_file, "rb")

        def update_signal(length_f, length_r, signal_f, signal_r, read, p1, p2):
            if not read.is_reverse:
                length_f.append(abs(read.template_length))
                cut_site = read.pos + self.forward_shift
                if p1 <= cut_site < p2:
                    signal_f[cut_site - p1] += 1.0
            else:
                length_r.append(abs(read.template_length))
                cut_site = read.aend + self.reverse_shift - 1
                if p1 <= cut_site < p2:
                    signal_r[cut_site - p1] += 1.0

        def output2file(output_loc, output_prefix, outputs):
            fname = os.path.join(output_loc, "{}.txt".format(output_prefix))
            with open(fname, "w") as f:
                f.write("values" + "\t" + output_prefix + "\n")
                for out in outputs:
                    f.write(str(out) + "\t" + output_prefix + "\n")

        length_f_Nfr_I = list()
        length_r_Nfr_I = list()
        length_f_Nfr_II = list()
        length_r_Nfr_II = list()

        length_f_1N_I = list()
        length_r_1N_I = list()
        length_f_1N_II = list()
        length_r_1N_II = list()
        length_f_1N_III = list()
        length_r_1N_III = list()

        length_f_2N_I = list()
        length_r_2N_I = list()
        length_f_2N_II = list()
        length_r_2N_II = list()
        length_f_2N_III = list()
        length_r_2N_III = list()

        signal_f_Nfr_I = np.zeros(self.window_size)
        signal_r_Nfr_I = np.zeros(self.window_size)
        signal_f_Nfr_II = np.zeros(self.window_size)
        signal_r_Nfr_II = np.zeros(self.window_size)

        signal_f_1N_I = np.zeros(self.window_size)
        signal_r_1N_I = np.zeros(self.window_size)
        signal_f_1N_II = np.zeros(self.window_size)
        signal_r_1N_II = np.zeros(self.window_size)
        signal_f_1N_III = np.zeros(self.window_size)
        signal_r_1N_III = np.zeros(self.window_size)

        signal_f_2N_I = np.zeros(self.window_size)
        signal_r_2N_I = np.zeros(self.window_size)
        signal_f_2N_II = np.zeros(self.window_size)
        signal_r_2N_II = np.zeros(self.window_size)
        signal_f_2N_III = np.zeros(self.window_size)
        signal_r_2N_III = np.zeros(self.window_size)

        mpbs_regions = GenomicRegionSet("Motif Predicted Binding Sites")
        mpbs_regions.read(self.motif_file)

        for region in mpbs_regions:
            if str(region.name).split(":")[-1] == "Y":
                mid = (region.initial + region.final) / 2
                p1 = mid - (self.window_size / 2)
                p2 = mid + (self.window_size / 2)

                for read in bam.fetch(region.chrom, p1, p2):
                    if not read.is_read1: continue
                    if not read.is_reverse:
                        start = read.pos - mid
                        end = read.pos + abs(read.template_length) - mid
                    else:
                        start = read.aend - abs(read.template_length) - mid
                        end = read.aend - mid

                    # Nucleosome free I
                    if abs(read.template_length) <= 146:
                        if (start in range(-142, 0) and end in range(-142, 0)) or (
                                start in range(0, 142) and end in range(0, 142)):
                            update_signal(length_f_Nfr_I, length_r_Nfr_I, signal_f_Nfr_I, signal_r_Nfr_I, read, p1, p2)

                        # Nucleosome free II
                        elif start in range(-142, 0) and end in range(0, 142):
                            update_signal(length_f_Nfr_II, length_r_Nfr_II, signal_f_Nfr_II, signal_r_Nfr_II, read, p1,
                                          p2)

                    elif 146 < abs(read.template_length) <= 307:
                        # 1 Nucleosome I
                        if (start in range(-322, -142) and end in range(-142, 0)) or (
                                start in range(0, 142) and end in range(142, 322)):
                            update_signal(length_f_1N_I, length_r_1N_I, signal_f_1N_I, signal_r_1N_I, read, p1, p2)

                        # 1 Nucleosome II
                        elif (start in range(-322, -142) and end in range(0, 142)) or (
                                start in range(-142, 0) and end in range(142, 322)):
                            update_signal(length_f_1N_II, length_r_1N_II, signal_f_1N_II, signal_r_1N_II, read, p1, p2)

                        # 1 Nucleosome III
                        elif (start in range(-500, -322) and end in range(-322, -142)) or (
                                start in range(142, 322) and end in range(322, 500)):
                            update_signal(length_f_1N_III, length_r_1N_III, signal_f_1N_III, signal_r_1N_III, read, p1,
                                          p2)

                    elif 307 < abs(read.template_length):
                        # 2 Nucleosome I
                        if (start in range(-500, -322) and end in range(-142, 0)) or (
                                start in range(0, 142) and end in range(322, 500)):
                            update_signal(length_f_2N_I, length_r_2N_I, signal_f_2N_I, signal_r_2N_I, read, p1, p2)

                        # 2 Nucleosome II
                        elif (start in range(-500, -322) and end in range(0, 142)) or (
                                start in range(-142, 0) and end in range(322, 500)):
                            update_signal(length_f_2N_II, length_r_2N_II, signal_f_2N_II, signal_r_2N_II, read, p1, p2)

                        # 2 Nucleosome III
                        elif start in range(-322, -142) and end in range(142, 322):
                            update_signal(length_f_2N_III, length_r_2N_III, signal_f_2N_III, signal_r_2N_III, read, p1,
                                          p2)

        fig, ax = plt.subplots(8, 2, figsize=(10, 15))
        p1, p2, p3, p4 = 178, 358, 642, 822
        start = -(self.window_size / 2)
        end = (self.window_size / 2) - 1
        x = np.linspace(start, end, num=self.window_size)
        x_ticks = [start, p1 - 500, p2 - 500, 0, p3 - 500, p4 - 500, end]

        self.update_axes_for_density(ax[0, 0], length_f_Nfr_I, length_r_Nfr_I, "Nucleosome free I")
        self.update_axes_for_density(ax[1, 0], length_f_Nfr_II, length_r_Nfr_II, "Nucleosome free II")
        self.update_axes_for_density(ax[2, 0], length_f_1N_I, length_r_1N_I, "1 Nucleosome I")
        self.update_axes_for_density(ax[3, 0], length_f_1N_II, length_r_1N_II, "1 Nucleosome II")
        self.update_axes_for_density(ax[4, 0], length_f_1N_III, length_r_1N_III, "1 Nucleosome III")
        self.update_axes_for_density(ax[5, 0], length_f_2N_I, length_r_2N_I, "2 Nucleosome I")
        self.update_axes_for_density(ax[6, 0], length_f_2N_II, length_r_2N_II, "2 Nucleosome II")
        self.update_axes_for_density(ax[7, 0], length_f_2N_III, length_r_2N_III, "2 Nucleosome III")
        self.update_axes_for_signal_by_position(ax[0, 1], x, x_ticks, start, end, signal_f_Nfr_I, signal_r_Nfr_I, p1,
                                                p2,
                                                p3, p4, "Nucleosome free I")
        self.update_axes_for_signal_by_position(ax[1, 1], x, x_ticks, start, end, signal_f_Nfr_II, signal_r_Nfr_II, p1,
                                                p2,
                                                p3, p4, "Nucleosome free II")
        self.update_axes_for_signal_by_position(ax[2, 1], x, x_ticks, start, end, signal_f_1N_I, signal_r_1N_I, p1, p2,
                                                p3, p4, "1 Nucleosome I")
        self.update_axes_for_signal_by_position(ax[3, 1], x, x_ticks, start, end, signal_f_1N_II, signal_r_1N_II, p1,
                                                p2,
                                                p3, p4, "1 Nucleosome II")
        self.update_axes_for_signal_by_position(ax[4, 1], x, x_ticks, start, end, signal_f_1N_III, signal_r_1N_III, p1,
                                                p2,
                                                p3, p4, "1 Nucleosome III")
        self.update_axes_for_signal_by_position(ax[5, 1], x, x_ticks, start, end, signal_f_2N_I, signal_r_2N_I, p1, p2,
                                                p3, p4, "2 Nucleosome I")
        self.update_axes_for_signal_by_position(ax[6, 1], x, x_ticks, start, end, signal_f_2N_II, signal_r_2N_II, p1,
                                                p2,
                                                p3, p4, "2 Nucleosome II")
        self.update_axes_for_signal_by_position(ax[7, 1], x, x_ticks, start, end, signal_f_2N_III, signal_r_2N_III, p1,
                                                p2,
                                                p3, p4, "2 Nucleosome III")
        figure_fname = os.path.join(self.output_loc, "{}.pdf".format(self.output_prefix))
        fig.subplots_adjust(bottom=.2, hspace=.5)
        fig.tight_layout()
        fig.savefig(figure_fname, format="pdf", dpi=300)

        output2file(self.output_loc, "{}_length_f_Nfr_I".format(self.output_prefix), length_f_Nfr_I)
        output2file(self.output_loc, "{}_length_r_Nfr_I".format(self.output_prefix), length_r_Nfr_I)
        output2file(self.output_loc, "{}_length_f_Nfr_II".format(self.output_prefix), length_f_Nfr_II)
        output2file(self.output_loc, "{}_length_r_Nfr_II".format(self.output_prefix), length_r_Nfr_II)
        output2file(self.output_loc, "{}_length_f_1N_I".format(self.output_prefix), length_f_1N_I)
        output2file(self.output_loc, "{}_length_r_1N_I".format(self.output_prefix), length_r_1N_I)
        output2file(self.output_loc, "{}_length_f_1N_II".format(self.output_prefix), length_f_1N_II)
        output2file(self.output_loc, "{}_length_r_1N_II".format(self.output_prefix), length_r_1N_II)
        output2file(self.output_loc, "{}_length_f_1N_III".format(self.output_prefix), length_f_1N_III)
        output2file(self.output_loc, "{}_length_r_1N_III".format(self.output_prefix), length_r_1N_III)
        output2file(self.output_loc, "{}_length_f_2N_I".format(self.output_prefix), length_f_2N_I)
        output2file(self.output_loc, "{}_length_r_2N_I".format(self.output_prefix), length_r_2N_I)
        output2file(self.output_loc, "{}_length_f_2N_II".format(self.output_prefix), length_f_2N_II)
        output2file(self.output_loc, "{}_length_r_2N_II".format(self.output_prefix), length_r_2N_II)
        output2file(self.output_loc, "{}_length_f_2N_III".format(self.output_prefix), length_f_2N_III)
        output2file(self.output_loc, "{}_length_r_2N_III".format(self.output_prefix), length_r_2N_III)

        output2file(self.output_loc, "{}_signal_f_Nfr_I".format(self.output_prefix), np.roll(signal_f_Nfr_I, -4))
        output2file(self.output_loc, "{}_signal_r_Nfr_I".format(self.output_prefix), np.roll(signal_r_Nfr_I, 4))
        output2file(self.output_loc, "{}_signal_f_Nfr_II".format(self.output_prefix), np.roll(signal_f_Nfr_II, -4))
        output2file(self.output_loc, "{}_signal_r_Nfr_II".format(self.output_prefix), np.roll(signal_r_Nfr_II, 4))
        output2file(self.output_loc, "{}_signal_f_1N_I".format(self.output_prefix), np.roll(signal_f_1N_I, -4))
        output2file(self.output_loc, "{}_signal_r_1N_I".format(self.output_prefix), np.roll(signal_r_1N_I, 4))
        output2file(self.output_loc, "{}_signal_f_1N_II".format(self.output_prefix), np.roll(signal_f_1N_II, -4))
        output2file(self.output_loc, "{}_signal_r_1N_II".format(self.output_prefix), np.roll(signal_r_1N_II, 4))
        output2file(self.output_loc, "{}_signal_f_1N_III".format(self.output_prefix), np.roll(signal_f_1N_III, -4))
        output2file(self.output_loc, "{}_signal_r_1N_III".format(self.output_prefix), np.roll(signal_r_1N_III, 4))
        output2file(self.output_loc, "{}_signal_f_2N_I".format(self.output_prefix), np.roll(signal_f_2N_I, -4))
        output2file(self.output_loc, "{}_signal_r_2N_I".format(self.output_prefix), np.roll(signal_r_2N_I, 4))
        output2file(self.output_loc, "{}_signal_f_2N_II".format(self.output_prefix), np.roll(signal_f_2N_II, -4))
        output2file(self.output_loc, "{}_signal_r_2N_II".format(self.output_prefix), np.roll(signal_r_2N_II, 4))
        output2file(self.output_loc, "{}_signal_f_2N_III".format(self.output_prefix), np.roll(signal_f_2N_III, -4))
        output2file(self.output_loc, "{}_signal_r_2N_III".format(self.output_prefix), np.roll(signal_r_2N_III, 4))

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

                    if not read.is_reverse:
                        cut_site = read.pos + self.forward_shift
                        if p1 <= cut_site < p2:
                            signal_raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + self.reverse_shift - 1
                        if p1 <= cut_site < p2:
                            signal_raw_r[cut_site - p1] += 1.0

                    if length <= 145:
                        # nucleosome free
                        if not read.is_reverse:
                            if mid - 142 < read.pos < mid < read.pos + length < mid + 142:
                                lengths_class_1.append(length)
                                lengths_class_1_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if p1 <= cut_site < p2:
                                    signal_raw_class_1_f[cut_site - p1] += 1.0
                        else:
                            # if mid - 142 < read.pos < mid < read.pos + length < mid + 142:
                            if mid - 142 < read.aend - length < mid < read.aend < mid + 142:
                                lengths_class_1.append(length)
                                lengths_class_1_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if p1 <= cut_site < p2:
                                    signal_raw_class_1_r[cut_site - p1] += 1.0

                    elif 145 < length <= 307:
                        # nucleosome only
                        if not read.is_reverse:
                            if (mid - 248 < read.pos < mid - 142 <= read.pos + length <= mid) or \
                                    (mid < read.pos <= mid + 142 < read.pos + length < mid + 284):
                                lengths_class_2.append(length)
                                lengths_class_2_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if p1 <= cut_site < p2:
                                    signal_raw_class_2_f[cut_site - p1] += 1.0
                        else:
                            if (mid - 248 < read.aend - length < mid - 142 <= read.aend <= mid) or \
                                    (mid < read.aend - length <= mid + 142 < read.aend < mid + 284):
                                lengths_class_2.append(length)
                                lengths_class_2_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if p1 <= cut_site < p2:
                                    signal_raw_class_2_r[cut_site - p1] += 1.0

                        # nucleosome and TF
                        if not read.is_reverse:
                            if read.pos < mid - 142 <= read.pos + length > mid or \
                                    (mid - 142 <= read.pos <= mid and read.pos + length > mid + 142):
                                lengths_class_3.append(length)
                                lengths_class_3_f.append(length)
                                cut_site = read.pos + self.forward_shift
                                if p1 <= cut_site < p2:
                                    signal_raw_class_3_f[cut_site - p1] += 1.0
                        else:
                            if (read.aend - length < mid - 142 and mid <= read.aend < mid + 142) or \
                                    (read.aend > mid + 142 and mid - 142 <= read.aend - length < mid):
                                lengths_class_3.append(length)
                                lengths_class_3_r.append(length)
                                cut_site = read.aend + self.reverse_shift - 1
                                if p1 <= cut_site < p2:
                                    signal_raw_class_3_r[cut_site - p1] += 1.0

                    elif length > 307:
                        lengths_class_4.append(length)
                        if not read.is_reverse:
                            lengths_class_4_f.append(length)
                            cut_site = read.pos + self.forward_shift
                            if p1 <= cut_site < p2:
                                signal_raw_class_4_f[cut_site - p1] += 1.0
                        else:
                            lengths_class_4_r.append(length)
                            cut_site = read.aend + self.reverse_shift - 1
                            if p1 <= cut_site < p2:
                                signal_raw_class_4_r[cut_site - p1] += 1.0

        # output the plot
        unique, counts = np.unique(lengths, return_counts=True)
        count_list = list()
        length_dict = dict(list(zip(unique, counts)))
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

        ax[1, 0].set_title("Histogram of Nucleosome free", fontweight='bold')
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

        ax[2, 0].set_title("Histogram of Nucleosome I", fontweight='bold')
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

        ax[3, 0].set_title("Histogram of Nucleosome II", fontweight='bold')
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

        ax[4, 0].set_title("Histogram of Nucleosome III", fontweight='bold')
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

    def get_linkers_position(self, signal_f_146_307, signal_r_146_307, signal_f_min_307, signal_r_min_307):
        smooth_signal_f_146_307 = savgol_filter(signal_f_146_307, window_length=51, polyorder=2)
        smooth_signal_r_146_307 = savgol_filter(signal_r_146_307, window_length=51, polyorder=2)
        smooth_signal_f_min_307 = savgol_filter(signal_f_min_307, window_length=51, polyorder=2)
        smooth_signal_r_min_307 = savgol_filter(signal_r_min_307, window_length=51, polyorder=2)

        position_f_1 = np.argmax(smooth_signal_f_146_307[:400])
        position_f_2 = np.argmax(smooth_signal_f_min_307[:position_f_1])

        position_r_1 = np.argmax(smooth_signal_r_146_307[600:]) + 600
        position_r_2 = np.argmax(smooth_signal_r_min_307[position_r_1:]) + position_r_1

        return position_f_1, position_r_1, position_f_2, position_r_2

    def update_axes_for_fragment_size_line(self, ax, x, x_ticks, start, end, signal_f, signal_r, p1, p2,
                                           p3, p4):
        max_signal = max(max(signal_f), max(signal_r))
        min_signal = min(min(signal_f), min(signal_r))
        ax.plot(x, signal_f, color='red', label='Forward')
        ax.plot(x, signal_r, color='green', label='Reverse')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 15))
        ax.tick_params(direction='out')
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(list(map(str, x_ticks)))
        ax.set_xlim(start, end)
        ax.set_yticks([min_signal, max_signal])
        ax.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax.set_ylim([min_signal, max_signal])
        ax.legend().set_visible(False)
        f_1, r_1 = sum(signal_f[:p1]) / sum(signal_r), sum(signal_r[:p1]) / sum(signal_r)
        f_2, r_2 = sum(signal_f[p1:p2]) / sum(signal_r), sum(signal_r[p1:p2]) / sum(signal_r)
        f_3, r_3 = sum(signal_f[p2:500]) / sum(signal_r), sum(signal_r[p2:500]) / sum(signal_r)
        f_4, r_4 = sum(signal_f[500:p3]) / sum(signal_r), sum(signal_r[500:p3]) / sum(signal_r)
        f_5, r_5 = sum(signal_f[p3:p4]) / sum(signal_r), sum(signal_r[p3:p4]) / sum(signal_r)
        f_6, r_6 = sum(signal_f[p4:]) / sum(signal_r), sum(signal_r[p4:]) / sum(signal_r)
        text_x_1 = ((p1 - 0) / 2.0 + 0) / 1000
        text_x_2 = ((p2 - p1) / 2.0 + p1) / 1000
        text_x_3 = ((500 - p2) / 2.0 + p2) / 1000
        text_x_4 = ((p3 - 500) / 2.0 + 500) / 1000
        text_x_5 = ((p4 - p3) / 2.0 + p3) / 1000
        text_x_6 = ((1000 - p4) / 2.0 + p4) / 1000
        ax.text(text_x_1, 1.0, str(round(f_1, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_1, 0.9, str(round(r_1, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_2, 1.0, str(round(f_2, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_2, 0.9, str(round(r_2, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_3, 1.0, str(round(f_3, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_3, 0.9, str(round(r_3, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_4, 1.0, str(round(f_4, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_4, 0.9, str(round(r_4, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_5, 1.0, str(round(f_5, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_5, 0.9, str(round(r_5, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_6, 1.0, str(round(f_6, 2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)
        ax.text(text_x_6, 0.9, str(round(r_6, 2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=12)

    def update_axes_for_signal_by_position(self, ax, x, x_ticks, start, end, signal_f_1, signal_r_1, p1, p2,
                                           p3, p4, title):

        signal_f = np.roll(signal_f_1, -4)
        signal_r = np.roll(signal_r_1, 4)
        max_signal = max(max(signal_f), max(signal_r))
        min_signal = min(min(signal_f), min(signal_r))
        ax.plot(x, signal_f, color='red', label='Forward')
        ax.plot(x, signal_r, color='green', label='Reverse')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 15))
        ax.tick_params(direction='out')
        ax.set_xticks(x_ticks)
        ax.set_xticklabels(list(map(str, x_ticks)))
        ax.set_xlim(start, end)
        ax.set_yticks([min_signal, max_signal])
        ax.set_yticklabels([str(int(min_signal)), str(int(max_signal))], rotation=90)
        ax.set_ylim([min_signal, max_signal])
        ax.legend().set_visible(False)
        ax.set_title(title, fontweight='bold')
        f_1, r_1 = sum(signal_f[:p1]), sum(signal_r[:p1])
        f_2, r_2 = sum(signal_f[p1:p2]), sum(signal_r[p1:p2])
        f_3, r_3 = sum(signal_f[p2:500]), sum(signal_r[p2:500])
        f_4, r_4 = sum(signal_f[500:p3]), sum(signal_r[500:p3])
        f_5, r_5 = sum(signal_f[p3:p4]), sum(signal_r[p3:p4])
        f_6, r_6 = sum(signal_f[p4:]), sum(signal_r[p4:])
        text_x_1 = ((p1 - 0) / 2.0 + 0) / 1000
        text_x_2 = ((p2 - p1) / 2.0 + p1) / 1000
        text_x_3 = ((500 - p2) / 2.0 + p2) / 1000
        text_x_4 = ((p3 - 500) / 2.0 + 500) / 1000
        text_x_5 = ((p4 - p3) / 2.0 + p3) / 1000
        text_x_6 = ((1000 - p4) / 2.0 + p4) / 1000
        ax.text(text_x_1, 1.0, str(int(f_1)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_1, 0.9, str(int(r_1)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_2, 1.0, str(int(f_2)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_2, 0.9, str(int(r_2)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_3, 1.0, str(int(f_3)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_3, 0.9, str(int(r_3)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_4, 1.0, str(int(f_4)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_4, 0.9, str(int(r_4)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_5, 1.0, str(int(f_5)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_5, 0.9, str(int(r_5)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_6, 1.0, str(int(f_6)), verticalalignment='center', color='red',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)
        ax.text(text_x_6, 0.9, str(int(r_6)), verticalalignment='center', color='green',
                horizontalalignment='center', transform=ax.transAxes, fontsize=8)

    def update_axes_for_density(self, ax, length_f, length_r, title):
        unique, counts = np.unique(length_f, return_counts=True)
        ax.hist(length_f, bins=unique, color='red', label='Forward', alpha=1)
        unique, counts = np.unique(length_r, return_counts=True)
        ax.hist(length_r, bins=unique, color='green', label='Reverse', alpha=0.5)
        ax.legend(loc="upper right", frameon=False)

    def rescaling(self, vector):
        maxN = max(vector)
        minN = min(vector)
        return [(e - minN) / (maxN - minN) for e in vector]
