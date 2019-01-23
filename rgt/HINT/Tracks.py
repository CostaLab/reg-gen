import os
from argparse import SUPPRESS
import numpy as np
from pysam import Samfile, Fastafile
from scipy.stats import scoreatpercentile

# Internal
from rgt.Util import GenomeData, HmmData, ErrorHandler
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.HINT.biasTable import BiasTable
from rgt.HINT.signalProcessing import GenomicSignal


def tracks_args(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19, hg38. mm9, and mm10. DEFAULT: hg19")
    parser.add_argument("--bias-table", type=str, metavar="FILE1_F,FILE1_R", default=None,
                        help="Bias table files used to generate bias corrected tracks. DEFAULT: None")

    # Hidden Options
    parser.add_argument("--initial-clip", type=int, metavar="INT", default=50, help=SUPPRESS)
    parser.add_argument("--downstream-ext", type=int, metavar="INT", default=1, help=SUPPRESS)
    parser.add_argument("--upstream-ext", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--forward-shift", type=int, metavar="INT", default=5, help=SUPPRESS)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=-4, help=SUPPRESS)
    parser.add_argument("--k-nb", type=int, metavar="INT", default=6, help=SUPPRESS)

    # Output Options
    parser.add_argument("--raw", action="store_true", default=False,
                        help="If set, the raw signals from DNase-seq or ATAC-seq data will be generated. DEFAULT: False")
    parser.add_argument("--bc", action="store_true", default=False,
                        help="If set, the bias corrected signals from DNase-seq or ATAC-seq data will be generated. "
                             "DEFAULT: False")
    parser.add_argument("--norm", action="store_true", default=False,
                        help="If set, the normalised signals from DNase-seq or ATAC-seq data will be generated. "
                             "DEFAULT: False")
    parser.add_argument("--bigWig", action="store_true", default=False,
                        help="If set, all .wig files will be converted to .bw files. DEFAULT: False")
    parser.add_argument("--strand-specific", action="store_true", default=False,
                        help="If set, the tracks will be splitted into two files, one for forward and another for "
                             "reverse strand. DEFAULT: False")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="tracks",
                        help="The prefix for results files. DEFAULT: tracks")

    parser.add_argument('input_files', metavar='reads.bam regions.bed', type=str, nargs='*',
                        help='BAM file of reads and BED files of interesting regions')


def tracks_run(args):
    if args.raw:
        get_raw_tracks(args)
    if args.bc:
        get_bc_tracks(args)


def get_raw_tracks(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    output_fname = os.path.join(args.output_location, "{}.wig".format(args.output_prefix))

    bam = Samfile(args.input_files[0], "rb")
    regions = GenomicRegionSet("Interested regions")
    regions.read(args.input_files[1])
    regions.merge()
    reads_file = GenomicSignal()

    with open(output_fname, "a") as output_f:
        for region in regions:
            # Raw counts
            signal = [0.0] * (region.final - region.initial)
            for read in bam.fetch(region.chrom, region.initial, region.final):
                if not read.is_reverse:
                    cut_site = read.pos + args.forward_shift
                    if region.initial <= cut_site < region.final:
                        signal[cut_site - region.initial] += 1.0
                else:
                    cut_site = read.aend + args.reverse_shift - 1
                    if region.initial <= cut_site < region.final:
                        signal[cut_site - region.initial] += 1.0

            if args.norm:
                signal = reads_file.boyle_norm(signal)
                perc = scoreatpercentile(signal, 98)
                std = np.std(signal)
                signal = reads_file.hon_norm_atac(signal, perc, std)

            output_f.write("fixedStep chrom=" + region.chrom + " start=" + str(region.initial + 1) + " step=1\n" +
                           "\n".join([str(e) for e in np.nan_to_num(signal)]) + "\n")
    output_f.close()

    if args.bigWig:
        genome_data = GenomeData(args.organism)
        chrom_sizes_file = genome_data.get_chromosome_sizes()
        bw_filename = os.path.join(args.output_location, "{}.bw".format(args.output_prefix))
        os.system(" ".join(["wigToBigWig", output_fname, chrom_sizes_file, bw_filename, "-verbose=0"]))
        os.remove(output_fname)


def get_bc_tracks(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    regions = GenomicRegionSet("Interested regions")
    regions.read(args.input_files[1])
    regions.merge()

    reads_file = GenomicSignal()

    bam = Samfile(args.input_files[0], "rb")
    genome_data = GenomeData(args.organism)
    fasta = Fastafile(genome_data.get_genome())

    hmm_data = HmmData()
    if args.bias_table:
        bias_table_list = args.bias_table.split(",")
        bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                            table_file_name_R=bias_table_list[1])
    else:
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table = BiasTable().load_table(table_file_name_F=table_F,
                                            table_file_name_R=table_R)

    if args.strand_specific:
        fname_forward = os.path.join(args.output_location, "{}_forward.wig".format(args.output_prefix))
        fname_reverse = os.path.join(args.output_location, "{}_reverse.wig".format(args.output_prefix))

        f_forward = open(fname_forward, "a")
        f_reverse = open(fname_reverse, "a")
        for region in regions:
            signal_f, signal_r = reads_file.get_bc_signal_by_fragment_length(
                ref=region.chrom, start=region.initial, end=region.final, bam=bam, fasta=fasta, bias_table=bias_table,
                forward_shift=args.forward_shift, reverse_shift=args.reverse_shift, min_length=None, max_length=None,
                strand=True)

            if args.norm:
                signal_f = reads_file.boyle_norm(signal_f)
                perc = scoreatpercentile(signal_f, 98)
                std = np.std(signal_f)
                signal_f = reads_file.hon_norm_atac(signal_f, perc, std)

                signal_r = reads_file.boyle_norm(signal_r)
                perc = scoreatpercentile(signal_r, 98)
                std = np.std(signal_r)
                signal_r = reads_file.hon_norm_atac(signal_r, perc, std)

            f_forward.write("fixedStep chrom=" + region.chrom + " start=" + str(region.initial + 1) + " step=1\n" +
                            "\n".join([str(e) for e in np.nan_to_num(signal_f)]) + "\n")

            f_reverse.write("fixedStep chrom=" + region.chrom + " start=" + str(region.initial + 1) + " step=1\n" +
                            "\n".join([str(-e) for e in np.nan_to_num(signal_r)]) + "\n")

        f_forward.close()
        f_reverse.close()

        if args.bigWig:
            genome_data = GenomeData(args.organism)
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            bw_filename = os.path.join(args.output_location, "{}_forward.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", fname_forward, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(fname_forward)

            bw_filename = os.path.join(args.output_location, "{}_reverse.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", fname_reverse, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(fname_reverse)

    else:
        output_fname = os.path.join(args.output_location, "{}.wig".format(args.output_prefix))
        with open(output_fname, "a") as output_f:
            for region in regions:
                signal = reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=region.initial,
                                                                     end=region.final,
                                                                     bam=bam, fasta=fasta, bias_table=bias_table,
                                                                     forward_shift=args.forward_shift,
                                                                     reverse_shift=args.reverse_shift,
                                                                     min_length=None, max_length=None, strand=False)

                if args.norm:
                    signal = reads_file.boyle_norm(signal)
                    perc = scoreatpercentile(signal, 98)
                    std = np.std(signal)
                    signal = reads_file.hon_norm_atac(signal, perc, std)

                output_f.write("fixedStep chrom=" + region.chrom + " start=" + str(region.initial + 1) + " step=1\n" +
                               "\n".join([str(e) for e in np.nan_to_num(signal)]) + "\n")
        output_f.close()

        if args.bigWig:
            genome_data = GenomeData(args.organism)
            chrom_sizes_file = genome_data.get_chromosome_sizes()
            bw_filename = os.path.join(args.output_location, "{}.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", output_fname, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(output_fname)
