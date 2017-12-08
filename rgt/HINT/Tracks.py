import os
from argparse import SUPPRESS
from numpy import int

# Internal
from ..Util import GenomeData, HmmData
from ..GenomicRegionSet import GenomicRegionSet
from biasTable import BiasTable
from signalProcessing import GenomicSignal


def tracks_args(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19, hg38. mm9, and mm10. DEFAULT: hg19")
    parser.add_argument("--reads-file", type=str, metavar="STRING", default=None,
                        help="A bam file containing all the DNase-seq or ATAC-seq reads. DEFAULT: None")
    parser.add_argument("--regions-file", type=str, metavar="STRING", default=None,
                        help="A bed file containing all the interested regions. DEFAULT: None")
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
    parser.add_argument("--raw-signal", action="store_true", default=False,
                        help="If set, the raw signals from DNase-seq or ATAC-seq data will be generated. DEFAULT: False")
    parser.add_argument("--bc-signal", action="store_true", default=False,
                        help="If set, the bias corrected signals from DNase-seq or ATAC-seq data will be generated. DEFAULT: False")
    parser.add_argument("--norm-signal", action="store_true", default=False,
                        help="If set, the normalised signals from DNase-seq or ATAC-seq data will be generated. DEFAULT: False")
    parser.add_argument("--bigWig", action="store_true", default=False,
                        help="If set, all .wig files will be converted to .bw files. DEFAULT: False")
    parser.add_argument("--strand-specific", action="store_true", default=False,
                        help=(
                        "If set, the tracks will be splitted into two files, one for forward and another for reverse strand. DEFAULT: False"))

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="tracks",
                        help="The prefix for results files. DEFAULT: tracks")


def tracks_run(args):
    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

    raw_signal_file = None
    bc_signal_file = None
    norm_signal_file = None
    # Output wig signal
    if args.raw_signal:
        raw_signal_file = os.path.join(args.output_location, "{}.raw.wig".format(args.output_prefix))
        open(raw_signal_file, "a").close()
    if args.bc_signal:
        bc_signal_file = os.path.join(args.output_location, "{}.bc.wig".format(args.output_prefix))
        open(bc_signal_file, "a").close()
    if args.norm_signal:
        norm_signal_file = os.path.join(args.output_location, "{}.norm.wig".format(args.output_prefix))
        open(norm_signal_file, "a").close()


    signal = GenomicSignal(args.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    regions = GenomicRegionSet("Interested regions")
    regions.read(args.regions_file)

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    if args.bias_table:
        bias_table_list = args.bias_table.split(",")
        bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                            table_file_name_R=bias_table_list[1])
    else:
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table = BiasTable().load_table(table_file_name_F=table_F,
                                            table_file_name_R=table_R)

    for region in regions:
        signal.print_signal(ref=region.chrom, start=region.initial, end=region.final,
                            downstream_ext=args.downstream_ext,
                            bias_table=bias_table,
                            upstream_ext=args.upstream_ext,
                            forward_shift=args.forward_shift,
                            reverse_shift=args.reverse_shift,
                            genome_file_name=genome_data.get_genome(),
                            raw_signal_file=raw_signal_file,
                            bc_signal_file=bc_signal_file,
                            norm_signal_file=norm_signal_file,
                            strand_specific=args.strand_specific)

    chrom_sizes_file = genome_data.get_chromosome_sizes()
    if args.bigWig:
        if args.raw_signal:
            bw_filename = os.path.join(args.output_location, "{}.raw.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", raw_signal_file, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(raw_signal_file)

        if args.bc_signal:
            bw_filename = os.path.join(args.output_location, "{}.bc.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", bc_signal_file, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(bc_signal_file)

            if args.strand_specific:
                bw_filename = os.path.join(args.output_location, "{}_Forward.bc.bw".format(args.output_prefix))
                wig_filename = os.path.join(args.output_location, "{}_Forward.bc.wig".format(args.output_prefix))
                os.system(" ".join(["wigToBigWig", wig_filename, chrom_sizes_file, bw_filename, "-verbose=0"]))
                os.remove(wig_filename)

                bw_filename = os.path.join(args.output_location, "{}_Reverse.bc.bw".format(args.output_prefix))
                wig_filename = os.path.join(args.output_location, "{}_Reverse.bc.wig".format(args.output_prefix))
                os.system(" ".join(["wigToBigWig", wig_filename, chrom_sizes_file, bw_filename, "-verbose=0"]))
                os.remove(wig_filename)

        if args.norm_signal:
            bw_filename = os.path.join(args.output_location, "{}.norm.bw".format(args.output_prefix))
            os.system(" ".join(["wigToBigWig", norm_signal_file, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(norm_signal_file)

            if args.strand_specific:
                bw_filename = os.path.join(args.output_location, "{}_Forward.norm.bw".format(args.output_prefix))
                wig_filename = os.path.join(args.output_location, "{}_Forward.norm.wig".format(args.output_prefix))
                os.system(" ".join(["wigToBigWig", wig_filename, chrom_sizes_file, bw_filename, "-verbose=0"]))
                os.remove(wig_filename)

                bw_filename = os.path.join(args.output_location, "{}_Reverse.norm.bw".format(args.output_prefix))
                wig_filename = os.path.join(args.output_location, "{}_Reverse.norm.wig".format(args.output_prefix))
                os.system(" ".join(["wigToBigWig", wig_filename, chrom_sizes_file, bw_filename, "-verbose=0"]))
                os.remove(wig_filename)
