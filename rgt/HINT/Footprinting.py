

import os
from copy import deepcopy
from argparse import SUPPRESS

# Internal
from rgt.Util import ErrorHandler, HmmData, GenomeData, OverlapType
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.HINT.signalProcessing import GenomicSignal
from rgt.HINT.hmm import HMM, _compute_log_likelihood
from rgt.HINT.biasTable import BiasTable

# External
import types
import pysam
from numpy import array, sum, isnan
from hmmlearn.hmm import GaussianHMM
import joblib

# Test
import numpy as np
from pysam import Fastafile, Samfile
from scipy.stats import scoreatpercentile


def footprinting_args(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19, hg38. mm9, and mm10. DEFAULT: hg19")
    parser.add_argument("--hmm-file", type=str, metavar="FILE", default=None,
                        help="If the argument is not given, then a default HMM will be used.")
    parser.add_argument("--bias-table", type=str, metavar="FILE_F,FILE_R", default=None,
                        help=("List of files with all possible k-mers (for any k) and their bias estimates. "
                              "Each line should contain a kmer and the bias estimate separated by tab."))

    parser.add_argument("--paired-end", action="store_true", default=False,
                        help="Set it if your ATAC-seq data is paired-end sequenced. "
                             "Note that this option is only applied to ATAC-seq data. DEFAULT: False")
    parser.add_argument("--bias-correction", action="store_true", default=False,
                        help="If set, footprint calling will based on bias corrected DNase-seq signal. "
                             "This option is only applied to DNase-seq. DEFAULT: False")
    parser.add_argument("--bias-type", dest="bias_type", type=str, metavar="STRING", default="SH",
                        help=("Type of protocol used to generate the DNase-seq. "
                              "Available options are: 'SH' (DNase-seq single-hit protocol), 'DH' "
                              "(DNase-seq double-hit protocol). DEFAULT: SH"))

    #Hidden Options
    parser.add_argument("--initial-clip", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--sg-window-size", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--norm-per", type=float, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--slope-per", type=float, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--downstream-ext", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--upstream-ext", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--forward-shift", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=None, help=SUPPRESS)

    parser.add_argument("--region-total-ext", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-max-size", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-min-size", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-ext", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--tc-ext", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-limit", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-state", type=int, metavar="INT", default=None, help=SUPPRESS)
    parser.add_argument("--fp-bed-fname", type=str, metavar="STRING", default=None, help=SUPPRESS)

    parser.add_argument("--dnase-initial-clip", type=int, metavar="INT", default=1000, help=SUPPRESS)
    parser.add_argument("--dnase-sg-window-size", type=int, metavar="INT", default=9, help=SUPPRESS)
    parser.add_argument("--dnase-norm-per", type=float, metavar="INT", default=98, help=SUPPRESS)
    parser.add_argument("--dnase-slope-per", type=float, metavar="INT", default=98, help=SUPPRESS)
    parser.add_argument("--dnase-downstream-ext", type=int, metavar="INT", default=1, help=SUPPRESS)
    parser.add_argument("--dnase-upstream-ext", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--dnase-forward-shift", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--dnase-reverse-shift", type=int, metavar="INT", default=0, help=SUPPRESS)

    parser.add_argument("--histone-initial-clip", type=int, metavar="INT", default=1000, help=SUPPRESS)
    parser.add_argument("--histone-sg-window-size", type=int, metavar="INT", default=201, help=SUPPRESS)
    parser.add_argument("--histone-norm-per", type=float, metavar="INT", default=98, help=SUPPRESS)
    parser.add_argument("--histone-slope-per", type=float, metavar="INT", default=98, help=SUPPRESS)
    parser.add_argument("--histone-downstream-ext", type=int, metavar="INT", default=200, help=SUPPRESS)
    parser.add_argument("--histone-upstream-ext", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--histone-forward-shift", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--histone-reverse-shift", type=int, metavar="INT", default=0, help=SUPPRESS)

    parser.add_argument("--model", type=str, metavar="STRING", default=None, help=SUPPRESS)
    parser.add_argument("--unstrand-specific", action="store_true", default=False, help=SUPPRESS)

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="footprints",
                        help="The prefix for results files. DEFAULT: footprints")

    parser.add_argument("--atac-seq", default=False, action='store_true',
                        help="If set, footprint calling will based on ATAC-seq model. DEFAULT: False")
    parser.add_argument("--dnase-seq", default=False, action='store_true',
                        help="If set, footprint calling will based on DNase-seq model. DEFAULT: False")
    parser.add_argument("--histone", default=False, action='store_true',
                        help="If set, footprint calling will based on histone modification model. DEFAULT: False")
    parser.add_argument("--dnase-histone", default=False, action='store_true', help=SUPPRESS)

    parser.add_argument('input_files', metavar='reads.bam regions.bed', type=str, nargs='*',
                        help='BAM file of reads and BED files of interesting regions')


def footprinting_run(args):
    if args.atac_seq:
        atac_seq(args)
    if args.dnase_seq:
        dnase_seq(args)
    if args.histone:
        histone(args)
    if args.dnase_histone:
        dnase_histone(args)


def atac_seq(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    # Check if the index file exists
    base_name = "{}.bai".format(args.input_files[0])
    if not os.path.exists(base_name):
        print("The index file of {} doesn't exist, generating with Pysam".format(args.input_files[0]))
        pysam.index(args.input_files[0])

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Load parameters for ATAC-seq
    ###################################################################################################

    if args.hmm_file:
        hmm_file = args.hmm_file
    else:
        if args.paired_end:
            hmm_file = hmm_data.get_default_hmm_atac_paired()
        else:
            hmm_file = hmm_data.get_default_hmm_atac_single()

    # FIXME: make a clear-text file, like for dnase, and remove this pickle business
    import pickle
    with open(hmm_file, "rb") as f:
        hmm = pickle.load(f, encoding="latin1")

    hmm._compute_log_likelihood = types.MethodType(_compute_log_likelihood, hmm)

    if args.bias_table:
        bias_table_list = args.bias_table.split(",")
        bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                            table_file_name_R=bias_table_list[1])
    else:
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table = BiasTable().load_table(table_file_name_F=table_F,
                                            table_file_name_R=table_R)

    if args.initial_clip is None:
        initial_clip = 50
    else:
        initial_clip = args.initial_clip

    if args.sg_window_size is None:
        sg_window_size = 9
    else:
        sg_window_size = args.sg_window_size

    if args.norm_per is None:
        norm_per = 98
    else:
        norm_per = args.norm_per

    if args.slope_per is None:
        slope_per = 98
    else:
        slope_per = args.slope_per

    if args.downstream_ext is None:
        downstream_ext = 1
    else:
        downstream_ext = args.downstream_ext

    if args.upstream_ext is None:
        upstream_ext = 0
    else:
        upstream_ext = args.upstream_ext

    if args.forward_shift is None:
        forward_shift = 5
    else:
        forward_shift = args.forward_shift

    if args.reverse_shift is None:
        reverse_shift = -5
    else:
        reverse_shift = args.reverse_shift

    if args.region_total_ext is None:
        region_total_ext = 0
    else:
        region_total_ext = args.region_total_ext

    if args.fp_max_size is None:
        fp_max_size = 50
    else:
        fp_max_size = args.fp_max_size

    if args.fp_min_size is None:
        fp_min_size = 5
    else:
        fp_min_size = args.fp_min_size

    if args.fp_ext is None:
        fp_ext = 5
    else:
        fp_ext = args.fp_ext

    if args.tc_ext is None:
        tc_ext = 100
    else:
        tc_ext = args.tc_ext

    if args.paired_end:
        if args.fp_state is None:
            fp_state = 8
        else:
            fp_state = args.fp_state
    else:
        if args.fp_state is None:
            fp_state = 6
        else:
            fp_state = args.fp_state

    # Initializing result set
    footprints = GenomicRegionSet(args.output_prefix)

    reads_file = GenomicSignal(args.input_files[0])
    reads_file.load_sg_coefs(sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read(args.input_files[1])

    regions = deepcopy(original_regions)
    regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
    regions.merge()

    bam = Samfile(args.input_files[0], "rb")
    fasta = Fastafile(genome_data.get_genome())

    if args.paired_end:
        for region in original_regions:
            input_sequence = list()

            try:
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=region.initial, end=region.final,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_min_145, signal_bc_r_min_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=region.initial, end=region.final,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=None, strand=True)


            except Exception:
                err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([region.chrom, str(region.initial), str(
                    region.final)]) + "). This iteration will be skipped.")
                continue

            signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
            perc = scoreatpercentile(signal_bc_f_max_145, 98)
            std = np.array(signal_bc_f_max_145).std()
            signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
            signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

            signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
            perc = scoreatpercentile(signal_bc_r_max_145, 98)
            std = np.array(signal_bc_r_max_145).std()
            signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
            signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

            signal_bc_f_min_145 = reads_file.boyle_norm(signal_bc_f_min_145)
            perc = scoreatpercentile(signal_bc_f_min_145, 98)
            std = np.array(signal_bc_f_min_145).std()
            signal_bc_f_min_145 = reads_file.hon_norm_atac(signal_bc_f_min_145, perc, std)
            signal_bc_f_min_145_slope = reads_file.slope(signal_bc_f_min_145, reads_file.sg_coefs)

            signal_bc_r_min_145 = reads_file.boyle_norm(signal_bc_r_min_145)
            perc = scoreatpercentile(signal_bc_r_min_145, 98)
            std = np.array(signal_bc_r_min_145).std()
            signal_bc_r_min_145 = reads_file.hon_norm_atac(signal_bc_r_min_145, perc, std)
            signal_bc_r_min_145_slope = reads_file.slope(signal_bc_r_min_145, reads_file.sg_coefs)

            input_sequence.append(signal_bc_f_max_145)
            input_sequence.append(signal_bc_f_max_145_slope)
            input_sequence.append(signal_bc_r_max_145)
            input_sequence.append(signal_bc_r_max_145_slope)

            input_sequence.append(signal_bc_f_min_145)
            input_sequence.append(signal_bc_f_min_145_slope)
            input_sequence.append(signal_bc_r_min_145)
            input_sequence.append(signal_bc_r_min_145_slope)

            # Applying HMM
            try:
                posterior_list = hmm.predict(np.array(input_sequence).T)
            except Exception:
                err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([region.chrom, str(region.initial), str(
                    region.final)]) + "). This iteration will be skipped.")
                continue

            if args.fp_bed_fname is not None:
                output_bed_file(region.chrom, region.initial, region.final, posterior_list, args.fp_bed_fname, fp_state)

            # Formatting results
            start_pos = 0
            flag_start = False
            for k in range(region.initial, region.initial + len(posterior_list)):
                curr_index = k - region.initial
                if flag_start:
                    if posterior_list[curr_index] != fp_state:
                        if k - start_pos < fp_max_size:
                            fp = GenomicRegion(region.chrom, start_pos, k)
                            footprints.add(fp)
                        flag_start = False
                else:
                    if posterior_list[curr_index] == fp_state:
                        flag_start = True
                        start_pos = k
            if flag_start:
                if region.initial + len(posterior_list) - start_pos < fp_max_size:
                    fp = GenomicRegion(region.chrom, start_pos, region.final)
                    footprints.add(fp)
    else:
        for region in original_regions:

            input_sequence = list()
            atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r = \
                reads_file.get_signal_atac(region.chrom, region.initial, region.final, downstream_ext,
                                           upstream_ext, forward_shift, reverse_shift,
                                           initial_clip, norm_per, slope_per,
                                           bias_table, genome_data.get_genome())

            input_sequence.append(atac_norm_f)
            input_sequence.append(atac_slope_f)
            input_sequence.append(atac_norm_r)
            input_sequence.append(atac_slope_r)

            try:
                posterior_list = hmm.predict(np.array(input_sequence).T)
            except Exception:
                err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([region.chrom, str(region.initial), str(
                    region.final)]) + "). This iteration will be skipped.")
                continue

            if args.fp_bed_fname is not None:
                output_bed_file(region.chrom, region.initial, region.final, posterior_list, args.fp_bed_fname, fp_state)
            # Formatting results
            start_pos = 0
            flag_start = False
            for k in range(region.initial, region.initial + len(posterior_list)):
                curr_index = k - region.initial
                if flag_start:
                    if posterior_list[curr_index] != fp_state:
                        if k - start_pos < fp_max_size:
                            fp = GenomicRegion(region.chrom, start_pos, k)
                            footprints.add(fp)
                        flag_start = False
                else:
                    if posterior_list[curr_index] == fp_state:
                        flag_start = True
                        start_pos = k
            if flag_start:
                if region.initial + len(posterior_list) - start_pos < fp_max_size:
                    fp = GenomicRegion(region.chrom, start_pos, region.final)
                    footprints.add(fp)

    ###################################################################################################
    # Post-processing
    ###################################################################################################

    post_processing(footprints=footprints, original_regions=original_regions, fp_min_size=fp_min_size,
                    fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                    reads_file=reads_file, downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                    forward_shift=forward_shift, reverse_shift=reverse_shift,
                    initial_clip=initial_clip, output_location=args.output_location,
                    output_prefix=args.output_prefix)


def dnase_seq(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    bias_table = None
    if args.bias_correction:
        if args.bias_table:
            bias_table_list = args.bias_table.split(",")
            bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                                table_file_name_R=bias_table_list[1])
        else:
            if args.bias_type == 'SH':
                table_F = hmm_data.get_default_bias_table_F_SH()
                table_R = hmm_data.get_default_bias_table_R_SH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)
            elif args.bias_type == 'DH':
                table_F = hmm_data.get_default_bias_table_F_DH()
                table_R = hmm_data.get_default_bias_table_R_DH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    hmm_file = None
    if args.hmm_file:
        hmm_file = args.hmm_file
    else:
        if args.bias_correction:
            hmm_file = hmm_data.get_default_hmm_dnase_bc()
        else:
            hmm_file = hmm_data.get_default_hmm_dnase()

    scikit_hmm = None
    try:
        hmm_scaffold = HMM()
        hmm_scaffold.load_hmm(hmm_file)
        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
        scikit_hmm._compute_log_likelihood = types.MethodType(_compute_log_likelihood, scikit_hmm)
        scikit_hmm.startprob_ = array(hmm_scaffold.pi)
        scikit_hmm.transmat_ = array(hmm_scaffold.A)
        scikit_hmm.means_ = array(hmm_scaffold.means)
        scikit_hmm.covars_ = array(hmm_scaffold.covs)
    except Exception:
        err.throw_error("FP_HMM_FILES")

    initial_clip = 1000 if not args.initial_clip else args.initial_clip
    sg_window_size = 9 if not args.sg_window_size else args.sg_window_size
    norm_per = 98 if not args.norm_per else args.norm_per
    slope_per = 98 if not args.slope_per else args.slope_per
    downstream_ext = 1 if not args.downstream_ext else args.downstream_ext
    upstream_ext = 0 if not args.upstream_ext else args.upstream_ext
    region_total_ext = 10000 if not args.region_total_ext else args.region_total_ext
    forward_shift = 0 if not args.forward_shift else args.forward_shift
    reverse_shift = 0 if not args.reverse_shift else args.reverse_shift
    fp_max_size = 50 if not args.fp_max_size else args.fp_max_size
    fp_min_size = 10 if not args.fp_min_size else args.fp_min_size
    fp_ext = 5 if not args.fp_ext else args.fp_ext
    tc_ext = 100 if not args.tc_ext else args.tc_ext

    # Initializing result set
    footprints = GenomicRegionSet(args.output_prefix)

    reads_file = GenomicSignal(args.input_files[0])
    reads_file.load_sg_coefs(sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read(args.input_files[1])

    regions = deepcopy(original_regions)
    regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        dnase_norm, dnase_slope = reads_file.get_signal(r.chrom, r.initial, r.final, downstream_ext,
                                                        upstream_ext, forward_shift,
                                                        reverse_shift, initial_clip, norm_per,
                                                        slope_per, bias_table, genome_data.get_genome())

        try:
            input_sequence = array([dnase_norm, dnase_slope]).T
        except Exception:
            err.throw_warning("FP_SEQ_FORMAT", add_msg="for region (" + ",".join([r.chrom, str(r.initial), str(
                r.final)]) + "). This iteration will be skipped.")
            continue

        # Applying HMM
        try:
            posterior_list = scikit_hmm.predict(input_sequence)
        except Exception:
            err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([r.chrom, str(r.initial), str(
                r.final)]) + "). This iteration will be skipped.")
            continue

        # Formatting results
        start_pos = 0
        flag_start = False
        fp_state_nb = 4
        for k in range(r.initial, r.initial + len(posterior_list)):
            curr_index = k - r.initial
            if flag_start:
                if posterior_list[curr_index] != fp_state_nb:
                    if k - start_pos < fp_max_size:
                        fp = GenomicRegion(r.chrom, start_pos, k)
                        footprints.add(fp)
                    flag_start = False
            else:
                if posterior_list[curr_index] == fp_state_nb:
                    flag_start = True
                    start_pos = k
        if flag_start:
            fp = GenomicRegion(r.chrom, start_pos, r.final)
            footprints.add(fp)

    ###################################################################################################
    # Post-processing
    ###################################################################################################
    post_processing(footprints=footprints, original_regions=original_regions, fp_min_size=fp_min_size,
                    fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                    reads_file=reads_file, downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                    forward_shift=forward_shift, reverse_shift=reverse_shift,
                    initial_clip=initial_clip, output_location=args.output_location,
                    output_prefix=args.output_prefix)


def histone(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    hmm_file = None
    if args.hmm_file:
        hmm_file = args.hmm_file
    else:
        hmm_file = hmm_data.get_default_hmm_histone()

    scikit_hmm = None
    try:
        hmm_scaffold = HMM()
        hmm_scaffold.load_hmm(hmm_file)
        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
        scikit_hmm.startprob_ = array(hmm_scaffold.pi)
        scikit_hmm.transmat_ = array(hmm_scaffold.A)
        scikit_hmm.means_ = array(hmm_scaffold.means)
        scikit_hmm.covars_ = array(hmm_scaffold.covs)
        scikit_hmm._compute_log_likelihood = types.MethodType(_compute_log_likelihood, scikit_hmm)
    except Exception:
        err.throw_error("FP_HMM_FILES")

    initial_clip = 1000 if not args.initial_clip else args.initial_clip
    sg_window_size = 201 if not args.sg_window_size else args.sg_window_size
    norm_per = 98 if not args.norm_per else args.norm_per
    slope_per = 98 if not args.slope_per else args.slope_per
    downstream_ext = 200 if not args.downstream_ext else args.downstream_ext
    upstream_ext = 0 if not args.upstream_ext else args.upstream_ext
    region_total_ext = 10000 if not args.region_total_ext else args.region_total_ext
    forward_shift = 0 if not args.forward_shift else args.forward_shift
    reverse_shift = 0 if not args.reverse_shift else args.reverse_shift
    fp_max_size = 2000 if not args.fp_max_size else args.fp_max_size
    fp_min_size = 200 if not args.fp_min_size else args.fp_min_size
    fp_ext = 50 if not args.fp_ext else args.fp_ext
    tc_ext = 500 if not args.tc_ext else args.tc_ext

    # Initializing result set
    footprints = GenomicRegionSet(args.output_prefix)

    reads_file = GenomicSignal(args.input_files[0])
    reads_file.load_sg_coefs(sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read(args.input_files[1])

    regions = deepcopy(original_regions)
    regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        histone_norm, histone_slope = reads_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                            downstream_ext=downstream_ext,
                                                            upstream_ext=upstream_ext,
                                                            forward_shift=forward_shift,
                                                            reverse_shift=reverse_shift,
                                                            initial_clip=initial_clip,
                                                            per_norm=norm_per,
                                                            per_slope=slope_per,
                                                            genome_file_name=genome_data.get_genome())

        try:
            input_sequence = array([histone_norm, histone_slope]).T
        except Exception:
            err.throw_warning("FP_SEQ_FORMAT", add_msg="for region (" + ",".join([r.chrom, str(r.initial), str(
                r.final)]) + "). This iteration will be skipped.")
            continue

        # Applying HMM
        try:
            posterior_list = scikit_hmm.predict(input_sequence)
        except Exception:
            err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([r.chrom, str(r.initial), str(
                r.final)]) + "). This iteration will be skipped.")
            continue

        # Formatting results
        start_pos = 0
        flag_start = False
        fp_state_nb = 4
        for k in range(r.initial, r.initial + len(posterior_list)):
            curr_index = k - r.initial
            if flag_start:
                if posterior_list[curr_index] != fp_state_nb:
                    if k - start_pos < fp_max_size:
                        fp = GenomicRegion(r.chrom, start_pos, k)
                        footprints.add(fp)
                    flag_start = False
            else:
                if posterior_list[curr_index] == fp_state_nb:
                    flag_start = True
                    start_pos = k
        if flag_start:
            fp = GenomicRegion(r.chrom, start_pos, r.final)
            footprints.add(fp)

    ###################################################################################################
    # Post-processing
    ###################################################################################################
    post_processing(footprints=footprints, original_regions=original_regions, fp_min_size=fp_min_size,
                    fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                    reads_file=reads_file, downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                    forward_shift=forward_shift, reverse_shift=reverse_shift,
                    initial_clip=initial_clip, output_location=args.output_location,
                    output_prefix=args.output_prefix)


def dnase_histone(args):
    # Initializing Error Handler
    err = ErrorHandler()

    if len(args.input_files) < 3:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify DNase and histone reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    bias_table = None
    if args.bias_correction:
        if args.bias_table:
            bias_table_list = args.bias_table.split(",")
            bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                                table_file_name_R=bias_table_list[1])
        else:
            if args.bias_type == 'SH':
                table_F = hmm_data.get_default_bias_table_F_SH()
                table_R = hmm_data.get_default_bias_table_R_SH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)
            elif args.bias_type == 'DH':
                table_F = hmm_data.get_default_bias_table_F_DH()
                table_R = hmm_data.get_default_bias_table_R_DH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    if args.hmm_file:
        hmm_file = args.hmm_file
    else:
        if args.bias_correction:
            hmm_file = hmm_data.get_default_hmm_dnase_histone_bc()
        else:
            hmm_file = hmm_data.get_default_hmm_dnase_histone()
    scikit_hmm = None
    try:
        hmm_scaffold = HMM()
        hmm_scaffold.load_hmm(hmm_file)
        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
        scikit_hmm.startprob_ = array(hmm_scaffold.pi)
        scikit_hmm.transmat_ = array(hmm_scaffold.A)
        scikit_hmm.means_ = array(hmm_scaffold.means)
        scikit_hmm.covars_ = array(hmm_scaffold.covs)
        scikit_hmm._compute_log_likelihood = types.MethodType(_compute_log_likelihood, scikit_hmm)
    except Exception:
        err.throw_error("FP_HMM_FILES")


    region_total_ext = 10000 if not args.region_total_ext else args.region_total_ext
    fp_max_size = 50 if not args.fp_max_size else args.fp_max_size
    fp_min_size = 10 if not args.fp_min_size else args.fp_min_size
    fp_ext = 5 if not args.fp_ext else args.fp_ext
    tc_ext = 100 if not args.tc_ext else args.tc_ext

    # Initializing result set
    footprints = GenomicRegionSet(args.output_prefix)

    dnase_reads_file = GenomicSignal(args.input_files[0])
    dnase_reads_file.load_sg_coefs(args.dnase_sg_window_size)

    histone_reads_file_list = list()
    for f in args.input_files[1:-1]:
        histone_reads_file_list.append(GenomicSignal(f))

    original_regions = GenomicRegionSet("regions")
    original_regions.read(args.input_files[-1])

    regions = deepcopy(original_regions)
    regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        dnase_norm, dnase_slope = dnase_reads_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                              downstream_ext=args.dnase_downstream_ext,
                                                              upstream_ext=args.dnase_upstream_ext,
                                                              forward_shift=args.dnase_forward_shift,
                                                              reverse_shift=args.dnase_reverse_shift,
                                                              initial_clip=args.dnase_initial_clip,
                                                              per_norm=args.dnase_norm_per,
                                                              per_slope=args.dnase_slope_per,
                                                              bias_table=bias_table,
                                                              genome_file_name=genome_data.get_genome())
        # Iterating over histone modifications
        for i in range(0, len(histone_reads_file_list)):
            # Fetching histone signal
            histone_file = None
            try:
                histone_file = histone_reads_file_list[i]
                histone_file.load_sg_coefs(args.histone_sg_window_size)
                histone_norm, histone_slope = histone_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                                      downstream_ext=args.histone_downstream_ext,
                                                                      upstream_ext=args.histone_upstream_ext,
                                                                      forward_shift=args.histone_forward_shift,
                                                                      reverse_shift=args.histone_reverse_shift,
                                                                      initial_clip=args.histone_initial_clip,
                                                                      per_norm=args.histone_norm_per,
                                                                      per_slope=args.histone_slope_per,
                                                                      genome_file_name=genome_data.get_genome())
            except Exception:
                err.throw_warning("FP_HISTONE_PROC", add_msg="for region (" + ",".join([r.chrom, str(r.initial),
                                str(r.final)]) + ") and histone modification " + histone_file.file_name +
                                ". This iteration will be skipped for this histone.")
                continue

            # Formatting sequence
            input_sequence = array([dnase_norm, dnase_slope, histone_norm, histone_slope]).T
            if isnan(sum(input_sequence)): continue  # Handling NAN's in signal / hmmlearn throws error TODO ERROR
            try:
                posterior_list = scikit_hmm.predict(input_sequence)
            except Exception:
                err.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([r.chrom, str(r.initial), str(
                    r.final)]) + ") and histone modification " + histone_file.file_name + ". This iteration will be skipped.")
                continue

            # Formatting results
            start_pos = 0
            flag_start = False
            fp_state_nb = 7
            for k in range(r.initial, r.initial + len(posterior_list)):
                curr_index = k - r.initial
                if flag_start:
                    if posterior_list[curr_index] != fp_state_nb:
                        if k - start_pos < fp_max_size:
                            fp = GenomicRegion(r.chrom, start_pos, k)
                            footprints.add(fp)
                        flag_start = False
                else:
                    if posterior_list[curr_index] == fp_state_nb:
                        flag_start = True
                        start_pos = k
            if flag_start:
                fp = GenomicRegion(r.chrom, start_pos, r.final)
                footprints.add(fp)

    ###################################################################################################
    # Post-processing
    ###################################################################################################
    post_processing(footprints=footprints, original_regions=original_regions, fp_min_size=fp_min_size,
                    fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                    reads_file=dnase_reads_file, downstream_ext=args.dnase_downstream_ext,
                    upstream_ext=args.dnase_upstream_ext,
                    forward_shift=args.dnase_forward_shift, reverse_shift=args.dnase_reverse_shift,
                    initial_clip=args.dnase_initial_clip, output_location=args.output_location,
                    output_prefix=args.output_prefix)


def post_processing(footprints, original_regions, fp_min_size, fp_ext, genome_data, tc_ext, reads_file,
                    downstream_ext, upstream_ext, forward_shift, reverse_shift, initial_clip, output_location,
                    output_prefix):
    # Overlapping results with original regions
    footprints_overlap = footprints.intersect(original_regions, mode=OverlapType.ORIGINAL)

    # Sorting and Merging
    footprints_overlap.merge()

    # Extending footprints
    for f in footprints_overlap.sequences:
        if f.final - f.initial < fp_min_size:
            f.initial = max(0, f.initial - fp_ext)
            f.final = f.final + fp_ext

    # Sorting and Merging
    footprints_overlap.merge()

    # Fetching chromosome sizes
    chrom_sizes_file_name = genome_data.get_chromosome_sizes()
    chrom_sizes_file = open(chrom_sizes_file_name, "r")
    chrom_sizes_dict = dict()
    for chrom_sizes_entry_line in chrom_sizes_file:
        chrom_sizes_entry_vec = chrom_sizes_entry_line.strip().split("\t")
        chrom_sizes_dict[chrom_sizes_entry_vec[0]] = int(chrom_sizes_entry_vec[1])
    chrom_sizes_file.close()

    # Evaluating TC
    for f in footprints_overlap.sequences:
        mid = (f.initial + f.final) / 2
        p1 = max(mid - tc_ext, 0)
        p2 = min(mid + tc_ext, chrom_sizes_dict[f.chrom])
        try:
            tag_count = reads_file.get_tag_count(ref=f.chrom, start=p1, end=p2,
                                                 downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                                                 forward_shift=forward_shift, reverse_shift=reverse_shift,
                                                 initial_clip=initial_clip)
        except Exception:
            tag_count = 0
        f.data = str(int(tag_count))

    ###################################################################################################
    # Writing output
    ###################################################################################################
    output_file_name = os.path.join(output_location, "{}.bed".format(output_prefix))
    footprints_overlap.write(output_file_name)

    # the number of reads
    num_reads = pysam.AlignmentFile(reads_file.file_name).count(until_eof=True)

    # the number of peaks and tag count within peaks
    num_peaks = 0
    num_tc = 0
    for r in original_regions:
        num_peaks += 1
        try:
            tag_count = reads_file.get_tag_count(r.chrom, r.initial, r.final, downstream_ext, upstream_ext,
                                                 forward_shift, reverse_shift, initial_clip)
        except Exception:
            tag_count = 0
        num_tc += tag_count

    # the number of footprints
    num_fp = len(footprints_overlap)

    output_file_name = os.path.join(output_location, "{}.info".format(output_prefix))
    with open(output_file_name, "w") as f:
        f.write("Number of reads: " + str(num_reads) + "\n")
        f.write("Number of peaks: " + str(num_peaks) + "\n")
        f.write("Number of tag counts within peaks: " + str(num_tc) + "\n")
        f.write("Number of footprints: " + str(num_fp) + "\n")


def output_bed_file(chrom, start, end, states, output_fname, fp_state):
    state_dict = dict([(0, "0"), (1, "1"), (2, "2"), (3, "3"), (4, "4"), (5, "5"), (6, "6"), (7, "7"), (8, "8")])
    color_dict = dict([(0, "0,0,0"), (1, "102,0,51"), (2, "153,0,153"), (3, "102,0,204"), (4, "0,0,255"),
                       (5, "51,153,255"), (6, "102,255,255"), (7, "102,204,255"), (8, "204,255,255")])
    state_dict[int(fp_state)] = "FP"

    current_state = states[0]
    start_postion = start
    is_print = False
    with open(output_fname, "a") as bed_file:
        for i in range(len(states)):
            if states[i] != current_state:
                end_position = start + i
                is_print = True
            elif i == len(states) - 1:
                end_position = end
                is_print = True

            if is_print:
                bed_file.write(chrom + " " + str(start_postion) + " " + str(end_position) + " "
                               + state_dict[current_state] + " " + str(1000) + " . "
                               + str(start_postion) + " " + str(end_position) + " "
                               + color_dict[current_state] + "\n")
                start_postion = end_position
                current_state = states[i]
                is_print = False


def atac_test(args):
    # Initializing Error Handler
    err = ErrorHandler()


    if len(args.input_files) != 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(args.organism)
    hmm_data = HmmData()

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

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    hmm = joblib.load(args.hmm_file)
    hmm._compute_log_likelihood = types.MethodType(_compute_log_likelihood, hmm)

    initial_clip = 50 if not args.initial_clip else args.initial_clip
    sg_window_size = 9 if not args.sg_window_size else args.sg_window_size
    norm_per = 98 if not args.norm_per else args.norm_per
    slope_per = 98 if not args.slope_per else args.slope_per
    downstream_ext = 1 if not args.downstream_ext else args.downstream_ext
    upstream_ext = 0 if not args.upstream_ext else args.upstream_ext
    region_total_ext = 0 if not args.region_total_ext else args.region_total_ext
    forward_shift = 5 if not args.forward_shift else args.forward_shift
    reverse_shift = -4 if not args.reverse_shift else args.reverse_shift
    fp_max_size = 50 if not args.fp_max_size else args.fp_max_size
    fp_min_size = 5 if not args.fp_min_size else args.fp_min_size
    fp_ext = 5 if not args.fp_ext else args.fp_ext
    tc_ext = 100 if not args.tc_ext else args.tc_ext
    fp_state = 8 if not args.fp_state else args.fp_state




    # Initializing result set
    footprints = GenomicRegionSet(args.output_prefix)

    reads_file = GenomicSignal(args.input_files[0])
    reads_file.load_sg_coefs(sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read(args.input_files[1])

    regions = deepcopy(original_regions)
    regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
    regions.merge()

    bam = Samfile(args.input_files[0], "rb")
    fasta = Fastafile(genome_data.get_genome())

    if not args.unstrand_specific:
        for region in regions:
            p1 = region.initial
            p2 = region.final

            input_sequence = list()
            if args.model == "a":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)

                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)
            elif args.model == "b":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_min_145, signal_bc_r_min_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=None, strand=True)

                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                signal_bc_f_min_145 = reads_file.boyle_norm(signal_bc_f_min_145)
                perc = scoreatpercentile(signal_bc_f_min_145, 98)
                std = np.array(signal_bc_f_min_145).std()
                signal_bc_f_min_145 = reads_file.hon_norm_atac(signal_bc_f_min_145, perc, std)
                signal_bc_f_min_145_slope = reads_file.slope(signal_bc_f_min_145, reads_file.sg_coefs)

                signal_bc_r_min_145 = reads_file.boyle_norm(signal_bc_r_min_145)
                perc = scoreatpercentile(signal_bc_r_min_145, 98)
                std = np.array(signal_bc_r_min_145).std()
                signal_bc_r_min_145 = reads_file.hon_norm_atac(signal_bc_r_min_145, perc, std)
                signal_bc_r_min_145_slope = reads_file.slope(signal_bc_r_min_145, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)

                input_sequence.append(signal_bc_f_min_145)
                input_sequence.append(signal_bc_f_min_145_slope)
                input_sequence.append(signal_bc_r_min_145)
                input_sequence.append(signal_bc_r_min_145_slope)
            elif args.model == "c":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_145_307, signal_bc_r_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=True)

                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                signal_bc_f_145_307 = reads_file.boyle_norm(signal_bc_f_145_307)
                perc = scoreatpercentile(signal_bc_f_145_307, 98)
                std = np.array(signal_bc_f_145_307).std()
                signal_bc_f_145_307 = reads_file.hon_norm_atac(signal_bc_f_145_307, perc, std)
                signal_bc_f_145_307_slope = reads_file.slope(signal_bc_f_145_307, reads_file.sg_coefs)

                signal_bc_r_145_307 = reads_file.boyle_norm(signal_bc_r_145_307)
                perc = scoreatpercentile(signal_bc_r_145_307, 98)
                std = np.array(signal_bc_r_145_307).std()
                signal_bc_r_145_307 = reads_file.hon_norm_atac(signal_bc_r_145_307, perc, std)
                signal_bc_r_145_307_slope = reads_file.slope(signal_bc_r_145_307, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)

                input_sequence.append(signal_bc_f_145_307)
                input_sequence.append(signal_bc_f_145_307_slope)
                input_sequence.append(signal_bc_r_145_307)
                input_sequence.append(signal_bc_r_145_307_slope)
            elif args.model == "d":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_145_307, signal_bc_r_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=True)
                signal_bc_f_min_307, signal_bc_r_min_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=None, strand=True)
                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                signal_bc_f_145_307 = reads_file.boyle_norm(signal_bc_f_145_307)
                perc = scoreatpercentile(signal_bc_f_145_307, 98)
                std = np.array(signal_bc_f_145_307).std()
                signal_bc_f_145_307 = reads_file.hon_norm_atac(signal_bc_f_145_307, perc, std)
                signal_bc_f_145_307_slope = reads_file.slope(signal_bc_f_145_307, reads_file.sg_coefs)

                signal_bc_r_145_307 = reads_file.boyle_norm(signal_bc_r_145_307)
                perc = scoreatpercentile(signal_bc_r_145_307, 98)
                std = np.array(signal_bc_r_145_307).std()
                signal_bc_r_145_307 = reads_file.hon_norm_atac(signal_bc_r_145_307, perc, std)
                signal_bc_r_145_307_slope = reads_file.slope(signal_bc_r_145_307, reads_file.sg_coefs)

                signal_bc_f_min_307 = reads_file.boyle_norm(signal_bc_f_min_307)
                perc = scoreatpercentile(signal_bc_f_min_307, 98)
                std = np.array(signal_bc_f_min_307).std()
                signal_bc_f_min_307 = reads_file.hon_norm_atac(signal_bc_f_min_307, perc, std)
                signal_bc_f_min_307_slope = reads_file.slope(signal_bc_f_min_307, reads_file.sg_coefs)

                signal_bc_r_min_307 = reads_file.boyle_norm(signal_bc_r_min_307)
                perc = scoreatpercentile(signal_bc_r_min_307, 98)
                std = np.array(signal_bc_r_min_307).std()
                signal_bc_r_min_307 = reads_file.hon_norm_atac(signal_bc_r_min_307, perc, std)
                signal_bc_r_min_307_slope = reads_file.slope(signal_bc_r_min_307, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)

                input_sequence.append(signal_bc_f_145_307)
                input_sequence.append(signal_bc_f_145_307_slope)
                input_sequence.append(signal_bc_r_145_307)
                input_sequence.append(signal_bc_r_145_307_slope)

                input_sequence.append(signal_bc_f_min_307)
                input_sequence.append(signal_bc_f_min_307_slope)
                input_sequence.append(signal_bc_r_min_307)
                input_sequence.append(signal_bc_r_min_307_slope)
            elif args.model == "e":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_145_307, signal_bc_r_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=True)
                signal_bc_f_307_500, signal_bc_r_307_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=500, strand=True)
                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                signal_bc_f_145_307 = reads_file.boyle_norm(signal_bc_f_145_307)
                perc = scoreatpercentile(signal_bc_f_145_307, 98)
                std = np.array(signal_bc_f_145_307).std()
                signal_bc_f_145_307 = reads_file.hon_norm_atac(signal_bc_f_145_307, perc, std)
                signal_bc_f_145_307_slope = reads_file.slope(signal_bc_f_145_307, reads_file.sg_coefs)

                signal_bc_r_145_307 = reads_file.boyle_norm(signal_bc_r_145_307)
                perc = scoreatpercentile(signal_bc_r_145_307, 98)
                std = np.array(signal_bc_r_145_307).std()
                signal_bc_r_145_307 = reads_file.hon_norm_atac(signal_bc_r_145_307, perc, std)
                signal_bc_r_145_307_slope = reads_file.slope(signal_bc_r_145_307, reads_file.sg_coefs)

                signal_bc_f_307_500 = reads_file.boyle_norm(signal_bc_f_307_500)
                perc = scoreatpercentile(signal_bc_f_307_500, 98)
                std = np.array(signal_bc_f_307_500).std()
                signal_bc_f_307_500 = reads_file.hon_norm_atac(signal_bc_f_307_500, perc, std)
                signal_bc_f_307_500_slope = reads_file.slope(signal_bc_f_307_500, reads_file.sg_coefs)

                signal_bc_r_307_500 = reads_file.boyle_norm(signal_bc_r_307_500)
                perc = scoreatpercentile(signal_bc_r_307_500, 98)
                std = np.array(signal_bc_r_307_500).std()
                signal_bc_r_307_500 = reads_file.hon_norm_atac(signal_bc_r_307_500, perc, std)
                signal_bc_r_307_500_slope = reads_file.slope(signal_bc_r_307_500, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)

                input_sequence.append(signal_bc_f_145_307)
                input_sequence.append(signal_bc_f_145_307_slope)
                input_sequence.append(signal_bc_r_145_307)
                input_sequence.append(signal_bc_r_145_307_slope)

                input_sequence.append(signal_bc_f_307_500)
                input_sequence.append(signal_bc_f_307_500_slope)
                input_sequence.append(signal_bc_r_307_500)
                input_sequence.append(signal_bc_r_307_500_slope)
            elif args.model == "f":
                signal_bc_f_max_145, signal_bc_r_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=True)
                signal_bc_f_145_307, signal_bc_r_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=True)
                signal_bc_f_307_500, signal_bc_r_307_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=500, strand=True)
                signal_bc_f_min_500, signal_bc_r_min_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=500, max_length=None, strand=True)
                signal_bc_f_max_145 = reads_file.boyle_norm(signal_bc_f_max_145)
                perc = scoreatpercentile(signal_bc_f_max_145, 98)
                std = np.array(signal_bc_f_max_145).std()
                signal_bc_f_max_145 = reads_file.hon_norm_atac(signal_bc_f_max_145, perc, std)
                signal_bc_f_max_145_slope = reads_file.slope(signal_bc_f_max_145, reads_file.sg_coefs)

                signal_bc_r_max_145 = reads_file.boyle_norm(signal_bc_r_max_145)
                perc = scoreatpercentile(signal_bc_r_max_145, 98)
                std = np.array(signal_bc_r_max_145).std()
                signal_bc_r_max_145 = reads_file.hon_norm_atac(signal_bc_r_max_145, perc, std)
                signal_bc_r_max_145_slope = reads_file.slope(signal_bc_r_max_145, reads_file.sg_coefs)

                signal_bc_f_145_307 = reads_file.boyle_norm(signal_bc_f_145_307)
                perc = scoreatpercentile(signal_bc_f_145_307, 98)
                std = np.array(signal_bc_f_145_307).std()
                signal_bc_f_145_307 = reads_file.hon_norm_atac(signal_bc_f_145_307, perc, std)
                signal_bc_f_145_307_slope = reads_file.slope(signal_bc_f_145_307, reads_file.sg_coefs)

                signal_bc_r_145_307 = reads_file.boyle_norm(signal_bc_r_145_307)
                perc = scoreatpercentile(signal_bc_r_145_307, 98)
                std = np.array(signal_bc_r_145_307).std()
                signal_bc_r_145_307 = reads_file.hon_norm_atac(signal_bc_r_145_307, perc, std)
                signal_bc_r_145_307_slope = reads_file.slope(signal_bc_r_145_307, reads_file.sg_coefs)

                signal_bc_f_307_500 = reads_file.boyle_norm(signal_bc_f_307_500)
                perc = scoreatpercentile(signal_bc_f_307_500, 98)
                std = np.array(signal_bc_f_307_500).std()
                signal_bc_f_307_500 = reads_file.hon_norm_atac(signal_bc_f_307_500, perc, std)
                signal_bc_f_307_500_slope = reads_file.slope(signal_bc_f_307_500, reads_file.sg_coefs)

                signal_bc_r_307_500 = reads_file.boyle_norm(signal_bc_r_307_500)
                perc = scoreatpercentile(signal_bc_r_307_500, 98)
                std = np.array(signal_bc_r_307_500).std()
                signal_bc_r_307_500 = reads_file.hon_norm_atac(signal_bc_r_307_500, perc, std)
                signal_bc_r_307_500_slope = reads_file.slope(signal_bc_r_307_500, reads_file.sg_coefs)

                signal_bc_f_min_500 = reads_file.boyle_norm(signal_bc_f_min_500)
                perc = scoreatpercentile(signal_bc_f_min_500, 98)
                std = np.array(signal_bc_f_min_500).std()
                signal_bc_f_min_500 = reads_file.hon_norm_atac(signal_bc_f_min_500, perc, std)
                signal_bc_f_min_500_slope = reads_file.slope(signal_bc_f_min_500, reads_file.sg_coefs)

                signal_bc_r_min_500 = reads_file.boyle_norm(signal_bc_r_min_500)
                perc = scoreatpercentile(signal_bc_r_min_500, 98)
                std = np.array(signal_bc_r_min_500).std()
                signal_bc_r_min_500 = reads_file.hon_norm_atac(signal_bc_r_min_500, perc, std)
                signal_bc_r_min_500_slope = reads_file.slope(signal_bc_r_min_500, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f_max_145)
                input_sequence.append(signal_bc_f_max_145_slope)
                input_sequence.append(signal_bc_r_max_145)
                input_sequence.append(signal_bc_r_max_145_slope)

                input_sequence.append(signal_bc_f_145_307)
                input_sequence.append(signal_bc_f_145_307_slope)
                input_sequence.append(signal_bc_r_145_307)
                input_sequence.append(signal_bc_r_145_307_slope)

                input_sequence.append(signal_bc_f_307_500)
                input_sequence.append(signal_bc_f_307_500_slope)
                input_sequence.append(signal_bc_r_307_500)
                input_sequence.append(signal_bc_r_307_500_slope)

                input_sequence.append(signal_bc_f_min_500)
                input_sequence.append(signal_bc_f_min_500_slope)
                input_sequence.append(signal_bc_r_min_500)
                input_sequence.append(signal_bc_r_min_500_slope)
            elif args.model == "all":
                signal_bc_f, signal_bc_r = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=None, strand=True)

                signal_bc_f = reads_file.boyle_norm(signal_bc_f)
                perc = scoreatpercentile(signal_bc_f, 98)
                std = np.array(signal_bc_f).std()
                signal_bc_f = reads_file.hon_norm_atac(signal_bc_f, perc, std)
                signal_bc_f_slope = reads_file.slope(signal_bc_f, reads_file.sg_coefs)

                signal_bc_r = reads_file.boyle_norm(signal_bc_r)
                perc = scoreatpercentile(signal_bc_r, 98)
                std = np.array(signal_bc_r).std()
                signal_bc_r = reads_file.hon_norm_atac(signal_bc_r, perc, std)
                signal_bc_r_slope = reads_file.slope(signal_bc_r, reads_file.sg_coefs)

                input_sequence.append(signal_bc_f)
                input_sequence.append(signal_bc_f_slope)
                input_sequence.append(signal_bc_r)
                input_sequence.append(signal_bc_r_slope)


            posterior_list = hmm.predict(np.array(input_sequence).T)

            if args.fp_bed_fname is not None:
                output_bed_file(region.chrom, region.initial, region.final, posterior_list,
                                args.fp_bed_fname, fp_state)

            # Formatting results
            start_pos = 0
            flag_start = False
            fp_state_nb = fp_state
            for k in range(p1, p1 + len(posterior_list)):
                curr_index = k - p1
                if flag_start:
                    if posterior_list[curr_index] != fp_state_nb:
                        if k - start_pos < fp_max_size:
                            fp = GenomicRegion(region.chrom, start_pos, k)
                            footprints.add(fp)
                        flag_start = False
                else:
                    if posterior_list[curr_index] == fp_state_nb:
                        flag_start = True
                        start_pos = k
            if flag_start:
                if region.initial + len(posterior_list) - start_pos < fp_max_size:
                    fp = GenomicRegion(region.chrom, start_pos, region.final)
                    footprints.add(fp)
    else:
        for region in regions:
            p1 = region.initial
            p2 = region.final

            input_sequence = list()
            if args.model == "a":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)

                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
            elif args.model == "b":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_min_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=None, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)

                signal_bc_min_145 = reads_file.boyle_norm(signal_bc_min_145)
                perc = scoreatpercentile(signal_bc_min_145, 98)
                std = np.array(signal_bc_min_145).std()
                signal_bc_min_145 = reads_file.hon_norm_atac(signal_bc_min_145, perc, std)
                signal_bc_min_145_slope = reads_file.slope(signal_bc_min_145, reads_file.sg_coefs)

                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
                input_sequence.append(signal_bc_min_145)
                input_sequence.append(signal_bc_min_145_slope)
            elif args.model == "c":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)

                signal_bc_145_307 = reads_file.boyle_norm(signal_bc_145_307)
                perc = scoreatpercentile(signal_bc_145_307, 98)
                std = np.array(signal_bc_145_307).std()
                signal_bc_145_307 = reads_file.hon_norm_atac(signal_bc_145_307, perc, std)
                signal_bc_145_307_slope = reads_file.slope(signal_bc_145_307, reads_file.sg_coefs)

                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
                input_sequence.append(signal_bc_145_307)
                input_sequence.append(signal_bc_145_307_slope)
            elif args.model == "d":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=False)

                signal_bc_min_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=None, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)

                signal_bc_145_307 = reads_file.boyle_norm(signal_bc_145_307)
                perc = scoreatpercentile(signal_bc_145_307, 98)
                std = np.array(signal_bc_145_307).std()
                signal_bc_145_307 = reads_file.hon_norm_atac(signal_bc_145_307, perc, std)
                signal_bc_145_307_slope = reads_file.slope(signal_bc_145_307, reads_file.sg_coefs)

                signal_bc_min_307 = reads_file.boyle_norm(signal_bc_min_307)
                perc = scoreatpercentile(signal_bc_min_307, 98)
                std = np.array(signal_bc_min_307).std()
                signal_bc_min_307 = reads_file.hon_norm_atac(signal_bc_min_307, perc, std)
                signal_bc_min_307_slope = reads_file.slope(signal_bc_min_307, reads_file.sg_coefs)

                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
                input_sequence.append(signal_bc_145_307)
                input_sequence.append(signal_bc_145_307_slope)
                input_sequence.append(signal_bc_min_307)
                input_sequence.append(signal_bc_min_307_slope)
            elif args.model == "e":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=False)

                signal_bc_307_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=500, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)


                signal_bc_145_307 = reads_file.boyle_norm(signal_bc_145_307)
                perc = scoreatpercentile(signal_bc_145_307, 98)
                std = np.array(signal_bc_145_307).std()
                signal_bc_145_307 = reads_file.hon_norm_atac(signal_bc_145_307, perc, std)
                signal_bc_145_307_slope = reads_file.slope(signal_bc_145_307, reads_file.sg_coefs)


                signal_bc_307_500 = reads_file.boyle_norm(signal_bc_307_500)
                perc = scoreatpercentile(signal_bc_307_500, 98)
                std = np.array(signal_bc_307_500).std()
                signal_bc_307_500 = reads_file.hon_norm_atac(signal_bc_307_500, perc, std)
                signal_bc_307_500_slope = reads_file.slope(signal_bc_307_500, reads_file.sg_coefs)


                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
                input_sequence.append(signal_bc_145_307)
                input_sequence.append(signal_bc_145_307_slope)
                input_sequence.append(signal_bc_307_500)
                input_sequence.append(signal_bc_307_500_slope)
            elif args.model == "f":
                signal_bc_max_145 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=145, strand=False)

                signal_bc_145_307 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=145, max_length=307, strand=False)

                signal_bc_307_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=307, max_length=500, strand=False)

                signal_bc_min_500 = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=500, max_length=None, strand=False)

                signal_bc_max_145 = reads_file.boyle_norm(signal_bc_max_145)
                perc = scoreatpercentile(signal_bc_max_145, 98)
                std = np.array(signal_bc_max_145).std()
                signal_bc_max_145 = reads_file.hon_norm_atac(signal_bc_max_145, perc, std)
                signal_bc_max_145_slope = reads_file.slope(signal_bc_max_145, reads_file.sg_coefs)


                signal_bc_145_307 = reads_file.boyle_norm(signal_bc_145_307)
                perc = scoreatpercentile(signal_bc_145_307, 98)
                std = np.array(signal_bc_145_307).std()
                signal_bc_145_307 = reads_file.hon_norm_atac(signal_bc_145_307, perc, std)
                signal_bc_145_307_slope = reads_file.slope(signal_bc_145_307, reads_file.sg_coefs)


                signal_bc_307_500 = reads_file.boyle_norm(signal_bc_307_500)
                perc = scoreatpercentile(signal_bc_307_500, 98)
                std = np.array(signal_bc_307_500).std()
                signal_bc_307_500 = reads_file.hon_norm_atac(signal_bc_307_500, perc, std)
                signal_bc_307_500_slope = reads_file.slope(signal_bc_307_500, reads_file.sg_coefs)


                signal_bc_min_500 = reads_file.boyle_norm(signal_bc_min_500)
                perc = scoreatpercentile(signal_bc_min_500, 98)
                std = np.array(signal_bc_min_500).std()
                signal_bc_min_500 = reads_file.hon_norm_atac(signal_bc_min_500, perc, std)
                signal_bc_min_500_slope = reads_file.slope(signal_bc_min_500, reads_file.sg_coefs)

                input_sequence.append(signal_bc_max_145)
                input_sequence.append(signal_bc_max_145_slope)
                input_sequence.append(signal_bc_145_307)
                input_sequence.append(signal_bc_145_307_slope)
                input_sequence.append(signal_bc_307_500)
                input_sequence.append(signal_bc_307_500_slope)
                input_sequence.append(signal_bc_min_500)
                input_sequence.append(signal_bc_min_500_slope)
            elif args.model == "all":
                signal_bc = \
                    reads_file.get_bc_signal_by_fragment_length(ref=region.chrom, start=p1, end=p2,
                                                                bam=bam, fasta=fasta, bias_table=bias_table,
                                                                forward_shift=forward_shift,
                                                                reverse_shift=reverse_shift,
                                                                min_length=None, max_length=None, strand=False)

                signal_bc = reads_file.boyle_norm(signal_bc)
                perc = scoreatpercentile(signal_bc, 98)
                std = np.array(signal_bc).std()
                signal_bc = reads_file.hon_norm_atac(signal_bc, perc, std)
                signal_bc_slope = reads_file.slope(signal_bc, reads_file.sg_coefs)

                input_sequence.append(signal_bc)
                input_sequence.append(signal_bc_slope)


            posterior_list = hmm.predict(np.array(input_sequence).T)

            # Formatting results
            start_pos = 0
            flag_start = False
            fp_state_nb = fp_state
            for k in range(p1, p1 + len(posterior_list)):
                curr_index = k - p1
                if flag_start:
                    if posterior_list[curr_index] != fp_state_nb:
                        if k - start_pos < fp_max_size:
                            fp = GenomicRegion(region.chrom, start_pos, k)
                            footprints.add(fp)
                        flag_start = False
                else:
                    if posterior_list[curr_index] == fp_state_nb:
                        flag_start = True
                        start_pos = k
            if flag_start:
                if region.initial + len(posterior_list) - start_pos < fp_max_size:
                    fp = GenomicRegion(region.chrom, start_pos, region.final)
                    footprints.add(fp)
    ###################################################################################################
    # Post-processing
    ###################################################################################################

    post_processing(footprints=footprints, original_regions=original_regions, fp_min_size=fp_min_size,
                    fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                    reads_file=reads_file, downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                    forward_shift=forward_shift, reverse_shift=reverse_shift,
                    initial_clip=initial_clip, output_location=args.output_location,
                    output_prefix=args.output_prefix)
