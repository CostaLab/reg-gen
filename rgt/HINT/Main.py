###################################################################################################
# Libraries
###################################################################################################

# Python
from os import remove, system, getcwd
from sys import exit
from copy import deepcopy
from optparse import SUPPRESS_HELP
import warnings

warnings.filterwarnings("ignore")

# Internal
from rgt import __version__
from ..Util import PassThroughOptionParser, ErrorHandler, HmmData, GenomeData, OverlapType
from ..ExperimentalMatrix import ExperimentalMatrix
from ..GenomicRegion import GenomicRegion
from ..GenomicRegionSet import GenomicRegionSet
from signalProcessing import GenomicSignal
from hmm import HMM
from biasTable import BiasTable
from evaluation import Evaluation
from train import TrainHMM
from evidence import Evidence
from plot import Plot

# External
import os
from numpy import array, sum, isnan
from hmmlearn.hmm import GaussianHMM
from hmmlearn import __version__ as hmm_ver

"""
HINT - HMM-based Identification of TF Footprints.
Finds transcription factor footprints given open chromatin data.

Basic Input:
- Regions (bed) in which to find footprints (i.e. enriched regions or hypersensitivity regions).
- Reads (bam) containing the open chromatin signal for DNase/ATAC and 0 <= N <= 3 histone modifications.

Dependencies:
- python >= 2.7
- numpy >= 1.4.0
- scipy >= 0.7.0
- scikit >= 0.14
- hmmlearn >= 0.1.1
- pysam >= 0.7.5
- matplotlib >= 1.6
Authors: Eduardo G. Gusmao, Manuel Allhoff, Joseph Kuo and Ivan G. Costa.
"""


def main():
    """
    Main function that performs footprint analysis.

    Keyword arguments: None

    Return: None
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing ErrorHandler
    error_handler = ErrorHandler()

    # Parameters
    usage_message = ("\n--------------------------------------------------\n"
                     "The 'hint' program predicts TFBSs given open chromatin data.\n"
                     "In order to use this tools, please type: \n\n"
                     "%prog [options] <experiment_matrix>\n\n"
                     "The minimal <experiment matrix> should contain:\n"
                     "- One region file representing the regions in which the HMM\n"
                     "  will be applied. It should contain 'regions' in the type field\n"
                     "  and 'HS' in the data field\n"
                     "- One DNase-seq or ATAC-seq aligned reads file (bam) file with\n"
                     "  'reads' in the type field and 'DNASE' or 'ATAC' in the data field.\n"
                     "- Zero to Three histone modification aligned reads file (bam)\n"
                     "  with 'reads' in the type field and 'HISTONE' in the data field.\n\n"

                     "For more information, please refer to:\n"
                     "http://www.regulatory-genomics.org/hint/introduction/\n\n"

                     "For further questions or comments please refer to our group:\n"
                     "https://groups.google.com/forum/#!forum/rgtusers\n"
                     "--------------------------------------------------")
    version_message = "HINT - Regulatory Analysis Toolbox (RGT). Version: " + str(__version__)

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message, version=version_message)

    # Optional Input Options
    parser.add_option("--hmm-file", dest="hmm_file", type="string",
                      metavar="FILE_1_1[[,...,FILE_N_1];...;FILE_1_M[,...,FILE_N_M]]", default=None,
                      help=("List of HMM files separated by comma. Rules: "
                            "- If DNase/ATAC-only analysis, provide only one HMM file per group. "
                            "- If also using histone modifications: "
                            "  * If only one file is provided it will be applied for all histone data, "
                            "  * Otherwise, the list must have the same number of input histone data. "
                            "    The order of the list should be the order of the histones in the "
                            "    input_matrix file. "
                            "- In case multiple input groups are used, then multiple file lists can "
                            "  be passed using semicolon. The number of groups of lists should "
                            "  equals the number of input groups. "
                            "- If the argument is not given, then a default HMM will be used."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE1_F,FILE1_R[;...;FILEM_F,FILEM_R]", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates. Each input group"
                            "should have two files: one for the forward and one for the reverse strand. "
                            "Each line should contain a kmer and the bias estimate separated by tab. "
                            "Leave an empty set for histone-only analysis groups. Eg. FILE1;;FILE3."))

    # Parameters Options
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--estimate-bias-correction", dest="estimate_bias_correction",
                      action="store_true", default=False,
                      help=("Applies DNase-seq sequence cleavage bias correction with k-mer bias estimated "
                            "from the given DNase-seq data (SLOW HINT-BC)."))
    parser.add_option("--original-regions", dest="original_regions", type="string",
                      metavar="STRING", default=None,
                      help=("The regions that used to estimate the bias table "
                            "should be bed file containing HS regions or FASTA file containing naked DNA"))
    parser.add_option("--estimate-bias-type", dest="estimate_bias_type", type="string",
                      metavar="STRING", default=None,
                      help=("The methods that used to estimate the bias table "
                            "Available options are: 'FRE' (the bias estimation is computed as the ratio "
                            "between the observed and background cleavage frequency) "
                            "and 'PWM' (frequency is replaced by pwm score)."))
    parser.add_option("--default-bias-correction", dest="default_bias_correction",
                      action="store_true", default=False,
                      help=("Applies DNase-seq cleavage bias correction with default k-mer bias "
                            "estimates (FAST HINT-BC). Please set the correct --default-bias-type "
                            "option that matches your experimental settings."))
    parser.add_option("--default-bias-type", dest="default_bias_type", type="string",
                      metavar="STRING", default="DH",
                      help=("Type of protocol used to generate the DNase-seq or ATAC-seq data. "
                            "Available options are: 'SH' (DNase-seq single-hit protocol), 'DH' "
                            "(DNase-seq double-hit protocol) and 'ATAC' (ATAC-seq data)."))

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH",
                      default=getcwd(),
                      help=("Path where the output files will be written."))
    parser.add_option("--output-fname", dest="output_fname", type="string", metavar="STRING",
                      default=None)
    parser.add_option("--print-raw-signal", dest="print_raw_signal", type="string", metavar="STRING",
                      default=None,
                      help=("If used, it will print the base overlap (raw) signals from DNase-seq "
                            " or ATAC-seq data. The option should equal the file name."
                            "The extension must be (.wig)."))
    parser.add_option("--print-bc-signal", dest="print_bc_signal", type="string", metavar="STRING",
                      default=None,
                      help=("If used, it will print the DNase-seq or ATAC-seq bias-corrected signal. "
                            "The option should equal the file name. The extension must be (.wig)."))
    parser.add_option("--print-norm-signal", dest="print_norm_signal", type="string", metavar="STRING",
                      default=None,
                      help=("If used, it will print the normalized signals from DNase-seq "
                            " or ATAC-seq data. The option should equal the file name."
                            "The extension must be (.wig)."))
    parser.add_option("--print-slope-signal", dest="print_slope_signal", type="string", metavar="STRING",
                      default=None,
                      help=("If used, it will print the slope signals from DNase-seq "
                            " or ATAC-seq data. The option should equal the file name."
                            "The extension must be (.wig)."))
    parser.add_option("--print-line-plot", dest="print_line_plot",
                      action="store_true", default=False,
                      help=("If used, it will print the line plot of raw signal and bias corrected"
                            "signal for the particular motif."))
    parser.add_option("--window-size", dest="window_size", type="int", metavar="INT", default=50)
    parser.add_option("--motif-file", dest="motif_file", type="string", metavar="STRING",
                      default=None,
                      help=("A bed file containing all motif-predicted binding sites (MPBSs)."
                            "The extension must be (.bed)."))
    parser.add_option("--motif-name", dest="motif_name", type="string", metavar="STRING",
                      default=None)
    parser.add_option("--protection-score", dest="protection_score",
                      action="store_true", default=False,
                      help=("If used, it will print the protection score"))
    parser.add_option("--strands-specific", dest="strands_specific",
                      action="store_true", default=False,
                      help=("A boolean indicating if DNaseI digestion site "
                            "counts should be aggregated across strands. "
                            "default: False"))

    # Train Options
    parser.add_option("--train-hmm", dest="train_hmm",
                      action="store_true", default=False,
                      help=("If used, HINT will train a hidden Markov model (HMM) based on "
                            "the annotation data"))
    parser.add_option("--bam-file", dest="bam_file", type="string", metavar="STRING",
                      default=None,
                      help=("A bam file containing all the DNase-seq reads."))
    parser.add_option("--annotate-file", dest="annotate_file", type="string", metavar="STRING",
                      default=None,
                      help=("A annotate file containing all the states."))
    parser.add_option("--print-bed-file", dest="print_bed_file",
                      action="store_true", default=False,
                      help=("If used, HINT will output the bed file containing "
                            "the annotated regions, so that you can visualize "
                            "you HMM annotation for potential errors"))
    parser.add_option("--model-fname", dest="model_fname", type="string", metavar="STRING",
                      default="model",
                      help=("The output file name"))

    # Evaluation Options
    parser.add_option("--evaluate-footprints", dest="evaluate_footprints",
                      action="store_true", default=False,
                      help=("If used, HINT will evaluate the footprints prediction."))
    parser.add_option("--tf-name", dest="tf_name", type="string", metavar="STRING",
                      default=None,
                      help=("The name of transcription factor"))
    parser.add_option("--tfbs-file", dest="tfbs_file", type="string", metavar="STRING",
                      default=None,
                      help=("A bed file containing all motif-predicted binding sites (MPBSs)."
                            "The values in the bed SCORE field will be used to rank the MPBSs."
                            "The extension must be (.bed)."))
    parser.add_option("--footprint-file", dest="footprint_file", type="string",
                      metavar="FILE1,FILE2,FILE3,FILE4...",
                      default=None,
                      help=("The bed files containing the footprint predictions."
                            "The extension must be (.bed)."))
    parser.add_option("--footprint-name", dest="footprint_name", type="string",
                      metavar="NAME1,NAME2,NAME3,NAME4...",
                      default=None,
                      help=("The methods name used to predicted the footprint."
                            "The number of methods name must be consistent with that of footprint file"))
    parser.add_option("--alignment-file", dest="alignment_file", type="string",
                      metavar="FILE",
                      default=None,
                      help=("A bam file containing all the DNase-seq reads."
                            "Used to fetch the tag counts for SEG (segmentation approach)"))
    parser.add_option("--footprint-type", dest="footprint_type", type="string",
                      metavar="TYPE1,TYPE2,TYPE3,TYPE4...",
                      default=None,
                      help=("The methods type used to predicted the footprint."
                            "Must be SC (site centric) or SEG (segmentation approach)"))
    parser.add_option("--print-roc-curve", dest="print_roc_curve",
                      action="store_true", default=False,
                      help=("If used, HINT will print the receiver operating characteristic curve."))
    parser.add_option("--print-pr-curve", dest="print_pr_curve",
                      action="store_true", default=False,
                      help=("If used, HINT will print the precision recall curve."))

    # Evidence Options
    parser.add_option("--create-evidence", dest="create_evidence",
                      action="store_true", default=False,
                      help=("If used, HINT will create a bed file with MPBSs with and "
                            "without evidence. Also, the name of the instances will be Y "
                            "for evidence, N for non evidence."))
    parser.add_option("--peak-ext", dest="peak_ext", type="int", metavar="INT", default=100,
                      help=("The number used to extend the ChIP-seq summit"))
    parser.add_option("--mpbs-name", dest="mpbs_name", type="string",
                      default=None,
                      help=("motif predicted binding sites name"))
    parser.add_option("--tfbs-summit-fname", dest="tfbs_summit_fname", type="string",
                      default=None,
                      help=("the ChIP-seq peak files"))
    parser.add_option("--mpbs-fname", dest="mpbs_fname", type="string",
                      default=None,
                      help=("motif predicted binding sites file"))

    # GENERAL Hidden Options
    parser.add_option("--region-total-ext", dest="region_total_ext", type="int", metavar="INT", default=10000,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size", dest="fp_limit_size", type="int", metavar="INT", default=50,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size-histone", dest="fp_limit_size_histone", type="int",
                      metavar="INT", default=2000, help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size-ext", dest="fp_limit_size_ext", type="int", metavar="INT", default=10,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size-ext-histone", dest="fp_limit_size_ext_histone",
                      type="int", metavar="INT", default=200, help=SUPPRESS_HELP)
    parser.add_option("--fp-ext", dest="fp_ext", type="int", metavar="INT", default=5,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-ext-histone", dest="fp_ext_histone", type="int", metavar="INT", default=50,
                      help=SUPPRESS_HELP)
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=50,
                      help=SUPPRESS_HELP)
    parser.add_option("--tc-ext-histone", dest="tc_ext_histone", type="int", metavar="INT", default=500,
                      help=SUPPRESS_HELP)

    # DNASE Hidden Options
    parser.add_option("--dnase-initial-clip", dest="dnase_initial_clip", type="int",
                      metavar="INT", default=1000, help=SUPPRESS_HELP)
    parser.add_option("--dnase-sg-window-size", dest="dnase_sg_window_size", type="int",
                      metavar="INT", default=9, help=SUPPRESS_HELP)
    parser.add_option("--dnase-norm-per", dest="dnase_norm_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--dnase-slope-per", dest="dnase_slope_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--dnase-downstream-ext", dest="dnase_downstream_ext", type="int",
                      metavar="INT", default=1, help=SUPPRESS_HELP)
    parser.add_option("--dnase-upstream-ext", dest="dnase_upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--dnase-forward-shift", dest="dnase_forward_shift", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--dnase-reverse-shift", dest="dnase_reverse_shift", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--dnase-bias-correction-k", dest="dnase_bias_correction_k", type="int",
                      metavar="INT", default=6, help=SUPPRESS_HELP)

    # ATAC Hidden Options
    parser.add_option("--atac-initial-clip", dest="atac_initial_clip", type="int",
                      metavar="INT", default=1000, help=SUPPRESS_HELP)
    parser.add_option("--atac-sg-window-size", dest="atac_sg_window_size", type="int",
                      metavar="INT", default=9, help=SUPPRESS_HELP)
    parser.add_option("--atac-norm-per", dest="atac_norm_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--atac-slope-per", dest="atac_slope_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--atac-downstream-ext", dest="atac_downstream_ext", type="int",
                      metavar="INT", default=1, help=SUPPRESS_HELP)
    parser.add_option("--atac-upstream-ext", dest="atac_upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--atac-forward-shift", dest="atac_forward_shift", type="int",
                      metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--atac-reverse-shift", dest="atac_reverse_shift", type="int",
                      metavar="INT", default=-4, help=SUPPRESS_HELP)
    parser.add_option("--atac-bias-correction-k", dest="atac_bias_correction_k", type="int",
                      metavar="INT", default=6, help=SUPPRESS_HELP)

    # HISTONE Hidden Options
    parser.add_option("--histone-initial-clip", dest="histone_initial_clip", type="int",
                      metavar="INT", default=1000, help=SUPPRESS_HELP)
    parser.add_option("--histone-sg-window-size", dest="histone_sg_window_size", type="int",
                      metavar="INT", default=201, help=SUPPRESS_HELP)
    parser.add_option("--histone-norm-per", dest="histone_norm_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--histone-slope-per", dest="histone_slope_per", type="float",
                      metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--histone-downstream-ext", dest="histone_downstream_ext", type="int",
                      metavar="INT", default=200, help=SUPPRESS_HELP)
    parser.add_option("--histone-upstream-ext", dest="histone_upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--histone-forward-shift", dest="histone_forward_shift", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--histone-reverse-shift", dest="histone_reverse_shift", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)

    # Processing Options
    options, arguments = parser.parse_args()
    if len(arguments) < 1:
        print(usage_message)
        exit(1)
    # if(not arguments or len(arguments) > 1): error_handler.throw_error("FP_WRONG_ARGUMENT")

    # General hidden options ###############################################################
    region_total_ext = options.region_total_ext
    fp_limit_size = options.fp_limit_size
    fp_limit_size_histone = options.fp_limit_size_histone
    fp_limit_size_ext = options.fp_limit_size_ext
    fp_limit_size_ext_histone = options.fp_limit_size_ext_histone
    fp_ext = options.fp_ext
    fp_ext_histone = options.fp_ext_histone
    tc_ext = options.tc_ext
    tc_ext_histone = options.tc_ext_histone
    # DNASE Hidden Options
    dnase_initial_clip = options.dnase_initial_clip
    dnase_sg_window_size = options.dnase_sg_window_size
    dnase_norm_per = options.dnase_norm_per
    dnase_slope_per = options.dnase_slope_per
    dnase_downstream_ext = options.dnase_downstream_ext
    dnase_upstream_ext = options.dnase_upstream_ext
    dnase_forward_shift = options.dnase_forward_shift
    dnase_reverse_shift = options.dnase_reverse_shift
    dnase_bias_correction_k = options.dnase_bias_correction_k
    # ATAC Hidden Options
    atac_initial_clip = options.atac_initial_clip
    atac_sg_window_size = options.atac_sg_window_size
    atac_norm_per = options.atac_norm_per
    atac_slope_per = options.atac_slope_per
    atac_downstream_ext = options.atac_downstream_ext
    atac_upstream_ext = options.atac_upstream_ext
    atac_forward_shift = options.atac_forward_shift
    atac_reverse_shift = options.atac_reverse_shift
    atac_bias_correction_k = options.atac_bias_correction_k
    # HISTONE Hidden Options
    histone_initial_clip = options.histone_initial_clip
    histone_sg_window_size = options.histone_initial_clip
    histone_norm_per = options.histone_norm_per
    histone_slope_per = options.histone_slope_per
    histone_downstream_ext = options.histone_downstream_ext
    histone_upstream_ext = options.histone_upstream_ext
    histone_forward_shift = options.histone_forward_shift
    histone_reverse_shift = options.histone_reverse_shift

    ########################################################################################

    ##########################################################################################

    # IF HINT is required to output the line plot and motif logo
    if options.print_line_plot:
        plot = Plot(options.bam_file, options.motif_file, options.motif_name, options.window_size,
                    atac_downstream_ext, atac_upstream_ext, atac_forward_shift, atac_reverse_shift,
                    atac_initial_clip, options.organism, options.bias_table,
                    atac_bias_correction_k, options.protection_score, options.strands_specific,
                    options.output_location)
        plot.line()
        return

    # If HINT is required to create the validation data set
    if options.create_evidence:
        evidence = Evidence(options.peak_ext, options.mpbs_name, options.tfbs_summit_fname,
                            options.mpbs_fname, options.output_location)
        evidence.create_file()
        return

    # If HINT is required to evaluate the existing footprint predictions
    if options.evaluate_footprints:
        evaluation = Evaluation(options.tf_name, options.tfbs_file, options.footprint_file,
                                options.footprint_name, options.footprint_type,
                                options.print_roc_curve, options.print_roc_curve,
                                options.output_location, options.alignment_file, options.organism)
        evaluation.chip_evaluate()
        return

    # If HINT is required to train a hidden Markov model (HMM)
    if options.train_hmm:
        train_hmm_model = TrainHMM(options.bam_file, options.annotate_file, options.print_bed_file,
                                   options.output_location, options.model_fname,
                                   options.print_raw_signal, options.print_bc_signal,
                                   options.print_norm_signal, options.print_slope_signal,
                                   atac_initial_clip, atac_downstream_ext, atac_upstream_ext,
                                   atac_forward_shift, atac_reverse_shift,
                                   options.estimate_bias_correction, options.estimate_bias_type,
                                   options.bias_table,
                                   options.original_regions, options.organism,
                                   atac_bias_correction_k)
        train_hmm_model.train()
        return

    ########################################################################################################
    # Output wig signal
    if (options.print_raw_signal):
        system("touch " + options.print_raw_signal + " | echo -n "" > " + options.print_raw_signal)
    if (options.print_bc_signal):
        system("touch " + options.print_bc_signal + " | echo -n "" > " + options.print_bc_signal)
    if (options.print_norm_signal):
        system("touch " + options.print_norm_signal + " | echo -n "" > " + options.print_norm_signal)
    if (options.print_slope_signal):
        system("touch " + options.print_raw_signal + " | echo -n "" > " + options.print_slope_signal)

    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading input argument
    input_matrix = arguments[0]

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception:
        error_handler.throw_error("FP_WRONG_EXPMAT")

    ###################################################################################################
    # Reading Input
    ###################################################################################################

    # Group class
    class Group:
        def __init__(self):
            self.name = None
            self.original_regions = None
            self.regions = None
            self.dnase_file = None
            self.histone_file_list = []
            self.dnase_only = True
            self.histone_only = True
            self.hmm = []
            self.flag_multiple_hmms = False
            self.bias_table = None
            self.is_atac = False

    # Initialization
    name_list = exp_matrix.names
    type_list = exp_matrix.types
    file_dict = exp_matrix.files
    fields_dict = exp_matrix.fieldsDict
    objects_dict = exp_matrix.objectsDict

    # Populating fields dict data
    for e in ["HS", "DNASE", "ATAC", "HISTONE"]:
        try:
            fields_dict["data"][e]
        except Exception:
            fields_dict["data"][e] = []

    # Fetching chromosome sizes
    chrom_sizes_file_name = genome_data.get_chromosome_sizes()
    chrom_sizes_file = open(chrom_sizes_file_name, "r")
    chrom_sizes_dict = dict()
    for chrom_sizes_entry_line in chrom_sizes_file:
        chrom_sizes_entry_vec = chrom_sizes_entry_line.strip().split("\t")
        chrom_sizes_dict[chrom_sizes_entry_vec[0]] = int(chrom_sizes_entry_vec[1])
    chrom_sizes_file.close()

    # Fetching files per group
    group_list = []
    for g in fields_dict["group"].keys():
        group = Group()
        group.name = g
        for i in range(0, len(fields_dict["group"][g])):
            if (name_list[i] in fields_dict["data"]["HS"]):
                group.original_regions = objects_dict[name_list[i]]
                group.regions = deepcopy(group.original_regions)
                group.regions.extend(int(region_total_ext / 2), int(region_total_ext / 2))  # Extending
                group.regions.merge()  # Sort & Merge
            elif (name_list[i] in fields_dict["data"]["DNASE"]):
                group.dnase_file = GenomicSignal(file_dict[name_list[i]])
                group.dnase_file.load_sg_coefs(dnase_sg_window_size)
                group.is_atac = False
            elif (name_list[i] in fields_dict["data"]["ATAC"]):
                group.dnase_file = GenomicSignal(file_dict[name_list[i]])
                group.dnase_file.load_sg_coefs(atac_sg_window_size)
                group.is_atac = True
            elif (name_list[i] in fields_dict["data"]["HISTONE"]):
                group.histone_file_list.append(GenomicSignal(file_dict[name_list[i]]))
                group.histone_file_list[-1].load_sg_coefs(histone_sg_window_size)
            else:
                pass  # TODO ERROR (Category of data outside "HS, DNASE, ATAC, HISTONE")
        if (group.dnase_file): group.histone_only = False
        if (group.histone_file_list): group.dnase_only = False
        if (group.histone_only and group.dnase_only): pass  # TODO ERROR (There is no DNase or histone data)
        if (not group.original_regions): pass  # TODO ERROR (There is no HS regions)
        group_list.append(group)

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################

    bias_correction = False
    if (options.bias_table):

        bias_table_group_list = options.bias_table.split(";")
        if (len(bias_table_group_list) != len(group_list)): pass  # TODO ERROR
        for g in range(0, len(group_list)):
            group = group_list[g]
            bias_table_list = bias_table_group_list[g].split(",")
            if (group.histone_only): continue
            group.bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                                      table_file_name_R=bias_table_list[1])
        bias_correction = True

    elif (options.estimate_bias_correction):
        for group in group_list:
            if (group.histone_only): continue
            if (group.is_atac):
                my_k_nb = atac_bias_correction_k
                my_shift = atac_downstream_ext
                group.bias_table = BiasTable().estimate_table(regions=group.original_regions,
                                                              dnase_file_name=group.dnase_file.file_name,
                                                              genome_file_name=genome_data.get_genome(),
                                                              k_nb=my_k_nb,
                                                              forward_shift=atac_forward_shift,
                                                              reverse_shift=atac_reverse_shift)
            else:
                my_k_nb = dnase_bias_correction_k
                my_shift = dnase_downstream_ext
                group.bias_table = BiasTable().estimate_table(regions=group.original_regions,
                                                              dnase_file_name=group.dnase_file.file_name,
                                                              genome_file_name=genome_data.get_genome(),
                                                              k_nb=my_k_nb,
                                                              forward_shift=dnase_forward_shift,
                                                              reverse_shift=dnase_reverse_shift)
        bias_correction = True

    elif (options.default_bias_correction):

        for group in group_list:
            if (group.histone_only): continue

            if (options.default_bias_type == "SH"):
                my_table_file_F = hmm_data.get_default_bias_table_F_SH()
                my_table_file_R = hmm_data.get_default_bias_table_R_SH()
            elif (options.default_bias_type == "DH"):
                my_table_file_F = hmm_data.get_default_bias_table_F_DH()
                my_table_file_R = hmm_data.get_default_bias_table_R_DH()
            else:
                my_table_file_F = hmm_data.get_default_bias_table_F_ATAC()
                my_table_file_R = hmm_data.get_default_bias_table_R_ATAC()

            group.bias_table = BiasTable().load_table(table_file_name_F=my_table_file_F,
                                                      table_file_name_R=my_table_file_R)
        bias_correction = True
    ###################################################################################################
    # Creating HMMs
    ###################################################################################################

    # Fetching HMM input
    flag_multiple_hmms = False
    if (options.hmm_file):  # Argument is passed

        hmm_group_list = options.hmm_file.split(";")
        if (len(hmm_group_list) != len(group_list)): pass  # TODO ERROR
        for g in range(0, len(group_list)):

            group = group_list[g]

            # Fetching list of HMM files
            group.hmm = hmm_group_list[g].split(",")

            # Verifying HMM application mode (one HMM or multiple HMM files)
            if (len(group.hmm) == 1):
                group.flag_multiple_hmms = False
                group.hmm = group.hmm[0]
            elif (len(group.hmm) == len(histone_file_name_list)):
                flag_multiple_hmms = True
            else:
                error_handler.throw_error("FP_NB_HMMS")

    else:  # Argument was not passed

        for group in group_list:

            group.flag_multiple_hmms = False
            if (group.dnase_only):
                if (bias_correction):
                    if (group.is_atac):
                        group.hmm = hmm_data.get_default_hmm_atac_bc()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_bc()
                else:
                    if (group.is_atac):
                        group.hmm = hmm_data.get_default_hmm_atac()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase()
            elif (group.histone_only):
                group.hmm = hmm_data.get_default_hmm_histone()
            else:
                if (bias_correction):
                    if (group.is_atac):
                        group.hmm = hmm_data.get_default_hmm_atac_histone_bc()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_histone_bc()
                else:
                    if (group.is_atac):
                        group.hmm = hmm_data.get_default_hmm_atac_histone()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_histone()

    # Creating scikit HMM list
    for group in group_list:

        if (group.flag_multiple_hmms):

            hmm_list = []
            for hmm_file_name in group.hmm:

                try:
                    hmm_scaffold = HMM()
                    hmm_scaffold.load_hmm(hmm_file_name)
                    if (int(hmm_ver.split(".")[0]) <= 0 and int(hmm_ver.split(".")[1]) <= 1):
                        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full",
                                                 transmat=array(hmm_scaffold.A), startprob=array(hmm_scaffold.pi))
                        scikit_hmm.means_ = array(hmm_scaffold.means)
                        scikit_hmm.covars_ = array(hmm_scaffold.covs)
                    else:
                        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
                        scikit_hmm.startprob_ = array(hmm_scaffold.pi)
                        scikit_hmm.transmat_ = array(hmm_scaffold.A)
                        scikit_hmm.means_ = array(hmm_scaffold.means)
                        scikit_hmm.covars_ = array(hmm_scaffold.covs)

                except Exception:
                    error_handler.throw_error("FP_HMM_FILES")
                hmm_list.append(scikit_hmm)

            group.hmm = hmm_list

        else:

            scikit_hmm = None
            try:
                hmm_scaffold = HMM()
                hmm_scaffold.load_hmm(group.hmm)
                if (int(hmm_ver.split(".")[0]) <= 0 and int(hmm_ver.split(".")[1]) <= 1):
                    scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full",
                                             transmat=array(hmm_scaffold.A), startprob=array(hmm_scaffold.pi))
                    scikit_hmm.means_ = array(hmm_scaffold.means)
                    scikit_hmm.covars_ = array(hmm_scaffold.covs)
                else:
                    scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
                    scikit_hmm.startprob_ = array(hmm_scaffold.pi)
                    scikit_hmm.transmat_ = array(hmm_scaffold.A)
                    scikit_hmm.means_ = array(hmm_scaffold.means)
                    scikit_hmm.covars_ = array(hmm_scaffold.covs)


            except Exception:
                error_handler.throw_error("FP_HMM_FILES")
            group.hmm = scikit_hmm

    ###################################################################################################
    # Main Pipeline
    ###################################################################################################

    # Iterating over groups
    for group in group_list:

        # Initializing result set
        footprints = GenomicRegionSet(group.name)

        # Iterating over regions
        for r in group.regions.sequences:

            ###################################################################################################
            # DNASE ONLY
            ###################################################################################################

            if (group.dnase_only):

                # Fetching DNase signal
                try:
                    if (group.is_atac):
                        dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                              atac_downstream_ext, atac_upstream_ext,
                                                                              atac_forward_shift, atac_reverse_shift,
                                                                              atac_initial_clip, atac_norm_per,
                                                                              atac_slope_per,
                                                                              group.bias_table,
                                                                              genome_data.get_genome(),
                                                                              options.print_raw_signal,
                                                                              options.print_bc_signal,
                                                                              options.print_norm_signal,
                                                                              options.print_slope_signal)
                    else:
                        dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                              dnase_downstream_ext, dnase_upstream_ext,
                                                                              dnase_forward_shift, dnase_reverse_shift,
                                                                              dnase_initial_clip, dnase_norm_per,
                                                                              dnase_slope_per,
                                                                              group.bias_table,
                                                                              genome_data.get_genome(),
                                                                              options.print_raw_signal,
                                                                              options.print_bc_signal,
                                                                              options.print_norm_signal,
                                                                              options.print_slope_signal)
                except Exception:
                    raise
                    error_handler.throw_warning("FP_DNASE_PROC", add_msg="for region (" + ",".join([r.chrom,
                                                                                                    str(r.initial), str(
                            r.final)]) + "). This iteration will be skipped.")
                    continue

                # Formatting sequence
                try:
                    input_sequence = array([dnase_norm, dnase_slope]).T
                except Exception:
                    raise
                    error_handler.throw_warning("FP_SEQ_FORMAT", add_msg="for region (" + ",".join([r.chrom,
                                                                                                    str(r.initial), str(
                            r.final)]) + "). This iteration will be skipped.")
                    continue

                # Applying HMM
                if (isinstance(group.hmm, list)): continue  # TODO ERROR
                if (isnan(sum(input_sequence))): continue  # Handling NAN's in signal / hmmlearn throws error TODO ERROR
                try:
                    posterior_list = group.hmm.predict(input_sequence)
                except Exception:
                    raise
                    error_handler.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([r.chrom,
                                                                                                   str(r.initial), str(
                            r.final)]) + "). This iteration will be skipped.")
                    continue

                # Formatting results
                start_pos = 0
                flag_start = False
                fp_state_nb = 4
                for k in range(r.initial, r.initial + len(posterior_list)):
                    curr_index = k - r.initial
                    if (flag_start):
                        if (posterior_list[curr_index] != fp_state_nb):
                            if (k - start_pos < fp_limit_size):
                                fp = GenomicRegion(r.chrom, start_pos, k)
                                footprints.add(fp)
                            flag_start = False
                    else:
                        if (posterior_list[curr_index] == fp_state_nb):
                            flag_start = True
                            start_pos = k
                if (flag_start):
                    fp = GenomicRegion(r.chrom, start_pos, r.final)
                    footprints.add(fp)

            ###################################################################################################
            # HISTONES
            ###################################################################################################

            else:

                # Fetching DNase signal
                if (not group.histone_only):
                    try:
                        if (group.is_atac):
                            dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                                  atac_downstream_ext,
                                                                                  atac_upstream_ext,
                                                                                  atac_forward_shift,
                                                                                  atac_reverse_shift,
                                                                                  dnase_initial_clip, dnase_norm_per,
                                                                                  dnase_slope_per,
                                                                                  group.bias_table,
                                                                                  genome_data.get_genome(),
                                                                                  options.print_raw_signal,
                                                                                  options.print_bc_signal,
                                                                                  options.print_norm_signal,
                                                                                  options.print_slope_signal)
                        else:
                            dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                                  dnase_downstream_ext,
                                                                                  dnase_upstream_ext,
                                                                                  dnase_forward_shift,
                                                                                  dnase_reverse_shift,
                                                                                  dnase_initial_clip, dnase_norm_per,
                                                                                  dnase_slope_per,
                                                                                  group.bias_table,
                                                                                  genome_data.get_genome(),
                                                                                  options.print_raw_signal,
                                                                                  options.print_bc_signal,
                                                                                  options.print_norm_signal,
                                                                                  options.print_slope_signal)
                    except Exception:
                        raise
                        error_handler.throw_warning("FP_DNASE_PROC", add_msg="for region (" + ",".join([r.chrom,
                                                                                                        str(r.initial),
                                                                                                        str(
                                                                                                            r.final)]) + "). This iteration will be skipped.")
                        continue

                # Iterating over histone modifications
                for i in range(0, len(group.histone_file_list)):

                    # Fetching histone signal
                    try:
                        histone_file = group.histone_file_list[i]
                        histone_norm, histone_slope = histone_file.get_signal(r.chrom, r.initial, r.final,
                                                                              histone_downstream_ext,
                                                                              histone_upstream_ext,
                                                                              histone_forward_shift,
                                                                              histone_reverse_shift,
                                                                              histone_initial_clip, histone_norm_per,
                                                                              histone_slope_per, False, False, False,
                                                                              False, False)
                    except Exception:
                        raise
                        error_handler.throw_warning("FP_HISTONE_PROC", add_msg="for region (" + ",".join([r.chrom,
                                                                                                          str(
                                                                                                              r.initial),
                                                                                                          str(
                                                                                                              r.final)]) + ") and histone modification " + histone_file.file_name + ". This iteration will be skipped for this histone.")
                        continue

                    # Formatting sequence
                    try:
                        if (group.histone_only):
                            input_sequence = array([histone_norm, histone_slope]).T
                        else:
                            input_sequence = array([dnase_norm, dnase_slope, histone_norm, histone_slope]).T
                    except Exception:
                        raise
                        error_handler.throw_warning("FP_SEQ_FORMAT", add_msg="for region (" + ",".join(
                            [r.chrom, str(r.initial), str(
                                r.final)]) + ") and histone modification " + histone_file.file_name + ". This iteration will be skipped.")
                        continue

                    # Applying HMM
                    if (flag_multiple_hmms):
                        current_hmm = group.hmm[i]
                    else:
                        current_hmm = group.hmm
                    if (
                            isnan(sum(
                                input_sequence))): continue  # Handling NAN's in signal / hmmlearn throws error TODO ERROR
                    try:
                        posterior_list = current_hmm.predict(input_sequence)
                    except Exception:
                        raise
                        error_handler.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join(
                            [r.chrom, str(r.initial), str(
                                r.final)]) + ") and histone modification " + histone_file.file_name + ". This iteration will be skipped.")
                        continue

                    # Histone-only limit size
                    if (group.histone_only):
                        fp_limit_size = fp_limit_size_histone
                        fp_state_nb = 4
                    else:
                        fp_state_nb = 7

                    # Formatting results
                    start_pos = 0
                    flag_start = False
                    for k in range(r.initial, r.initial + len(posterior_list)):
                        curr_index = k - r.initial
                        if (flag_start):
                            if (posterior_list[curr_index] != fp_state_nb):
                                if (k - start_pos < fp_limit_size):
                                    fp = GenomicRegion(r.chrom, start_pos, k)
                                    footprints.add(fp)
                                flag_start = False
                        else:
                            if (posterior_list[curr_index] == fp_state_nb):
                                flag_start = True
                                start_pos = k
                    if (flag_start):
                        fp = GenomicRegion(r.chrom, start_pos, r.final)
                        footprints.add(fp)

        ###################################################################################################
        # Post-processing
        ###################################################################################################

        # Parameters
        if (group.histone_only):
            fp_limit = fp_limit_size_ext_histone
            fp_ext = fp_ext_histone
            tc_ext = tc_ext_histone
            tcsignal = group.histone_file_list[0]
            tcext1 = histone_downstream_ext
            tcext2 = histone_upstream_ext
            tcshift1 = histone_forward_shift
            tcshift2 = histone_reverse_shift
            tcinitialclip = histone_initial_clip
        else:
            fp_limit = fp_limit_size_ext
            fp_ext = fp_ext
            tc_ext = tc_ext
            tcsignal = group.dnase_file
            if (group.is_atac):
                tcext1 = atac_downstream_ext
                tcext2 = atac_upstream_ext
                tcshift1 = atac_forward_shift
                tcshift2 = atac_reverse_shift
                tcinitialclip = atac_initial_clip
            else:
                tcext1 = dnase_downstream_ext
                tcext2 = dnase_upstream_ext
                tcshift1 = dnase_forward_shift
                tcshift2 = dnase_reverse_shift
                tcinitialclip = dnase_initial_clip

        # Sorting and Merging
        footprints.merge()

        # Overlapping results with original regions
        footprints = footprints.intersect(group.original_regions, mode=OverlapType.ORIGINAL)

        # Extending footprints
        for f in footprints.sequences:
            if (f.final - f.initial < fp_limit):
                f.initial = max(0, f.initial - fp_ext)
                f.final = f.final + fp_ext
            if (f.final - f.initial > 2 * fp_limit):
                mid = (f.initial + f.final) / 2
                f.initial = max(mid - fp_limit, 0)
                f.final = f.final + fp_limit

        # Evaluating TC
        for f in footprints.sequences:
            mid = (f.initial + f.final) / 2
            p1 = max(mid - tc_ext, 0)
            p2 = min(mid + tc_ext, chrom_sizes_dict[f.chrom])
            try:
                tag_count = tcsignal.get_tag_count(f.chrom, p1, p2, tcext1, tcext2, tcshift1, tcshift2, tcinitialclip)
            except Exception:
                tag_count = 0
            f.data = str(int(tag_count))

        ###################################################################################################
        # Writing output
        ###################################################################################################

        # Creating output file
        # output_file_name = options.output_location + options.output_fname + ".bed"
        output_file_name = os.path.join(options.output_location, "{}.bed".format(options.output_fname))
        footprints.write_bed(output_file_name)