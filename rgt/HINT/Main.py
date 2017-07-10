###################################################################################################
# Libraries
###################################################################################################

# Python
from os import system, getcwd
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
from diff_footprints import DiffFootprints

# External
import os
import sys
import pysam
from numpy import array, sum, isnan, subtract, absolute
from hmmlearn.hmm import GaussianHMM
from sklearn.externals import joblib

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

    if sys.argv[1] == "--estimate-bias":
        estimate_bias()
    elif sys.argv[1] == "--train-hmm":
        train_hmm()
    elif sys.argv[1] == "--evaluate-footprints":
        evaluate_footprints()
    elif sys.argv[1] == "--print-lines":
        print_lines()
    elif sys.argv[1] == "--print-signal":
        print_signal()
    elif sys.argv[1] == "--create-evidence":
        create_evidence()
    elif sys.argv[1] == "--diff-footprints":
        diff_footprints()
    elif sys.argv[1] == "--atac-footprints":
        atac_footprints()
    elif sys.argv[1] == "--dnase-footprints":
        dnase_footprints()
    elif sys.argv[1] == "--histone-footprints":
        histone_footprints()
    elif sys.argv[1] == "--dnase-histone-footprints":
        dnase_histone_footprints()


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
    parser.add_option("--original-regions", dest="original_regions", type="string",
                      metavar="STRING", default=None,
                      help=("The regions that used to estimate the bias table "
                            "should be bed file containing HS regions or FASTA file containing naked DNA"))
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
    parser.add_option("--output-prefix", dest="output_prefix", type="string", metavar="STRING",
                      default=None)

    # GENERAL Hidden Options
    parser.add_option("--region-total-ext", dest="region_total_ext", type="int", metavar="INT", default=0,
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
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=100,
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
                      metavar="INT", default=50, help=SUPPRESS_HELP)
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

    ########################################################################################################
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
                if bias_correction:
                    if (group.is_atac):
                        #group.hmm = hmm_data.get_default_hmm_atac_bc()
                        hmm_file = hmm_data.get_default_hmm_atac_bc()
                        group.hmm = joblib.load(hmm_file)
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_bc()
                else:
                    if (group.is_atac):
                        #group.hmm = hmm_data.get_default_hmm_atac()
                        hmm_file = hmm_data.get_default_hmm_atac_bc()
                        group.hmm = joblib.load(hmm_file)
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase()
            elif group.histone_only:
                group.hmm = hmm_data.get_default_hmm_histone()
            else:
                if bias_correction:
                    if group.is_atac:
                        group.hmm = hmm_data.get_default_hmm_atac_histone_bc()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_histone_bc()
                else:
                    if group.is_atac:
                        group.hmm = hmm_data.get_default_hmm_atac_histone()
                    else:
                        group.hmm = hmm_data.get_default_hmm_dnase_histone()

    # Creating scikit HMM list
    for group in group_list:

        if (group.flag_multiple_hmms):

            hmm_list = []
            for hmm_file_name in group.hmm:
                scikit_hmm = None
                try:
                    hmm_scaffold = HMM()
                    hmm_scaffold.load_hmm(hmm_file_name)
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

            if group.dnase_only:

                # Fetching DNase signal
                try:
                    if group.is_atac:
                        atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r = \
                            group.dnase_file.get_signal_atac(r.chrom, r.initial, r.final, atac_downstream_ext,
                                                             atac_upstream_ext, atac_forward_shift, atac_reverse_shift,
                                                             atac_initial_clip, atac_norm_per, atac_slope_per,
                                                             group.bias_table, genome_data.get_genome())

                        input_sequence = array([atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r]).T
                    else:
                        dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                              dnase_downstream_ext, dnase_upstream_ext,
                                                                              dnase_forward_shift, dnase_reverse_shift,
                                                                              dnase_initial_clip, dnase_norm_per,
                                                                              dnase_slope_per,
                                                                              group.bias_table,
                                                                              genome_data.get_genome())
                        input_sequence = array([dnase_norm, dnase_slope]).T
                except Exception:
                    error_handler.throw_warning("FP_DNASE_PROC", add_msg="for region (" + ",".join([r.chrom,
                                                                                                    str(r.initial), str(
                            r.final)]) + "). This iteration will be skipped.")
                    continue

                # Applying HMM
                if (isinstance(group.hmm, list)): continue  # TODO ERROR
                try:
                    posterior_list = group.hmm.predict(input_sequence)
                except Exception:
                    error_handler.throw_warning("FP_HMM_APPLIC", add_msg="in region (" + ",".join([r.chrom,
                                                                                                   str(r.initial), str(
                            r.final)]) + "). This iteration will be skipped.")
                    continue

                # Formatting results
                start_pos = 0
                flag_start = False
                if group.is_atac:
                    fp_state_nb = 1
                else:
                    fp_state_nb = 4
                for k in range(r.initial, r.initial + len(posterior_list)):
                    curr_index = k - r.initial
                    if (flag_start):
                        if (posterior_list[curr_index] != fp_state_nb):
                            if (k - start_pos < options.fp_limit_size):
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
                                                                                  genome_data.get_genome())
                        else:
                            dnase_norm, dnase_slope = group.dnase_file.get_signal(r.chrom, r.initial, r.final,
                                                                                  dnase_downstream_ext,
                                                                                  dnase_upstream_ext,
                                                                                  dnase_forward_shift,
                                                                                  dnase_reverse_shift,
                                                                                  dnase_initial_clip, dnase_norm_per,
                                                                                  dnase_slope_per,
                                                                                  group.bias_table,
                                                                                  genome_data.get_genome())
                    except Exception:
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
                    isnan(sum(input_sequence))): continue  # Handling NAN's in signal / hmmlearn throws error TODO ERROR
                    try:
                        posterior_list = current_hmm.predict(input_sequence)
                    except Exception:
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
            reads_file = group.histone_file_list[0]
            downstream_ext = histone_downstream_ext
            upstream_ext = histone_upstream_ext
            forward_shift = histone_forward_shift
            reverse_shift = histone_reverse_shift
            initial_clip = histone_initial_clip
        else:
            fp_limit = fp_limit_size_ext
            fp_ext = fp_ext
            tc_ext = tc_ext
            reads_file = group.dnase_file
            if (group.is_atac):
                downstream_ext = atac_downstream_ext
                upstream_ext = atac_upstream_ext
                forward_shift = atac_forward_shift
                reverse_shift = atac_reverse_shift
                initial_clip = atac_initial_clip
            else:
                downstream_ext = dnase_downstream_ext
                upstream_ext = dnase_upstream_ext
                forward_shift = dnase_forward_shift
                reverse_shift = dnase_reverse_shift
                initial_clip = dnase_initial_clip

        ###################################################################################################
        # Post-processing
        ###################################################################################################
        post_processing(footprints=footprints, original_regions=group.original_regions, fp_limit=fp_limit,
                        fp_ext=fp_ext, genome_data=genome_data, tc_ext=tc_ext,
                        reads_file=reads_file, downstream_ext=downstream_ext, upstream_ext=upstream_ext,
                        forward_shift=forward_shift, reverse_shift=reverse_shift,
                        initial_clip=initial_clip, output_location=options.output_location,
                        output_prefix=options.output_prefix)

def atac_footprints():
    """
    Performs footprints identification from ATAC-seq.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --atac-footprints [options] [reads.bam regions.bed]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    parser.add_option("--hmm-file", dest="hmm_file", type="string",
                      metavar="FILE", default=None,
                      help=("If the argument is not given, then a default HMM will be used."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE_F,FILE_R", default=None,
                      help=("List of files with all possible k-mers (for any k) and their bias estimates. "
                            "Each line should contain a kmer and the bias estimate separated by tab."))
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))

    # ATAC Hidden Options
    parser.add_option("--initial-clip", dest="initial_clip", type="int", metavar="INT", default=50, help=SUPPRESS_HELP)
    parser.add_option("--sg-window-size", dest="sg_window_size", type="int", metavar="INT", default=9,
                      help=SUPPRESS_HELP)
    parser.add_option("--norm-per", dest="norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--slope-per", dest="slope_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int", metavar="INT", default=1,
                      help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int", metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int", metavar="INT", default=-4,
                      help=SUPPRESS_HELP)
    parser.add_option("--bias-correction-k", dest="bias_correction_k", type="int", metavar="INT", default=8,
                      help=SUPPRESS_HELP)

    parser.add_option("--fp-limit-size", dest="fp_limit_size", type="int", metavar="INT", default=30,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-ext", dest="fp_ext", type="int", metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=100, help=SUPPRESS_HELP)
    parser.add_option("--fp-limit", dest="fp_limit", type="int", metavar="INT", default=10, help=SUPPRESS_HELP)

    parser.add_option("--fp-bed-fname", dest="fp_bed_fname", type="string", metavar="STRING", default=None)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string", metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    # Processing Options
    options, arguments = parser.parse_args()

    if len(arguments) < 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    if options.bias_table:
        bias_table_list = options.bias_table.split(",")
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
    hmm = None
    if options.hmm_file:
        hmm = joblib.load(options.hmm_file)
    else:
        hmm_file = hmm_data.get_default_hmm_atac_bc()
        hmm = joblib.load(hmm_file)

    # Initializing result set
    footprints = GenomicRegionSet(options.output_prefix)

    reads_file = GenomicSignal(arguments[0])
    reads_file.load_sg_coefs(options.sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read_bed(arguments[1])

    for r in original_regions:
        atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r = \
            reads_file.get_signal_atac(r.chrom, r.initial, r.final, options.downstream_ext,
                                       options.upstream_ext, options.forward_shift, options.reverse_shift,
                                       options.initial_clip, options.norm_per, options.slope_per,
                                       bias_table, genome_data.get_genome())
        try:
            input_sequence = array([atac_norm_f, atac_slope_f, atac_norm_r, atac_slope_r]).T
        except Exception:
            err.throw_warning("FP_SEQ_FORMAT", add_msg="for region (" + ",".join([r.chrom, str(r.initial), str(
                r.final)]) + "). This iteration will be skipped.")
            continue

        # Applying HMM
        try:
            posterior_list = hmm.predict(input_sequence)
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
            if (flag_start):
                if (posterior_list[curr_index] != fp_state_nb):
                    if (2 < k - start_pos < options.fp_limit_size):
                        fp = GenomicRegion(r.chrom, start_pos, k)
                        footprints.add(fp)
                    flag_start = False
            else:
                if (posterior_list[curr_index] == fp_state_nb):
                    flag_start = True
                    start_pos = k
        if (flag_start):
            if (2 < r.initial + len(posterior_list) - start_pos < options.fp_limit_size):
                fp = GenomicRegion(r.chrom, start_pos, r.final)
                footprints.add(fp)

    ###################################################################################################
    # Post-processing
    ###################################################################################################

    post_processing(footprints=footprints, original_regions=original_regions, fp_limit=options.fp_limit,
                    fp_ext=options.fp_ext, genome_data=genome_data, tc_ext=options.tc_ext,
                    reads_file=reads_file, downstream_ext=options.downstream_ext, upstream_ext=options.upstream_ext,
                    forward_shift=options.forward_shift, reverse_shift=options.reverse_shift,
                    initial_clip=options.initial_clip, output_location=options.output_location,
                    output_prefix=options.output_prefix)
    # TODO
    exit(0)


def output_bed_file(chrom, start, end, states, output_fname, fp_state):
    #state_dict = dict([(0, "0"), (1, "1"), (2, "2"), (3, "3"), (4, "4")])
    #color_dict = dict([(0, "0,0,0"), (1, "102,0,51"), (2, "153,0,153"), (3, "102,0,204"), (4, "0,0,255")])
    state_dict = dict([(0, "0"), (1, "1"), (2, "2"), (3, "3"), (4, "4"), (5, "5"), (6, "6")])
    color_dict = dict([(0, "0,0,0"), (1, "102,0,51"), (2, "153,0,153"), (3, "102,0,204"), (4, "0,0,255"),
                       (5, "51,153,255"), (6, "102,255,255")])
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


def dnase_footprints():
    """
    Performs footprints identification from DNase-seq.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --dnase-footprints [options] [reads.bam regions.bed]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    parser.add_option("--hmm-file", dest="hmm_file", type="string",
                      metavar="FILE", default=None,
                      help=("If the argument is not given, then a default HMM will be used."))
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE_F,FILE_R", default=None,
                      help=("List of files with all possible k-mers (for any k) and their bias estimates. "
                            "Each line should contain a kmer and the bias estimate separated by tab."))
    parser.add_option("--bias-correction", dest="bias_correction",
                      action="store_true", default=False,
                      help=("Applies DNase-seq cleavage bias correction with default k-mer bias "
                            "estimates (FAST HINT-BC). Please set the correct --bias-type "
                            "option that matches your experimental settings."))
    parser.add_option("--bias-type", dest="bias_type", type="string",
                      metavar="STRING", default="SH",
                      help=("Type of protocol used to generate the DNase-seq. "
                            "Available options are: 'SH' (DNase-seq single-hit protocol), 'DH' "
                            "(DNase-seq double-hit protocol). "))

    # DNASE Hidden Options
    parser.add_option("--initial-clip", dest="initial_clip", type="int", metavar="INT", default=1000,
                      help=SUPPRESS_HELP)
    parser.add_option("--sg-window-size", dest="sg_window_size", type="int", metavar="INT", default=9,
                      help=SUPPRESS_HELP)
    parser.add_option("--norm-per", dest="norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--slope-per", dest="slope_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int", metavar="INT", default=1,
                      help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--bias-correction-k", dest="bias_correction_k", type="int", metavar="INT", default=6,
                      help=SUPPRESS_HELP)

    parser.add_option("--region-total-ext", dest="region_total_ext", type="int", metavar="INT", default=10000,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size", dest="fp_limit_size", type="int", metavar="INT", default=50,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-ext", dest="fp_ext", type="int", metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=100, help=SUPPRESS_HELP)
    parser.add_option("--fp-limit", dest="fp_limit", type="int", metavar="INT", default=10, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    # Processing Options
    options, arguments = parser.parse_args()

    if len(arguments) < 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    bias_table = None
    if options.bias_correction:
        if options.bias_table:
            bias_table_list = options.bias_table.split(",")
            bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                                table_file_name_R=bias_table_list[1])
        else:
            if options.bias_type == 'SH':
                table_F = hmm_data.get_default_bias_table_F_SH()
                table_R = hmm_data.get_default_bias_table_R_SH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)
            elif options.bias_type == 'DH':
                table_F = hmm_data.get_default_bias_table_F_DH()
                table_R = hmm_data.get_default_bias_table_R_DH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    hmm_file = None
    if options.hmm_file:
        hmm_file = options.hmm_file
    else:
        if options.bias_correction:
            hmm_file = hmm_data.get_default_hmm_dnase_bc()
        else:
            hmm_file = hmm_data.get_default_hmm_dnase()

    scikit_hmm = None
    try:
        hmm_scaffold = HMM()
        hmm_scaffold.load_hmm(hmm_file)
        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full")
        scikit_hmm.startprob_ = array(hmm_scaffold.pi)
        scikit_hmm.transmat_ = array(hmm_scaffold.A)
        scikit_hmm.means_ = array(hmm_scaffold.means)
        scikit_hmm.covars_ = array(hmm_scaffold.covs)
    except Exception:
        err.throw_error("FP_HMM_FILES")

    # Initializing result set
    footprints = GenomicRegionSet(options.output_prefix)

    reads_file = GenomicSignal(arguments[0])
    reads_file.load_sg_coefs(options.sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read_bed(arguments[1])

    regions = deepcopy(original_regions)
    regions.extend(int(options.region_total_ext / 2), int(options.region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        dnase_norm, dnase_slope = reads_file.get_signal(r.chrom, r.initial, r.final, options.downstream_ext,
                                                        options.upstream_ext, options.forward_shift,
                                                        options.reverse_shift, options.initial_clip, options.norm_per,
                                                        options.slope_per, bias_table, genome_data.get_genome())

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
            if (flag_start):
                if (posterior_list[curr_index] != fp_state_nb):
                    if (k - start_pos < options.fp_limit_size):
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
    post_processing(footprints=footprints, original_regions=original_regions, fp_limit=options.fp_limit,
                    fp_ext=options.fp_ext, genome_data=genome_data, tc_ext=options.tc_ext,
                    reads_file=reads_file, downstream_ext=options.downstream_ext, upstream_ext=options.upstream_ext,
                    forward_shift=options.forward_shift, reverse_shift=options.reverse_shift,
                    initial_clip=options.initial_clip, output_location=options.output_location,
                    output_prefix=options.output_prefix)
    # TODO
    exit(0)


def histone_footprints():
    """
    Performs footprints identification from Histone data.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --histone-footprints [options] [reads.bam regions.bed]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    parser.add_option("--hmm-file", dest="hmm_file", type="string",
                      metavar="FILE", default=None,
                      help=("If the argument is not given, then a default HMM will be used."))
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))

    parser.add_option("--initial-clip", dest="initial_clip", type="int", metavar="INT", default=1000,
                      help=SUPPRESS_HELP)
    parser.add_option("--sg-window-size", dest="sg_window_size", type="int", metavar="INT", default=201,
                      help=SUPPRESS_HELP)
    parser.add_option("--norm-per", dest="norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--slope-per", dest="slope_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int", metavar="INT", default=200,
                      help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)

    parser.add_option("--region-total-ext", dest="region_total_ext", type="int", metavar="INT", default=10000,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size", dest="fp_limit_size", type="int", metavar="INT", default=2000,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-ext", dest="fp_ext", type="int", metavar="INT", default=50, help=SUPPRESS_HELP)
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=500, help=SUPPRESS_HELP)
    parser.add_option("--fp-limit", dest="fp_limit", type="int", metavar="INT", default=200, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    # Processing Options
    options, arguments = parser.parse_args()

    if len(arguments) < 2:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    hmm_file = None
    if options.hmm_file:
        hmm_file = options.hmm_file
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
    except Exception:
        err.throw_error("FP_HMM_FILES")

    # Initializing result set
    footprints = GenomicRegionSet(options.output_prefix)

    reads_file = GenomicSignal(arguments[0])
    reads_file.load_sg_coefs(options.sg_window_size)

    original_regions = GenomicRegionSet("regions")
    original_regions.read_bed(arguments[1])

    regions = deepcopy(original_regions)
    regions.extend(int(options.region_total_ext / 2), int(options.region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        histone_norm, histone_slope = reads_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                        downstream_ext=options.downstream_ext,
                                                        upstream_ext=options.upstream_ext,
                                                        forward_shift=options.forward_shift,
                                                        reverse_shift=options.reverse_shift,
                                                        initial_clip=options.initial_clip,
                                                        per_norm=options.norm_per,
                                                        per_slope=options.slope_per,
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
            if (flag_start):
                if (posterior_list[curr_index] != fp_state_nb):
                    if (k - start_pos < options.fp_limit_size):
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
    post_processing(footprints=footprints, original_regions=original_regions, fp_limit=options.fp_limit,
                    fp_ext=options.fp_ext, genome_data=genome_data, tc_ext=options.tc_ext,
                    reads_file=reads_file, downstream_ext=options.downstream_ext, upstream_ext=options.upstream_ext,
                    forward_shift=options.forward_shift, reverse_shift=options.reverse_shift,
                    initial_clip=options.initial_clip, output_location=options.output_location,
                    output_prefix=options.output_prefix)
    # TODO
    exit(0)


def dnase_histone_footprints():
    """
    Performs footprints identification from DNase-seq combining with histone modification data.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --dnase-histone-footprints [options] [DNase.bam histone.bam regions.bed]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    parser.add_option("--hmm-file", dest="hmm_file", type="string",
                      metavar="FILE", default=None,
                      help=("If the argument is not given, then a default HMM will be used."))
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE_F,FILE_R", default=None,
                      help=("List of files with all possible k-mers (for any k) and their bias estimates. "
                            "Each line should contain a kmer and the bias estimate separated by tab."))
    parser.add_option("--bias-correction", dest="bias_correction",
                      action="store_true", default=False,
                      help=("Applies DNase-seq cleavage bias correction with default k-mer bias "
                            "estimates (FAST HINT-BC). Please set the correct --bias-type "
                            "option that matches your experimental settings."))
    parser.add_option("--bias-type", dest="bias_type", type="string",
                      metavar="STRING", default="SH",
                      help=("Type of protocol used to generate the DNase-seq. "
                            "Available options are: 'SH' (DNase-seq single-hit protocol), 'DH' "
                            "(DNase-seq double-hit protocol). "))


    parser.add_option("--dnase-initial-clip", dest="dnase_initial_clip", type="int", metavar="INT", default=1000,
                      help=SUPPRESS_HELP)
    parser.add_option("--dnase-sg-window-size", dest="dnase_sg_window_size", type="int", metavar="INT", default=9,
                      help=SUPPRESS_HELP)
    parser.add_option("--dnase-norm-per", dest="dnase_norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--dnase-slope-per", dest="dnase_slope_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--dnase-downstream-ext", dest="dnase_downstream_ext", type="int", metavar="INT", default=1,
                      help=SUPPRESS_HELP)
    parser.add_option("--dnase-upstream-ext", dest="dnase_upstream_ext", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--dnase-forward-shift", dest="dnase_forward_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--dnase-reverse-shift", dest="dnase_reverse_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)


    parser.add_option("--histone-initial-clip", dest="histone_initial_clip", type="int", metavar="INT", default=1000,
                      help=SUPPRESS_HELP)
    parser.add_option("--histone-sg-window-size", dest="histone_sg_window_size", type="int", metavar="INT", default=201,
                      help=SUPPRESS_HELP)
    parser.add_option("--histone-norm-per", dest="histone_norm_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--histone-slope-per", dest="histone_slope_per", type="float", metavar="INT", default=98, help=SUPPRESS_HELP)
    parser.add_option("--histone-downstream-ext", dest="histone_downstream_ext", type="int", metavar="INT", default=200,
                      help=SUPPRESS_HELP)
    parser.add_option("--histone-upstream-ext", dest="histone_upstream_ext", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--histone-forward-shift", dest="histone_forward_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--histone-reverse-shift", dest="histone_reverse_shift", type="int", metavar="INT", default=0, help=SUPPRESS_HELP)


    parser.add_option("--region-total-ext", dest="region_total_ext", type="int", metavar="INT", default=10000,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-limit-size", dest="fp_limit_size", type="int", metavar="INT", default=50,
                      help=SUPPRESS_HELP)
    parser.add_option("--fp-ext", dest="fp_ext", type="int", metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--tc-ext", dest="tc_ext", type="int", metavar="INT", default=100, help=SUPPRESS_HELP)
    parser.add_option("--fp-limit", dest="fp_limit", type="int", metavar="INT", default=10, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))


    # Processing Options
    options, arguments = parser.parse_args()

    if len(arguments) < 3:
        err.throw_error("ME_FEW_ARG", add_msg="You must specify DNase and histone reads and regions file.")

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    bias_table = None
    if options.bias_correction:
        if options.bias_table:
            bias_table_list = options.bias_table.split(",")
            bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                                table_file_name_R=bias_table_list[1])
        else:
            if options.bias_type == 'SH':
                table_F = hmm_data.get_default_bias_table_F_SH()
                table_R = hmm_data.get_default_bias_table_R_SH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)
            elif options.bias_type == 'DH':
                table_F = hmm_data.get_default_bias_table_F_DH()
                table_R = hmm_data.get_default_bias_table_R_DH()
                bias_table = BiasTable().load_table(table_file_name_F=table_F, table_file_name_R=table_R)

    ###################################################################################################
    # Creating HMMs
    ###################################################################################################
    if options.hmm_file:
        hmm_file = options.hmm_file
    else:
        if options.bias_correction:
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
    except Exception:
        err.throw_error("FP_HMM_FILES")


    # Initializing result set
    footprints = GenomicRegionSet(options.output_prefix)

    dnase_reads_file = GenomicSignal(arguments[0])
    dnase_reads_file.load_sg_coefs(options.dnase_sg_window_size)

    histone_reads_file_list = list()
    for f in arguments[1:-1]:
        histone_reads_file_list.append(GenomicSignal(f))

    original_regions = GenomicRegionSet("regions")
    original_regions.read_bed(arguments[-1])

    regions = deepcopy(original_regions)
    regions.extend(int(options.region_total_ext / 2), int(options.region_total_ext / 2))  # Extending
    regions.merge()

    for r in regions:
        dnase_norm, dnase_slope = dnase_reads_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                              downstream_ext=options.dnase_downstream_ext,
                                                              upstream_ext=options.dnase_upstream_ext,
                                                              forward_shift=options.dnase_forward_shift,
                                                              reverse_shift=options.dnase_reverse_shift,
                                                              initial_clip=options.dnase_initial_clip,
                                                              per_norm=options.dnase_norm_per,
                                                              per_slope=options.dnase_slope_per,
                                                              bias_table=bias_table,
                                                              genome_file_name=genome_data.get_genome())
        # Iterating over histone modifications
        for i in range(0, len(histone_reads_file_list)):
            # Fetching histone signal
            histone_file = None
            try:
                histone_file = histone_reads_file_list[i]
                histone_file.load_sg_coefs(options.histone_sg_window_size)
                histone_norm, histone_slope = histone_file.get_signal(ref=r.chrom, start=r.initial, end=r.final,
                                                                      downstream_ext=options.histone_downstream_ext,
                                                                      upstream_ext=options.histone_upstream_ext,
                                                                      forward_shift=options.histone_forward_shift,
                                                                      reverse_shift=options.histone_reverse_shift,
                                                                      initial_clip=options.histone_initial_clip,
                                                                      per_norm=options.histone_norm_per,
                                                                      per_slope=options.histone_slope_per,
                                                                      genome_file_name=genome_data.get_genome())
            except Exception:
                err.throw_warning("FP_HISTONE_PROC", add_msg="for region (" + ",".join([r.chrom,str(r.initial),
                                str(r.final)]) + ") and histone modification " + histone_file.file_name +
                                ". This iteration will be skipped for this histone.")
                continue

            # Formatting sequence
            input_sequence = array([dnase_norm, dnase_slope, histone_norm, histone_slope]).T
            if (isnan(sum(input_sequence))): continue  # Handling NAN's in signal / hmmlearn throws error TODO ERROR
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
                if (flag_start):
                    if (posterior_list[curr_index] != fp_state_nb):
                        if (k - start_pos < options.fp_limit_size):
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
    post_processing(footprints=footprints, original_regions=original_regions, fp_limit=options.fp_limit,
                    fp_ext=options.fp_ext, genome_data=genome_data, tc_ext=options.tc_ext,
                    reads_file=dnase_reads_file, downstream_ext=options.dnase_downstream_ext,
                    upstream_ext=options.dnase_upstream_ext,
                    forward_shift=options.dnase_forward_shift, reverse_shift=options.dnase_reverse_shift,
                    initial_clip=options.dnase_initial_clip, output_location=options.output_location,
                    output_prefix=options.output_prefix)

    # TODO
    exit(0)


def post_processing(footprints, original_regions, fp_limit, fp_ext, genome_data, tc_ext, reads_file,
                    downstream_ext, upstream_ext, forward_shift, reverse_shift, initial_clip, output_location, output_prefix):
    # Overlapping results with original regions
    footprints_overlap = footprints.intersect(original_regions, mode=OverlapType.ORIGINAL)

    # Sorting and Merging
    footprints_overlap.merge()

    # Extending footprints
    for f in footprints_overlap.sequences:
        if (f.final - f.initial < fp_limit):
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
                                                 forward_shift=forward_shift,reverse_shift=reverse_shift,
                                                 initial_clip=initial_clip)
        except Exception:
            tag_count = 0
        f.data = str(int(tag_count))

    ###################################################################################################
    # Writing output
    ###################################################################################################
    output_file_name = os.path.join(output_location, "{}.bed".format(output_prefix))
    footprints_overlap.write_bed(output_file_name)

    # the number of reads
    lines = pysam.idxstats(reads_file.file_name).splitlines()
    num_reads = sum(map(int, [x.split("\t")[2] for x in lines if not x.startswith("#")]))

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

def estimate_bias():
    """
    Performs bias table estimation.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --estimate-bias [options] [reads.bam regions.bed]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--bias-type", dest="bias_type", type="string",
                      metavar="STRING", default=None,
                      help=("The methods that used to estimate the bias table "
                            "Available options are: 'KMER', 'PPM' and 'LSlim"))
    parser.add_option("--reads-file", dest="reads_file", type="string",
                      metavar="FILE", default=None,
                      help=("The BAM file containing the DNase-seq or ATAC-seq reads"))
    parser.add_option("--regions-file", dest="regions_file", type="string",
                      metavar="FILE", default=None,
                      help=("BED file of the regions you want to estimate the bias"))
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int",
                      metavar="INT", default=1,
                      help=("Number of bps to extend towards the downstream region"))
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to extend towards the upstream region"))
    parser.add_option("--forward-shift", dest="forward_shift", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to shift the reads aligned to the forward strand"))
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to shift the reads aligned to the reverse strand"))
    parser.add_option("--k-nb", dest="k_nb", type="int",
                      metavar="INT", default=6,
                      help=("K-mers size used to estimate the bias"))

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    options, arguments = parser.parse_args()

    genome_data = GenomeData(options.organism)

    bias_table = BiasTable(reads_file=options.reads_file, regions_file=options.regions_file,
                           genome_file=genome_data.get_genome(), k_nb=options.k_nb,
                           forward_shift=options.forward_shift, reverse_shift=options.reverse_shift,
                           bias_type=options.bias_type, output_location=options.output_location)

    table = bias_table.estimate()
    bias_table.write(options.output_prefix, table)

    # TODO
    exit(0)


def train_hmm():
    """
    Performs Hidden Markov Model training.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --train-hmm [options] [file1 file2]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--organism", dest="organism", type="string",
                      metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--reads-file", dest="reads_file", type="string",
                      metavar="FILE", default=None,
                      help=("The BAM file containing the DNase-seq or ATAC-seq reads"))
    parser.add_option("--chrom", dest="chrom", type="string",
                      metavar="STRING", default="chr1",
                      help=("The name of chromosome used to train HMM"))
    parser.add_option("--start", dest="start", type="int",
                      metavar="INT", default=211428000)
    parser.add_option("--end", dest="end", type="int",
                      metavar="INT", default=211438000)
    parser.add_option("--annotate-file", dest="annotate_file", type="string", metavar="STRING",
                      default=None,
                      help=("A annotate file containing all the states."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE1_F,FILE1_R[;...;FILEM_F,FILEM_R]", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates. Each input group"
                            "should have two files: one for the forward and one for the reverse strand. "
                            "Each line should contain a kmer and the bias estimate separated by tab. "
                            "Leave an empty set for histone-only analysis groups. Eg. FILE1;;FILE3."))
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int",
                      metavar="INT", default=1,
                      help=("Number of bps to extend towards the downstream region"))
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to extend towards the upstream region"))
    parser.add_option("--forward-shift", dest="forward_shift", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to shift the reads aligned to the forward strand"))
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int",
                      metavar="INT", default=0,
                      help=("Number of bps to shift the reads aligned to the reverse strand"))
    parser.add_option("--k-nb", dest="k_nb", type="int",
                      metavar="INT", default=6,
                      help=("K-mers size used to estimate the bias"))

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))
    parser.add_option("--print-bed-file", dest="print_bed_file",
                      action="store_true", default=False,
                      help=("If used, HINT will output the bed file containing "
                            "the annotated regions, so that you can visualize "
                            "you HMM annotation for potential errors"))

    options, arguments = parser.parse_args()

    # If HINT is required to train a hidden Markov model (HMM)
    hmm_model = TrainHMM(reads_file=options.reads_file, annotate_file=options.annotate_file,
                         print_bed_file=options.print_bed_file, output_location=options.output_location,
                         output_prefix=options.output_prefix, bias_table=options.bias_table,
                         organism=options.organism, k_nb=options.k_nb, chrom=options.chrom,
                         start=options.start, end=options.end,
                         downstream_ext=options.downstream_ext, upstream_ext=options.upstream_ext,
                         forward_shift=options.forward_shift, reverse_shift=options.reverse_shift)
    hmm_model.train()

    # TODO
    exit(0)


def evaluate_footprints():
    """
    Performs prediction evaluation.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --evaluate-footprints [options] [file1 file2]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--tfbs-file", dest="tfbs_file", type="string",
                      metavar="FILE", default=None,
                      help=("A bed file containing all motif-predicted binding sites (MPBSs)."
                            "The values in the bed SCORE field will be used to rank the MPBSs."
                            "The extension must be (.bed)."))
    parser.add_option("--footprint-file", dest="footprint_file", type="string",
                      metavar="FILE1,FILE2,FILE3,FILE4...", default=None,
                      help=("The bed files containing the footprint predictions."
                            "The extension must be (.bed)."))
    parser.add_option("--footprint-name", dest="footprint_name", type="string",
                      metavar="NAME1,NAME2,NAME3,NAME4...", default=None,
                      help=("The name used to distinguish different footprint files."
                            "The number of methods name must be consistent with that of footprint files"))
    parser.add_option("--footprint-type", dest="footprint_type", type="string",
                      metavar="TYPE1,TYPE2,TYPE3,TYPE4...", default=None,
                      help=("The methods type used to predicted the footprint."
                            "Must be SC (site centric) or SEG (segmentation approach)"))

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))
    parser.add_option("--print-roc-curve", dest="print_roc_curve",
                      action="store_true", default=False,
                      help=("If used, HINT will print the receiver operating characteristic curve."))
    parser.add_option("--print-pr-curve", dest="print_pr_curve",
                      action="store_true", default=False,
                      help=("If used, HINT will print the precision recall curve."))

    options, arguments = parser.parse_args()

    evaluation = Evaluation(options.tfbs_file, options.footprint_file,
                            options.footprint_name, options.footprint_type,
                            options.output_location, options.output_prefix,
                            options.print_roc_curve, options.print_pr_curve)
    evaluation.chip_evaluate()

    # TODO
    exit(0)


def print_lines():
    """
    Performs lines plotting.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --print-lines [options] [file1 file2]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--organism", dest="organism", type="string",
                      metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--reads-file", dest="reads_file", type="string",
                      metavar="FILE", default=None,
                      help=("A bam file containing all the DNase-seq or ATAC-seq reads."))
    parser.add_option("--motif-file", dest="motif_file", type="string",
                      metavar="FILE", default=None,
                      help=("A bed file containing all motif-predicted binding sites (MPBSs)."
                            "The extension must be (.bed)."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE1_F,FILE1_R[;...;FILEM_F,FILEM_R]", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates. Each input group"
                            "should have two files: one for the forward and one for the reverse strand. "
                            "Each line should contain a kmer and the bias estimate separated by tab. "
                            "Leave an empty set for histone-only analysis groups. Eg. FILE1;;FILE3."))
    parser.add_option("--window-size", dest="window_size", type="int",
                      metavar="INT", default=100)

    # Hidden Options
    parser.add_option("--initial-clip", dest="initial_clip", type="int",
                      metavar="INT", default=50, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int",
                      metavar="INT", default=1, help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int",
                      metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int",
                      metavar="INT", default=-4, help=SUPPRESS_HELP)
    parser.add_option("--k-nb", dest="k_nb", type="int",
                      metavar="INT", default=6, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))
    parser.add_option("--print-corrected-plot", dest="print_corrected_plot",
                      action="store_true", default=False,
                      help=("If used, the plot containing bias signal for forward and reverse strand,"
                            "corrected and uncorredted signal will be printed."))
    parser.add_option("--print-strand-plot", dest="print_strand_plot",
                      action="store_true", default=False,
                      help=("If used, the plot containing the aggregated signal,"
                            "and strand-specific signal will be printed."))

    options, arguments = parser.parse_args()

    plot = Plot(options.organism, options.reads_file, options.motif_file, options.window_size,
                options.downstream_ext, options.upstream_ext, options.forward_shift, options.reverse_shift,
                options.initial_clip, options.bias_table, options.k_nb,
                options.output_location, options.output_prefix)

    if options.print_corrected_plot:
        plot.line()
    if options.print_strand_plot:
        plot.line1()

    # TODO
    exit(0)


def print_signal():
    """
    Performs signals output.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --print-signal [options] [file1 file2]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--organism", dest="organism", type="string",
                      metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--reads-file", dest="reads_file", type="string",
                      metavar="STRING", default=None,
                      help=("A bam file containing all the DNase-seq reads."))
    parser.add_option("--regions-file", dest="regions_file", type="string",
                      metavar="STRING", default=None,
                      help=("A bed file containing all the interested regions."))
    parser.add_option("--bias-table", dest="bias_table", type="string",
                      metavar="FILE1_F,FILE1_R[;...;FILEM_F,FILEM_R]", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates. Each input group"
                            "should have two files: one for the forward and one for the reverse strand. "
                            "Each line should contain a kmer and the bias estimate separated by tab. "
                            "Leave an empty set for histone-only analysis groups. Eg. FILE1;;FILE3."))

    # Hidden Options
    parser.add_option("--initial-clip", dest="initial_clip", type="int",
                      metavar="INT", default=50, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int",
                      metavar="INT", default=1, help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int",
                      metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int",
                      metavar="INT", default=-4, help=SUPPRESS_HELP)
    parser.add_option("--k-nb", dest="k_nb", type="int",
                      metavar="INT", default=6, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))
    parser.add_option("--raw-signal", dest="raw_signal",
                      action="store_true", default=False,
                      help=("If used, it will print the base overlap (raw) signals from DNase-seq "
                            " or ATAC-seq data. "))
    parser.add_option("--bc-signal", dest="bc_signal",
                      action="store_true", default=False,
                      help=("If used, it will print the base overlap (bias corrected) signals from DNase-seq "
                            " or ATAC-seq data. "))
    parser.add_option("--bigWig", dest="bigWig",
                      action="store_true", default=False,
                      help=("If used, all .wig files will be converted to .bw files"))

    options, arguments = parser.parse_args()

    ########################################################################################################
    # Global class initialization
    genome_data = GenomeData(options.organism)
    hmm_data = HmmData()

    raw_signal_file = None
    bc_signal_file = None
    # Output wig signal
    if options.raw_signal:
        raw_signal_file = os.path.join(options.output_location, "{}.raw.wig".format(options.output_prefix))
        #system("touch " + raw_signal_file + " | echo -n "" > " + raw_signal_file)
        system("touch " + raw_signal_file)
    if options.bc_signal:
        bc_signal_file = os.path.join(options.output_location, "{}.bc.wig".format(options.output_prefix))
        #system("touch " + bc_signal_file + " | echo -n "" > " + bc_signal_file)
        system("touch " + bc_signal_file)

    signal = GenomicSignal(options.reads_file)
    signal.load_sg_coefs(slope_window_size=9)
    regions = GenomicRegionSet("Interested regions")
    regions.read_bed(options.regions_file)

    ###################################################################################################
    # Fetching Bias Table
    ###################################################################################################
    if options.bias_table:
        bias_table_list = options.bias_table.split(",")
        bias_table = BiasTable().load_table(table_file_name_F=bias_table_list[0],
                                            table_file_name_R=bias_table_list[1])
    else:
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        bias_table = BiasTable().load_table(table_file_name_F=table_F,
                                            table_file_name_R=table_R)

    for region in regions:
        signal.print_signal(ref=region.chrom, start=region.initial, end=region.final,
                            downstream_ext=options.downstream_ext,
                            bias_table=bias_table,
                            upstream_ext=options.upstream_ext,
                            forward_shift=options.forward_shift,
                            reverse_shift=options.reverse_shift,
                            genome_file_name=genome_data.get_genome(),
                            raw_signal_file=raw_signal_file,
                            bc_signal_file=bc_signal_file)

    chrom_sizes_file = genome_data.get_chromosome_sizes()
    if options.bigWig:
        if options.raw_signal:
            bw_filename = os.path.join(options.output_location, "{}.raw.bw".format(options.output_prefix))
            system(" ".join(["wigToBigWig", raw_signal_file, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(raw_signal_file)

        if options.bc_signal:
            bw_filename = os.path.join(options.output_location, "{}.bc.bw".format(options.output_prefix))
            system(" ".join(["wigToBigWig", bc_signal_file, chrom_sizes_file, bw_filename, "-verbose=0"]))
            os.remove(bc_signal_file)
    # TODO
    exit(0)


def create_evidence():
    """
    Performs evidence creation.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --create-evidence [options] [file1 file2]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--mpbs-file", dest="mpbs_file", type="string",
                      metavar="FILE", default=None,
                      help=("motif predicted binding sites file"))
    parser.add_option("--tfbs-summit-file", dest="tfbs_summit_file", type="string",
                      metavar="FILE", default=None,
                      help=("the ChIP-seq peak files"))
    parser.add_option("--peak-ext", dest="peak_ext", type="int",
                      metavar="INT", default=100,
                      help=("The number used to extend the ChIP-seq summit"))

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    options, arguments = parser.parse_args()

    evidence = Evidence(options.mpbs_file, options.tfbs_summit_file, options.peak_ext,
                        options.output_location, options.output_prefix)
    evidence.create_file()

    # TODO
    exit(0)


def diff_footprints():
    """
    Performs footprints different analysis.

    Authors: Eduardo G. Gusmao, Zhijian Li.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --diff-footprints [options] [file1 file2 ...]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Input Options
    parser.add_option("--organism", dest="organism", type="string",
                      metavar="STRING", default="hg19",
                      help=("Organism considered on the analysis. Check our full documentation for all available "
                            "options. All default files such as genomes will be based on the chosen organism "
                            "and the data.config file."))
    parser.add_option("--mpbs-file", dest="mpbs_file", type="string",
                      metavar="FILE", default=None,
                      help=("motif predicted binding sites file, must be .bed file"))
    parser.add_option("--reads-file1", dest="reads_file1", type="string",
                      metavar="FILE", default=None,
                      help=("The BAM file containing the DNase-seq or ATAC-seq reads"))
    parser.add_option("--reads-file2", dest="reads_file2", type="string",
                      metavar="FILE", default=None,
                      help=("The BAM file containing the DNase-seq or ATAC-seq reads"))
    parser.add_option("--bias-table1", dest="bias_table1", type="string",
                      metavar="FILE1_F,FILE1_R", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates."))
    parser.add_option("--bias-table2", dest="bias_table2", type="string",
                      metavar="FILE2_F,FILE2_R", default=None,
                      help=("List of files (for each input group; separated by semicolon) with all "
                            "possible k-mers (for any k) and their bias estimates."))
    parser.add_option("--window-size", dest="window_size", type="int",
                      metavar="INT", default=200)
    parser.add_option("--motif-ext", dest="motif_ext", type="int",
                      metavar="INT", default=20)
    parser.add_option("--min-value", dest="min_value", type="int",
                      metavar="INT", default=1)
    parser.add_option("--tc1", dest="tc1", type="float",
                      metavar="INT", default=None)
    parser.add_option("--tc2", dest="tc2", type="float",
                      metavar="INT", default=None)

    # Hidden Options
    parser.add_option("--initial-clip", dest="initial_clip", type="int",
                      metavar="INT", default=1000, help=SUPPRESS_HELP)
    parser.add_option("--downstream-ext", dest="downstream_ext", type="int",
                      metavar="INT", default=1, help=SUPPRESS_HELP)
    parser.add_option("--upstream-ext", dest="upstream_ext", type="int",
                      metavar="INT", default=0, help=SUPPRESS_HELP)
    parser.add_option("--forward-shift", dest="forward_shift", type="int",
                      metavar="INT", default=5, help=SUPPRESS_HELP)
    parser.add_option("--reverse-shift", dest="reverse_shift", type="int",
                      metavar="INT", default=-4, help=SUPPRESS_HELP)
    parser.add_option("--k-nb", dest="k_nb", type="int",
                      metavar="INT", default=6, help=SUPPRESS_HELP)

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string",
                      metavar="PATH", default=getcwd(),
                      help=("Path where the output bias table files will be written."))
    parser.add_option("--output-prefix", dest="output_prefix", type="string",
                      metavar="STRING", default=getcwd(),
                      help=("The prefix for results files."))

    options, arguments = parser.parse_args()

    hmm_data = HmmData()

    if options.bias_table1 == None and options.bias_table2 == None:
        table_F = hmm_data.get_default_bias_table_F_ATAC()
        table_R = hmm_data.get_default_bias_table_R_ATAC()
        options.bias_table1 = table_F + "," + table_R
        options.bias_table2 = table_F + "," + table_R


    diff_footprints = DiffFootprints(options.organism, options.mpbs_file, options.reads_file1,
                                     options.reads_file2, options.bias_table1, options.bias_table2,
                                     options.window_size, options.motif_ext, options.min_value,
                                     options.initial_clip, options.downstream_ext, options.upstream_ext,
                                     options.forward_shift, options.reverse_shift, options.k_nb,
                                     options.output_location, options.output_prefix)

    diff_footprints.diff(options.tc1, options.tc2)

    # TODO
    exit(0)
