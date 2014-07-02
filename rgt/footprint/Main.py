
###################################################################################################
# Libraries
###################################################################################################

# Python
from os import remove, system, getcwd
from sys import exit

# Internal
from .. Util import PassThroughOptionParser, ErrorHandler, HmmData, GenomeData
from .. ExperimentalMatrix import ExperimentalMatrix
from .. GenomicRegionSet import GenomicRegionSet
from signalProcessing import BamFile
from hmm import HMM

# External
from numpy import array
from sklearn.hmm import GaussianHMM

"""
Finds footprints in open chromatin data.

Basic Input:
- Regions (bed) in which to find footprints (i.e. enriched regions or hypersensitivity regions).
- Reads (bam) containing the open chromatin signal for DNase and 1 <= N <= 3 histone modifications.

Dependencies:
- numpy
- scipy
- scikit
- pysam

Authors: Eduardo G. Gusmao.
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

    # Parameters
    current_version = "0.0.1"
    usage_message = ("\n--------------------------------------------------\n"
                     "The 'footprint' program predicts TFBSs given open chromatin data.\n"
                     "In order to use this tools, please type: \n\n"
                     "%prog [options] <experiment_matrix>\n\n"
                     "The <experiment matrix> should contain:\n"
                     "- One region file representing the regions in which the HMM\n"
                     "  will be applied. It should contain 'regions' in the type field\n"
                     "- One DNase aligned reads file (bam) file with 'DNASE' in the name field.\n"
                     "- One to Three histone modification aligned reads file (bam).\n\n"

                     "For more information, please refer to our wiki:\n"
                     "https://code.google.com/p/reg-gen/wiki/RegGen\n"
                     "--------------------------------------------------")
    version_message = "Footprint - Regulatory Analysis Toolbox (RGT). Version: "+str(current_version)

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage = usage_message, version = version_message)

    # Optional Input Options
    parser.add_option("--hmm-file", dest = "hmm_file", type = "string", metavar="FILE_1[,FILE_2,...,FILE_N]", default = None,
                      help = ("List of HMM files separated by comma. If one file only, then this HMM will be "
                              "applied for all histone signals, otherwise, the list must have the same number."
                              "of histone files given. The order of the list should be the order of the"
                              "histones in the input_matrix file. If the argument is not given, then an HMM"
                              "trained with H3K4me3 in K562 will be used."))

    # Parameters Options
    parser.add_option("--organism", dest = "organism", type = "string", metavar="STRING", default = "hg19",
                      help = ("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file. This option is used only if a bigbed output is asked."))

    # Output Options
    parser.add_option("--output-location", dest = "output_location", type = "string", metavar="PATH", default = getcwd(),
                      help = ("Path where the output files will be written."))
    parser.add_option("--footprint-name", dest = "footprint_name", type = "string", metavar="STRING", default = "footprints",
                      help = ("Name of the footprint file (without extension)."))
    parser.add_option("--print-bb", dest = "print_bb", action = "store_true", default = False,
                      help = ("If used, the output will be a bigbed (.bb) file."))

    # Processing Options
    options, arguments = parser.parse_args()

    # Fixed Parameters ################
    region_total_ext = 10000
    fp_state_nb = 7
    fp_limit_size = 50
    ###
    dnase_initial_clip = 1000
    dnase_sg_window_size = 9
    dnase_norm_per = 98
    dnase_slope_per = 98
    dnase_frag_ext = 1
    ###
    histone_initial_clip = 1000
    histone_sg_window_size = 201
    histone_norm_per = 98
    histone_slope_per = 98
    histone_frag_ext = 200
    ###################################
    
    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading input argument
    input_matrix = arguments[0]
    if(not input_matrix): # TODO ERROR
        exit(0)

    # Create experimental matrix
    exp_matrix = ExperimentalMatrix()
    exp_matrix.read(input_matrix)

    ###################################################################################################
    # Reading Regions
    ###################################################################################################

    # Fetching region file
    region_set_list = exp_matrix.get_regionsets()
    if(len(region_set_list) != 1): # TODO ERROR
        exit(0)
    regions = region_set_list[0]
    
    # Sorting + Extending + Merging + Unique
    regions.sort()
    # TODO Extend, merge and unique 
    # Use region_total_ext

    ###################################################################################################
    # Reading Signals
    ###################################################################################################

    # Initialization
    name_list = exp_matrix.names
    type_list = exp_matrix.types
    file_dict = exp_matrix.files
    dnase_label = "DNASE"

    # Fetching signal files
    dnase_file = None
    histone_file_list = []
    for i in range(0,len(name_list)):
        if(type_list[i] == "regions"): continue
        if(name_list[i].upper() == dnase_label): # DNase signal
            dnase_file = BamFile(file_dict[name_list[i]])
            dnase_file.load_sg_coefs(dnase_sg_window_size)
        else: # Histone signal
            histone_file = BamFile(file_dict[name_list[i]])
            histone_file.load_sg_coefs(histone_sg_window_size)
            histone_file_list.append(histone_file)

    ###################################################################################################
    # Creating HMM list
    ###################################################################################################

    # Fetching HMM input
    flag_multiple_hmms = False
    if(options.hmm_file): # Argument is passed

        # Fetching list of HMM files
        hmm_file_list = options.hmm_file.split(",")

        # Verifying HMM application mode (one HMM or multiple HMM files)
        if(len(hmm_file_list) == 1): flag_multiple_hmms = False # One HMM file only
        elif(len(hmm_file_list) == len(histone_file_name_list)): flag_multiple_hmms = True # One HMM file for each histone
        else: # TODO ERROR - Incorrect number of HMMs passed
            exit(0)

    else: # Argument was not passed
        flag_multiple_hmms = False
        hmm_data = HmmData()
        hmm_file_list = [hmm_data.get_default_hmm()]

    # Creating scikit HMM list
    hmm_list = []
    for hmm_file_name in hmm_file_list:

        hmm_scaffold = HMM()
        hmm_scaffold.load_hmm(hmm_file_name)
        scikit_hmm = GaussianHMM(n_components=hmm_scaffold.states, covariance_type="full", 
                                 transmat=array(hmm_scaffold.A), startprob=array(hmm_scaffold.pi))
        scikit_hmm.means_ = array(hmm_scaffold.means)
        scikit_hmm.covars_ = array(hmm_scaffold.covs)
        hmm_list.append(scikit_hmm)

    ###################################################################################################
    # Main Pipeline
    ###################################################################################################

    # Creating output file
    output_file_name = options.output_location+options.footprint_name+".bed"
    output_file = open(output_file_name,"w")

    # Iterating over regions
    for r in regions.sequences:

        # Fetching DNase signal
        try:
            dnase_norm, dnase_slope = dnase_file.get_signal(r.chrom, r.initial, r.final, 
                                      dnase_frag_ext, dnase_initial_clip, dnase_norm_per, dnase_slope_per)
        except Exception: continue # TODO ERROR

        # Iterating over histone modifications
        for i in range(0,len(histone_file_list)):

            # Fetching histone signal
            try:
                histone_file = histone_file_list[i]
                histone_norm, histone_slope = histone_file.get_signal(r.chrom, r.initial, r.final, 
                                              histone_frag_ext, histone_initial_clip, histone_norm_per, histone_slope_per)
            except Exception: continue # TODO ERROR

            # Formatting sequence
            try:
                input_sequence = array([dnase_norm,dnase_slope,histone_norm,histone_slope]).T
            except Exception: continue # TODO ERROR

            # Applying HMM
            if(flag_multiple_hmms): current_hmm = hmm_list[i]
            else: current_hmm = hmm_list[0]
            try:
                posterior_list = current_hmm.predict(input_sequence)
            except Exception: continue # TODO ERROR

            # Writing results
            start_pos = 0
            flag_start = False
            for k in range(r.initial, r.final):
                curr_index = k - r.initial
                if(flag_start):
                    if(posterior_list[curr_index] != fp_state_nb):
                        if(k-start_pos < fp_limit_size): output_file.write("\t".join([r.chrom,str(start_pos),str(k)])+"\n")
                        flag_start = False
                else:
                    if(posterior_list[curr_index] == fp_state_nb):
                        flag_start = True
                        start_pos = k
            if(flag_start): output_file.write("\t".join([r.chrom,str(start_pos),str(r.final)])+"\n")

    # Sorting + Merging + Unique on footprint results
    # TODO

    ###################################################################################################
    # Converting to big bed
    ###################################################################################################

    # Verifying condition to write bb
    if(options.print_bb):

        # Fetching file with chromosome sizes
        genome_data = GenomeData(options.organism)
        chrom_sizes_file = genome_data.get_chromosome_sizes()

        # Converting to big bed
        output_bb_name = options.output_location+options.footprint_name+".bb"
        try:
            system(" ".join("bedToBigBed",output_file_name,chrom_sizes_file,output_bb_name))
            #remove(output_file_name)
        except Exception: pass # TODO ERROR
        

