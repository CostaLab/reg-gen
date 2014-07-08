
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
from glob import glob

# Internal
from Motif import Motif
from .. Util import PassThroughOptionParser, ErrorHandler, MotifData
from .. ExperimentalMatrix import ExperimentalMatrix
from .. GeneSet import GeneSet
from .. GenomicRegionSet import GenomicRegionSet

# External
from Bio import motifs
from Bio.Seq import Seq

"""
Contains functions to common motif analyses.

Basic Input:
- XXX TODO

Dependencies:
- XXX TODO
- bedToBigBed, XXXX scripts in $PATH (if the option is used)

Authors: Eduardo G. Gusmao.
"""

def main():
    """
    Main function that redirects tool usage.

    Keyword arguments: None
        
    Return: None
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    current_version = "0.0.1"
    usage_message = ("\n--------------------------------------------------\n"
                     "The motif analysis program performs various motif-based analyses. "
                     "In order to use these tools, please type: \n\n"
                     "%prog [analysis type] [options]\n\n"
                     "Where [analysis type] refers to the type of the motif analysis performed "
                     "and [options] are the analysis-specific arguments.\n\n"
                     "Bellow you can find all current available analysis types. "
                     "To check the analyses specific options, please use:\n\n"
                     "%prog [analysis type] -h\n\n"
                     "For more information, please refer to our wiki:\n\n"
                     "https://code.google.com/p/reg-gen/wiki/RegGen\n\n"
                     "--------------------------------------------------\n\n"
                     "Options:\n"
                     "--version     show program's version number and exit.\n"
                     "-h, --help    show this help message and exit.\n"
                     "--matching    Performs motif matching analysis.\n"
                     "--enrichment  Performs motif enrichment analysis.\n")
    version_message = "Motif Analysis - Regulatory Analysis Toolbox (RGT). Version: "+str(current_version)

    # Processing Help/Version Options
    if(len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print(usage_message)
        sys.exit(0)
    elif(sys.argv[1] == "--version"):
        print(version_message)
        sys.exit(0)

    # Initializing Error Handler
    main_error_handler = ErrorHandler()

    ###################################################################################################
    # Redirecting to Specific Functions
    ###################################################################################################

    # Redirecting Loop
    if(sys.argv[1] == "--matching"):
        main_matching()
    elif(sys.argv[1] == "--enrichment"):
        main_enrichment()
    else:
        main_error_handler.throw_error("MOTIF_ANALYSIS_OPTION_ERROR")

def main_matching():
    """
    Performs motif matching.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    usage_message = "%prog --matching [options] <experiment_matrix>"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage = usage_message)

    # Parameters Options
    parser.add_option("--organism", dest = "organism", type = "string", metavar="STRING", default = "hg19",
                      help = ("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file."))
    parser.add_option("--fpr", dest = "fpr", type = "float", metavar="FLOAT", default = 0.0001,
                      help = ("False positive rate cutoff for motif matching."))
    parser.add_option("--precision", dest = "precision", type = "int", metavar="INT", default = 10000,
                      help = ("Score distribution precision for motif matching."))
    parser.add_option("--pseudocounts", dest = "pseudocounts", type = "float", metavar="FLOAT", default = 0.1,
                      help = ("Pseudocounts to be added to raw counts of each PWM."))
    parser.add_option("--prom-len", dest = "prom_len", type = "int", metavar="INT", default = 1000,
                      help = ("Length of the promoter region (in bp) considered on the creation of the regions-gene association."))
    parser.add_option("--max-len", dest = "max_len", type = "int", metavar="INT", default = 50000,
                      help = ("Maximum distance between a coordinate and a gene (in bp) in order for the former to "
                              "be considered associated with the latter."))
    parser.add_option("--rand-proportion", dest = "rand_proportion", type = "float", metavar="FLOAT", default = 10.0,
                      help = ("If random coordinates need to be created, then it will be created a number of "
                              "coordinates that equals this parameter x the number of input regions. If zero (0) "
                              "is passed, then no random coordinates are created."))

    # Output Options
    parser.add_option("--output-location", dest = "output_location", type = "string", metavar="PATH", default = os.getcwd(),
                      help = ("Path where the output files will be written."))
    parser.add_option("--bigbed", dest = "bigbed", action = "store_true", default = False,
                      help = ("If this option is used, all bed files will be written as bigbed."))

    # Processing Options
    options, arguments = parser.parse_args()

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading arguments
    try:
        input_matrix = arguments[0]
        if(len(arguments) > 1): pass # TODO WARNING
    except Exception: pass # TODO ERROR

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception: pass # TODO ERROR

    ###################################################################################################
    # Reading Regions & Gene Lists
    ###################################################################################################

    # Auxiliary object to store input
    class InputRegions:
        def __init__(self):
            self.gene_set = None
            self.region_list = []
            self.label_list = []

    # Initialization
    max_region_len = 0
    max_region = None
    input_regions_list = []
    try:
        exp_matrix_genelist_dict = exp_matrix.fieldsDict["genelist"]
    except Exception: pass # TODO ERROR
    try:
        exp_matrix_objects_dict = exp_matrix.objectsDict
    except Exception: pass # TODO ERROR

    # Iterating on experimental matrix dictionary
    for g in exp_matrix_genelist_dict.keys():

        # Initialization
        new_input_regions = InputRegions()
        flagGeneSet = False

        # Iterating on experimental matrix objects
        for k in exp_matrix_genelist_dict[g]:

            # If the object is a GenomicRegionSet
            if(isinstance(exp_matrix_objects_dict[k],GenomicRegionSet)):

                # Append label and GenomicRegionSet
                new_input_regions.label_list.append(k)
                new_input_regions.region_list.append(exp_matrix_objects_dict[k])
                curr_len = len(exp_matrix_objects_dict[k])
                if(curr_len > max_region_len):
                    max_region_len = curr_len
                    max_region = exp_matrix_objects_dict[k]

            # If the object is a GeneSet
            elif(isinstance(exp_matrix_objects_dict[k],GeneSet)):

                # Initialize GeneSet and gives warning if there are more than one genesets
                new_input_regions.gene_set = exp_matrix_objects_dict[k]
                if(not flagGeneSet): flagGeneSet = True
                else: pass # TODO WARNING

            # Warning if any additional file type
            else: pass # TODO WARNING

        # Append new_input_regions in input_regions_list
        input_regions_list.append(new_input_regions)

    ###################################################################################################
    # Creating random region
    ###################################################################################################

    # Create random coordinates
    rand_region = None
    if(options.rand_proportion > 0):
        rand_region = max_region.random_regions(options.organism, multiply_factor=float(options.rand_proportion), chrom_X=True)

    ###################################################################################################
    # Creating PWMs
    ###################################################################################################

    # TODO XXX Alterar setup.py para colocar os repositorios corretos de volta.

    # Initialization
    motif_data = MotifData()
    motif_list = []
    
    # Iterating on motif repositories
    for motif_repository in motif_data.get_pwm_list():

        # Iterating on motifs
        for motif_file_name in glob(os.path.join(motif_repository,"*.pwm")):
            motif = Motif(motif_file_name, float(options.pseudocounts), int(options.precision), float(options.fpr))
            motif_list.append(motif)

    ###################################################################################################
    # Motif Matching
    ###################################################################################################



    ###################################################################################################
    # Writing output
    ###################################################################################################



def main_enrichment():
    """
    Performs motif enrichment.

    Authors: Eduardo G. Gusmao, Manuel Allhoff, Joseph Kuo, Ivan G. Costa.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    usage_message = "%prog --enrichment --input-file <FILE> [options]"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage = usage_message)

    # Required Input Options
    parser.add_option("--input-file", dest = "input_file", type = "string", metavar="FILE", default = None,
                      help = ("Experimental matrix input file. This program will process only 'regions' "
                              "and 'genes' and must contain at least one 'regions' file. In case of multiple "
                              "'regions', the enrichment will be performed in each file separately. In case "
                              "of multiple 'regions' and 'genes', the enrichment will be performed in every "
                              "pairwise combination of these file types."))

    # Optional Input Options
    parser.add_option("--random-regions", dest = "random_regions", type = "string", metavar="FILE", default = None,
                      help = ("File containing the random regions for fisher's test. If None, create random "
                              "regions. The number of regions equals the size of the input regions file x "
                              "--random_proportion. A unique set of random regions will be created even if "
                              "multiple files are given in the experimental matrix. In this case, the input "
                              "file with the greatest number of regions will be used. The length of each genomic "
                              "region will follow the lengths of such input's regions."))

    # Parameters Options
    parser.add_option("--organism", dest = "organism", type = "string", metavar="STRING", default = "hg19",
                      help = ("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file."))
    parser.add_option("--motif-match-fpr", dest = "motif_match_fpr", type = "float", metavar="FLOAT", default = 0.0001,
                      help = ("False positive rate cutoff for motif matching."))
    parser.add_option("--motif-match-precision", dest = "motif_match_precision", type = "int", metavar="INT", default = 10000,
                      help = ("Score distribution precision for motif matching."))
    parser.add_option("--motif-match-pseudocounts", dest = "motif_match_pseudocounts", type = "float", metavar="FLOAT", default = 0.0,
                      help = ("Pseudocounts to be added to raw counts of each PWM."))
    parser.add_option("--multiple-test-alpha", dest = "multiple_test_alpha", type = "float", metavar="FLOAT", default = 0.05,
                      help = ("Alpha value for multiple test."))
    parser.add_option("--promoter-length", dest = "promoter_length", type = "int", metavar="INT", default = 1000,
                      help = ("Length of the promoter region (in bp) considered on the creation of the regions-gene association."))
    parser.add_option("--maximum-association-length", dest = "maximum_association_length", type = "int", metavar="INT", default = 50000,
                      help = ("Maximum distance between a coordinate and a gene (in bp) in order for the former to "
                              "be considered associated with the latter."))
    parser.add_option("--rand-proportion", dest = "rand_proportion", type = "float", metavar="FLOAT", default = 1.0,
                      help = ("If random coordinates need to be created, then it will be created a number of "
                              "coordinates that equals this parameter x the number of input regions. In case of "
                              "multiple input regions, the one with the greater number of regions will be "
                              "considered during this calculation."))
    parser.add_option("--all-regions", dest = "all_regions", action = "store_true", default = False,
                      help = ("If this option is used then all input regions will be considered as the evidence "
                              "set, i.e. it will be performed only the coordinate vs. background analysis."))

    # Output Options
    parser.add_option("--output-location", dest = "output_location", type = "string", metavar="PATH", default = os.getcwd(),
                      help = ("Path where the output files will be written."))
    parser.add_option("--print-association", dest = "print_association", type = "string", metavar="<Y|N>", default = "Y",
                      help = ("Whether to output a bigbed file containing the coordinate-gene association + MPBSs "
                              "that occured inside all coordinates."))
    parser.add_option("--print-mpbs", dest = "print_mpbs", type = "string", metavar="<Y|N>", default = "Y",
                      help = ("Whether to output a bigbed file containing all MPBSs found on input and random coordinates."))
    parser.add_option("--print-results-text", dest = "print_results_text", type = "string", metavar="<Y|N>", default = "Y",
                      help = ("Whether to output the fisher test results in text format."))
    parser.add_option("--print-results-html", dest = "print_results_html", type = "string", metavar="<Y|N>", default = "Y",
                      help = ("Whether to output the fisher test results in html format."))
    parser.add_option("--print-rand-coordinates", dest = "print_rand_coordinates", type = "string", metavar="<Y|N>", default = "Y",
                      help = ("Whether to output a bigbed file containing the random coordinates."))
    parser.add_option("--print-graph-mmscore", dest = "print_graph_mmscore", type = "string", metavar="<Y|N>", default = "N",
                      help = ("Whether to output graphs containing the motif matching score distribution for the "
                              "MPBSs found on the input and random coordinates."))
    parser.add_option("--print-graph-heatmap", dest = "print_graph_heatmap", type = "string", metavar="<Y|N>", default = "N",
                      help = ("Whether to output graphs containing heatmaps created based on the corrected p-values "
                              "of the multiple testing."))

    # Processing Options
    options, arguments = parser.parse_args()

    ###################################################################################################
    # Gene-Coordinate Association
    ###################################################################################################

    ###################################################################################################
    # Random Coordinates
    ###################################################################################################

    ###################################################################################################
    # Motif Matching
    ###################################################################################################

    ###################################################################################################
    # Statistics
    ###################################################################################################

    ###################################################################################################
    # Graphs
    ###################################################################################################




