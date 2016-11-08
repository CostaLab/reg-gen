###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import os
import sys
import warnings
from glob import glob
import time
from random import seed

warnings.filterwarnings("ignore")

# Internal
from rgt import __version__
from rgt.Util import PassThroughOptionParser, ErrorHandler, MotifData, GenomeData, ImageData, Html
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.GeneSet import GeneSet
from rgt.GenomicRegionSet import GenomicRegionSet
from Motif import Motif, Thresholds
from Match import match_single
from Statistics import multiple_test_correction, get_fisher_dict
from Util import Input, Result

# External
from pysam import Fastafile
from fisher import pvalue

"""
Contains functions to common motif analyses.

Dependencies:
- python >= 2.7
- numpy >= 1.4.0
- scipy >= 0.7.0
- biopython >= 1.64
- pysam >= 0.7.5
- fisher >= 0.1.4
- MOODS >= 1.0.1
- bedToBigBed and bigbedToBed scripts in $PATH (if the option is used)
- bedTools (deprecate this option)

Authors: Eduardo G. Gusmao.
"""


def main():
    start = time.time()
    """
    Main function that redirects tool usage.

    Keyword arguments: None

    Return: None
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    seed(42)
    usage_message = ("\n--------------------------------------------------\n"
                     "The motif analysis program performs various motif-based analyses. "
                     "In order to use these tools, please type: \n\n"
                     "%prog [analysis type] [options]\n\n"
                     "Where [analysis type] refers to the type of the motif analysis performed "
                     "and [options] are the analysis-specific arguments.\n\n"
                     "Below you can find all current available analysis types. "
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
    version_message = "Motif Analysis - Regulatory Analysis Toolbox (RGT). Version: " + str(__version__)

    # Processing Help/Version Options
    if len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print(usage_message)
        sys.exit(0)
    elif sys.argv[1] == "--version":
        print(version_message)
        sys.exit(0)

    # Initializing Error Handler
    main_error_handler = ErrorHandler()

    ###################################################################################################
    # Redirecting to Specific Functions
    ###################################################################################################

    # Redirecting Loop
    if sys.argv[1] == "--matching":
        main_matching()
    elif sys.argv[1] == "--enrichment":
        main_enrichment()
    else:
        main_error_handler.throw_error("MOTIF_ANALYSIS_OPTION_ERROR")

    print("Completed in", time.time() - start, "seconds")

def main_matching():
    """
    Performs motif matching.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    main_error_handler = ErrorHandler()

    # Parameters
    usage_message = "%prog --matching [options] <experiment_matrix>"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Parameters Options
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help="Organism considered on the analysis. Check our full documentation for all available "
                           "options. All default files such as genomes will be based on the chosen organism "
                           "and the data.config file.")
    parser.add_option("--fpr", dest="fpr", type="float", metavar="FLOAT", default=0.0001,
                      help="False positive rate cutoff for motif matching.")
    parser.add_option("--precision", dest="precision", type="int", metavar="INT", default=10000,
                      help="Score distribution precision for determining false positive rate cutoff.")
    parser.add_option("--pseudocounts", dest="pseudocounts", type="float", metavar="FLOAT", default=0.1,
                      help="Pseudocounts to be added to raw counts of each PFM.")
    parser.add_option("--rand-proportion", dest="rand_proportion", type="float", metavar="FLOAT", default=10.0,
                      help="If random coordinates need to be created (for further motif enrichment),"
                           "then it will be created a number of coordinates that equals this"
                           "parameter x the number of input regions (in case of multiple regions, the"
                           "larger is considered). If zero (0) is passed, then no random coordinates are created.")
    parser.add_option("--norm-threshold", dest="norm_threshold", action="store_true", default=False,
                      help="If this option is used, the thresholds for all PWMs will be normalized by their length."
                           "In this scheme, the threshold cutoff is evaluated in the regular way by the given fpr."
                           "Then, all thresholds are divided by the length of the motif. The final threshold"
                           "consists of the average between all normalized motif thresholds. This single threshold"
                           "will be applied to all motifs.")
    parser.add_option("--background-bed", dest="background_bed", type="string", metavar="STRING",
                      help="Bed file containing the genomic regions to use as background."
                           "If set, --rand-proportion has no effect.")
    parser.add_option("--selected-motifs", dest="selected_motifs_filename", type="string", metavar="STRING",
                      help="Only use the motifs contained within this file (one for each line).")
    # TODO: allow a more complex file type, specifying both motif name and database maybe?

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH",
                      default=os.getcwd(), help="Path where the output files will be written.")
    parser.add_option("--bigbed", dest="bigbed", action="store_true", default=False,
                      help="If this option is used, all bed files will be written as bigbed.")
    parser.add_option("--normalize-bitscore", dest="normalize_bitscore", action="store_true", default=False,
                      help="In order to print bigbed files the scores need to be normalized between 0 and 1000."
                           "Don't use this option if real bitscores should be printed in the resulting bed file."
                           "Without this option, bigbed files will never be created.")

    # Processing Options
    options, arguments = parser.parse_args()

    # Additional Parameters
    matching_folder_name = "Match"
    random_region_name = "random_regions"
    background_regions_name = "background_regions"

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    matching_output_location = os.path.join(options.output_location, matching_folder_name)
    try:
        if not os.path.isdir(matching_output_location):
            os.makedirs(matching_output_location)
    except Exception:
        main_error_handler.throw_error("MM_OUT_FOLDER_CREATION")

    # Default genomic data
    genome_data = GenomeData(options.organism)

    # Default motif data
    motif_data = MotifData()

    # Reading motif file
    selected_motifs = []

    if options.selected_motifs_filename:
        try:
            with open(options.selected_motifs_filename) as f:
                selected_motifs = f.read().splitlines()
                selected_motifs = filter(None, selected_motifs)
        except Exception:
            main_error_handler.throw_error("MM_MOTIFS_NOTFOUND", add_msg=options.selected_motifs_filename)

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading arguments
    try:
        input_matrix = arguments[0]
        if len(arguments) > 1:
            main_error_handler.throw_warning("MM_MANY_ARG")
    except Exception:
        main_error_handler.throw_error("MM_NO_ARGUMENT")

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception:
        main_error_handler.throw_error("MM_WRONG_EXPMAT")

    ###################################################################################################
    # Reading Regions
    ###################################################################################################

    # Initialization
    max_region_len = 0
    max_region = None
    regions_to_match = []

    try:
        exp_matrix_objects_dict = exp_matrix.objectsDict
    except Exception:
        main_error_handler.throw_error("MM_WRONG_EXPMAT")

    # Iterating on experimental matrix objects
    for k in exp_matrix_objects_dict.keys():

        curr_genomic_region = exp_matrix_objects_dict[k]

        # If the object is a GenomicRegionSet
        if isinstance(curr_genomic_region, GenomicRegionSet):

            # Sorting input region
            curr_genomic_region.sort()

            # Append label and GenomicRegionSet
            regions_to_match.append(curr_genomic_region)

            # Verifying max_region_len for random region generation
            curr_len = len(curr_genomic_region)
            if curr_len > max_region_len:
                max_region_len = curr_len
                max_region = exp_matrix_objects_dict[k]

    ###################################################################################################
    # Creating background regions
    ###################################################################################################

    # if a background file is provided, use it
    if options.background_bed:

        # TODO: should also accept big bed

        bg_regions = GenomicRegionSet(background_regions_name)

        try:
            bg_regions.read_bed(os.path.abspath(options.background_bed))
        except Exception:
            # FIXME: add another error type, requires some fixing around
            main_error_handler.throw_error("DEFAULT_ERROR")

        regions_to_match.append(bg_regions)

    # otherwise if a random proportion is set, create random regions
    elif options.rand_proportion > 0:

        # Create random coordinates and name it random_regions
        rand_region = max_region.random_regions(options.organism, multiply_factor=options.rand_proportion, chrom_X=True)
        rand_region.sort()
        rand_region.name = random_region_name

        # Add random regions to the list of regions to perform matching on
        regions_to_match.append(rand_region)

        # Writing random regions
        output_file_name = os.path.join(matching_output_location, random_region_name)
        rand_bed_file_name = output_file_name + ".bed"
        rand_region.write_bed(rand_bed_file_name)

        # Verifying condition to write bb
        if options.bigbed:

            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            rand_bb_file_name = output_file_name + ".bb"
            try:
                os.system(
                    " ".join(["bedToBigBed", rand_bed_file_name, chrom_sizes_file, rand_bb_file_name, "-verbose=0"]))
                os.remove(rand_bed_file_name)
            except Exception:
                main_error_handler.throw_warning("DEFAULT_WARNING")  # FIXME: maybe error instead?
    else:
        # FIXME: not very neat
        main_error_handler.throw_warning("DEFAULT_WARNING", add_msg="No background was selected")

    # if neither case was selected, we just do the motif matching over the input regions, without background

    ###################################################################################################
    # Creating PWMs
    ###################################################################################################

    # Initialization
    motif_list = []

    # Creating thresholds object
    thresholds = Thresholds(motif_data)

    # Fetching list with all motif file names
    motif_file_names = []
    for motif_repository in motif_data.get_pwm_list():
        for motif_file_name in glob(os.path.join(motif_repository, "*.pwm")):
            motif_name = os.path.basename(os.path.splitext(motif_file_name)[0])
            # if the user has given a list of motifs to use, we only
            # add those to our list
            if not selected_motifs or motif_name in selected_motifs:
                motif_file_names.append(motif_file_name)

    # Iterating on grouped file name list
    for motif_file_name in motif_file_names:
        # Append motif motif_list
        motif_list.append(Motif(motif_file_name, options.pseudocounts, options.precision, options.fpr, thresholds))

    # Performing normalized threshold strategy if requested
    if options.norm_threshold:
        threshold_list = [motif.threshold / motif.len for motif in motif_list]
        unique_threshold = sum(threshold_list) / len(threshold_list)
    else:
        unique_threshold = None

    ###################################################################################################
    # Motif Matching
    ###################################################################################################

    # Creating genome file
    genome_file = Fastafile(genome_data.get_genome())

    # Iterating on list of genomic regions
    for genomic_region_set in regions_to_match:

        # Initializing output bed file
        output_file_name = os.path.join(matching_output_location, genomic_region_set.name + "_mpbs")
        bed_file_name = output_file_name + ".bed"
        output_file = open(bed_file_name, "w")

        # Iterating on genomic regions
        for genomic_region in genomic_region_set.sequences:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            # Splitting the sequence in smaller sequences to remove the "N" regions
            sequence_list = filter(None, sequence.split("N"))

            # Perform motif matching for each motif in each sequence
            for seq in sequence_list:
                for motif in motif_list:
                    match_single(motif, seq, genomic_region, output_file, unique_threshold, options.normalize_bitscore)

        # Closing file
        output_file.close()

        # Verifying condition to write bb
        if options.bigbed and options.normalize_bitscore:
            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            sort_file_name = output_file_name + "_sort.bed"
            bb_file_name = output_file_name + ".bb"
            os.system("sort -k1,1 -k2,2n " + bed_file_name + " > " + sort_file_name)
            os.system(" ".join(["bedToBigBed", sort_file_name, chrom_sizes_file, bb_file_name, "-verbose=0"]))
            os.remove(bed_file_name)
            os.remove(sort_file_name)


def main_enrichment():
    """
    Performs motif enrichment.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    main_error_handler = ErrorHandler()

    # Parameters
    usage_message = "%prog --enrichment [options] <experiment_matrix> <input_path>"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage=usage_message)

    # Parameters Options
    parser.add_option("--organism", dest="organism", type="string", metavar="STRING", default="hg19",
                      help="Organism considered on the analysis. Check our full documentation for all available "
                           "options. All default files such as genomes will be based on the chosen organism "
                           "and the data.config file.")
    parser.add_option("--promoter-length", dest="promoter_length", type="int", metavar="INT", default=1000,
                      help="Length of the promoter region (in bp) considered on the creation of the "
                           "regions-gene association.")
    parser.add_option("--maximum-association-length", dest="maximum_association_length", type="int", metavar="INT",
                      default=50000,
                      help="Maximum distance between a coordinate and a gene (in bp) in order for the former to "
                           "be considered associated with the latter.")
    parser.add_option("--multiple-test-alpha", dest="multiple_test_alpha", type="float", metavar="FLOAT", default=0.05,
                      help="Alpha value for multiple test.")
    parser.add_option("--processes", dest="processes", type="int", metavar="INT", default=1,
                      help="Number of processes for multi-CPU based machines.")
    parser.add_option("--background-bed", dest="background_bed", type="string", metavar="STRING",
                      help="Bed file containing the genomic regions to use as background."
                           "If set, the random regions files will be ignored.")
    parser.add_option("--selected-motifs", dest="selected_motifs_filename", type="string", metavar="STRING",
                      help="Only use the motifs contained within this file (one for each line).")
    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH", default=None,
                      help="Path where the output files will be written. Default is the input PATH.")
    parser.add_option("--print-thresh", dest="print_thresh", type="float", metavar="FLOAT", default=0.05,
                      help="Only MPBSs whose factor's enrichment corrected p-value are less than equal "
                           "this option are print. Use 1.0 to print all MPBSs.")
    parser.add_option("--bigbed", dest="bigbed", action="store_true", default=False,
                      help="If this option is used, all bed files will be written as bigbed.")

    # Processing Options
    options, arguments = parser.parse_args()

    # Additional Parameters
    matching_folder_name = "Match"
    random_region_name = "random_regions"
    background_region_name = "background_regions"
    gene_column_name = "genegroup"
    output_association_name = "coord_association"
    # output_mpbs_filtered = "mpbs"
    output_mpbs_filtered_ev = "mpbs_ev"
    output_mpbs_filtered_nev = "mpbs_nev"
    output_stat_genetest = "genetest_statistics"
    output_stat_fulltest = "fulltest_statistics"
    ev_color = "0,130,0"
    nev_color = "130,0,0"
    results_header_text = "\t".join(
        ["FACTOR", "P-VALUE", "CORR.P-VALUE", "A", "B", "C", "D", "FREQ", "BACK.FREQ.", "GENES"])
    html_header = ["FACTOR", "MOTIF", "P-VALUE", "CORRECTED P-VALUE", "A", "B", "C", "D", "FREQUENCY",
                   "BACKGROUND FREQUENCY", "GO"]
    html_type_list = "sissssssssl"
    logo_width = 200
    if "hg" in options.organism:
        gprofiler_link = "http://biit.cs.ut.ee/gprofiler/index.cgi?significant=1&sort_by_structure=1&ordered_query=0&organism=hsapiens&query="
    else:
        gprofiler_link = "http://biit.cs.ut.ee/gprofiler/index.cgi?significant=1&sort_by_structure=1&ordered_query=0&organism=mmusculus&query="
    html_col_size = [300, logo_width, 100, 100, 50, 50, 50, 50, 100, 100, 50]

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if options.output_location:
        output_location = options.output_location  # Output location was given
    else:
        output_location = arguments[1]  # Output location is the same as the match folder location
    try:
        matrix_name_without_ext = ".".join(arguments[0].split(".")[:-1])
        output_location_results = os.path.join(output_location, os.path.basename(matrix_name_without_ext))
        if not os.path.isdir(output_location_results):
            os.makedirs(output_location_results)
    except Exception:
        main_error_handler.throw_error("ME_OUT_FOLDER_CREATION")

    # Default genomic data
    genome_data = GenomeData(options.organism)

    # Default motif data
    motif_data = MotifData()

    # Reading motif file
    selected_motifs = []

    if options.selected_motifs_filename:
        try:
            with open(options.selected_motifs_filename) as f:
                selected_motifs = f.read().splitlines()
                selected_motifs = filter(None, selected_motifs)
        except Exception:
            main_error_handler.throw_error("MM_MOTIFS_NOTFOUND", add_msg=options.selected_motifs_filename)

    # Default image data
    image_data = ImageData()

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading arguments
    try:
        input_matrix = arguments[0]
        input_location = arguments[1]
        if len(arguments) > 2:
            main_error_handler.throw_warning("ME_MANY_ARG")
    except Exception:
        main_error_handler.throw_error("ME_FEW_ARG")

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception:
        main_error_handler.throw_error("ME_WRONG_EXPMAT")

    ###################################################################################################
    # Reading Regions & Gene Lists
    ###################################################################################################

    # Initializations
    input_list = []

    # Reading dictionary grouped by fields
    flag_gene = True
    try:
        exp_matrix_fields_dict = exp_matrix.fieldsDict[gene_column_name]
    except Exception:
        flag_gene = False

    # Reading dictionary of objects
    try:
        exp_matrix_objects_dict = exp_matrix.objectsDict
    except Exception:
        main_error_handler.throw_error("ME_WRONG_EXPMAT")

    if flag_gene:  # Genelist and full site analysis will be performed

        # Iterating on experimental matrix fields
        for g in exp_matrix_fields_dict.keys():

            # Create input which will contain all regions associated with such gene group
            curr_input = Input(None, [])

            # This flag will be used to see if there are two gene files associated with the same gene label on genegroup column
            flag_foundgeneset = False

            # Iterating on experimental matrix objects
            for k in exp_matrix_fields_dict[g]:

                curr_object = exp_matrix_objects_dict[k]

                # If the current object is a GenomicRegionSet
                if isinstance(curr_object, GenomicRegionSet):
                    # Sorting input region
                    curr_object.sort()

                    # Updating Input object
                    curr_input.region_list.append(curr_object)

                # If the current object is a GeneSet
                if isinstance(curr_object, GeneSet):

                    # Updating Input object
                    curr_object.name = g  # The name in gene_group column will be used. The 'name' column for genes are not used.
                    if not flag_foundgeneset:
                        curr_input.gene_set = curr_object
                        flag_foundgeneset = True
                    else:
                        main_error_handler.throw_warning("ME_MANY_GENESETS")

            if not flag_foundgeneset:
                main_error_handler.throw_warning("ME_FEW_GENESETS")

            # Update input list
            input_list.append(curr_input)

    else:  # Only full site analysis will be performed

        # Create single input which will contain all regions
        single_input = Input(None, [])

        # Iterating on experimental matrix objects
        for k in exp_matrix_objects_dict.keys():

            curr_object = exp_matrix_objects_dict[k]

            # If the current object is a GenomicRegionSet
            if isinstance(curr_object, GenomicRegionSet):
                # Sorting input region
                curr_object.sort()

                # Updating Input object
                single_input.region_list.append(curr_object)

        # Updating input list with single input (only full site analysis will be performed)
        input_list = [single_input]

    ###################################################################################################
    # Fetching Motif List
    ###################################################################################################

    # Fetching list with all motif names
    motif_names = []
    for motif_repository in motif_data.get_pwm_list():
        for motif_file_name in glob(os.path.join(motif_repository, "*.pwm")):
            motif_name = os.path.basename(os.path.splitext(motif_file_name)[0])
            # if the user has given a list of motifs to use, we only
            # add those to our list
            if not selected_motifs or motif_name in selected_motifs:
                motif_names.append(motif_name)
    motif_names = sorted(motif_names)

    # Grouping motif file names by the number of processes requested
    if options.processes <= 0:
        main_error_handler.throw_error("ME_LOW_NPROC")
    elif options.processes == 1:
        motif_names_grouped = [[e] for e in motif_names]
    else:
        motif_names_grouped = map(None, *(iter(motif_names),) * options.processes)
    motif_names_grouped = [[e2 for e2 in e1 if e2 is not None] for e1 in motif_names_grouped]

    ###################################################################################################
    # Background Statistics
    ###################################################################################################

    # FIXME: this is brittle. The files to use as background, both normal regions and mpbs, should be
    # passed as arguments, not globbed - in both the "given" and the "random" cases

    # if the background is not defined, we expect random regions instead
    if not options.background_bed:
        background_region_name = random_region_name

        # Verifying background region file exists
        background_region_glob = glob(os.path.join(input_location, matching_folder_name, background_region_name + ".*"))
        try:
            background_region_file_name = background_region_glob[0]
        except Exception:
            main_error_handler.throw_error("ME_RAND_NOTFOUND")
    else:
        background_region_file_name = os.path.abspath(options.background_bed)

    # Verifying background region MPBS file exists
    background_region_glob = glob(os.path.join(input_location, matching_folder_name, background_region_name + "_mpbs.*"))
    try:
        background_mpbs_file_name = background_region_glob[0]
    except Exception:
        pass  # XXX TODO main_error_handler.throw_error("ME_RANDMPBS_NOTFOUND")

    # Converting regions bigbed file
    background_region_bed_name = ".".join(background_region_file_name.split(".")[:-1]) + ".bed"
    if background_region_file_name.split(".")[-1] == "bb":
        background_region_bed_name = os.path.join(output_location_results, background_region_name + ".bed")
        os.system(" ".join(["bigBedToBed", background_region_file_name, background_region_bed_name]))
    elif background_region_file_name.split(".")[-1] != "bed":
        main_error_handler.throw_error("ME_RAND_NOT_BED_BB")

    # Converting mpbs bigbed file
    background_mpbs_bed_name = ".".join(background_mpbs_file_name.split(".")[:-1]) + ".bed"
    if background_mpbs_file_name.split(".")[-1] == "bb":
        background_mpbs_bed_name = os.path.join(output_location_results, background_region_name + "_mpbs.bed")
        os.system(" ".join(["bigBedToBed", background_mpbs_file_name, background_mpbs_bed_name]))
    elif background_mpbs_file_name.split(".")[-1] != "bed":
        pass  # XXX TODO main_error_handler.throw_error("ME_RAND_NOT_BED_BB")

    # Evaluating background statistics
    bg_c_dict, bg_d_dict = get_fisher_dict(motif_names_grouped, background_region_bed_name, background_mpbs_bed_name,
                                           output_location_results, return_geneset=False)

    # Removing bed files if bb exist
    if background_region_file_name.split(".")[-1] == "bb":
        os.remove(background_region_bed_name)
    if background_mpbs_file_name.split(".")[-1] == "bb":
        os.remove(background_mpbs_bed_name)

    ###################################################################################################
    # Enrichment Statistics
    ###################################################################################################

    # Creating link dictionary for HTML file
    genetest_link_dict = dict()
    sitetest_link_dict = dict()
    link_location = "file://" + os.path.abspath(output_location_results)
    for curr_input in input_list:
        for grs in curr_input.region_list:
            if curr_input.gene_set:
                link_name = grs.name + " (" + curr_input.gene_set.name + ")"
                genetest_link_dict[link_name] = os.path.join(link_location, grs.name + "__" + curr_input.gene_set.name,
                                                             output_stat_genetest + ".html")
                sitetest_link_dict[link_name] = os.path.join(link_location, grs.name + "__" + curr_input.gene_set.name,
                                                             output_stat_fulltest + ".html")
            else:
                link_name = grs.name
                sitetest_link_dict[link_name] = os.path.join(link_location, link_name, output_stat_fulltest + ".html")

    # Iterating on each input object
    for curr_input in input_list:

        # Iterating on each input genomic region set
        for grs in curr_input.region_list:

            # Initialization
            original_name = grs.name
            to_remove_list = []

            # Creating output folder
            if curr_input.gene_set:
                curr_output_folder_name = os.path.join(output_location_results,
                                                       grs.name + "__" + curr_input.gene_set.name)
            else:
                curr_output_folder_name = os.path.join(output_location_results, grs.name)
            if not os.path.isdir(curr_output_folder_name):
                os.makedirs(curr_output_folder_name)

            # Verifying if MPBS file exists
            curr_mpbs_glob = glob(os.path.join(input_location, matching_folder_name, original_name + "_mpbs.*"))
            try:
                curr_mpbs_file_name = curr_mpbs_glob[0]
            except Exception:
                pass  # TODO main_error_handler.throw_error("ME_RAND_NOTFOUND")

            # Converting ev mpbs bigbed file
            curr_mpbs_bed_name = ".".join(curr_mpbs_file_name.split(".")[:-1]) + ".bed"
            if curr_mpbs_file_name.split(".")[-1] == "bb":
                curr_mpbs_bed_name = os.path.join(curr_output_folder_name, original_name + "_mpbs.bed")
                os.system(" ".join(["bigBedToBed", curr_mpbs_file_name, curr_mpbs_bed_name]))
                to_remove_list.append(curr_mpbs_bed_name)
            elif curr_mpbs_file_name.split(".")[-1] != "bed":
                pass  # XXX TODO main_error_handler.throw_error("ME_RAND_NOT_BED_BB")

            ###################################################################################################
            # Gene Evidence Statistics
            ###################################################################################################

            if curr_input.gene_set:

                # Performing association of input region with gene_set
                grs = grs.gene_association(curr_input.gene_set, options.organism, options.promoter_length,
                                           options.maximum_association_length)

                # Writing gene-coordinate association
                output_file_name = os.path.join(curr_output_folder_name, output_association_name + ".bed")
                output_file = open(output_file_name, "w")
                for gr in grs:
                    if gr.name == ".":
                        curr_name = "."
                    else:
                        curr_gene_list = [e if e[0] != "." else e[1:] for e in gr.name.split(":")]
                        curr_prox_list = gr.proximity.split(":")
                        curr_name = ":".join([e[0] + "_" + e[1] for e in zip(curr_gene_list, curr_prox_list)])
                    output_file.write("\t".join([str(e) for e in [gr.chrom, gr.initial, gr.final, curr_name]]) + "\n")
                output_file.close()
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    bb_file_name = output_file_name + ".bb"
                    try:
                        os.system(
                            " ".join(["bedToBigBed", output_file_name, chrom_sizes_file, bb_file_name, "-verbose=0"]))
                        os.remove(output_file_name)
                    except Exception:
                        pass  # WARNING

                # Writing ev and nev regions to temporary bed files in order to evaluate statistics
                ev_regions_file_name = os.path.join(curr_output_folder_name, "ev_regions.bed")
                nev_regions_file_name = os.path.join(curr_output_folder_name, "nev_regions.bed")
                output_file_ev = open(ev_regions_file_name, "w")
                output_file_nev = open(nev_regions_file_name, "w")
                for gr in grs:
                    if len([e for e in gr.name.split(":") if e[0] != "."]) > 0:
                        output_file_ev.write("\t".join([str(e) for e in
                                                        [gr.chrom, gr.initial, gr.final, gr.name, gr.data,
                                                         gr.orientation]]) + "\n")
                    else:
                        output_file_nev.write("\t".join([str(e) for e in
                                                         [gr.chrom, gr.initial, gr.final, gr.name, gr.data,
                                                          gr.orientation]]) + "\n")
                output_file_ev.close()
                output_file_nev.close()
                to_remove_list.append(ev_regions_file_name)
                to_remove_list.append(nev_regions_file_name)

                # Calculating statistics
                ev_mpbs_file_name_temp = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + "_temp.bed")
                nev_mpbs_file_name_temp = os.path.join(curr_output_folder_name, output_mpbs_filtered_nev + "_temp.bed")
                ev_mpbs_file = open(ev_mpbs_file_name_temp, "w")
                nev_mpbs_file = open(nev_mpbs_file_name_temp, "w")
                curr_a_dict, curr_b_dict, ev_genelist_dict = get_fisher_dict(motif_names_grouped, ev_regions_file_name,
                                                                             curr_mpbs_bed_name,
                                                                             curr_output_folder_name,
                                                                             return_geneset=True,
                                                                             output_mpbs_file=ev_mpbs_file,
                                                                             color=ev_color)
                curr_c_dict, curr_d_dict, nev_genelist_dict = get_fisher_dict(motif_names_grouped,
                                                                              nev_regions_file_name, curr_mpbs_bed_name,
                                                                              curr_output_folder_name,
                                                                              return_geneset=True,
                                                                              output_mpbs_file=nev_mpbs_file,
                                                                              color=nev_color)
                ev_mpbs_file.close()
                nev_mpbs_file.close()
                to_remove_list.append(ev_mpbs_file_name_temp)
                to_remove_list.append(nev_mpbs_file_name_temp)

                # Performing fisher test
                result_list = []
                for k in motif_names:
                    r = Result()
                    r.name = k
                    r.a = curr_a_dict[k]
                    r.b = curr_b_dict[k]
                    r.c = curr_c_dict[k]
                    r.d = curr_d_dict[k]
                    r.percent = float(r.a) / float(r.a + r.b)
                    r.back_percent = float(r.c) / float(r.c + r.d)
                    r.genes = ev_genelist_dict[k]
                    try:
                        p = pvalue(r.a, r.b, r.c, r.d)
                        r.p_value = p.right_tail
                    except Exception:
                        r.p_value = 1.0
                    result_list.append(r)

                # Performing multiple test correction
                multuple_corr_rej, multiple_corr_list = multiple_test_correction([e.p_value for e in result_list],
                                                                                 alpha=options.multiple_test_alpha,
                                                                                 method='indep')
                corr_pvalue_dict = dict()  # Needed to filter the mpbs in a fast way
                for i in range(0, len(multiple_corr_list)):
                    result_list[i].corr_p_value = multiple_corr_list[i]
                    corr_pvalue_dict[result_list[i].name] = result_list[i].corr_p_value

                # Sorting result list
                result_list = sorted(result_list, key=lambda x: x.name)
                result_list = sorted(result_list, key=lambda x: x.percent, reverse=True)
                result_list = sorted(result_list, key=lambda x: x.p_value)
                result_list = sorted(result_list, key=lambda x: x.corr_p_value)

                # Preparing results for printing
                for r in result_list:
                    r.p_value = "%.4e" % r.p_value
                    r.corr_p_value = "%.4e" % r.corr_p_value
                    r.percent = str(round(r.percent, 4) * 100) + "%"
                    r.back_percent = str(round(r.back_percent, 4) * 100) + "%"

                # Printing ev and nev mpbs
                ev_mpbs_file = open(ev_mpbs_file_name_temp, "r")
                nev_mpbs_file = open(nev_mpbs_file_name_temp, "r")
                ev_mpbs_file_name_thresh = os.path.join(curr_output_folder_name,
                                                        output_mpbs_filtered_ev + "_thresh.bed")
                nev_mpbs_file_name_thresh = os.path.join(curr_output_folder_name,
                                                         output_mpbs_filtered_nev + "_thresh.bed")
                output_file_ev = open(ev_mpbs_file_name_thresh, "w")
                output_file_nev = open(nev_mpbs_file_name_thresh, "w")
                for line in ev_mpbs_file:
                    ll = line.strip().split("\t")
                    if corr_pvalue_dict[ll[3]] > options.print_thresh:
                        continue
                    output_file_ev.write(line)
                for line in nev_mpbs_file:
                    ll = line.strip().split("\t")
                    if corr_pvalue_dict[ll[3]] > options.print_thresh:
                        continue
                    output_file_nev.write(line)
                output_file_ev.close()
                output_file_nev.close()
                to_remove_list.append(ev_mpbs_file_name_thresh)
                to_remove_list.append(nev_mpbs_file_name_thresh)

                # Sorting ev and nev mpbs
                output_file_name_ev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bed")
                output_file_name_nev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_nev + ".bed")
                os.system(
                    "sort -k1,1 -k2,2n " + ev_mpbs_file_name_thresh + " > " + output_file_name_ev_bed)  # Sorting ev file
                os.system(
                    "sort -k1,1 -k2,2n " + nev_mpbs_file_name_thresh + " > " + output_file_name_nev_bed)  # Sorting nev file

                # Converting ev and nev mpbs to bigbed
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    output_file_name_ev_bb = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bb")
                    output_file_name_nev_bb = os.path.join(curr_output_folder_name, output_mpbs_filtered_nev + ".bb")
                    try:
                        os.system(" ".join(
                            ["bedToBigBed", output_file_name_ev_bed, chrom_sizes_file, output_file_name_ev_bb,
                             "-verbose=0"]))
                        os.system(" ".join(
                            ["bedToBigBed", output_file_name_nev_bed, chrom_sizes_file, output_file_name_nev_bb,
                             "-verbose=0"]))
                        to_remove_list.append(output_file_name_ev_bed)
                        to_remove_list.append(output_file_name_nev_bed)
                    except Exception:
                        pass  # WARNING

                # Printing statistics text
                output_file_name_stat_text = os.path.join(curr_output_folder_name, output_stat_genetest + ".txt")
                output_file = open(output_file_name_stat_text, "w")
                output_file.write(results_header_text + "\n")
                for r in result_list:
                    output_file.write(str(r) + "\n")
                output_file.close()

                # Printing statistics html - Creating data table
                data_table = []
                for r in result_list:
                    curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                    for rep in motif_data.get_logo_list():
                        logo_file_name = os.path.join(rep, r.name + ".png")
                        if os.path.isfile(logo_file_name):
                            curr_motif_tuple = [logo_file_name, logo_width]
                            break
                    curr_gene_tuple = ["View", gprofiler_link + ",".join(r.genes.genes)]
                    data_table.append(
                        [r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                         str(r.c), str(r.d), str(r.percent), str(r.back_percent), curr_gene_tuple])

                # Printing statistics html - Writing to HTML
                output_file_name_html = os.path.join(curr_output_folder_name, output_stat_genetest + ".html")
                html = Html("Motif Enrichment Analysis", genetest_link_dict)
                html.add_heading(
                    "Results for <b>" + original_name + "</b> region <b>Gene Test*</b> using genes from <b>" + curr_input.gene_set.name + "</b>",
                    align="center", bold=False)
                html.add_heading(
                    "* This gene test considered regions associated to the gene list given against regions not associated to the gene list",
                    align="center", bold=False, size=3)
                html.add_zebra_table(html_header, html_col_size, html_type_list, data_table, align="center")
                html.write(output_file_name_html)

            else:

                # Association still needs to be done with all genes in order to print gene list
                grs = grs.gene_association(None, options.organism, options.promoter_length,
                                           options.maximum_association_length)

                # If there is no gene list, then the current evidence set consists of all coordinates
                ev_regions_file_name = os.path.join(curr_output_folder_name, "ev_regions.bed")
                output_file_ev = open(ev_regions_file_name, "w")
                for gr in grs:
                    output_file_ev.write("\t".join([str(e) for e in [gr.chrom, gr.initial, gr.final,
                                                                     gr.name, gr.data, gr.orientation]]) + "\n")
                output_file_ev.close()
                to_remove_list.append(ev_regions_file_name)

                # Calculating statistics
                ev_mpbs_file_name_temp = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + "_temp.bed")
                ev_mpbs_file = open(ev_mpbs_file_name_temp, "w")
                curr_a_dict, curr_b_dict, ev_genelist_dict = get_fisher_dict(motif_names_grouped, ev_regions_file_name,
                                                                             curr_mpbs_bed_name,
                                                                             curr_output_folder_name,
                                                                             return_geneset=True,
                                                                             output_mpbs_file=ev_mpbs_file,
                                                                             color=ev_color)
                ev_mpbs_file.close()
                to_remove_list.append(ev_mpbs_file_name_temp)

            ###################################################################################################
            # Final wrap-up
            ###################################################################################################

            # Performing fisher test
            result_list = []
            for k in motif_names:
                r = Result()
                r.name = k
                r.a = curr_a_dict[k]
                r.b = curr_b_dict[k]
                r.c = bg_c_dict[k]
                r.d = bg_d_dict[k]
                r.percent = float(r.a) / float(r.a + r.b)
                r.back_percent = float(r.c) / float(r.c + r.d)
                r.genes = ev_genelist_dict[k]
                try:
                    p = pvalue(r.a, r.b, r.c, r.d)
                    r.p_value = p.right_tail
                except Exception:
                    r.p_value = 1.0
                result_list.append(r)

            # Performing multiple test correction
            multuple_corr_rej, multiple_corr_list = multiple_test_correction([e.p_value for e in result_list],
                                                                             alpha=options.multiple_test_alpha,
                                                                             method='indep')
            corr_pvalue_dict = dict()  # Needed to filter the mpbs in a fast way
            for i in range(0, len(multiple_corr_list)):
                result_list[i].corr_p_value = multiple_corr_list[i]
                corr_pvalue_dict[result_list[i].name] = result_list[i].corr_p_value

            # Sorting result list
            result_list = sorted(result_list, key=lambda x: x.name)
            result_list = sorted(result_list, key=lambda x: x.percent, reverse=True)
            result_list = sorted(result_list, key=lambda x: x.p_value)
            result_list = sorted(result_list, key=lambda x: x.corr_p_value)

            # Preparing results for printing
            for r in result_list:
                r.p_value = "%.4e" % r.p_value
                r.corr_p_value = "%.4e" % r.corr_p_value
                r.percent = str(round(r.percent, 4) * 100) + "%"
                r.back_percent = str(round(r.back_percent, 4) * 100) + "%"

            # Printing ev if it was not already print in geneset
            if not curr_input.gene_set:

                # Printing ev and nev mpbs
                ev_mpbs_file = open(ev_mpbs_file_name_temp, "r")
                ev_mpbs_file_name_thresh = os.path.join(curr_output_folder_name,
                                                        output_mpbs_filtered_ev + "_thresh.bed")
                output_file_ev = open(ev_mpbs_file_name_thresh, "w")
                for line in ev_mpbs_file:
                    ll = line.strip().split("\t")
                    if corr_pvalue_dict[ll[3]] > options.print_thresh:
                        continue
                    output_file_ev.write(line)
                output_file_ev.close()
                to_remove_list.append(ev_mpbs_file_name_thresh)

                # Sorting ev mpbs
                output_file_name_ev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bed")
                os.system(
                    "sort -k1,1 -k2,2n " + ev_mpbs_file_name_thresh + " > " + output_file_name_ev_bed)  # Sorting ev file

                # Converting ev and nev mpbs to bigbed
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    output_file_name_ev_bb = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bb")
                    try:
                        os.system(" ".join(
                            ["bedToBigBed", output_file_name_ev_bed, chrom_sizes_file, output_file_name_ev_bb,
                             "-verbose=0"]))
                        to_remove_list.append(output_file_name_ev_bed)
                    except Exception:
                        pass  # WARNING

            # Printing statistics text
            output_file_name_stat_text = os.path.join(curr_output_folder_name, output_stat_fulltest + ".txt")
            output_file = open(output_file_name_stat_text, "w")
            output_file.write(results_header_text + "\n")
            for r in result_list:
                output_file.write(str(r) + "\n")
            output_file.close()

            # Printing statistics html - Creating data table
            data_table = []
            for r in result_list:
                curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                for rep in motif_data.get_logo_list():
                    logo_file_name = os.path.join(rep, r.name + ".png")
                    if os.path.isfile(logo_file_name):
                        curr_motif_tuple = [logo_file_name, logo_width]
                        break
                curr_gene_tuple = ["View", gprofiler_link + ",".join(r.genes.genes)]
                data_table.append([r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                                   str(r.c), str(r.d), str(r.percent), str(r.back_percent), curr_gene_tuple])

            # Printing statistics html
            output_file_name_html = os.path.join(curr_output_folder_name, output_stat_fulltest + ".html")
            html = Html("Motif Enrichment Analysis", sitetest_link_dict)
            if curr_input.gene_set:
                html.add_heading(
                    "Results for <b>" + original_name + "</b> region <b>Site Test*</b> using genes from <b>" + curr_input.gene_set.name + "</b>",
                    align="center", bold=False)
                html.add_heading(
                    "* This test considered regions associated to the gene list given against background regions",
                    align="center", bold=False, size=3)
            else:
                html.add_heading(
                    "Results for <b>" + original_name + "</b> region <b>Site Test*</b> using all input regions",
                    align="center", bold=False)
                html.add_heading("* This test considered all regions against background regions",
                                 align="center", bold=False, size=3)

            html.add_zebra_table(html_header, html_col_size, html_type_list, data_table, align="center")
            html.write(output_file_name_html)

            # Removing files
            for e in to_remove_list:
                os.remove(e)

    ###################################################################################################
    # Heatmap
    ###################################################################################################

    # TODO

    ###################################################################################################
    # Network
    ###################################################################################################

    # TODO

    if __name__ == "__main__":
        main()
