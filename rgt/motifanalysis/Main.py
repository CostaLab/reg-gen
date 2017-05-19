###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import os
import sys
from glob import glob
import time
from random import seed
from optparse import OptionGroup
from shutil import copy

# Internal
from rgt import __version__
from rgt.Util import PassThroughOptionParser, ErrorHandler, MotifData, GenomeData, ImageData, Html, npath
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.GeneSet import GeneSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion
from Motif import Motif, Thresholds
from Match import match_single
from Statistics import multiple_test_correction, get_fisher_dict
from Util import Input, Result
from rgt.AnnotationSet import AnnotationSet

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

Authors: Eduardo G. Gusmao, Fabio Ticconi
"""


def is_bed(filename):
    _, ext = os.path.splitext(filename)

    if ext.lower() == ".bed":
        return True

    return False


def is_bb(filename):
    _, ext = os.path.splitext(filename)

    if ext.lower() == ".bb":
        return True

    return False


def bb_to_bed(filename):
    path, ext = os.path.splitext(filename)

    if ext.lower() == ".bb":
        # convert BB to BED
        bed_filename = os.path.join(path + ".bed")
        os.system(" ".join(["bigBedToBed", filename, bed_filename]))

        return bed_filename
    elif ext.lower() == ".bed":
        return filename
    else:
        raise ValueError("{} is neither a BED nor a BB".format(filename))


def bed_to_bb(filename, chrom_sizes_filename):
    path, ext = os.path.splitext(filename)

    if ext.lower() == ".bed":
        # convert BED to BB
        bb_filename = os.path.join(path + ".bb")
        os.system(" ".join(["bedToBigBed", filename, chrom_sizes_filename, bb_filename, "-verbose=0"]))

        return bb_filename
    elif ext.lower() == ".bb":
        return filename
    else:
        raise ValueError("{} is neither a BED nor a BB".format(filename))


def write_bed_color(region_set, filename, color):
    with open(filename, 'w') as f:
        for s in region_set:
            append = "\t".join([str(s.initial), str(s.final), color])
            print(str(s) + append, file=f)


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
    err = ErrorHandler()

    ###################################################################################################
    # Redirecting to Specific Functions
    ###################################################################################################

    # Redirecting Loop
    if sys.argv[1] == "--matching":
        main_matching()
    elif sys.argv[1] == "--enrichment":
        main_enrichment()
    else:
        err.throw_error("MOTIF_ANALYSIS_OPTION_ERROR")

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
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --matching [options] [input1.bed input2.bed ..]"

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
    parser.add_option("--rand-proportion", dest="rand_proportion", type="float", metavar="FLOAT",
                      help="If set, a random regions file will be created (eg, for later enrichment analysis). "
                           "The number of coordinates will be equal to this value times the size of the input regions. "
                           "We advise you use a value of at least 10.")
    parser.add_option("--norm-threshold", dest="norm_threshold", action="store_true", default=False,
                      help="If this option is used, the thresholds for all PWMs will be normalized by their length. "
                           "In this scheme, the threshold cutoff is evaluated in the regular way by the given fpr. "
                           "Then, all thresholds are divided by the length of the motif. The final threshold "
                           "consists of the average between all normalized motif thresholds. This single threshold "
                           "will be applied to all motifs.")
    parser.add_option("--use-only-motifs", dest="selected_motifs_filename", type="string", metavar="PATH",
                      help="Only use the motifs contained within this file (one for each line).")
    parser.add_option("--input-matrix", dest="input_matrix", type="string", metavar="PATH",
                      help="If an experimental matrix is provided, the input arguments will be ignored.")

    # Promoter-matching options
    group = OptionGroup(parser, "Promoter-regions matching options",
                        "Takes a list of genes, extracts their promoter regions and performs motif matching on these. "
                        "If a genes file is provided, the input files and experimental matrix will be ignored.")
    group.add_option("--gene-list", dest="promoter_genes_filename", type="string", metavar="PATH",
                     help="List of genes (one per line) to get the promoter regions from.")
    group.add_option("--make-background", dest="promoter_make_background", action="store_true", default=False,
                     help="If set, it will perform motif matching on the 'background regions', composed of "
                          "the promoters of all available genes. It doesn't require --gene-list.")
    group.add_option("--promoter-length", dest="promoter_length", type="int", metavar="INT", default=1000,
                     help="Length of the promoter region (in bp) to be extracted from each gene.")
    parser.add_option_group(group)

    # Output options
    group = OptionGroup(parser, "Output options",
                        "Where to put the output files and how to post-process them.")
    group.add_option("--output-location", dest="output_location", type="string", metavar="PATH",
                     help="Path where the output MPBS files will be written. Defaults to 'match' in the "
                          "current directory.")
    group.add_option("--bigbed", dest="bigbed", action="store_true", default=False,
                     help="If this option is used, all bed files will be written as bigbed.")
    group.add_option("--normalize-bitscore", dest="normalize_bitscore", action="store_true", default=False,
                     help="In order to print bigbed files the scores need to be normalized between 0 and 1000. "
                          "Don't use this option if real bitscores should be printed in the resulting bed file. "
                          "Without this option, bigbed files will never be created.")
    parser.add_option_group(group)

    # Processing Options
    options, arguments = parser.parse_args()

    # Additional Parameters
    matching_folder_name = "match"
    random_region_name = "random_regions"

    # we take care of conflicting parameters before going into the core of the method
    if options.promoter_genes_filename:
        # disable random regions and input matrix
        options.rand_proportion = None
        options.input_matrix = None

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if options.output_location:
        output_location = options.output_location
    else:
        output_location = npath(matching_folder_name)

    try:
        if not os.path.isdir(output_location):
            os.makedirs(output_location)
    except Exception:
        err.throw_error("MM_OUT_FOLDER_CREATION")

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
            err.throw_error("MM_MOTIFS_NOTFOUND", add_msg=options.selected_motifs_filename)

    ###################################################################################################
    # Reading Input Regions
    ###################################################################################################

    genomic_regions_dict = {}

    # get promoter regions from list of genes (both target and background)
    # TODO: should be more clever, allow precomputed regions etc
    if options.promoter_genes_filename:
        annotation = AnnotationSet(options.organism, alias_source=options.organism,
                                   protein_coding=True, known_only=True)

        target_genes = GeneSet("target_genes")
        target_genes.read(options.promoter_genes_filename)

        # TODO what do we do with unmapped genes? maybe just print them out
        target_regions = annotation.get_promoters(gene_set=target_genes, promoterLength=options.promoter_length)
        target_regions.name = "target_regions"
        target_regions.sort()
        output_file_name = npath(os.path.join(output_location, target_regions.name + ".bed"))
        target_regions.write_bed(output_file_name)

        genomic_regions_dict[target_regions.name] = target_regions

        if options.promoter_make_background:
            # background is made of all genes minus the target genes
            background_genes = GeneSet("background_genes")
            background_genes.get_all_genes(organism=options.organism)
            background_genes.subtract(target_genes)

            background_regions = annotation.get_promoters(gene_set=background_genes,
                                                          promoterLength=options.promoter_length)
            background_regions.name = "background_regions"
            background_regions.sort()
            output_file_name = npath(os.path.join(output_location, background_regions.name + ".bed"))
            background_regions.write_bed(output_file_name)

            genomic_regions_dict[background_regions.name] = background_regions

    # get experimental matrix, if available
    if options.input_matrix:
        try:
            exp_matrix = ExperimentalMatrix()
            exp_matrix.read(options.input_matrix)

            # if the matrix is present, the (empty) dictionary is overwritten
            genomic_regions_dict = exp_matrix.objectsDict
        except Exception:
            err.throw_error("MM_WRONG_EXPMAT")
    elif arguments:
        # get input files, if available
        for input_filename in arguments:
            name, _ = os.path.splitext(os.path.basename(input_filename))

            regions = GenomicRegionSet(name)
            regions.read_bed(npath(input_filename))

            genomic_regions_dict[name] = regions

    if not genomic_regions_dict:
        err.throw_error("DEFAULT_ERROR", add_msg="You must either specify an experimental matrix, or at least a "
                                                 "valid input file, or one of the 'promoter test' options.")

    max_region_len = 0
    max_region = None
    regions_to_match = []

    # Iterating on experimental matrix objects
    for k in genomic_regions_dict.keys():

        curr_genomic_region = genomic_regions_dict[k]

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
                max_region = curr_genomic_region

    ###################################################################################################
    # Creating random regions
    ###################################################################################################

    # if a random proportion is set, create random regions
    if options.rand_proportion:

        # Create random coordinates and name it random_regions
        rand_region = max_region.random_regions(options.organism, multiply_factor=options.rand_proportion, chrom_X=True)
        rand_region.sort()
        rand_region.name = random_region_name

        # Add random regions to the list of regions to perform matching on
        regions_to_match.append(rand_region)

        # Writing random regions
        output_file_name = npath(os.path.join(output_location, random_region_name))
        rand_bed_file_name = output_file_name + ".bed"
        rand_region.write_bed(rand_bed_file_name)

        # Verifying condition to write bb
        if options.bigbed:

            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            try:
                # Converting to big bed
                bed_to_bb(rand_bed_file_name, chrom_sizes_file)

                # removing previously-created BED file
                os.remove(rand_bed_file_name)
            except Exception:
                err.throw_warning("DEFAULT_WARNING")  # FIXME: maybe error instead?

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
        for motif_file_name in glob(npath(os.path.join(motif_repository, "*.pwm"))):
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
        output_bed_file = os.path.join(output_location, genomic_region_set.name + "_mpbs.bed")

        # GenomicRegionSet where all found MPBS regions are added
        output_grs = GenomicRegionSet("output")

        # Iterating on genomic regions
        for genomic_region in genomic_region_set.sequences:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            for motif in motif_list:
                grs = match_single(motif, sequence, genomic_region, unique_threshold, options.normalize_bitscore,
                                   # supposedly, python sort implementation works best with partially sorted sets
                                   sort=True)
                output_grs.combine(grs, change_name=False)

        output_grs.sort()

        # writing sorted regions to BED file
        output_grs.write_bed(output_bed_file)

        # Verifying condition to write bb
        if options.bigbed and options.normalize_bitscore:
            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            bed_to_bb(output_bed_file, chrom_sizes_file)

            # removing BED file
            os.remove(output_bed_file)


def main_enrichment():
    """
    Performs motif enrichment.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Parameters
    usage_message = "%prog --enrichment [options] <background_bed_file> [input1.bed input2.bed ..]"

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
    parser.add_option("--use-only-motifs", dest="selected_motifs_filename", type="string", metavar="PATH",
                      help="Only use the motifs contained within this file (one for each line).")
    parser.add_option("--matching-location", dest="match_location", type="string", metavar="PATH",
                      help="Directory where the matching output containing the MPBS files resides. "
                           "Defaults to 'match' in the current directory.")
    parser.add_option("--input-matrix", dest="input_matrix", type="string", metavar="PATH",
                      help="If an experimental matrix is provided, the input arguments will be ignored.")

    # Output Options
    parser.add_option("--output-location", dest="output_location", type="string", metavar="PATH",
                      help="Path where the output MPBS files will be written. Defaults to 'enrichment' in the "
                           "current directory.")
    parser.add_option("--print-thresh", dest="print_thresh", type="float", metavar="FLOAT", default=0.05,
                      help="Only MPBSs whose factor's enrichment corrected p-value are less than equal "
                           "this option are printed. Use 1.0 to print all MPBSs.")
    parser.add_option("--bigbed", dest="bigbed", action="store_true", default=False,
                      help="If this option is used, all bed files will be written as bigbed.")
    parser.add_option("--no-copy-logos", dest="no_copy_logos", action="store_true", default=False,
                      help="If set, the logos to be showed on the enrichment statistics page will NOT be "
                           "copied to a local directory; instead, the HTML result file will contain absolute "
                           "paths to the logos in your rgtdata folder.")

    # Processing Options
    options, arguments = parser.parse_args()

    if len(arguments) < 1:
        print(usage_message)
        sys.exit(1)

    background_filename = arguments[0]
    input_files = arguments[1:]  # empty list if no input files

    # Additional Parameters
    matching_folder_name = "match"
    enrichment_folder_name = "enrichment"
    gene_column_name = "genegroup"
    output_association_name = "coord_association"
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
        output_location = options.output_location
    else:
        output_location = os.path.join(os.getcwd(), enrichment_folder_name)

    # if the output folder doesn't exist, we create it
    try:
        if not os.path.isdir(output_location):
            os.makedirs(output_location)
    except Exception:
        err.throw_error("ME_OUT_FOLDER_CREATION")

    # Matching folder
    if options.match_location:
        match_location = options.match_location
    else:
        match_location = os.path.join(os.getcwd(), matching_folder_name)

    # the matching directory must exist, for obvious reasons
    if not os.path.isdir(match_location):
        err.throw_error("ME_MATCH_NOTFOUND")

    # Background file must exist
    background_original_filename = background_filename
    if not os.path.isfile(background_filename):
        err.throw_error("DEFAULT_ERROR", add_msg="Background file does not exist or is not readable.")
    elif is_bb(background_filename):
        background_filename = bb_to_bed(background_filename)
    elif is_bed(background_filename):
        pass
    else:
        err.throw_error("DEFAULT_ERROR", add_msg="Background file must be in either BED or BigBed format.")

    # Background MPBS file must exist
    path, ext = os.path.splitext(background_filename)

    # first we check at matching folder location
    background_mpbs_filename = os.path.join(match_location, os.path.basename(path) + "_mpbs" + ext)

    if not os.path.isfile(background_mpbs_filename):
        # if not found, we search at background file location
        background_mpbs_filename = os.path.join(path + "_mpbs" + ext)

        if not os.path.isfile(background_mpbs_filename):
            err.throw_error("DEFAULT_ERROR", add_msg="Background MPBS file does not exist or is not readable. "
                                                     "It must be located at either the matching location, or in the "
                                                     "same directory of the Background BED/BigBed file. "
                                                     "Note: it must be consistent with the background file, ie "
                                                     "if the background is BED, the MPBS must be a BED file too. Same "
                                                     "for BigBed.")

    background_mpbs_original_filename = background_mpbs_filename

    if is_bb(background_mpbs_filename):
        background_mpbs_filename = bb_to_bed(background_mpbs_filename)
    elif is_bed(background_mpbs_filename):
        pass
    else:
        err.throw_error("DEFAULT_ERROR", add_msg="Background MPBS file must be in either BED or BigBed format.")

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
            err.throw_error("MM_MOTIFS_NOTFOUND", add_msg=options.selected_motifs_filename)

    # Default image data
    image_data = ImageData()

    genomic_regions_dict = {}
    exp_matrix_fields_dict = {}

    # will be set if genelists are used in the experimental matrix
    flag_gene = False

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # get experimental matrix, if available
    if options.input_matrix:
        try:
            exp_matrix = ExperimentalMatrix()
            exp_matrix.read(options.input_matrix)

            # if the matrix is present, the (empty) dictionary is overwritten
            genomic_regions_dict = exp_matrix.objectsDict

            # Reading dictionary grouped by fields (only for gene association)
            try:
                exp_matrix_fields_dict = exp_matrix.fieldsDict[gene_column_name]
                flag_gene = True
            except KeyError:
                flag_gene = False

            del exp_matrix

        except Exception:
            err.throw_error("MM_WRONG_EXPMAT")
    elif input_files:
        # get input files, if available
        for input_filename in input_files:
            name, _ = os.path.splitext(os.path.basename(input_filename))

            regions = GenomicRegionSet(name)

            try:
                regions.read_bed(os.path.abspath(input_filename))
            except:
                err.throw_error("DEFAULT_ERROR", add_msg="Input file {} could not be loaded.".format(input_filename))
            genomic_regions_dict[name] = regions

    if not genomic_regions_dict:
        err.throw_error("DEFAULT_ERROR", add_msg="You must either specify an experimental matrix, "
                                                 "or at least a valid input file.")

    ###################################################################################################
    # Reading Regions & Gene Lists
    ###################################################################################################

    # Initializations
    input_list = []

    if flag_gene:  # Genelist and full site analysis will be performed

        # Iterating on experimental matrix fields
        for g in exp_matrix_fields_dict.keys():

            # Create input which will contain all regions associated with such gene group
            curr_input = Input(None, [])

            # This flag will be used to see if there are two gene files associated with
            # the same gene label on genegroup column
            flag_foundgeneset = False

            # Iterating over the genomic regions
            for k in exp_matrix_fields_dict[g]:

                curr_object = genomic_regions_dict[k]

                # If the current object is a GenomicRegionSet
                if isinstance(curr_object, GenomicRegionSet):
                    # Sorting input region
                    curr_object.sort()

                    # Updating Input object
                    curr_input.region_list.append(curr_object)

                # If the current object is a GeneSet
                if isinstance(curr_object, GeneSet):

                    # Updating Input object
                    # The name in gene_group column will be used. The 'name' column for genes are not used.
                    curr_object.name = g
                    if not flag_foundgeneset:
                        curr_input.gene_set = curr_object
                        flag_foundgeneset = True
                    else:
                        err.throw_warning("ME_MANY_GENESETS")

            if not flag_foundgeneset:
                err.throw_warning("ME_FEW_GENESETS")

            # Update input list
            input_list.append(curr_input)

    else:  # Only full site analysis will be performed

        # Create single input which will contain all regions
        single_input = Input(None, [])

        # Iterating on experimental matrix objects
        for k in genomic_regions_dict.keys():

            curr_object = genomic_regions_dict[k]

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
    motif_names.sort()

    ###################################################################################################
    # Background Statistics
    ###################################################################################################

    background = GenomicRegionSet("background")
    background.read_bed(background_filename)
    background_mpbs = GenomicRegionSet("background_mpbs")
    background_mpbs.read_bed(background_mpbs_filename)

    # Evaluating background statistics
    bg_c_dict, bg_d_dict, _, _ = get_fisher_dict(motif_names, background, background_mpbs)

    # removing temporary BED files if the originals were BBs
    if is_bb(background_original_filename):
        os.remove(background_filename)
    if is_bb(background_mpbs_original_filename):
        os.remove(background_mpbs_filename)

    # scheduling region sets for garbage collection
    del background
    del background_mpbs

    ###################################################################################################
    # Enrichment Statistics
    ###################################################################################################

    # Creating link dictionary for HTML file
    genetest_link_dict = dict()
    sitetest_link_dict = dict()
    link_location = "../"
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
                curr_output_folder_name = os.path.join(output_location, grs.name + "__" + curr_input.gene_set.name)
            else:
                curr_output_folder_name = os.path.join(output_location, grs.name)

            curr_output_folder_name = npath(curr_output_folder_name)

            if not os.path.isdir(curr_output_folder_name):
                os.makedirs(curr_output_folder_name)

            # Verifying if MPBS file exists
            curr_mpbs_glob = glob(os.path.join(match_location, original_name + "_mpbs.*"))
            try:
                curr_mpbs_file_name = npath(curr_mpbs_glob[0])
            except Exception:
                err.throw_warning("DEFAULT_ERROR", add_msg="File {} does not have a matching MPBS file. "
                                                           "Ignoring.".format(original_name))
                # skip to next genomic region set
                continue

            if is_bb(curr_mpbs_file_name):
                curr_mpbs_bed_name = bb_to_bed(curr_mpbs_file_name)

                # at the end of calculations, we'll remove all the temporary bed files we created
                to_remove_list.append(curr_mpbs_bed_name)
            elif is_bed(curr_mpbs_file_name):
                curr_mpbs_bed_name = curr_mpbs_file_name
            else:
                err.throw_warning("DEFAULT_ERROR", add_msg="The matching MPBS file for {} is neither in BED nor BigBed "
                                                           "format. Ignoring.".format(original_name))
                continue

            curr_mpbs = GenomicRegionSet("curr_mpbs")
            curr_mpbs.read_bed(curr_mpbs_bed_name)
            curr_mpbs.sort()

            ###################################################################################################
            # Gene Evidence Statistics
            ###################################################################################################

            if curr_input.gene_set:

                # Performing association of input region with gene_set
                grs = grs.gene_association(curr_input.gene_set, options.organism, options.promoter_length,
                                           options.maximum_association_length)

                # Writing gene-coordinate association
                output_file_name = npath(os.path.join(curr_output_folder_name, output_association_name + ".bed"))
                output_file = open(output_file_name, "w")
                for gr in grs:
                    if gr.name == ".":
                        curr_name = "."
                    elif not gr.proximity:
                        curr_name = ":".join([e if e[0] != "." else e[1:] for e in gr.name.split(":")])
                    else:
                        curr_gene_list = [e if e[0] != "." else e[1:] for e in gr.name.split(":")]
                        curr_prox_list = gr.proximity.split(":")
                        curr_name = ":".join([e[0] + "_" + e[1] for e in zip(curr_gene_list, curr_prox_list)])
                    output_file.write("\t".join([str(e) for e in [gr.chrom, gr.initial, gr.final, curr_name]]) + "\n")
                output_file.close()
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    try:
                        bed_to_bb(output_file_name, chrom_sizes_file)
                        os.remove(output_file_name)
                    except Exception:
                        pass  # TODO: warning?

                # Writing ev and nev regions
                ev_regions = GenomicRegionSet("ev_regions")
                nev_regions = GenomicRegionSet("nev_regions")

                for gr in grs:
                    if len([e for e in gr.name.split(":") if e[0] != "."]) > 0:
                        ev_regions.add(GenomicRegion(gr.chrom, gr.initial, gr.final,
                                                     name=gr.name, orientation=gr.orientation, data=gr.data))
                    else:
                        nev_regions.add(GenomicRegion(gr.chrom, gr.initial, gr.final,
                                                      name=gr.name, orientation=gr.orientation, data=gr.data))

                # Calculating statistics
                a_dict, b_dict, ev_genes_dict, ev_mpbs_dict = get_fisher_dict(motif_names, ev_regions, curr_mpbs,
                                                                              gene_set=True, mpbs_set=True)

                c_dict, d_dict, _, nev_mpbs_dict = get_fisher_dict(motif_names, nev_regions, curr_mpbs,
                                                                   gene_set=True, mpbs_set=True)

                # Performing fisher test
                result_list = []
                for k in motif_names:
                    r = Result()
                    r.name = k
                    r.a = a_dict[k]
                    r.b = b_dict[k]
                    r.c = c_dict[k]
                    r.d = d_dict[k]
                    r.percent = float(r.a) / float(r.a + r.b)
                    r.back_percent = float(r.c) / float(r.c + r.d)
                    r.genes = ev_genes_dict[k]
                    try:
                        p = pvalue(r.a, r.b, r.c, r.d)
                        r.p_value = p.right_tail
                    except Exception:
                        r.p_value = 1.0
                    result_list.append(r)

                # Performing multiple test correction
                multiple_corr_rej, multiple_corr_list = multiple_test_correction([e.p_value for e in result_list],
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

                # filtering out MPBS with low corr. p-value
                ev_mpbs_grs_filtered = GenomicRegionSet("ev_mpbs_filtered")
                for m, ev_mpbs_grs in ev_mpbs_dict.items():
                    for region in ev_mpbs_grs:
                        if corr_pvalue_dict[m] <= options.print_thresh:
                            ev_mpbs_grs_filtered.add(region)
                del ev_mpbs_dict

                nev_mpbs_grs_filtered = GenomicRegionSet("nev_mpbs_filtered")
                for m, nev_mpbs_grs in nev_mpbs_dict.items():
                    for region in nev_mpbs_grs:
                        if corr_pvalue_dict[m] <= options.print_thresh:
                            nev_mpbs_grs_filtered.add(region)
                del nev_mpbs_dict

                output_mpbs_filtered_ev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bed")
                output_mpbs_filtered_nev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_nev + ".bed")

                # sorting and saving to BED
                ev_mpbs_grs_filtered.sort()
                write_bed_color(ev_mpbs_grs_filtered, npath(output_mpbs_filtered_ev_bed), ev_color)
                del ev_mpbs_grs_filtered

                nev_mpbs_grs_filtered.sort()
                write_bed_color(nev_mpbs_grs_filtered, npath(output_mpbs_filtered_nev_bed), nev_color)
                del nev_mpbs_grs_filtered

                # Converting ev and nev mpbs to bigbed
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    try:
                        # create big bed files
                        bed_to_bb(output_mpbs_filtered_ev_bed, chrom_sizes_file)
                        bed_to_bb(output_mpbs_filtered_nev_bed, chrom_sizes_file)

                        # temporary BED files to be removed later
                        to_remove_list.append(output_mpbs_filtered_ev_bed)
                        to_remove_list.append(output_mpbs_filtered_nev_bed)
                    except Exception:
                        pass  # TODO: warning?

                # Printing statistics text
                output_file_name_stat_text = os.path.join(curr_output_folder_name, output_stat_genetest + ".txt")
                output_file = open(npath(output_file_name_stat_text), "w")
                output_file.write(results_header_text + "\n")
                for r in result_list:
                    output_file.write(str(r) + "\n")
                output_file.close()

                # unless explicitly forbidden, we copy the logo images locally
                if not options.no_copy_logos:
                    logo_dir_path = npath(os.path.join(output_location, "logos"))
                    try:
                        os.stat(logo_dir_path)
                    except:
                        os.mkdir(logo_dir_path)

                # Printing statistics html - Creating data table
                data_table = []
                for r in result_list:
                    curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                    for rep in motif_data.get_logo_list():
                        logo_file_name = npath(os.path.join(rep, r.name + ".png"))

                        if os.path.isfile(logo_file_name):
                            if not options.no_copy_logos:
                                copy(logo_file_name, npath(os.path.join(logo_dir_path, r.name + ".png")))

                                # use relative paths in the html
                                # FIXME can we do it in a better way? (inside the Html class)
                                logo_file_name = os.path.join("..", "logos", r.name + ".png")

                            curr_motif_tuple = [logo_file_name, logo_width]
                            break
                    curr_gene_tuple = ["View", gprofiler_link + "+".join(r.genes.genes)]
                    data_table.append(
                        [r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                         str(r.c), str(r.d), str(r.percent), str(r.back_percent), curr_gene_tuple])

                # Printing statistics html - Writing to HTML
                output_file_name_html = os.path.join(curr_output_folder_name, output_stat_genetest + ".html")
                fig_path = npath(os.path.join(output_location, "fig"))
                html = Html("Motif Enrichment Analysis", genetest_link_dict, fig_dir=fig_path)
                html.add_heading("Results for <b>" + original_name + "</b> "
                                 "region <b>Gene Test*</b> using genes from <b>" + curr_input.gene_set.name +
                                 "</b>",
                                 align="center", bold=False)
                html.add_heading("* This gene test considered regions associated to the given "
                                 "gene list against regions not associated to the gene list",
                                 align="center", bold=False, size=3)
                html.add_zebra_table(html_header, html_col_size, html_type_list, data_table, align="center")
                html.write(npath(output_file_name_html))

            else:

                # Association still needs to be done with all genes in order to print gene list
                grs = grs.gene_association(None, options.organism, options.promoter_length,
                                           options.maximum_association_length)

                # Calculating statistics
                a_dict, b_dict, ev_genes_dict, ev_mpbs_dict = get_fisher_dict(motif_names, grs, curr_mpbs,
                                                                              gene_set=True, mpbs_set=True)

            ###################################################################################################
            # Final wrap-up
            ###################################################################################################

            # Performing fisher test
            result_list = []
            for k in motif_names:
                r = Result()
                r.name = k
                r.a = a_dict[k]
                r.b = b_dict[k]
                r.c = bg_c_dict[k]
                r.d = bg_d_dict[k]
                r.percent = float(r.a) / float(r.a + r.b)
                r.back_percent = float(r.c) / float(r.c + r.d)
                r.genes = ev_genes_dict[k]
                try:
                    p = pvalue(r.a, r.b, r.c, r.d)
                    r.p_value = p.right_tail
                except Exception:
                    r.p_value = 1.0
                result_list.append(r)

            # Performing multiple test correction
            multiple_corr_rej, multiple_corr_list = multiple_test_correction([e.p_value for e in result_list],
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

                # filtering out MPBS with low corr. p-value,
                ev_mpbs_grs_filtered = GenomicRegionSet("ev_mpbs_filtered")

                for m, ev_mpbs_grs in ev_mpbs_dict.items():
                    for region in ev_mpbs_grs:
                        if corr_pvalue_dict[m] <= options.print_thresh:
                            ev_mpbs_grs_filtered.add(region)
                del ev_mpbs_dict

                output_mpbs_filtered_ev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bed")
                output_mpbs_filtered_ev_bed = npath(output_mpbs_filtered_ev_bed)

                # sorting and saving to BED
                ev_mpbs_grs_filtered.sort()
                write_bed_color(ev_mpbs_grs_filtered, output_mpbs_filtered_ev_bed, ev_color)
                del ev_mpbs_grs_filtered

                # Converting ev mpbs to bigbed
                if options.bigbed:
                    chrom_sizes_file = genome_data.get_chromosome_sizes()

                    try:
                        # create big bed file
                        bed_to_bb(output_mpbs_filtered_ev_bed, chrom_sizes_file)

                        # temporary BED file to be removed later
                        to_remove_list.append(output_mpbs_filtered_ev_bed)
                    except Exception:
                        pass  # WARNING

            # Printing statistics text
            output_file_name_stat_text = os.path.join(curr_output_folder_name, output_stat_fulltest + ".txt")
            output_file = open(npath(output_file_name_stat_text), "w")
            output_file.write(results_header_text + "\n")
            for r in result_list:
                output_file.write(str(r) + "\n")
            output_file.close()

            # unless explicitly forbidden, we copy the logo images locally
            if not options.no_copy_logos:
                logo_dir_path = npath(os.path.join(output_location, "logos"))
                try:
                    os.stat(logo_dir_path)
                except:
                    os.mkdir(logo_dir_path)

            # Printing statistics html - Creating data table
            data_table = []
            for r in result_list:
                curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                for rep in motif_data.get_logo_list():
                    logo_file_name = npath(os.path.join(rep, r.name + ".png"))

                    if os.path.isfile(logo_file_name):
                        if not options.no_copy_logos:
                            copy(logo_file_name, npath(os.path.join(logo_dir_path, r.name + ".png")))

                            # use relative paths in the html
                            # FIXME can we do it in a better way? (inside the Html class)
                            logo_file_name = os.path.join("..", "logos", r.name + ".png")

                        curr_motif_tuple = [logo_file_name, logo_width]
                        break
                curr_gene_tuple = ["View", gprofiler_link + "+".join(r.genes.genes)]
                data_table.append([r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                                   str(r.c), str(r.d), str(r.percent), str(r.back_percent), curr_gene_tuple])

            # Printing statistics html
            output_file_name_html = os.path.join(curr_output_folder_name, output_stat_fulltest + ".html")
            fig_path = npath(os.path.join(output_location, "fig"))
            html = Html("Motif Enrichment Analysis", sitetest_link_dict, fig_dir=fig_path)
            if curr_input.gene_set:
                html.add_heading("Results for <b>" + original_name +
                                 "</b> region <b>Site Test*</b> using genes from <b>" + curr_input.gene_set.name +
                                 "</b>", align="center", bold=False)
                html.add_heading("* This test considered regions associated to the given gene list "
                                 "against background regions", align="center", bold=False, size=3)
            else:
                html.add_heading("Results for <b>" + original_name +
                                 "</b> region <b>Site Test*</b> using all input regions", align="center", bold=False)
                html.add_heading("* This test considered all input regions against background regions",
                                 align="center", bold=False, size=3)

            html.add_zebra_table(html_header, html_col_size, html_type_list, data_table, align="center")
            html.write(npath(output_file_name_html))

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
