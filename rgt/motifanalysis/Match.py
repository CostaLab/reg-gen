###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import os
from glob import glob
from optparse import OptionGroup

# Internal
from rgt.Util import PassThroughOptionParser, ErrorHandler, MotifData, GenomeData, npath
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.GeneSet import GeneSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion
from rgt.AnnotationSet import AnnotationSet
from Motif import Motif, Thresholds
from Util import bed_to_bb

# External
from pysam import Fastafile
from MOODS import tools, scan


###################################################################################################
# Functions
###################################################################################################

def main_matching():
    """
    Performs motif matching.
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
    parser.add_option("--pseudocounts", dest="pseudocounts", type="float", metavar="FLOAT", default=1.0,
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
        target_regions.write(output_file_name)

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
            background_regions.write(output_file_name)

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
            regions.read(npath(input_filename))

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
        rand_region.write(rand_bed_file_name)

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
        motif_list.append(Motif(motif_file_name, options.pseudocounts, options.fpr, thresholds))

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

        # must remove it because we append the MPBS
        if os.path.isfile(output_bed_file):
            os.remove(output_bed_file)

        # Iterating on genomic regions
        for genomic_region in genomic_region_set:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            grs = GenomicRegionSet("tmp")

            for motif in motif_list:
                match_single(motif, sequence, genomic_region, unique_threshold, options.normalize_bitscore, output=grs)

            # TODO: measure and document speed/memory improvements, if any.
            # grs = match_multiple(motif_list, sequence, genomic_region)

            grs.write(output_bed_file, mode="a")

        # Verifying condition to write bb
        if options.bigbed and options.normalize_bitscore:
            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            bed_to_bb(output_bed_file, chrom_sizes_file)

            # removing BED file
            os.remove(output_bed_file)


def match_single(motif, sequence, genomic_region, unique_threshold=None, normalize_bitscore=False, output=None):
    """
    Performs motif matching given sequence and the motif.pssm passed as parameter.
    The genomic_region is needed to evaluate the correct binding position.

    Keyword arguments:
    motif -- a Motif object to use for genome scanning.
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
    unique_threshold -- If this argument is provided, the motif search will be made using a threshold of 0 and
                        then accepting only the motif matches with bitscore/motif_length >= unique_threshold.
    normalize_bitscore -- If True, it normalises the scores between 0 and 1000. Necessary for bigbed conversion.
    output -- A GenomicRegionSet where all matching GenomicRegions will be appended.
        
    Return:
    Either the "output" GenomicRegionSet if provided, or a newly-instantiated one.
    """

    # Establishing threshold
    if unique_threshold:
        current_threshold = 0.0
        threshold = unique_threshold
        motif_max = motif.max / motif.len
    else:
        current_threshold = motif.threshold
        threshold = motif.threshold
        motif_max = motif.max

    # Performing motif matching
    results = scan.scan(sequence, [motif.pssm], motif.bg, [current_threshold], 7, motif.alphabet)

    if output is None:
        output = GenomicRegionSet("mpbs")

    pos_start = genomic_region.initial
    chrom = genomic_region.chrom

    for search_result in results:
        for r in search_result:
            position = r.pos
            score = r.score

            score_len = score / motif.len

            # Verifying unique threshold acceptance
            if unique_threshold and score_len < unique_threshold:
                continue

            # If match forward strand
            if position >= 0:
                p1 = pos_start + position
                strand = "+"
            # If match reverse strand
            elif not motif.is_palindrome:
                p1 = pos_start - position
                strand = "-"
            else:
                continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if normalize_bitscore:
                # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if motif_max > threshold:
                    score = int(((score - threshold) * 1000.0) / (motif_max - threshold))
                else:
                    score = 1000
            elif unique_threshold:
                score = score_len

            output.add(GenomicRegion(chrom, int(p1), int(p2), name=motif.name, orientation=strand, data=str(score)))

    return output


# TODO must double-check/fix the normalisation.
def match_multiple(motifs, sequence, genomic_region, unique_threshold=None, normalize_bitscore=False, output=None):
    """
    More efficient than calling match_single on every motif.

    Keyword arguments:
    motifs -- a list of Motif objects to use for genome scanning.
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
    unique_threshold -- If this argument is provided, the motif search will be made using a threshold of 0 and
                        then accepting only the motif matches with bitscore/motif_length >= unique_threshold.
    normalize_bitscore -- If True, it normalises the scores between 0 and 1000. Necessary for bigbed conversion.
    output -- A GenomicRegionSet where all matching GenomicRegions will be appended.

    Return:
    Either the "output" GenomicRegionSet if provided, or a newly-instantiated one.
    """

    pssm_list = []
    thresholds = []
    for motif in motifs:
        if unique_threshold:
            thresholds.append(0.0)
        else:
            thresholds.append(motif.threshold)
            pssm_list.append(motif.pssm)

    # Performing motif matching
    # TODO: we can expand this to use bg from sequence, for example,
    # or from organism.
    bg = tools.flat_bg(4)
    results = scan.scan_dna(sequence, pssm_list, bg, thresholds, 7)

    if output is None:
        output = GenomicRegionSet("mpbs")

    pos_start = genomic_region.initial
    chrom = genomic_region.chrom

    for i, search_result in enumerate(results):
        motif = motifs[i]

        if unique_threshold:
            motif_max = motif.max / motif.len
            threshold = unique_threshold
        else:
            motif_max = motif.max
            threshold = motif.threshold

        for r in search_result:
            position = r.pos
            score = r.score

            score_len = score / motif.len

            # Verifying unique threshold acceptance
            if unique_threshold and score_len < unique_threshold:
                continue

            # If match forward strand
            if position >= 0:
                p1 = pos_start + position
                strand = "+"
            # If match reverse strand
            elif not motif.is_palindrome:
                p1 = pos_start - position
                strand = "-"
            else:
                continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if normalize_bitscore:
                # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if motif_max > threshold:
                    score = int(((score - threshold) * 1000.0) / (motif_max - threshold))
                else:
                    score = 1000
            elif unique_threshold:
                score = score_len

            output.add(GenomicRegion(chrom, int(p1), int(p2), name=motif.name, orientation=strand, data=str(score)))

    return output
