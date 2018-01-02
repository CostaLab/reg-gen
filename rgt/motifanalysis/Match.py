###################################################################################################
# Libraries
###################################################################################################

# Python 3 compatibility
from __future__ import print_function
from __future__ import division

# Python
import os
from glob import glob
import time
import sys

# Internal
from ..Util import ErrorHandler, MotifData, GenomeData, npath
from ..ExperimentalMatrix import ExperimentalMatrix
from ..GeneSet import GeneSet
from ..GenomicRegionSet import GenomicRegionSet
from ..GenomicRegion import GenomicRegion
from ..AnnotationSet import AnnotationSet
from .Motif import Motif, ThresholdTable
from .Util import bed_to_bb

# External
from pysam import Fastafile
from MOODS import tools, scan


###################################################################################################
# Functions
###################################################################################################

def options(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19 or hg38.")
    parser.add_argument("--fpr", type=float, metavar="FLOAT", default=0.0001,
                        help="False positive rate cutoff.")
    parser.add_argument("--pseudocounts", type=float, metavar="FLOAT", default=1.0,
                        help="Pseudocounts to be added to raw counts of each PFM.")
    parser.add_argument("--rand-proportion", type=float, metavar="FLOAT",
                        help="If set, a random regions file will be created (eg, for later enrichment analysis). "
                             "The number of coordinates will be equal to this value times the size of the input"
                             "regions. We advise you use a value of at least 10.")
    parser.add_argument("--norm-threshold", action="store_true", default=False,
                        help="If this option is used, the thresholds for all PWMs will be normalized by their length. "
                             "In this scheme, the threshold cutoff is evaluated in the regular way by the given fpr. "
                             "Then, all thresholds are divided by the length of the motif. The final threshold "
                             "consists of the average between all normalized motif thresholds. This single threshold "
                             "will be applied to all motifs.")
    parser.add_argument("--use-only-motifs", dest="selected_motifs_filename", type=str, metavar="PATH",
                        help="Only use the motifs contained within this file (one for each line).")
    parser.add_argument("--motif-dbs", type=str, metavar="PATH", nargs="+",
                        help="New 'motif DB' folders to use instead of the ones within "
                             "the RGTDATA folder. Each folder must contain PWM files.")

    # Promoter-matching args
    group = parser.add_argument_group("Promoter-regions matching",
                                      "These arguments are only used with the --promoters-only option (for the "
                                      "purpose of matching only on the promoters of all or a subset of genes)")
    group.add_argument("--target-genes", dest="target_genes_filename", type=str, metavar="PATH",
                       help="List of genes (one per line) to get the promoter regions from.")
    group.add_argument("--make-background", dest="promoter_make_background", action="store_true", default=False,
                       help="If set, it will perform motif matching on the 'background regions', composed of "
                            "the promoters of all available genes (minus the target genes, if specified). "
                            "It doesn't require --target-genes.")
    group.add_argument("--promoter-length", type=int, metavar="INT", default=1000,
                       help="Length of the promoter region (in bp) to be extracted from each gene.")

    # Output args
    group = parser.add_argument_group("Output",
                                      "Where to put the output files and how to post-process them.")
    group.add_argument("--output-location", type=str, metavar="PATH",
                       help="Path where the output MPBS files will be written. Defaults to 'match' in the "
                            "current directory.")
    group.add_argument("--bigbed", action="store_true", default=False,
                       help="If this option is used, all bed files will be written as bigbed.")
    group.add_argument("--normalize-bitscore", action="store_true", default=False,
                       help="In order to print bigbed files the scores need to be normalized between 0 and 1000. "
                            "Don't use this option if real bitscores should be printed in the resulting bed file. "
                            "Without this option, bigbed files will never be created.")

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--input-matrix", type=str, metavar="matrix.txt",
                       help="The experimental matrix allows the specification of gene-association rules among "
                            "input files (see online documentation for details).")
    group.add_argument('--promoters-only', action="store_true",
                       help="If you ONLY want to perform promoter matching without providing any input file/matrix. "
                            "If --target-genes is not provided, then all available promoters will be matched against. "
                            "Note how this makes '--make-background' redundant.")
    group.add_argument('--input-files', metavar="regions.bed", nargs='+', type=str,
                       help='BED files to perform motif matching on.')


def main(args):
    """
    Performs motif matching.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Initializing Error Handler
    err = ErrorHandler()

    # Additional Parameters
    matching_folder_name = "match"
    random_region_name = "random_regions"

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if args.output_location:
        output_location = args.output_location
    else:
        output_location = npath(matching_folder_name)
    print(">> output location:", output_location)

    # Default genomic data
    genome_data = GenomeData(args.organism)

    print(">> genome:", genome_data.organism)
    print(">> pseudocounts:", args.pseudocounts)
    print(">> fpr threshold:", args.fpr)

    # Reading motif file
    selected_motifs = []

    if args.selected_motifs_filename:
        try:
            with open(args.selected_motifs_filename) as f:
                selected_motifs = f.read().splitlines()
                selected_motifs = filter(None, selected_motifs)
                print(">> motif file loaded:", len(selected_motifs), "motifs")
        except Exception:
            err.throw_error("MM_MOTIFS_NOTFOUND", add_msg=args.selected_motifs_filename)

    ###################################################################################################
    # Reading Input Regions
    ###################################################################################################

    genomic_regions_dict = {}

    # get experimental matrix, if available
    if args.input_matrix:
        try:
            exp_matrix = ExperimentalMatrix()
            exp_matrix.read(args.input_matrix)

            # if the matrix is present, the (empty) dictionary is overwritten
            genomic_regions_dict = exp_matrix.objectsDict

            print(">> experimental matrix loaded")

        except Exception:
            err.throw_error("MM_WRONG_EXPMAT")
    elif args.input_files:
        # get input files, if available
        for input_filename in args.input_files:
            name, _ = os.path.splitext(os.path.basename(input_filename))

            regions = GenomicRegionSet(name)
            regions.read(npath(input_filename))

            genomic_regions_dict[name] = regions

        print(">> input regions BED files loaded")

    # we put this here because we don't want to create the output directory unless we
    # are sure the initialisation (including loading input files) worked
    try:
        if not os.path.isdir(output_location):
            os.makedirs(output_location)
    except Exception:
        err.throw_error("MM_OUT_FOLDER_CREATION")

    annotation = None
    target_genes = None
    # get promoter regions from list of genes (both target and background)
    # TODO: should be more clever, allow precomputed regions etc
    if args.target_genes_filename:
        annotation = AnnotationSet(args.organism, alias_source=args.organism,
                                   protein_coding=True, known_only=True)

        target_genes = GeneSet("target_genes")
        target_genes.read(args.target_genes_filename)

        # TODO what do we do with unmapped genes? maybe just print them out
        target_regions = annotation.get_promoters(gene_set=target_genes, promoter_length=args.promoter_length)
        target_regions.name = "target_regions"
        target_regions.sort()
        output_file_name = npath(os.path.join(output_location, target_regions.name + ".bed"))
        target_regions.write(output_file_name)

        genomic_regions_dict[target_regions.name] = target_regions

        print(">> target promoter file created:", len(target_regions), "regions")

    # we make a background in case it's requested, but also in case a list of target genes has not been
    # provided
    if args.promoter_make_background or (args.promoters_only and not args.target_genes_filename):
        if not annotation:
            annotation = AnnotationSet(args.organism, alias_source=args.organism,
                                       protein_coding=True, known_only=True)

        # background is made of all known genes minus the target genes (if any)
        background_genes = GeneSet("background_genes")
        background_genes.get_all_genes(organism=args.organism)

        if target_genes:
            background_genes.subtract(target_genes)

        background_regions = annotation.get_promoters(gene_set=background_genes,
                                                      promoter_length=args.promoter_length)
        background_regions.name = "background_regions"
        background_regions.sort()
        output_file_name = npath(os.path.join(output_location, background_regions.name + ".bed"))
        background_regions.write(output_file_name)

        genomic_regions_dict[background_regions.name] = background_regions

        print(">> background promoter file created:", len(background_regions), "regions")

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
    if args.rand_proportion:

        # Create random coordinates and name it random_regions
        rand_region = max_region.random_regions(args.organism, multiply_factor=args.rand_proportion, chrom_X=True)
        rand_region.sort()
        rand_region.name = random_region_name

        # Add random regions to the list of regions to perform matching on
        regions_to_match.append(rand_region)

        # Writing random regions
        output_file_name = npath(os.path.join(output_location, random_region_name))
        rand_bed_file_name = output_file_name + ".bed"
        rand_region.write(rand_bed_file_name)

        # Verifying condition to write bb
        if args.bigbed:

            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            try:
                # Converting to big bed
                bed_to_bb(rand_bed_file_name, chrom_sizes_file)

                # removing previously-created BED file
                os.remove(rand_bed_file_name)
            except Exception:
                err.throw_warning("DEFAULT_WARNING")  # FIXME: maybe error instead?

        print(">> random regions file created:", len(rand_region), "regions")

    ###################################################################################################
    # Creating PWMs
    ###################################################################################################

    # Initialization
    motif_list = []

    # Default motif data
    motif_data = MotifData()
    if args.motif_dbs:
        # must overwrite the default DBs
        motif_data.set_custom(args.motif_dbs)
        print(">> custom motif repositories:", motif_data.repositories_list)
    else:
        print(">> motif repositories:", motif_data.repositories_list)

    # Creating thresholds object
    threshold_table = ThresholdTable(motif_data)

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
        motif_list.append(Motif(motif_file_name, args.pseudocounts, args.fpr, threshold_table))

    # Performing normalized threshold strategy if requested
    if args.norm_threshold:
        threshold_list = [motif.threshold / motif.len for motif in motif_list]
        unique_threshold = sum(threshold_list) / len(threshold_list)
    else:
        unique_threshold = None

    scanner = scan.Scanner(7)
    pssm_list = []
    thresholds = []
    for motif in motif_list:
        if unique_threshold:
            thresholds.append(0.0)
        else:
            thresholds.append(motif.threshold)
            pssm_list.append(motif.pssm)

    # Performing motif matching
    # TODO: we can expand this to use bg from sequence, for example,
    # or from organism.
    bg = tools.flat_bg(4)
    scanner.set_motifs(pssm_list, bg, thresholds)

    ###################################################################################################
    # Motif Matching
    ###################################################################################################

    # Creating genome file
    genome_file = Fastafile(genome_data.get_genome())

    print()

    # Iterating on list of genomic regions
    for grs in regions_to_match:

        start = time.time()
        print(">> matching [", grs.name, "], ", len(grs), " regions... ", end='', sep="")
        sys.stdout.flush()

        # Initializing output bed file
        output_bed_file = os.path.join(output_location, grs.name + "_mpbs.bed")

        # must remove it because we append the MPBS
        if os.path.isfile(output_bed_file):
            os.remove(output_bed_file)

        # Iterating on genomic regions
        for genomic_region in grs:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            grs = match_multiple(scanner, motif_list, sequence, genomic_region)
            grs.write(output_bed_file, mode="a")

        del grs.sequences[:]

        # Verifying condition to write bb
        if args.bigbed and args.normalize_bitscore:
            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            bed_to_bb(output_bed_file, chrom_sizes_file)

            # removing BED file
            os.remove(output_bed_file)

        secs = time.time() - start
        print("[", "%02.3f" % secs, " seconds]", sep="")


# TODO must double-check/fix the normalisation.
def match_multiple(scanner, motifs, sequence, genomic_region, unique_threshold=None, normalize_bitscore=False, output=None):
    """
    Efficient matching of a set of motifs against a single DNA sequence.

    Keyword arguments:
    scanner -- a MOODS.scan.Scanner instance with all properties already set.
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

    results = scanner.scan(sequence)

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
