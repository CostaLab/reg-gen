###################################################################################################
# Libraries
###################################################################################################

# Python 3 compatibility



# Python
import os
import time
import sys

# Internal
from ..Util import ErrorHandler, GenomeData, npath
from ..ExperimentalMatrix import ExperimentalMatrix
from ..GeneSet import GeneSet
from ..GenomicRegionSet import GenomicRegionSet
from ..GenomicRegion import GenomicRegion
from ..AnnotationSet import AnnotationSet
from ..MotifSet import MotifSet
from .Util import bed_to_bb, parse_filter

# External
from pysam import Fastafile
from MOODS import tools, scan


###################################################################################################
# Functions
###################################################################################################

def options(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", required=True,
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
    parser.add_argument("--motif-dbs", type=str, metavar="PATH", nargs="+",
                        help="New 'motif DB' folders to use instead of the ones within "
                             "the RGTDATA folder. Each folder must contain PWM files.")
    parser.add_argument("--remove-strand-duplicates", action="store_true", default=False,
                        help="Certain motifs are 'palindromic', or more specifically they have a palindromic "
                             "consensus sequence. When this happens, the output MPBS file will have duplicates: "
                             "same chromosome and initial and final position, but opposing strand. Select this option "
                             "to only retain the 'strand duplicate' with the highest score. Duplicates due to "
                             "overlapping input regions are NOT affected by this.")
    parser.add_argument("--rmdup", action="store_true", default=False,
                        help="Remove any duplicate region from the input BED files.")
    parser.add_argument("--filter", type=str, default="", metavar="KEY_VALUE_PATTERN",
                        help="List of key-value patterns to select a subset of TFs using the metadata (MTF files), "
                             "e.g. for Mouse and Human on Selex data use: \"species:sapiens,mus;data_source:selex\". "
                             "NB: the DATABASE values must be written in full - exact matching is always performed."
                             "Valid key types are \"name\", \"gene_names\", \"family\", \"uniprot_ids\", "
                             "\"data_source\", \"tax_group\", \"species\", \"database\", \"name_file\" "
                             "and \"gene_names_file\"")
    parser.add_argument("--filter-type", choices=("inexact", "exact", "regex"), default="inexact",
                        help="Only useful together with the --filter argument."
                             "Exact will only match perfect matching of the value for each key. "
                             "Inexact will match in case the value pattern is contained within the motif. "
                             "Regex allows for a more complex pattern use.")

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

    filter_values = parse_filter(args.filter)

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if args.output_location:
        output_location = args.output_location
    else:
        output_location = npath(matching_folder_name)
    print(">> output location:", output_location)
    print()

    # Default genomic data
    genome_data = GenomeData(args.organism)

    print(">> genome:", genome_data.organism)
    print(">> pseudocounts:", args.pseudocounts)
    print(">> fpr threshold:", args.fpr)
    print()

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
            print()

        except Exception:
            err.throw_error("MM_WRONG_EXPMAT")
    elif args.input_files:
        print(">> loading input files..")
        # get input files, if available
        for input_filename in args.input_files:
            name, _ = os.path.splitext(os.path.basename(input_filename))

            regions = GenomicRegionSet(name)
            regions.read(npath(input_filename))

            genomic_regions_dict[name] = regions

            print(">>> ", name, ", ", len(regions), " regions", sep="")
        print()

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
        print(">> creating target promoter file..")
        annotation = AnnotationSet(args.organism, alias_source=args.organism,
                                   protein_coding=True, known_only=True)

        target_genes = GeneSet("target_genes")
        target_genes.read(args.target_genes_filename)

        # TODO: what do we do with unmapped genes? maybe just print them out
        target_regions = annotation.get_promoters(gene_set=target_genes, promoter_length=args.promoter_length)
        target_regions.name = "target_regions"
        target_regions.sort()
        output_file_name = npath(os.path.join(output_location, target_regions.name + ".bed"))
        target_regions.write(output_file_name)

        genomic_regions_dict[target_regions.name] = target_regions

        print(">>> ", len(target_regions), " regions", sep="")
        print()

    # we make a background in case it's requested, but also in case a list of target genes has not been
    # provided
    if args.promoter_make_background or (args.promoters_only and not args.target_genes_filename):
        print(">> creating background promoter file..")
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

        print(">>> ", len(background_regions), " regions", sep="")
        print()

    if not genomic_regions_dict:
        err.throw_error("DEFAULT_ERROR", add_msg="You must either specify an experimental matrix, or at least a "
                                                 "valid input file, or one of the 'promoter test' options.")

    max_region_len = 0
    max_region = None
    regions_to_match = []

    # Iterating on experimental matrix objects
    for k in list(genomic_regions_dict.keys()):

        curr_genomic_region = genomic_regions_dict[k]

        # If the object is a GenomicRegionSet
        if isinstance(curr_genomic_region, GenomicRegionSet):

            if args.rmdup:
                # remove duplicates and sort regions
                curr_genomic_region.remove_duplicates(sort=True)
            else:
                # sort regions
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
        print(">> creating random regions file..")

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

        print(">>> ", len(rand_region), " regions", sep="")
        print()

    ###################################################################################################
    # Creating PWMs
    ###################################################################################################

    # Load motif_set (includes MotifData object), is either customized or default
    if args.motif_dbs:
        print(">> loading external motif databases..")
        # args.motif_dbs is a list of paths to pwm files
        motif_set = MotifSet(preload_motifs=args.motif_dbs, motif_dbs=True)

        # filter for dbs only if --motif_dbs is not set
        if 'database' in filter_values:
            del filter_values['database']
    else:
        print(">> loading motif databases..")
        if 'database' in filter_values:
            motif_set = MotifSet(preload_motifs=filter_values['database'])
        else:
            motif_set = MotifSet(preload_motifs="default")

    for db in motif_set.motif_data.get_repositories_list():
        print(">>>", db)
    print()

    # applying filtering pattern, taking a subset of the motif set
    if args.filter:
        print(">> applying motif filter..")
        motif_set = motif_set.filter(filter_values, search=args.filter_type)

    motif_list = motif_set.get_motif_list(args.pseudocounts, args.fpr)
    print(">> motifs loaded:", len(motif_list))
    print()

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
            thresholds.append(0.0)
        else:
            thresholds.append(motif.threshold)
            thresholds.append(motif.threshold)

        pssm_list.append(motif.pssm)
        pssm_list.append(motif.pssm_rc)

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

    # Iterating on list of genomic region sets
    for grs in regions_to_match:

        start = time.time()
        print(">> matching [", grs.name, "], ", len(grs), " regions... ", sep="", end='')
        sys.stdout.flush()

        # Initializing output bed file
        output_bed_file = os.path.join(output_location, grs.name + "_mpbs.bed")

        # must remove it because we append the MPBS
        if os.path.isfile(output_bed_file):
            os.remove(output_bed_file)

        # Iterating on genomic region set
        for genomic_region in grs:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            grs_tmp = match_multiple(scanner, motif_list, sequence, genomic_region)

            # post-processing: if required, remove duplicate regions on opposing strands (keep highest score)
            if len(grs_tmp) > 1 and args.remove_strand_duplicates:
                grs_tmp.sort()
                seqs = grs_tmp.sequences
                seqs_new = []
                cur_pos = 0
                end_pos = len(seqs) - 1
                while cur_pos < end_pos:
                    gr = seqs[cur_pos]

                    new_pos = cur_pos + 1
                    while new_pos < end_pos:
                        gr2 = seqs[new_pos]

                        # if this sequence is unrelated, we move on
                        if gr.name != gr2.name or gr.chrom != gr2.chrom or gr.initial != gr2.initial or gr.final != gr2.final or gr.orientation == gr2.orientation:
                            break

                        if float(gr.data) < float(gr2.data):
                            gr = gr2

                        new_pos = new_pos + 1

                    # adding the currently-selected genomic region
                    seqs_new.append(gr)

                    # at the next loop, we start from the next right-handed sequences
                    cur_pos = new_pos

                # edge case: the last element was not considered
                # (when it is, cur_pos == end_pos+1)
                if cur_pos == end_pos:
                    seqs_new.append(seqs[cur_pos])

                grs_tmp.sequences = seqs_new

            grs_tmp.write(output_bed_file, mode="a")

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
        # every second result will refer to the reverse complement of the previous
        # motif. We need to handle this appropriately.
        rc = i % 2 == 1

        motif = motifs[i//2]

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

            p1 = pos_start + position
            strand = "-" if rc else "+"

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
