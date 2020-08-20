###################################################################################################
# Libraries
###################################################################################################

# Python
import base64
import os
import sys
from glob import glob
from shutil import copy
import time

# Internal
from ..Util import ErrorHandler, GenomeData, ImageData, Html, npath
from ..ExperimentalMatrix import ExperimentalMatrix
from ..GeneSet import GeneSet
from ..GenomicRegionSet import GenomicRegionSet
from ..GenomicRegion import GenomicRegion
from ..MotifSet import MotifSet
from .Statistics import multiple_test_correction, get_fisher_dict
from .Util import *

# External
from fisher import pvalue


###################################################################################################
# Functions
###################################################################################################

def options(parser):
    parser.add_argument("--organism", type=str, metavar="STRING", required=True,
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19 or hg38.")
    parser.add_argument("--matching-location", type=str, metavar="PATH",
                        help="Directory where the matching output containing the MPBS files resides. "
                             "Defaults to 'match' in the current directory.")
    parser.add_argument("--use-only-motifs", dest="selected_motifs_filename", type=str, metavar="PATH",
                        help="Only use the motifs contained within this file (one for each line).")
    parser.add_argument("--input-matrix", type=str, metavar="PATH",
                        help="If an experimental matrix is provided, the input arguments will be ignored.")
    parser.add_argument("--multiple-test-alpha", type=float, metavar="FLOAT", default=0.05,
                        help="Alpha value for multiple test.")
    parser.add_argument("--motif-dbs", type=str, metavar="PATH", nargs="+",
                        help="New 'motif DB' folders to use instead of the ones within "
                             "the RGTDATA folder. Each folder must contain PWM files.")
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

    group = parser.add_argument_group("Promoter-regions enrichment",
                                      "Used both for gene set via experimental matrix (see documentation), "
                                      "and for reporting the gene names associated to each motif.")
    group.add_argument("--promoter-length", type=int, metavar="INT", default=1000,
                       help="Length of the promoter region (in bp) to be extracted from each gene.")
    group.add_argument("--maximum-association-length", type=int, metavar="INT", default=50000,
                       help="Maximum distance between a coordinate and a gene (in bp) in order for the former to "
                            "be considered associated with the latter.")
    group.add_argument("--exclude-target-genes", action="store_true", help="If set the specified target genes are"
                                                                           "excluded from background file")

    # Output Options
    group = parser.add_argument_group("Output",
                                      "Where to put the output files and how to post-process them.")
    group.add_argument("--output-location", type=str, metavar="PATH",
                       help="Path where the output MPBS files will be written. Defaults to 'enrichment' in the "
                            "current directory.")
    group.add_argument("--print-thresh", type=float, metavar="FLOAT", default=0.05,
                       help="Only MPBSs whose factor's enrichment corrected p-value are less than equal "
                            "this option are printed. Use 1.0 to print all MPBSs.")
    group.add_argument("--bigbed", action="store_true", default=False,
                       help="If this option is used, all bed files will be written as bigbed.")

    # Logo mutually exclusive:
    group = parser.add_mutually_exclusive_group(required=False)
    group.add_argument("--logo-copy", action="store_true", default=False,
                       help="The logos are copied to a local directory. The HTML report will contain relative "
                            "paths to this directory.")
    group.add_argument("--logo-embed", action="store_true", default=False,
                       help="The logos are embedded directly into the HTML report.")

    parser.add_argument('background_file', metavar='background.bed', type=str,
                        help='BED file containing background regions.')
    parser.add_argument('input_files', metavar='input.bed', type=str, nargs='*',
                        help='BED files to be enriched against the background.')


def main(args):
    """
    Performs motif enrichment.
    """

    # Initializing Error Handler
    err = ErrorHandler()

    background_filename = args.background_file
    input_files = args.input_files  # empty list if no input files

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
                   "BACKGROUND FREQUENCY"]
    html_type_list = "sissssssssl"
    logo_width = 200
    html_col_size = [300, logo_width, 100, 100, 50, 50, 50, 50, 100, 100, 50]

    filter_values = parse_filter(args.filter)

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if args.output_location:
        output_location = args.output_location
    else:
        output_location = os.path.join(os.getcwd(), enrichment_folder_name)
    print(">> output location:", output_location)

    # Matching folder
    if args.matching_location:
        match_location = args.matching_location
    else:
        match_location = os.path.join(os.getcwd(), matching_folder_name)
    print(">> match location:", match_location)
    print()

    # Default genomic data
    genome_data = GenomeData(args.organism)
    print(">> genome:", genome_data.organism)
    print()

    # the matching directory must exist, for obvious reasons
    if not os.path.isdir(match_location):
        err.throw_error("ME_MATCH_NOTFOUND")

    # Background file must exist
    print(">> loading background file..")
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
    background_name = os.path.basename(path)

    # first we check at matching folder location
    background_mpbs_filename = os.path.join(match_location, background_name + "_mpbs" + ext)

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

    background = GenomicRegionSet("background")
    background.read(background_filename)
    print(">>> ", background_name, ", ", len(background), " regions", sep="")
    print()

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

    motif_names = list(motif_set.motifs_map.keys())
    print(">> motifs loaded:", len(motif_names))
    print()

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
    if args.input_matrix:
        try:
            exp_matrix = ExperimentalMatrix()
            exp_matrix.read(args.input_matrix)

            # if the matrix is present, the (empty) dictionary is overwritten
            genomic_regions_dict = exp_matrix.objectsDict

            # Reading dictionary grouped by fields (only for gene association)
            try:
                exp_matrix_fields_dict = exp_matrix.fieldsDict[gene_column_name]
                flag_gene = True
            except KeyError:
                flag_gene = False

            del exp_matrix

            print(">> experimental matrix loaded")
            print()

        except Exception:
            err.throw_error("MM_WRONG_EXPMAT")
    elif input_files:
        print(">> loading input files..")
        # get input files, if available
        for input_filename in input_files:
            name, _ = os.path.splitext(os.path.basename(input_filename))

            regions = GenomicRegionSet(name)
            regions.read(os.path.abspath(input_filename))

            genomic_regions_dict[name] = regions

            print(">>> ", name, ", ", len(regions), " regions", sep="")
            print()

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
        for g in list(exp_matrix_fields_dict.keys()):

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
        for k in list(genomic_regions_dict.keys()):

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
    # Background Statistics
    ###################################################################################################

    # if the output folder doesn't exist, we create it
    # (we do this here to not create the output folder when the program fails early)
    try:
        if not os.path.isdir(output_location):
            os.makedirs(output_location)
    except Exception:
        err.throw_error("ME_OUT_FOLDER_CREATION")

    start = time.time()
    print(">> collecting background statistics...", sep="", end="")
    sys.stdout.flush()
    background_mpbs = GenomicRegionSet("background_mpbs")
    background_mpbs.read(background_mpbs_filename)

    # Evaluating background statistics
    bg_c_dict, bg_d_dict, _, _ = get_fisher_dict(motif_names, background, background_mpbs)

    # removing temporary BED files if the originals were BBs
    if is_bb(background_original_filename):
        os.remove(background_filename)
    if is_bb(background_mpbs_original_filename):
        os.remove(background_mpbs_filename)

    # # scheduling region sets for garbage collection
    # del background.sequences[:]
    # del background
    # del background_mpbs.sequences[:]
    # del background_mpbs

    secs = time.time() - start
    print("[", "%02.3f" % secs, " seconds]", sep="")

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

            start = time.time()
            print(">>> enriching [", grs.name, "], ", len(grs), " regions...", end="", sep="")
            sys.stdout.flush()

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
            curr_mpbs.read(curr_mpbs_bed_name)
            curr_mpbs.sort()

            ###################################################################################################
            # Gene Evidence Statistics
            ###################################################################################################

            if curr_input.gene_set:

                # Performing association of input region with gene_set
                grs = grs.gene_association(organism=args.organism, gene_set=curr_input.gene_set,
                                           promoter_length=args.promoter_length,
                                           thresh_dist=args.maximum_association_length)

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
                if args.bigbed:
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
                                                                                 alpha=args.multiple_test_alpha,
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
                for m, ev_mpbs_grs in list(ev_mpbs_dict.items()):
                    for region in ev_mpbs_grs:
                        if corr_pvalue_dict[m] <= args.print_thresh:
                            ev_mpbs_grs_filtered.add(region)
                del ev_mpbs_dict

                nev_mpbs_grs_filtered = GenomicRegionSet("nev_mpbs_filtered")
                for m, nev_mpbs_grs in list(nev_mpbs_dict.items()):
                    for region in nev_mpbs_grs:
                        if corr_pvalue_dict[m] <= args.print_thresh:
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
                if args.bigbed:
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

                # we copy the logo images locally
                if args.logo_copy:
                    logo_dir_path = npath(os.path.join(output_location, "logos"))
                    try:
                        os.stat(logo_dir_path)
                    except Exception:
                        os.mkdir(logo_dir_path)

                # Printing statistics html - Creating data table
                data_table = []
                for r in result_list:
                    curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                    for rep in motif_set.motif_data.get_logo_list():
                        logo_file_name = npath(os.path.join(rep, r.name + ".png"))

                        if os.path.isfile(logo_file_name):
                            if args.logo_copy:
                                copy(logo_file_name, npath(os.path.join(logo_dir_path, r.name + ".png")))

                                # use relative paths in the html
                                # FIXME can we do it in a better way? (inside the Html class)
                                logo_file_name = os.path.join("..", "logos", r.name + ".png")
                            elif args.logo_embed:
                                # encoding image with Base64 and adding HTML special URI to embed it
                                data_uri = base64.b64encode(open(logo_file_name, 'rb').read()).decode('utf-8')
                                logo_file_name = "data:image/png;base64,{0}".format(data_uri)

                            curr_motif_tuple = [logo_file_name, logo_width]
                            break
                    data_table.append(
                        [r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                         str(r.c), str(r.d), str(r.percent), str(r.back_percent)])

                # Printing statistics html - Writing to HTML
                output_file_name_html = os.path.join(curr_output_folder_name, output_stat_genetest + ".html")
                fig_path = npath(os.path.join(output_location, "fig"))
                html = Html("Motif Enrichment Analysis", genetest_link_dict, fig_dir=fig_path)
                html.add_heading("Results for <b>" + original_name + "</b> "
                                                                     "region <b>Gene Test*</b> using genes from <b>"
                                 + curr_input.gene_set.name + "</b>",
                                 align="center", bold=False)
                html.add_heading("* This gene test considered regions associated to the given "
                                 "gene list against regions not associated to the gene list",
                                 align="center", bold=False, size=3)
                html.add_zebra_table(html_header, html_col_size, html_type_list, data_table, align="center")
                html.write(npath(output_file_name_html))

            else:
                # Association still needs to be done with all genes in order to print gene list
                grs = grs.gene_association(organism=args.organism, gene_set=None,
                                           promoter_length=args.promoter_length,
                                           thresh_dist=args.maximum_association_length)

                # Calculating statistics
                a_dict, b_dict, ev_genes_dict, ev_mpbs_dict = get_fisher_dict(motif_names, grs, curr_mpbs,
                                                                              gene_set=True, mpbs_set=True)

                if args.exclude_target_genes:
                    # subtract target_genes
                    background_tmp = background.subtract(grs, exact=True)

                    # fisher dict for new (smaller) background
                    bg_c_dict, bg_d_dict, _, _ = get_fisher_dict(motif_names, background_tmp, background_mpbs)

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
                                                                             alpha=args.multiple_test_alpha,
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

                for m, ev_mpbs_grs in list(ev_mpbs_dict.items()):
                    for region in ev_mpbs_grs:
                        if corr_pvalue_dict[m] <= args.print_thresh:
                            ev_mpbs_grs_filtered.add(region)
                del ev_mpbs_dict

                output_mpbs_filtered_ev_bed = os.path.join(curr_output_folder_name, output_mpbs_filtered_ev + ".bed")
                output_mpbs_filtered_ev_bed = npath(output_mpbs_filtered_ev_bed)

                # sorting and saving to BED
                ev_mpbs_grs_filtered.sort()
                write_bed_color(ev_mpbs_grs_filtered, output_mpbs_filtered_ev_bed, ev_color)
                del ev_mpbs_grs_filtered

                # Converting ev mpbs to bigbed
                if args.bigbed:
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

            # we copy the logo images locally
            if args.logo_copy:
                logo_dir_path = npath(os.path.join(output_location, "logos"))
                try:
                    os.stat(logo_dir_path)
                except Exception:
                    os.mkdir(logo_dir_path)

            # Printing statistics html - Creating data table
            data_table = []
            for r in result_list:
                curr_motif_tuple = [image_data.get_default_motif_logo(), logo_width]
                for rep in motif_set.motif_data.get_logo_list():
                    logo_file_name = npath(os.path.join(rep, r.name + ".png"))

                    if os.path.isfile(logo_file_name):
                        if args.logo_copy:
                            copy(logo_file_name, npath(os.path.join(logo_dir_path, r.name + ".png")))

                            # use relative paths in the html
                            # FIXME can we do it in a better way? (inside the Html class)
                            logo_file_name = os.path.join("..", "logos", r.name + ".png")
                        elif args.logo_embed:
                            # encoding image with Base64 and adding HTML special URI to embed it
                            data_uri = base64.b64encode(open(logo_file_name, 'rb').read()).decode('utf-8')
                            logo_file_name = "data:image/png;base64,{0}".format(data_uri)

                        curr_motif_tuple = [logo_file_name, logo_width]
                        break
                data_table.append([r.name, curr_motif_tuple, str(r.p_value), str(r.corr_p_value), str(r.a), str(r.b),
                                   str(r.c), str(r.d), str(r.percent), str(r.back_percent)])

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

            secs = time.time() - start
            print("[", "%02.3f" % secs, " seconds]", sep="")
