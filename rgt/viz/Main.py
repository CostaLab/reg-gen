# Python Libraries

import os
import sys
import time
import getpass
import argparse
import datetime
import matplotlib

matplotlib.use('Agg')

from .boxplot import Boxplot
from .lineplot import Lineplot
from .jaccard_test import Jaccard
from .projection_test import Projection
from .intersection_test import Intersect
from .bed_profile import BedProfile
from .shared_function import check_dir, print2, output_parameters, \
    copy_em, list_all_index, output
from .plotTools import Venn
from .. import __version__

current_dir = os.getcwd()
"""
Statistical analysis methods and plotting tools for ExperimentalMatrix

Author: Joseph C.C. Kuo
"""


def main():
    ###############################################################################
    ##### PARAMETERS ##############################################################
    ###############################################################################

    # Some general help descriptions
    ######### Some general plotting arguments descriptions ###############
    helpinput = 'The file name of the input Experimental Matrix file. Recommended to add more columns for more information for ploting. For example, cell type or factors. (default: %(default)s)'
    helpoutput = 'The directory name for the output files. For example, project name. (default: %(default)s)'
    helptitle = 'The title shown on the top of the plot and also the folder name. (default: %(default)s)'
    helpgroup = "Group the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None. (default: %(default)s)"
    helpgroupbb = "Group the data by any optional column (for example, 'cell') of experimental matrix, or None. (default: %(default)s)"
    helpsort = "Sort the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None. (default: %(default)s)"
    helpcolor = "Color the data by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None. (default: %(default)s)"
    helpcolorbb = "Color the data by any optional column (for example, 'cell') of experimental matrix, or None. (default: %(default)s)"
    help_define_color = 'Define the specific colors with the given column "color" in experimental matrix. The color should be in the format of matplotlib.colors. For example, "r" for red, "b" for blue, or "(100, 35, 138)" for RGB. (default: %(default)s)'
    helpreference = 'The file name of the reference Experimental Matrix. Multiple references are acceptable. (default: %(default)s)'
    helpquery = 'The file name of the query Experimental Matrix. Multiple queries are acceptable. (default: %(default)s)'
    helpcol = "Group the data in columns by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None. (default: %(default)s)"
    helprow = "Group the data in rows by reads(needs 'factor' column), regions(needs 'factor' column), another name of column (for example, 'cell')in the header of experimental matrix, or None. (default: %(default)s)"
    helpmp = "Define the number of cores for parallel computation. (default: %(default)s)"

    version_message = "viz - Regulatory Analysis Toolbox (RGT). Version: " + str(__version__)
    parser = argparse.ArgumentParser(description='Provides various Statistical analysis methods and plotting tools for ExperimentalMatrix.\
    \nAuthor: Joseph C.C. Kuo, Ivan Gesteira Costa Filho', formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     add_help=True)
    parser.add_argument('--version', action='version', version=version_message)

    subparsers = parser.add_subparsers(help='sub-command help', dest='mode')

    ################### BED profile ##########################################
    parser_bedprofile = subparsers.add_parser('bed_profile',
                                              help='BED profile analyzes the given BED file(s) by their length, distribution and composition of the sequences.')
    parser_bedprofile.add_argument('-i', metavar='  ',
                                   help="Input experimental matrix or Input BED file or Input directory which contains BED files")
    parser_bedprofile.add_argument('-o', metavar='  ', help=helpoutput)
    parser_bedprofile.add_argument('-t', metavar='  ', default='bed_profile', help=helptitle)
    parser_bedprofile.add_argument('-organism', metavar='  ', default=None,
                                   help='Define the organism. (default: %(default)s)')
    parser_bedprofile.add_argument('-biotype', metavar='  ', default=False,
                                   help='Define the directory for biotype BED files.')
    parser_bedprofile.add_argument('-repeats', metavar='  ', default=False,
                                   help='Define the directory for repeats BED files.')
    parser_bedprofile.add_argument('-genposi', metavar='  ', default=False,
                                   help='Define the directory for the generic position BED files. (exons, introns, and intergenic regions)')
    parser_bedprofile.add_argument('-labels', metavar='  ', default=None, help='Define the labels for more BED sets')
    parser_bedprofile.add_argument('-sources', metavar='  ', default=None,
                                   help='Define the directories for more BED sets corresponding to the labels')
    parser_bedprofile.add_argument('-strand', metavar='  ', default=None,
                                   help='Define whether to perform strand-specific comparison for each reference corresponding to the labels (T or F)')
    parser_bedprofile.add_argument('-other', metavar='  ', default=None,
                                   help='Define whether to count "else" for each reference corresponding to the labels (T or F)')
    parser_bedprofile.add_argument('-background', metavar='  ', default=None,
                                   help='Add the background to the first row of the figures (T or F)')
    parser_bedprofile.add_argument('-coverage', action="store_true", default=False,
                                   help='Calculate the overlapping region by coverage in bp instead of simple counting')
    parser_bedprofile.add_argument('-test', action="store_true", default=False,
                                   help='test script')

    ################### Projection test ##########################################
    parser_projection = subparsers.add_parser('projection',
                                              help='Projection test evaluates the association level by comparing to the random binomial model.')
    parser_projection.add_argument('-r', metavar='  ', help=helpreference)
    parser_projection.add_argument('-q', metavar='  ', help=helpquery)
    parser_projection.add_argument('-o', metavar='  ', help=helpoutput)
    parser_projection.add_argument('-t', metavar='  ', default='projection_test', help=helptitle)
    parser_projection.add_argument('-g', metavar='  ', default=None, help=helpgroupbb)
    parser_projection.add_argument('-c', metavar='  ', default="regions", help=helpcolorbb)
    parser_projection.add_argument('-bg', metavar='  ', type=str, default=None,
                                   help="Define a BED file as background. If not defined, the background is whole genome according to the given organism. (default: %(default)s)")
    parser_projection.add_argument('-union', action="store_true",
                                   help='Take the union of references as background for binominal test. (default: %(default)s)')
    parser_projection.add_argument('-organism', metavar='  ', default='hg19',
                                   help='Define the organism. (default: %(default)s)')
    parser_projection.add_argument('-log', action="store_true",
                                   help='Set y axis of the plot in log scale. (default: %(default)s)')
    parser_projection.add_argument('-color', action="store_true", help=help_define_color)
    parser_projection.add_argument('-show', action="store_true",
                                   help='Show the figure in the screen. (default: %(default)s)')
    parser_projection.add_argument('-table', action="store_true",
                                   help='Store the tables of the figure in text format. (default: %(default)s)')
    parser_projection.add_argument('-bed', action="store_true", default=False,
                                   help='Output BED files for the regions of query which overlap the reference. (default: %(default)s)')
    parser_projection.add_argument('-pw', metavar='  ', type=int, default=5,
                                   help='Define the width of single panel. (default: %(default)s)')
    parser_projection.add_argument('-ph', metavar='  ', type=int, default=3,
                                   help='Define the height of single panel. (default: %(default)s)')
    parser_projection.add_argument('-cfp', metavar='  ', type=float, default=0,
                                   help='Define the cutoff of the proportion. (default: %(default)s)')
    parser_projection.add_argument('-load', action="store_false", default=True,
                                   help='Load the BED files later during processing, which saves memory usage when dealing with large number of BED files.')

    ################### Intersect Test ##########################################
    parser_intersect = subparsers.add_parser('intersect',
                                             help='Intersection test provides various modes of intersection to test the association between references and queries.')
    parser_intersect.add_argument('-r', metavar='  ', help=helpreference)
    parser_intersect.add_argument('-q', metavar='  ', help=helpquery)
    parser_intersect.add_argument('-o', help=helpoutput)
    parser_intersect.add_argument('-t', metavar='  ', default='intersection_test', help=helptitle)
    parser_intersect.add_argument('-g', metavar='  ', default=None, help=helpgroupbb)
    parser_intersect.add_argument('-c', metavar='  ', default="regions", help=helpcolorbb)
    parser_intersect.add_argument('-organism', metavar='  ', default='hg19',
                                  help='Define the organism. (default: %(default)s)')
    parser_intersect.add_argument('-bg', metavar='  ',
                                  help="Define a BED file as background. If not defined, the background is whole genome according to the given organism. (default: %(default)s)")
    parser_intersect.add_argument('-m', metavar='  ', default="count", choices=['count', 'bp'],
                                  help="Define the mode of calculating intersection. 'count' outputs the number of overlapped regions.'bp' outputs the coverage(basepair) of intersection. (default: %(default)s)")
    parser_intersect.add_argument('-tc', metavar='  ', type=int, default=False,
                                  help="Define the threshold(in percentage) of reference length for intersection counting. For example, '20' means that the query which overlaps more than 20%% of reference is counted as intersection. (default: %(default)s)")
    parser_intersect.add_argument('-ex', metavar='  ', type=int, default=0,
                                  help="Define the extension(in bp) of reference length for intersection counting. For example, '20' means that each region of reference is extended by 20 bp in order to include proximal queries. (default: %(default)s)")
    parser_intersect.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
    parser_intersect.add_argument('-color', action="store_true", help=help_define_color)
    parser_intersect.add_argument('-show', action="store_true",
                                  help='Show the figure in the screen. (default: %(default)s)')
    parser_intersect.add_argument('-stest', metavar='  ', type=int, default=0,
                                  help='Define the repetition time of random subregion test between reference and query. (default: %(default)s)')
    parser_intersect.add_argument('-mp', metavar='  ', default=4, type=int, help=helpmp)
    parser_intersect.add_argument('-pw', metavar='  ', type=int, default=3,
                                  help='Define the width of single panel. (default: %(default)s)')
    parser_intersect.add_argument('-ph', metavar='  ', type=int, default=3,
                                  help='Define the height of single panel. (default: %(default)s)')

    ################### Jaccard test ##########################################

    parser_jaccard = subparsers.add_parser('jaccard',
                                           help='Jaccard test evaluates the association level by comparing with jaccard index from repeating randomization.')

    parser_jaccard.add_argument('-o', help=helpoutput)
    parser_jaccard.add_argument('-r', metavar='  ', help=helpreference)
    parser_jaccard.add_argument('-q', metavar='  ', help=helpquery)
    parser_jaccard.add_argument('-t', metavar='  ', default='jaccard_test', help=helptitle)
    parser_jaccard.add_argument('-rt', metavar='  ', type=int, default=500,
                                help='Define how many times to run the randomization. (default: %(default)s)')
    parser_jaccard.add_argument('-g', default=None, help=helpgroupbb)
    parser_jaccard.add_argument('-c', default="regions", help=helpcolorbb)
    parser_jaccard.add_argument('-organism', default='hg19', help='Define the organism. (default: %(default)s)')
    parser_jaccard.add_argument('-nlog', action="store_false",
                                help='Set y axis of the plot not in log scale. (default: %(default)s)')
    parser_jaccard.add_argument('-color', action="store_true", help=help_define_color)
    parser_jaccard.add_argument('-show', action="store_true",
                                help='Show the figure in the screen. (default: %(default)s)')
    parser_jaccard.add_argument('-table', action="store_true",
                                help='Store the tables of the figure in text format. (default: %(default)s)')
    parser_jaccard.add_argument('-pw', metavar='  ', type=int, default=3,
                                help='Define the width of single panel. (default: %(default)s)')
    parser_jaccard.add_argument('-ph', metavar='  ', type=int, default=3,
                                help='Define the height of single panel. (default: %(default)s)')

    ################### Combinatorial Test ##########################################
    parser_combinatorial = subparsers.add_parser('combinatorial',
                                                 help='Combinatorial test compare all combinatorial possibilities from reference to test the association between references and queries.')

    parser_combinatorial.add_argument('-o', help=helpoutput)
    parser_combinatorial.add_argument('-r', metavar='  ', help=helpreference)
    parser_combinatorial.add_argument('-q', metavar='  ', help=helpquery)
    parser_combinatorial.add_argument('-t', metavar='  ', default='combinatorial_test', help=helptitle)
    parser_combinatorial.add_argument('-g', default=None, help=helpgroupbb)
    parser_combinatorial.add_argument('-c', default="regions", help=helpcolorbb)
    parser_combinatorial.add_argument('-organism', default='hg19', help='Define the organism. (default: %(default)s)')
    parser_combinatorial.add_argument('-bg',
                                      help="Define a BED file as background. If not defined, the background is whole genome according to the given organism. (default: %(default)s)")
    parser_combinatorial.add_argument('-m', default="count", choices=['count', 'bp'],
                                      help="Define the mode of calculating intersection. 'count' outputs the number of overlapped regions.'bp' outputs the coverage(basepair) of intersection. (default: %(default)s)")
    parser_combinatorial.add_argument('-tc', type=int, default=False,
                                      help="Define the threshold(in percentage) of reference length for intersection counting. For example, '20' means that the query which overlaps more than 20%% of reference is counted as intersection. (default: %(default)s)")
    parser_combinatorial.add_argument('-ex', type=int, default=0,
                                      help="Define the extension(in percentage) of reference length for intersection counting. For example, '20' means that each region of reference is extended by 20%% in order to include proximal queries. (default: %(default)s)")
    parser_combinatorial.add_argument('-log', action="store_true",
                                      help='Set y axis of the plot in log scale. (default: %(default)s)')
    parser_combinatorial.add_argument('-color', action="store_true", help=help_define_color)
    parser_combinatorial.add_argument('-venn', action="store_true",
                                      help='Show the Venn diagram of the combinatorials of references. (default: %(default)s)')
    parser_combinatorial.add_argument('-show', action="store_true",
                                      help='Show the figure in the screen. (default: %(default)s)')
    parser_combinatorial.add_argument('-stest', type=int, default=0,
                                      help='Define the repetition time of random subregion test between reference and query. (default: %(default)s)')
    parser_combinatorial.add_argument('-pw', metavar='  ', type=int, default=3,
                                      help='Define the width of single panel. (default: %(default)s)')
    parser_combinatorial.add_argument('-ph', metavar='  ', type=int, default=3,
                                      help='Define the height of single panel. (default: %(default)s)')

    ################### Boxplot ##########################################

    parser_boxplot = subparsers.add_parser('boxplot',
                                           help='Boxplot based on the BAM and BED files for gene association analysis.')
    parser_boxplot.add_argument('input', help=helpinput)
    parser_boxplot.add_argument('-o', metavar='  ', help=helpoutput)
    parser_boxplot.add_argument('-t', metavar='  ', default='boxplot', help=helptitle)
    parser_boxplot.add_argument('-g', metavar='  ', default='reads', help=helpgroup)
    parser_boxplot.add_argument('-c', metavar='  ', default='regions', help=helpcolor)
    parser_boxplot.add_argument('-s', metavar='  ', default='None', help=helpsort)
    parser_boxplot.add_argument('-scol', action="store_true", help="Share y axis among columns. (default: %(default)s)")
    parser_boxplot.add_argument('-nlog', action="store_false",
                                help='Set y axis of the plot not in log scale. (default: %(default)s)')
    parser_boxplot.add_argument('-color', action="store_true", help=help_define_color)
    parser_boxplot.add_argument('-pw', metavar='  ', type=int, default=3,
                                help='Define the width of single panel. (default: %(default)s)')
    parser_boxplot.add_argument('-ph', metavar='  ', type=int, default=3,
                                help='Define the height of single panel. (default: %(default)s)')
    parser_boxplot.add_argument('-nqn', action="store_true",
                                help='No quantile normalization in calculation. (default: %(default)s)')
    parser_boxplot.add_argument('-df', action="store_true",
                                help="Show the difference of the two signals which share the same labels.The result is the subtraction of the first to the second. (default: %(default)s)")
    parser_boxplot.add_argument('-ylim', metavar='  ', type=int, default=None,
                                help="Define the limit of y axis. (default: %(default)s)")
    parser_boxplot.add_argument('-p', metavar='  ', type=float, default=0.05,
                                help='Define the significance level for multiple test.  (default: %(default)s)')
    parser_boxplot.add_argument('-show', action="store_true",
                                help='Show the figure in the screen. (default: %(default)s)')
    parser_boxplot.add_argument('-table', action="store_true",
                                help='Store the tables of the figure in text format. (default: %(default)s)')

    ################### Lineplot ##########################################
    parser_lineplot = subparsers.add_parser('lineplot', help='Generate lineplot with various modes.')

    choice_center = ['midpoint', 'bothends', 'upstream', 'downstream']
    # Be consist as the arguments of GenomicRegionSet.relocate_regions

    parser_lineplot.add_argument('input', help=helpinput)
    parser_lineplot.add_argument('-o', help=helpoutput)
    parser_lineplot.add_argument('-ga', action="store_true",
                                 help="Use genetic annotation data as input regions (e.g. TSS, TTS, exons and introns) instead of the BED files in the input matrix.")
    parser_lineplot.add_argument('-t', metavar='  ', default='lineplot', help=helptitle)
    parser_lineplot.add_argument('-center', metavar='  ', choices=choice_center, default='midpoint',
                                 help='Define the center to calculate coverage on the regions. Options are: ' + ', '.join(
                                     choice_center) + '. (default: %(default)s) The bothend mode will flap the right end region for calculation.')
    parser_lineplot.add_argument('-g', metavar='  ', default='None', help=helpgroup)
    parser_lineplot.add_argument('-row', metavar='  ', default='None', help=helprow)
    parser_lineplot.add_argument('-col', metavar='  ', default='regions', help=helpcol)
    parser_lineplot.add_argument('-c', metavar='  ', default='reads', help=helpcolor)
    parser_lineplot.add_argument('-e', metavar='  ', type=int, default=2000,
                                 help='Define the extend length of interested region for plotting. (default: %(default)s)')
    parser_lineplot.add_argument('-rs', metavar='  ', type=int, default=200,
                                 help='Define the readsize for calculating coverage. (default: %(default)s)')
    parser_lineplot.add_argument('-ss', metavar='  ', type=int, default=50,
                                 help='Define the stepsize for calculating coverage. (default: %(default)s)')
    parser_lineplot.add_argument('-bs', metavar='  ', type=int, default=100,
                                 help='Define the binsize for calculating coverage. (default: %(default)s)')
    parser_lineplot.add_argument('-log', action="store_true",
                                 help="Take log for the value before calculating average. (default: %(default)s)")
    parser_lineplot.add_argument('-scol', action="store_true",
                                 help="Share y axis among columns. (default: %(default)s)")
    parser_lineplot.add_argument('-srow', action="store_true", help="Share y axis among rows. (default: %(default)s)")
    parser_lineplot.add_argument('-organism', metavar='  ',
                                 help='Define the organism. (default: %(default)s)')
    parser_lineplot.add_argument('-color', action="store_true", help=help_define_color)
    parser_lineplot.add_argument('-pw', metavar='  ', type=int, default=3,
                                 help='Define the width of single panel. (default: %(default)s)')
    parser_lineplot.add_argument('-ph', metavar='  ', type=int, default=3,
                                 help='Define the height of single panel. (default: %(default)s)')
    parser_lineplot.add_argument('-test', action="store_true",
                                 help="Sample only the first 10 regions in all BED files for testing. (default: %(default)s)")
    parser_lineplot.add_argument('-mp', metavar='  ', type=int, default=0,
                                 help="Perform multiprocessing for faster computation. (default: %(default)s)")
    parser_lineplot.add_argument('-df', action="store_true",
                                 help="Show the difference of the two signals which share the same labels.The result is the subtraction of the first to the second. (default: %(default)s)")
    parser_lineplot.add_argument('-dft', metavar='  ', default=None,
                                 help="Add one more tag for calculating difference. (default: %(default)s)")
    parser_lineplot.add_argument('-show', action="store_true",
                                 help='Show the figure in the screen. (default: %(default)s)')
    parser_lineplot.add_argument('-table', action="store_true",
                                 help='Store the tables of the figure in text format. (default: %(default)s)')
    parser_lineplot.add_argument('-sense', action="store_true",
                                 help='Set the plot sense-specific. (default: %(default)s)')
    parser_lineplot.add_argument('-strand', action="store_true",
                                 help='Set the plot strand-specific. (default: %(default)s)')
    parser_lineplot.add_argument('-average', action="store_true",
                                 help='Show only the average of the replicates. (default: %(default)s)')
    parser_lineplot.add_argument('-flip_negative', action="store_true", default=False,
                                 help='Flip the negative strand (default: %(default)s)')
    parser_lineplot.add_argument('-extend_outside', action="store_true", default=False,
                                 help='Extend the window outside of the given regions and compress the given region into fixed internal. (default: %(default)s)')
    parser_lineplot.add_argument('-add_region_number', action="store_true", default=False,
                                 help="Add the number of regions in the axis label. (default: %(default)s)")

    ################### Heatmap ##########################################
    parser_heatmap = subparsers.add_parser('heatmap', help='Generate heatmap with various modes.')

    choice_center = ['midpoint', 'bothends', 'upstream', 'downstream']
    # Be consist as the arguments of GenomicRegionSet.relocate_regions

    parser_heatmap.add_argument('input', help=helpinput)
    parser_heatmap.add_argument('-o', metavar='  ', help=helpoutput)
    parser_heatmap.add_argument('-ga', action="store_true",
                                help="Use genetic annotation data as input regions (e.g. TSS, TTS, exons and introns) instead of the BED files in the input matrix. (default: %(default)s)")
    parser_heatmap.add_argument('-t', metavar='  ', default='heatmap', help=helptitle)
    parser_heatmap.add_argument('-center', metavar='  ', choices=choice_center, default='midpoint',
                                help='Define the center to calculate coverage on the regions. Options are: ' + ', '.join(
                                    choice_center) +
                                     '.(Default:midpoint) The bothend mode will flap the right end region for calculation. (default: %(default)s)')
    parser_heatmap.add_argument('-sort', metavar='  ', type=int, default=None,
                                help='Define the way to sort the signals.' +
                                     'Default is no sorting at all, the signals arrange in the order of their position; ' +
                                     '"0" is sorting by the average ranking of all signals; ' +
                                     '"1" is sorting by the ranking of 1st column; "2" is 2nd and so on...  (default: %(default)s)')
    parser_heatmap.add_argument('-col', metavar='  ', default='regions', help=helpcol)
    parser_heatmap.add_argument('-c', metavar='  ', default='reads', help=helpcolor)
    parser_heatmap.add_argument('-row', metavar='  ', default='None', help=helprow)
    parser_heatmap.add_argument('-e', metavar='  ', type=int, default=2000,
                                help='Define the extend length of interested region for plotting. (default: %(default)s)')
    parser_heatmap.add_argument('-rs', metavar='  ', type=int, default=200,
                                help='Define the readsize for calculating coverage. (default: %(default)s)')
    parser_heatmap.add_argument('-ss', metavar='  ', type=int, default=50,
                                help='Define the stepsize for calculating coverage. (default: %(default)s)')
    parser_heatmap.add_argument('-bs', metavar='  ', type=int, default=100,
                                help='Define the binsize for calculating coverage. (default: %(default)s)')
    parser_heatmap.add_argument('-organism', metavar='  ', default='hg19',
                                help='Define the organism. (default: %(default)s)')
    parser_heatmap.add_argument('-color', action="store_true", help=help_define_color)
    parser_heatmap.add_argument('-log', action="store_true", help='Set colorbar in log scale. (default: %(default)s)')
    parser_heatmap.add_argument('-mp', action="store_true",
                                help="Perform multiprocessing for faster computation. (default: %(default)s)")
    parser_heatmap.add_argument('-show', action="store_true",
                                help='Show the figure in the screen. (default: %(default)s)')
    parser_heatmap.add_argument('-table', action="store_true",
                                help='Store the tables of the figure in text format. (default: %(default)s)')

    ################### Venn Diagram ########################################
    parser_venn = subparsers.add_parser('venn', help='Generate Venn Diagram with peaks of gene list.')

    parser_venn.add_argument('-s1', metavar='  ', default=None,
                             help="Define the file for gene set 1 (BED or gene list)")
    parser_venn.add_argument('-s2', metavar='  ', default=None,
                             help="Define the file for gene set 2 (BED or gene list)")
    parser_venn.add_argument('-s3', metavar='  ', default=None,
                             help="Define the file for gene set 3 (BED or gene list)")
    parser_venn.add_argument('-s4', metavar='  ', default=None,
                             help="Define the file for gene set 3 (BED or gene list)")
    parser_venn.add_argument('-l1', metavar='  ', default=None, help="Define label on venn diagram for set 1")
    parser_venn.add_argument('-l2', metavar='  ', default=None, help="Define label on venn diagram for set 2")
    parser_venn.add_argument('-l3', metavar='  ', default=None, help="Define label on venn diagram for set 3")
    parser_venn.add_argument('-l4', metavar='  ', default=None, help="Define label on venn diagram for set 4")
    parser_venn.add_argument('-o', metavar='  ', help=helpoutput)
    parser_venn.add_argument('-t', metavar='  ', default='venn_diagram', help=helptitle)
    parser_venn.add_argument('-organism', metavar='  ', help='Define the organism. ')

    ################### Integration ##########################################
    parser_integration = subparsers.add_parser('integration',
                                               help='Provides some tools to deal with experimental matrix or other purposes.')
    parser_integration.add_argument('-ihtml', action="store_true",
                                    help='Integrate all the html files within the given directory and generate index.html for all plots.')
    parser_integration.add_argument('-l2m', help='Convert a given file list in txt format into a experimental matrix.')
    parser_integration.add_argument('-o', help='Define the folder of the output file.')
    ################### Parsing the arguments ################################
    # print(sys.argv)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    elif len(sys.argv) == 2:
        if sys.argv[1] == "-h" or sys.argv[1] == "--help":
            parser.print_help()
            sys.exit(0)
        elif sys.argv[1] == "-v" or sys.argv[1] == "--version":
            print(version_message)
            sys.exit(0)
        else:
            # retrieve subparsers from parser
            subparsers_actions = [action for action in parser._actions if
                                  isinstance(action, argparse._SubParsersAction)]
            # there will probably only be one subparser_action,but better save than sorry
            for subparsers_action in subparsers_actions:
                # get all subparsers and print help
                for choice, subparser in list(subparsers_action.choices.items()):
                    if choice == sys.argv[1]:
                        print("\nYou need more arguments.")
                        print("\nSubparser '{}'".format(choice))
                        subparser.print_help()
            sys.exit(1)
    else:
        args = parser.parse_args()
        if args.mode != 'integration':
            if not args.o:
                print("** Error: Please define the output directory (-o).")
                sys.exit(1)

            t0 = time.time()
            # Normalised output path
            args.o = os.path.normpath(os.path.join(current_dir, args.o))
            check_dir(args.o)
            check_dir(os.path.join(args.o, args.t))

            # Input parameters dictionary
            parameter = ["Time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                         "User: " + getpass.getuser(),
                         "\nCommand:\n\t$ " + " ".join(sys.argv)]

        #################################################################################################
        ##### Main #####################################################################################
        #################################################################################################

        if args.mode == 'bed_profile':
            ################### BED profile ##########################################
            print2(parameter, "\n############# BED profile #############")
            print2(parameter, "\tInput path:\t" + args.i)
            print2(parameter, "\tOutput path:\t" + os.path.join(args.o, args.t))
            if not args.organism:
                print("Please define organism...")
                sys.exit(1)
            else:
                print2(parameter, "\tOrganism:\t" + args.organism)

            if args.labels:
                args.labels = args.labels.split(",")
                args.sources = args.sources.split(",")
                if not args.sources:
                    print("Please define the sources files corresponding to the the labels.")
                    sys.exit(1)
                elif len(args.labels) != len(args.sources):
                    print("The number of labels doesn't match the number of sources.")
                    sys.exit(1)
                if args.strand:
                    strands = []
                    for i, bss in enumerate(args.strand.split(",")):
                        if bss == "T":
                            strands.append(True)
                            args.labels[i] += "(strand-specific)"
                        elif bss == "F":
                            strands.append(False)
                    args.strand = strands
                else:
                    args.strand = [True for i in args.labels]
                if args.other:
                    others = []
                    for i, bss in enumerate(args.other.split(",")):
                        if bss == "T":
                            others.append(True)
                        elif bss == "F":
                            others.append(False)
                    args.other = others
                else:
                    args.other = [True for i in args.labels]

            bed_profile = BedProfile(args.i, args.organism, args)
            bed_profile.cal_statistics()
            bed_profile.plot_distribution_length()
            bed_profile.plot_motif_composition()
            if args.biotype:
                bed_profile.plot_ref(ref_dir=args.biotype, tag="Biotype", other=True, strand=True, background=True)
            if args.repeats:
                bed_profile.plot_ref(ref_dir=args.repeats, tag="Repeats", other=True, background=True)
            if args.genposi:
                bed_profile.plot_ref(ref_dir=args.genposi, tag="Genetic position", other=False, strand=False)
            if args.labels:
                for i, label in enumerate(args.labels):
                    bed_profile.plot_ref(ref_dir=args.sources[i], tag=label, other=args.other[i], strand=args.strand[i], background=True)
            bed_profile.write_tables(args.o, args.t)
            bed_profile.save_fig(filename=os.path.join(args.o, args.t, "figure_" + args.t))
            bed_profile.gen_html(args.o, args.t)






        ################### Projection test ##########################################
        elif args.mode == 'projection':
            # Fetching reference and query EM
            print2(parameter, "\n############# Projection Test #############")
            print2(parameter, "\tReference:        " + args.r)
            print2(parameter, "\tQuery:            " + args.q)
            print2(parameter, "\tOutput directory: " + os.path.basename(args.o))
            print2(parameter, "\tExperiment title: " + args.t)

            projection = Projection(args.r, args.q, load_bed=args.load)
            projection.group_refque(args.g)
            projection.colors(args.c, args.color)

            if args.bg:
                print2(parameter, "\tBackground: " + args.bg)
                projection.set_background(bed_path=args.bg)
            if args.union:
                projection.ref_union()
                projection.projection_test(organism=args.organism)
                print2(parameter, "\tTaking union of references as the background. ")
            else:
                projection.projection_test(organism=args.organism)

            # generate pdf
            projection.plot(args.log, args.pw, args.ph)
            output(f=projection.fig, directory=args.o, folder=args.t, filename="projection_test",
                   extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)
            if args.bed:
                print2(parameter,
                       "\tOutput BED files: " + "/".join(os.path.join(args.o, args.t, "bed").split("/")[-3:]))
                projection.output_interq(directory=os.path.join(args.o, args.t, "bed"))
            # generate html 
            projection.gen_html(args.o, args.t, args=args)

            if args.table:
                projection.table(directory=args.o, folder=args.t)

            print("\nAll related files are saved in:  " + os.path.join(os.path.basename(args.o), args.t))
            t1 = time.time()
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1 - t0))))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ###########################################################################
        ################### Intersect Test ##########################################
        if args.mode == 'intersect':
            print2(parameter, "\n############ Intersection Test ############")
            print2(parameter, "\tReference:        " + args.r)
            print2(parameter, "\tQuery:            " + args.q)
            print2(parameter, "\tOutput directory: " + os.path.basename(args.o))
            print2(parameter, "\tExperiment title: " + args.t)

            # Fetching reference and query EM
            inter = Intersect(args.r, args.q, mode_count=args.m, organism=args.organism)

            # Grouping
            inter.group_refque(args.g)
            # Setting background
            inter.background(args.bg)

            # Extension
            if args.ex == 0:
                pass
            elif args.ex > 0:
                inter.extend_ref(args.ex)
            elif args.ex < 0:
                print("\n**** extension percentage(-ex) should be positive value, not negative.\n")
                sys.exit(1)

            inter.colors(args.c, args.color)
            print("\tProcessing data.", end="")
            sys.stdout.flush()
            inter.count_intersect(threshold=args.tc)

            # generate pdf
            print("\n\tGenerate graphics...")
            inter.barplot(logt=args.log)
            output(f=inter.bar, directory=args.o, folder=args.t, filename="intersection_bar",
                   extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)
            inter.stackedbar()
            output(f=inter.sbar, directory=args.o, folder=args.t, filename="intersection_stackedbar",
                   extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)
            inter.barplot(logt=args.log, percentage=True)
            output(f=inter.bar, directory=args.o, folder=args.t, filename="intersection_barp",
                   extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)

            if args.stest > 0:
                print("\tStatistical testing by randomizing the regions...")
                inter.stest(repeat=args.stest, threshold=args.tc, mp=args.mp)

            # generate html
            inter.gen_html(directory=args.o, title=args.t, align=50, args=args)

            t1 = time.time()
            print2(parameter, "\nAll related files are saved in:  " + os.path.join(os.path.basename(args.o), args.t))
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1 - t0))))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ###########################################################################
        ################### Jaccard test ##########################################
        if args.mode == "jaccard":
            """Return the jaccard test of every possible comparisons between two ExperimentalMatrix. 
            
            Method:
            The distribution of random jaccard index is calculated by randomizing query for given times. 
            Then, we compare the real jaccard index to the distribution and formulate p-value as 
            p-value = (# random jaccard > real jaccard)/(# random jaccard)
            
            """
            print("\n############## Jaccard Test ###############")
            jaccard = Jaccard(args.r, args.q)
            jaccard.group_refque(args.g)
            jaccard.colors(args.c, args.color)

            # jaccard test
            jaccard.jaccard_test(args.rt, args.organism)
            parameter = parameter + jaccard.parameter
            t1 = time.time()
            # ploting and generate pdf
            jaccard.plot(logT=args.nlog)
            for i, f in enumerate(jaccard.fig):
                output(f=f, directory=args.o, folder=args.t, filename="jaccard_test" + str(i + 1),
                       extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)
            # generate html
            jaccard.gen_html(args.o, args.t)

            if args.table:
                jaccard.table(directory=args.o, folder=args.t)

            print("\nAll related files are saved in:  " + os.path.join(current_dir, args.o, args.t))
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1 - t0))))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="Reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="Query_experimental_matrix.txt")
            list_all_index(path=args.o)

        ###########################################################################
        ################### Combinatorial Test ##########################################
        if args.mode == 'combinatorial':
            print("\n############ Combinatorial Test ############")
            # Fetching reference and query EM
            # comb = Combinatorial(args.r,args.q, mode_count=args.m, organism=args.organism)
            inter = Intersect(args.r, args.q, mode_count=args.m, organism=args.organism)
            # Setting background
            inter.background(args.bg)
            # Grouping
            inter.group_refque(args.g)

            # Extension
            if args.ex == 0:
                pass
            elif args.ex > 0:
                inter.extend_ref(args.ex)
            elif args.ex < 0:
                print("\n**** extension percentage(-ex) should be positive value, not negative.\n")
                sys.exit(1)
            # Combinatorial 
            print2(parameter, "Generating all combinatorial regions for further analysis...")
            inter.combinatorial()
            inter.count_intersect(threshold=args.tc, frequency=True)

            # generate pdf
            inter.colors_comb()
            # inter.barplot(args.log)
            # output(f=inter.bar, directory = args.output, folder = args.title, filename="intersection_bar",extra=matplotlib.pyplot.gci(),pdf=True,show=args.show)
            # if args.stackedbar:
            # inter.colors(args.c, args.color,ref_que = "ref")
            inter.comb_stacked_plot()

            output(f=inter.sbar, directory=args.o, folder=args.t, filename="intersection_stackedbar",
                   extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)
            if args.venn:
                inter.comb_venn(directory=os.path.join(args.o, args.t))

            # if args.lineplot:
            #    inter.comb_lineplot()
            if args.stest > 0:
                inter.stest(repeat=args.stest, threshold=args.tc, mp=args.mp)
            # generate html
            inter.gen_html_comb(directory=args.o, title=args.t, align=50, args=args)

            # parameter = parameter + inter.parameter
            t1 = time.time()
            print("\nAll related files are saved in:  " + os.path.join(current_dir, args.o, args.t))
            print2(parameter, "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t1 - t0))))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.r, directory=args.o, folder=args.t, filename="Reference_experimental_matrix.txt")
            copy_em(em=args.q, directory=args.o, folder=args.t, filename="Query_experimental_matrix.txt")
            # list_all_index(path=args.o)

        ###########################################################################
        ################### Boxplot ##########################################
        if args.mode == 'boxplot':
            print("\n################# Boxplot #################")
            boxplot = Boxplot(args.input, fields=[args.g, args.s, args.c], title=args.t, df=args.df)

            print2(parameter, "\nStep 1/5: Combining all regions")
            boxplot.combine_allregions()
            print2(parameter, "    " + str(len(boxplot.all_bed)) + " regions from all bed files are combined.")
            t1 = time.time()
            print2(parameter, "    --- finished in {0} secs\n".format(round(t1 - t0)))

            # Coverage of reads on all_bed
            print2(parameter, "Step 2/5: Calculating coverage of each bam file on all regions")
            boxplot.bedCoverage()
            t2 = time.time()
            print2(parameter, "    --- finished in {0} (H:M:S)\n".format(datetime.timedelta(seconds=round(t2 - t1))))

            # Quantile normalization
            print2(parameter, "Step 3/5: Quantile normalization of all coverage table")
            if args.nqn:
                print2(parameter, "    No quantile normalization.")
                boxplot.norm_table = boxplot.all_table
            else:
                boxplot.quantile_normalization()
            t3 = time.time()
            print2(parameter, "    --- finished in {0} secs\n".format(round(t3 - t2)))

            # Generate individual table for each bed
            print2(parameter, "Step 4/5: Constructing different tables for box plot")
            boxplot.tables_for_plot()
            # if args.table: boxplot.print_plot_table(directory = args.o, folder = args.t)
            t4 = time.time()
            print2(parameter, "    --- finished in {0} secs\n".format(round(t4 - t3)))

            # Plotting
            print2(parameter, "Step 5/5: Plotting")
            boxplot.group_tags(groupby=args.g, sortby=args.s, colorby=args.c)

            boxplot.group_data(directory=args.o, folder=args.t, log=args.nlog)
            boxplot.color_map(colorby=args.c, definedinEM=args.color)
            boxplot.plot(title=args.t, logT=args.nlog, scol=args.scol, ylim=args.ylim, pw=args.pw, ph=args.ph)
            if args.table:
                boxplot.print_plot_table(directory=args.o, folder=args.t)
            output(f=boxplot.fig, directory=args.o, folder=args.t, filename="boxplot", extra=matplotlib.pyplot.gci(),
                   pdf=True,
                   show=args.show)
            # HTML
            boxplot.gen_html(args.o, args.t, align=50)
            t5 = time.time()
            print2(parameter, "    --- finished in {0} secs\n".format(round(t5 - t4)))
            print2(parameter,
                   "Total running time is: " + str(datetime.timedelta(seconds=round(t5 - t0))) + " (H:M:S)\n")
            print("\nAll related files are saved in:  " + os.path.join(current_dir, args.o, args.t))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)

        ################### Lineplot #########################################
        if args.mode == 'lineplot':
            if args.scol and args.srow:
                print("** Err: -scol and -srow cannot be used simutaneously.")
                sys.exit(1)

            print("\n################ Lineplot #################")
            # Read experimental matrix
            t0 = time.time()
            if "reads" not in (args.g, args.col, args.c, args.row):
                print("Please add 'reads' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            # if "regions" not in (args.col, args.c, args.row):
            #     print("Please add 'regions' tag as one of grouping, sorting, or coloring argument.")
            #     sys.exit(1)

            if not os.path.isfile(args.input):
                print("Please check the input experimental matrix again. The given path is wrong.")
                sys.exit(1)
            print2(parameter, "Parameters:\tExtend length:\t" + str(args.e))
            print2(parameter, "\t\tRead size:\t" + str(args.rs))
            print2(parameter, "\t\tBin size:\t" + str(args.bs))
            print2(parameter, "\t\tStep size:\t" + str(args.ss))
            print2(parameter, "\t\tCenter mode:\t" + str(args.center + "\n"))

            lineplot = Lineplot(em_path=args.input, title=args.t, annotation=args.ga,
                                organism=args.organism, center=args.center, extend=args.e, rs=args.rs,
                                bs=args.bs, ss=args.ss, df=args.df, dft=args.dft,
                                fields=[args.g, args.col, args.row, args.c],
                                test=args.test, sense=args.sense, strand=args.strand, flipnegative=args.flip_negative,
                                outside=args.extend_outside, add_number=args.add_region_number)
            # Processing the regions by given parameters
            print2(parameter, "Step 1/3: Processing regions by given parameters")
            lineplot.relocate_bed()
            t1 = time.time()
            print2(parameter, "\t--- finished in {0} secs".format(str(round(t1 - t0))))

            if args.mp > 0:
                print2(parameter,
                       "\nStep 2/3: Calculating the coverage to all reads and averaging with multiprocessing ")
            else:
                print2(parameter, "\nStep 2/3: Calculating the coverage to all reads and averaging")
            lineplot.group_tags(groupby=args.g, rowby=args.row, columnby=args.col, colorby=args.c)
            lineplot.gen_cues()
            lineplot.coverage(sortby=args.row, mp=args.mp, log=args.log, average=args.average)
            t2 = time.time()
            print2(parameter, "\t--- finished in {0} (H:M:S)".format(str(datetime.timedelta(seconds=round(t2 - t1)))))

            # Plotting
            print2(parameter, "\nStep 3/3: Plotting the lineplots")
            lineplot.colormap(colorby=args.c, definedinEM=args.color)
            lineplot.plot(output=args.o, printtable=args.table, ylog=args.log,
                          scol=args.scol, srow=args.srow, w=args.pw, h=args.ph)
            for i, f in enumerate(lineplot.fig):
                output(f=f, directory=args.o, folder=args.t, filename="lineplot_" + lineplot.group_tags[i],
                       extra=matplotlib.pyplot.gci(), pdf=True, show=args.show)

            lineplot.gen_html(args.o, args.t)
            t3 = time.time()
            print2(parameter, "\t--- finished in {0} secs".format(str(round(t3 - t2))))
            print2(parameter,
                   "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t3 - t0))) + "(H:M:S)\n")
            print("\nAll related files are saved in:  " + os.path.join(current_dir, args.o, args.t))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)

        ################### Heatmap ##########################################
        if args.mode == 'heatmap':
            print("\n################# Heatmap #################")
            # Most part of heat map are the same as lineplot, so it share the same class as lineplot
            # Read experimental matrix
            t0 = time.time()

            if "reads" not in (args.g, args.col, args.c, args.row):
                print("Please add 'reads' tag as one of grouping, sorting, or coloring argument.")
                sys.exit(1)
            # if "regions" not in (args.g, args.col, args.c, args.row):
            #     print("Please add 'regions' tag as one of grouping, sorting, or coloring argument.")
            #     sys.exit(1)
            print2(parameter, "Parameters:\tExtend length:\t" + str(args.e))
            print2(parameter, "\t\tRead size:\t" + str(args.rs))
            print2(parameter, "\t\tBin size:\t" + str(args.bs))
            print2(parameter, "\t\tStep size:\t" + str(args.ss))
            print2(parameter, "\t\tCenter mode:\t" + str(args.center + "\n"))

            lineplot = Lineplot(em_path=args.input, title=args.t, annotation=args.ga,
                                organism=args.organism, center=args.center, extend=args.e, rs=args.rs,
                                bs=args.bs, ss=args.ss, df=False, fields=[args.col, args.row, args.c],
                                dft=args.dft, flipnegative=False, sense=False, strand=False, test=False)
            # Processing the regions by given parameters
            print2(parameter, "Step 1/4: Processing regions by given parameters")
            lineplot.relocate_bed()
            t1 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t1 - t0))))

            if args.mp:
                print2(parameter,
                       "\nStep 2/4: Calculating the coverage to all reads and averaging with multiprocessing ")
            else:
                print2(parameter, "\nStep 2/4: Calculating the coverage to all reads and averaging")
            lineplot.group_tags(groupby=args.col, sortby=args.row, colorby=args.c)
            lineplot.gen_cues()
            lineplot.coverage(sortby=args.s, heatmap=True, logt=args.log, mp=args.mp)
            t2 = time.time()
            print2(parameter, "    --- finished in {0} (h:m:s)".format(str(datetime.timedelta(seconds=round(t2 - t1)))))

            # Sorting 
            print2(parameter, "\nStep 3/4: Sorting the data for heatmap")
            lineplot.hmsort(sort=args.sort)
            t3 = time.time()
            print2(parameter, "    --- finished in {0} (h:m:s)".format(str(datetime.timedelta(seconds=round(t3 - t2)))))

            # Plotting
            print2(parameter, "\nStep 4/4: Plotting the heatmap")
            lineplot.hmcmlist(colorby=args.c, definedinEM=args.color)
            lineplot.heatmap(args.log)
            for i, name in enumerate(lineplot.hmfiles):
                output(f=lineplot.figs[i], directory=args.o, folder=args.t, filename=name, pdf=True, show=args.show)
            lineplot.gen_htmlhm(args.o, args.t)
            t4 = time.time()
            print2(parameter, "    --- finished in {0} secs".format(str(round(t4 - t3))))
            print2(parameter,
                   "\nTotal running time is : " + str(datetime.timedelta(seconds=round(t4 - t0))) + "(H:M:S)\n")
            print("\nAll related files are saved in:  " + os.path.join(current_dir, args.o, args.t))
            output_parameters(parameter, directory=args.o, folder=args.t, filename="parameters.txt")
            copy_em(em=args.input, directory=args.o, folder=args.t)
            list_all_index(path=args.o)

        ################### Venn Diagram ##########################################
        if args.mode == 'venn':
            print("\n################# Venn Diagram ###############")
            if not os.path.exists(os.path.join(args.o, args.t)):
                os.makedirs(os.path.join(args.o, args.t))
            sets = [s for s in [args.s1, args.s2, args.s3, args.s4] if s]
            venn = Venn(sets=sets, organism=args.organism)
            f = venn.venn_diagram(directory=args.o, title=args.t, labels=[args.l1, args.l2, args.l3, args.l4])
            output(f=f, directory=args.o, folder=args.t, filename="venn", pdf=True)

        ################### Integration ##########################################
        if args.mode == 'integration':
            print("\n################# Integration ###############")
            if args.ihtml:
                list_all_index(path=args.o)
