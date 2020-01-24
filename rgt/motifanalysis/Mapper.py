###################################################################################################
# Libraries
###################################################################################################

# Python
import os
from glob import glob
import time
import sys

# Internal
from rgt.Util import ErrorHandler, MotifData, GenomeData, npath
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.GeneSet import GeneSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion
from rgt.AnnotationSet import AnnotationSet
from .Motif import Motif
from .Util import bed_to_bb

# External
from pysam import Fastafile
from MOODS import tools, scan


###################################################################################################
# Functions
###################################################################################################

"""
The "Mapper" is meant to provide a way of doing the following:

    - converting Gene IDs to Motif IDs for specific TF databases
    - converting Motif IDs to Gene IDs
    - 
"""


def options(parser):
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--motif-to-genes", action="store_true",
                       help="Motif IDs will be converted to gene IDs")
    group.add_argument("--gene-to-motifs", action="store_true",
                       help="Gene IDs will be converted to motif IDs")

    parser.add_argument("--motif-db", choices=('jaspar_vertebrates', 'hocomoco',
                                               'uniprobe_primary', 'uniprobe_secondary'))
    parser.add_argument("--gene-db", choices=('ensembl', 'names'))

    parser.add_argument('input', metavar='PATH', type=str,
                        help='File containing either the genes or motifs to convert (one for each line)')


def main(args):
    """
    Performs motif mapping.
    """

    print(args)
