###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import sys
import time
from random import seed
import argparse

# Internal
from .. import __version__
from .Match import main_matching, matching_options
from .Enrichment import main_enrichment, enrichment_options

"""
Motif matching and enrichment based on motif PSSM. Can perform either Matching (finds all putative
Motif Predicted Binding Sites for a set of motifs and genomic regions) or Enrichment
(assigns statistical significance to each found MPBS using fisher test on target and background).

Authors: Eduardo G. Gusmao, Fabio Ticconi
"""


def main():
    start = time.time()

    seed(42)

    version_message = "Motif Analysis - Regulatory Analysis Toolbox (RGT) - v" + str(__version__)

    parser = argparse.ArgumentParser(prog='rgt-motifanalysis')
    parser.add_argument('--version', action='version', version=version_message)

    subparsers = parser.add_subparsers(help='Commands:')

    matching = subparsers.add_parser('matching', help='find all MPBS from input files')
    matching_options(matching)
    matching.set_defaults(func=main_matching)

    enrichment = subparsers.add_parser('enrichment', help='calculate statistics for MPBS files')
    enrichment_options(enrichment)
    enrichment.set_defaults(func=main_enrichment)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)

    secs = time.time() - start
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)
    print()
    print("[total time: ", "%dh %dm %ds" % (h, m, s), "]", sep="")


if __name__ == "__main__":
    main()
