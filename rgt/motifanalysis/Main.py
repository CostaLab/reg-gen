###################################################################################################
# Libraries
###################################################################################################

# Python 3 compatibility


# Python
import time
import sys
from random import seed
import argparse

# Internal
from .. import __version__
from . import Match
from . import Enrichment
from . import Mapper

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
    Match.options(matching)
    matching.set_defaults(func=Match.main)

    enrichment = subparsers.add_parser('enrichment', help='calculate statistics for MPBS files')
    Enrichment.options(enrichment)
    enrichment.set_defaults(func=Enrichment.main)

    # mapper = subparsers.add_parser('mapper', help='map genes names to TF IDs and viceversa')
    # Mapper.options(mapper)
    # mapper.set_defaults(func=Mapper.main)

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
