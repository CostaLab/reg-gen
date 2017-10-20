###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function
import sys
import time
from random import seed

# Internal
from rgt import __version__
from rgt.Util import ErrorHandler
from Match import main_matching
from Enrichment import main_enrichment

"""
Motif matching and enrichment based on motif PSSM. Can perform either Matching (finds all putative
Motif Predicted Binding Sites for a set of motifs and genomic regions) or Enrichment
(assigns statistical significance to each found MPBS using fisher test on target and background).

Authors: Eduardo G. Gusmao, Fabio Ticconi
"""


def main():
    start = time.time()

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
