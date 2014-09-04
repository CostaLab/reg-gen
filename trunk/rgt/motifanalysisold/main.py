#################################################################################################
# Motif Analysis:
#     1. enrichment: Fisher enrichment test.
#     2. matching: Motif matching.
#################################################################################################

# Libraries
import os
import sys
from enrichment import enrichment
from matching import matching

def main():
    
    #################################################################################################
    ##### PARAMETERS ################################################################################
    #################################################################################################

    # Parameters
    params = []

    params.append("\nRegulatoryGenomicsToolbox -- Version 0.0.1")
    params.append("Toolkit to perform common analysis of regulatory genomics data.")
    params.append("For full documentation visit the following webpage:")
    params.append("www.followingwebpage.com")

    params.append("\nAvailable Functions: ")
    params.append("  enrichment                Performs fisher exact test with multiple testing")
    params.append("                            correction on motif predicted binding sites.")
    params.append("  matching                  Searches for motif instances on the genome using")
    params.append("                            different algorithms.")

    # Print help manual if no arguments are given
    if(len(sys.argv) <= 1):
        for e in params: print e
        sys.exit(0)

    #################################################################################################
    ##### FUNCTIONS #################################################################################
    #################################################################################################

    # Motif Statistics
    if(sys.argv[1] == "enrichment"):
        enrichment.main(sys.argv[2:])

    # Motif Matching
    elif(sys.argv[1] == "matching"):
        matching.main(sys.argv[2:])

    # Function does not exist
    else:
        sys.stderr.write("The function does not exist. Run the tool without arguments to see all\n")
        sys.stderr.write("available options.\n")


#################################################################################################
##### MAIN ######################################################################################
#################################################################################################

if __name__ == '__main__':
    main()


