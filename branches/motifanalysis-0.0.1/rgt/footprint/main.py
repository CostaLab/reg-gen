#################################################################################################
# Footprinting with HMMs:
#     1. footprint: Predicts transcription factor binding sites using HMMs on DNase-seq and
#                   histone modification ChIP-seq data.
#################################################################################################

# Libraries
import os
import sys
import footprint

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
    params.append("  footprint                 Predicts transcription factor binding sites using")
    params.append("                            HMMs on DNase-seq and histone modification ChIP-seq")
    params.append("                            data.")

    # Print help manual if no arguments are given
    if(len(sys.argv) <= 1):
        for e in params: print e
        sys.exit(0)

    #################################################################################################
    ##### FUNCTIONS #################################################################################
    #################################################################################################

    # Footprint
    if(sys.argv[1] == "footprint"):
        footprint.main(sys.argv[2:])

    # Function does not exist
    else:
        sys.stderr.write("The function does not exist. Run the tool without arguments to see all\n")
        sys.stderr.write("available options.\n")


#################################################################################################
##### MAIN ######################################################################################
#################################################################################################

if __name__ == '__main__':
    main()


