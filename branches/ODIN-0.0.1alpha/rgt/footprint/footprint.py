#################################################################################################
# Predicts transcription factor binding sites using HMMs on DNase-seq and
# histone modification ChIP-seq data.
#################################################################################################

# Python Libraries
import os
import sys

# Local Libraries
import signalProcessing

def main(sysArg):

    #################################################################################################
    ##### PARAMETERS ################################################################################
    #################################################################################################

    # Parameters
    params = []

    params.append("\nfootprint")
    params.append("Predicts transcription factor binding sites using HMMs on DNase-seq and")
    params.append("histone modification ChIP-seq data.")

    params.append("\nRequired Input: ")
    params.append("  -dnase_file=<FILE>        File containing the DNase-seq signal. The file")
    params.append("                            consists of the genomic coverage (read overlap)")
    params.append("                            of the aligned reads. A high-resolution signal")
    params.append("                            should be given correspoding to the overlap of all")
    params.append("                            the 5' read base pairs.")
    params.append("                            Format: bigwig.")
    params.append("                            Default: None.")
    params.append("  -histone_file=<FILE>      File containing the ChIP-seq signal from histone")
    params.append("                            modifications. The file consists of the genomic")
    params.append("                            coverage of the aligned reads. Each read should be")
    params.append("                            extended to the average length of immunoprecipitated")
    params.append("                            fragments.")
    params.append("                            Format: bigwig.")
    params.append("                            Default: None.")
    params.append("  -hmm_file=<FILE>          File containing the HMM parameters. Check the manual")
    params.append("                            for more information on the format and how to obtain")
    params.append("                            these models.")
    params.append("                            Format: hmm.")
    params.append("                            Default: None.")

    params.append("\nOptional Input: ")
    params.append("  -coord_file=<FILE>        Coordinate file containing the regions of the genome")
    params.append("                            where the footprinting will be performed. If")
    params.append("                            None, the footprinting is applied to the complete")
    params.append("                            genome.")
    params.append("                            Format: bed.")
    params.append("                            Default: Complete genome tiled in 10000bp regions.")


    params.append("\nInput Parameters: ")
    params.append("  -organism=<STRING>        Organism considered on the analysis. Can be 'hg19'")
    params.append("                            or 'mm9'. All the default files are going to be")
    params.append("                            based on the chosen organism.")
    params.append("                            Default: hg19")


    params.append("\nOutput Options: ")
    params.append("  -output_location=<PATH>   Path where the output files will be written.")
    params.append("                            Default: current directory.")
    params.append("  -aaaaaaaaaaaaaaaa=<AAAA>  aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    params.append("                            aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    params.append("                            aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
    params.append("                            aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa.")
    params.append("                            Format: aaaaa.")
    params.append("                            Default: aaaaaaa.")

    params.append("")
    if(len(sysArg) < 1):
        for e in params: print e
        sys.exit(0)

    #################################################################################################
    ##### INPUT #####################################################################################
    #################################################################################################

    # Input parameters dictionary
    inputParameters = util.readInputParameters(sysArg)

    # Organism
    if("-organism" in inputParameters.keys()): pass
    else: inputParameters["-organism"] = "hg19"

    # Required Input
    flagD = False; flagH = False; flagM = False
    if("-dnase_file" in inputParameters.keys()): flagD = True
    else: inputParameters["-dnase_file"] = None
    if("-histone_file" in inputParameters.keys()): flagH = True
    else: inputParameters["-histone_file"] = None
    if("-hmm_file" in inputParameters.keys()): flagM = True
    else: inputParameters["-hmm_file"] = None

    # Required input verification
    if(not flagD or not flagH or not flagM): 
        print "ERROR: You must specify all of these files:\n       -dnase_file, -histone_file and -hmm_file."
        sys.exit(0)

    # Optional Input
    if("-coord_file" in inputParameters.keys()): pass
    else: inputParameters["-coord_file"] = None

    # Input Parameters

    # Output Options
    if("-output_location" in inputParameters.keys()):
        if(inputParameters["-output_location"][-1] != "/"): inputParameters["-output_location"] += "/"
    else: inputParameters["-output_location"] = "./"

    #################################################################################################
    ##### CREATING COORDINATES ######################################################################
    #################################################################################################
    
    # Creating coordinates
    if(inputParameters["-coord_file"]): # If coordinates were given by user
        pass
    else: # Create tiled genome-wide coordinates
        pass

    #################################################################################################
    ##### CREATING HMM ##############################################################################
    #################################################################################################

    

    #################################################################################################
    ##### NORMALIZING SIGNAL ########################################################################
    #################################################################################################

    

    #################################################################################################
    ##### APPLYING HMM ##############################################################################
    #################################################################################################

    

    #################################################################################################
    ##### PRINTING OUTPUT ###########################################################################
    #################################################################################################

    

    


