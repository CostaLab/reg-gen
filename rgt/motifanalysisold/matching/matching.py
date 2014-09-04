#################################################################################################
# Searches for motif instances on the genome using different algorithms.
#################################################################################################

# Python Libraries
import os
import sys
import math
import glob

# Local Libraries
import biopythonMM
import bitScoreMM
import fimoMM

# Distal Libraries
from .. util import *

def main(sysArg):

    #################################################################################################
    ##### PARAMETERS ################################################################################
    #################################################################################################

    # Parameters
    params = []

    params.append("\nmatching")
    params.append("Searches for motif instances on the genome using different algorithms.")

    params.append("\nRequired Input: ")
    params.append("  -pwm_list=<FILE1[,FILE2,...,FILEN]>")
    params.append("                            List of PWM files in jaspar format. If this")
    params.append("                            parameter is not used, then pwm_path need to be")
    params.append("                            used.")
    params.append("                            Format: jaspar.")
    params.append("                            Default: None.")
    params.append("  -pwm_path=<PATH>          A path containing many PWMs in jaspar format with")
    params.append("                            .pwm extension.")
    params.append("                            Default: None.")
    params.append("  -genome_list=<FILE1[,FILE2,...,FILEN]>")
    params.append("                            List of files containing the genomic sequences in")
    params.append("                            fasta format. Each chromosome must be put in the")
    params.append("                            header, as in UCSC format.")
    params.append("                            Format: fasta.")
    params.append("                            Default: None.")

    params.append("\nOptional Input: ")
    params.append("  -chrom_sizes_file=<FILE>  File containing the total length of each chromosome.")
    params.append("                            It is a plain text file containing the chromosome")
    params.append("                            name and the length in each line, separated by tab.")
    params.append("                            Format: plain text.")

    params.append("\nInput Parameters: ")
    params.append("  -organism=<STRING>        Organism considered on the analysis. Can be 'hg19'")
    params.append("                            or 'mm9'. All the default files are going to be")
    params.append("                            based on the chosen organism.")
    params.append("                            Default: hg19")
    params.append("  -search_method=<STRING>   Search method. The tool supports 'biopython', 'fimo'")
    params.append("                            and 'bitscore' (regular bit score-based search).")
    params.append("                            Default: bitscore.")
    params.append("  -scoring_method=<STRING>  Scoring method. The tool supports 'bitscore' (set")
    params.append("                            -bitscore parameter), 'fpr' (set -fpr and -precision")
    params.append("                            parameters), 'boyle' (set -high_cutoff and")
    params.append("                            -functional_depth parameters) and 'fimo' (set")
    params.append("                            -fimo_threshold parameter).")
    params.append("                            Default: bitscore.")
    params.append("  -pseudocounts=<FLOAT>     Pseudocounts to be added to raw counts of each PWM.")
    params.append("                            Default: 0.1.")
    params.append("  -bitscore=<FLOAT>         Bitscore cutoff. Used only when -scoring_method is")
    params.append("                            set to 'bitscore'.")
    params.append("                            Default: log_2(10000).")
    params.append("  -fpr=<FLOAT>              False positive rate cutoff. Used only when")
    params.append("                            -scoring_method is set to 'fpr'.")
    params.append("                            Default: 0.0001.")
    params.append("  -precision=<INT>          Score distribution precision. Used only when")
    params.append("                            -scoring_method is set to 'fpr'.")
    params.append("                            Default: 10000.")
    params.append("  -high_cutoff=<FLOAT>      High part of boyle's scoring method. Used only when")
    params.append("                            -scoring_method is set to 'boyle'.")
    params.append("                            Default: 0.7.")
    params.append("  -functional_depth=<FLOAT> Functional depth part of boyle's scoring method.")
    params.append("                            Used only when -scoring_method is set to 'boyle'.")
    params.append("                            Default: 0.9.")
    params.append("  -fimo_thresold=<FLOAT>    Threshold (p-value) for FIMO method only.")
    params.append("                            Default: 0.0001.")

    params.append("\nOutput Options: ")
    params.append("  -output_location=<PATH>   Path where the output files will be written.")
    params.append("                            Default: current directory.")
    params.append("  -bigbed=<Y|N>             Wether to output bed files as bigbed.")
    params.append("                            Default: Y.")
    params.append("  -print_mpbs=<Y|N>         Whether to output a bed file containing all MPBSs.")
    params.append("                            Default: Y.")
    params.append("  -print_graph_mmscore=<Y|N>")
    params.append("                            Whether to output graphs containing the motif")
    params.append("                            matching score distribution for the MPBSs.")
    params.append("                            Default: Y.")
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
    flagP = False; flagPP = False; flagG = False
    if("-pwm_list" in inputParameters.keys()):
        flagP = True
        inputParameters["-pwm_list"] = inputParameters["-pwm_list"].split(",")
    else: inputParameters["-pwm_list"] = None
    if("-pwm_path" in inputParameters.keys()):
        flagPP = True
        if(inputParameters["-pwm_path"][-1] != "/"): inputParameters["-pwm_path"] += "/"
    else: inputParameters["-pwm_path"] = None
    if("-genome_list" in inputParameters.keys()):
        flagG = True
        inputParameters["-genome_list"] = inputParameters["-genome_list"].split(",")
    else: inputParameters["-genome_list"] = None

    # Required input verification
    if(not flagP and not flagPP): 
        print "ERROR: You must specify at least one PWM file (-pwm_list) or path (-pwm_path)."
        sys.exit(0)
    if(not flagG): 
        print "ERROR: You must specify at least one genome file."
        sys.exit(0)

    # Optional Input
    if("-chrom_sizes_file" in inputParameters.keys()): pass
    else: 
        if(inputParameters["-organism"] == "hg19"): inputParameters["-chrom_sizes_file"] = constants.getChromSizes_HG19()
        elif(inputParameters["-organism"] == "mm9"): inputParameters["-chrom_sizes_file"] = constants.getChromSizes_MM9()
        else: inputParameters["-chrom_sizes_file"] = None

    # Input Parameters
    if("-search_method" in inputParameters.keys()): pass
    else: inputParameters["-search_method"] = "bitscore"
    if("-scoring_method" in inputParameters.keys()): pass
    else: inputParameters["-scoring_method"] = "bitscore"
    if("-pseudocounts" in inputParameters.keys()): inputParameters["-pseudocounts"] = float(inputParameters["-pseudocounts"])
    else: inputParameters["-pseudocounts"] = 0.1
    if("-bitscore" in inputParameters.keys()): inputParameters["-bitscore"] = float(inputParameters["-bitscore"])
    else: inputParameters["-bitscore"] = math.log(10000,2)
    if("-fpr" in inputParameters.keys()): inputParameters["-fpr"] = float(inputParameters["-fpr"])
    else: inputParameters["-fpr"] = 0.0001
    if("-precision" in inputParameters.keys()): inputParameters["-precision"] = int(inputParameters["-precision"])
    else: inputParameters["-precision"] = 10000
    if("-high_cutoff" in inputParameters.keys()): inputParameters["-high_cutoff"] = float(inputParameters["-high_cutoff"])
    else: inputParameters["-high_cutoff"] = 0.7
    if("-functional_depth" in inputParameters.keys()): inputParameters["-functional_depth"] = float(inputParameters["-functional_depth"])
    else: inputParameters["-functional_depth"] = 0.9
    if("-fimo_thresold" in inputParameters.keys()): inputParameters["-fimo_thresold"] = float(inputParameters["-fimo_thresold"])
    else: inputParameters["-fimo_thresold"] = 0.0001

    # Output Options
    if("-output_location" in inputParameters.keys()):
        if(inputParameters["-output_location"][-1] != "/"): inputParameters["-output_location"] += "/"
    else: inputParameters["-output_location"] = "./"
    if("-bigbed" in inputParameters.keys()):
        if(inputParameters["-bigbed"] == "Y"): inputParameters["-bigbed"] = True
        else: inputParameters["-bigbed"] = False
    else: inputParameters["-bigbed"] = True
    if("-print_mpbs" in inputParameters.keys()):
        if(inputParameters["-print_mpbs"] == "Y"): inputParameters["-print_mpbs"] = True
        else: inputParameters["-print_mpbs"] = False
    else: inputParameters["-print_mpbs"] = True
    if("-print_graph_mmscore" in inputParameters.keys()):
        if(inputParameters["-print_graph_mmscore"] == "Y"): inputParameters["-print_graph_mmscore"] = True
        else: inputParameters["-print_graph_mmscore"] = False
    else: inputParameters["-print_graph_mmscore"] = True

    #################################################################################################
    ##### MOTIF LIST ################################################################################
    #################################################################################################

    motifList = []
    if(flagP): motifList = inputParameters["-pwm_list"]
    else: motifList = glob.glob(inputParameters["-pwm_path"]+"*.pwm")
    
    #################################################################################################
    ##### MOTIF MATCHING ############################################################################
    #################################################################################################

    # Bit Score
    if(inputParameters["-search_method"] == "bitscore"):

        # Reading genome
        genomeDict = genome.readFastaFiles(inputParameters["-genome_list"])
        chrList = constants.getChromList(y=False,m=False,reference=[genomeDict]) # TODO Don't remove the Y and M like this. Has To remove in the genome.

        # Creating mpbsDict
        mpbsDict = dict([(e,[]) for e in chrList])        
        
        # Motif Matching
        for pwmFileName in motifList:
            bitScoreMM.bitScoreMM(pwmFileName,genomeDict,mpbsDict,inputParameters["-scoring_method"],inputParameters["-output_location"],inputParameters["-pseudocounts"],inputParameters["-bitscore"],inputParameters["-fpr"],inputParameters["-precision"],inputParameters["-high_cutoff"],inputParameters["-functional_depth"])

    # Biopython
    elif(inputParameters["-search_method"] == "biopython"):

        # Reading genome
        genomeDict = genome.readFastaFiles(inputParameters["-genome_list"])
        chrList = constants.getChromList(y=False,m=False,reference=[genomeDict]) # TODO Don't remove the Y and M like this. Has To remove in the genome.

        # Creating mpbsDict
        mpbsDict = dict([(e,[]) for e in chrList])        
        
        # Motif Matching
        for pwmFileName in motifList:
            biopythonMM.biopythonMM(pwmFileName,genomeDict,mpbsDict,inputParameters["-scoring_method"],inputParameters["-output_location"],inputParameters["-pseudocounts"],inputParameters["-bitscore"],inputParameters["-fpr"],inputParameters["-precision"],inputParameters["-high_cutoff"],inputParameters["-functional_depth"])

    # Fimo
    elif(inputParameters["-search_method"] == "fimo"):

        # Concatenate all genome files
        os.system("cat "+" ".join(inputParameters["-genome_list"])+" > "+inputParameters["-output_location"]+"genome.fa")

        # Creating mpbsDict
        mpbsDict = dict()

        # Motif Matching
        for pwmFileName in motifList:
            fimoMM.fimoMM(pwmFileName,inputParameters["-output_location"]+"genome.fa",mpbsDict,inputParameters["-scoring_method"],inputParameters["-output_location"],inputParameters["-pseudocounts"],inputParameters["-bitscore"],inputParameters["-fpr"],inputParameters["-precision"],inputParameters["-high_cutoff"],inputParameters["-functional_depth"],inputParameters["-fimo_thresold"])

        # Remove genome file
        os.system("rm "+inputParameters["-output_location"]+"genome.fa")

    # Method does not exist
    else:
        sys.stderr.write("The method does not exist. Run the tool without arguments to see all\n")
        sys.stderr.write("available options.\n")    

    #################################################################################################
    ##### PRINTING OUTPUT ###########################################################################
    #################################################################################################

    # Print MPBS
    if(inputParameters["-print_mpbs"]):
        
        # Fetching factorList
        factorList = []
        for k in mpbsDict.keys():
            for e in mpbsDict[k]:
                if(e[2] not in factorList): factorList.append(e[2])

        # Iterating on each factor
        for factorName in factorList:

            # Writing MPBSs
            outputFile = open(inputParameters["-output_location"]+factorName+".bed","w")
            for k in constants.getChromList(reference=[mpbsDict]):
                for e in mpbsDict[k]:
                    if(factorName == e[2]): outputFile.write("\t".join([k]+[str(b) for b in e])+"\n")
            outputFile.close()

            # Converting to bigBed
            if(inputParameters["-bigbed"]):
                os.system("sort -k1,1 -k2,2n "+inputParameters["-output_location"]+factorName+".bed > "+inputParameters["-output_location"]+factorName+"_sort.bed")
                bedFunctions.bedToBigBed(inputParameters["-output_location"]+factorName+"_sort.bed",inputParameters["-chrom_sizes_file"],inputParameters["-output_location"]+factorName+".bb",removeBed=True)
                os.system("rm "+inputParameters["-output_location"]+factorName+".bed")

    #################################################################################################
    ##### GRAPHS ####################################################################################
    #################################################################################################

    # Motif predicted binding sites score distribution
    if(inputParameters["-print_graph_mmscore"]):

        # Converting mpbsDict to format accepted by graph function
        mpbsDictGraph = bedFunctions.fitDictForBitscoreGraph(mpbsDict)

        # Creating graphs
        os.system("mkdir -p "+inputParameters["-output_location"]+"MM_score_distribution/")
        graphs.mpbsScoreHistogram([mpbsDictGraph],["black"],[1.0],["bitscore"],inputParameters["-output_location"]+"MM_score_distribution/",bins=100,outExt="png",useLegend=False)


