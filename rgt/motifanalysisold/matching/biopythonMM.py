#################################################################################################
# Performs biopython based motif matching.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Python Libraries
import os
import sys

# Distal Libraries
from .. util import *

# External Libraries
from Bio import Motif

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def biopythonMM(pwmFileName,genomeDict,mpbsDict,scoringMethod,tempLocation,pseudocounts=0.1,bitscore=12.0,fpr=0.01,precision=10**4,highCutoff=0.7,functionalDepth=0.9):
    """Performs Biopython based motif matching and writes the results to a dictionary indexed by chromosome.

    Keyword arguments:
    pwmFileName -- PWM file name.
    genomeDict -- Genome dictionary.
    mpbsDict -- Dictionary of MPBSs to insert the results.
    scoringMethod -- Method to evaluate which MPBSs are enriched.
    tempLocation -- Location to write temporary PWM files in order to help PWM creation and pseudocounting.
    pseudocounts -- Amount of pseudocounts to add in each PWM matrix's cell. (default 0.1)
    bitscore -- The cutoff bitscore value. (default 12.0)
    fpr -- False positive rate to determine the cutoff value. (default 0.01)
    precision -- Motif score distribution precision. (default 10**4)
    highCutoff -- High cutoff for Boyle's rule. (default 0.7)
    functionalDepth -- Functional depth for Boyle's rule. (default 0.9)

    Returns:
    mpbsDict -- This method inserts entries on the mpbsDict.
    """
    
    # Reading PWM
    pwm = readPwmFile(pwmFileName,tempLocation,pseudocounts)
    pwmName = pwmFileName.split("/")[-1].split(".")[0]
    pwmLen = len(pwm)

    # Evaluating threshold
    pwmThreshold = 0.0
    if(scoringMethod == "bitscore"):
        pwmThreshold = bitscore
    elif(scoringMethod == "fpr"):
        sd = Motif.ScoreDistribution(pwm,precision=precision)
        pwmThreshold = sd.threshold_fpr(fpr)
    elif(scoringMethod == "boyle"):
        maxScore = pwm.max_score()
        minScore = 0.0 # TODO Boyle's rule is not suited for negative values.
        pwmThreshold = min(highCutoff*maxScore,functionalDepth*(maxScore-minScore))
    else:
        sys.stderr.write("Choose a valid scoring method.\n")
        sys.exit(0)

    # Creating aditional parameters
    chrList = constants.getChromList(reference=[mpbsDict])
    tempMpbsDict = dict([(e,[]) for e in chrList])
    maxValue = -99.0

    # Iterating on chromosomes
    for chrName in chrList:

        # Reading genome
        sequence = genomeDict[chrName]

        # Performing biopython's motif matching
        for pos, score in pwm.search_pwm(sequence,threshold=pwmThreshold):
            if(score > maxValue): maxValue = score
            if(pos >= 0): tempMpbsDict[chrName].append([pos,pos+pwmLen,pwmName,score,"+"])
            else: tempMpbsDict[chrName].append([-pos,-pos+pwmLen,pwmName,score,"-"])

    # Update scores - new scores are within [0,1000]
    for chrName in chrList:
        for e in tempMpbsDict[chrName]:
            mpbsDict[chrName].append([e[0],e[1],e[2],int(1000*(e[3]-pwmThreshold)/(maxValue-pwmThreshold)),e[4]])
    
    return 0

def readPwmFile(pwmFileName, outputLocation, pseudocounts=0.0):
    """Reads a PWM file in Jaspar format and returns a Biopython PWM object.

    Keyword arguments:
    pwmFileName -- The path + name of the PWM file.
    outputLocation -- Path to write temporary pseudocount-pwm PWM.
    pseudocounts -- Amount of pseudocounts to add in each matrix cell. (default 0.0)

    Returns:
    pwm -- Biopython PWM object.
    """

    # Adding pseudocounts
    pwmFile = open(pwmFileName,"r");
    tempFileName = outputLocation+pwmFileName.split("/")[-1]+"temp"
    pwmFileT = open(tempFileName,"w")
    for line in pwmFile: pwmFileT.write(" ".join([str(float(e)+pseudocounts) for e in line.strip().split(" ")])+"\n")    
    pwmFile.close()
    pwmFileT.close()

    # Creating PWM from pseudocounted input
    pwmFile = open(tempFileName,"r")
    pwm = Motif.read(pwmFile,"jaspar-pfm")
    pwmFile.close()
    os.system("rm "+tempFileName)
    return pwm


