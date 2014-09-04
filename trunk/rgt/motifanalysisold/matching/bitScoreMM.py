#################################################################################################
# Performs biopython based motif matching.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Python Libraries
import os
import sys
import math

# Local Libraries
import biopythonMM

# Distal Libraries
from .. util import *

# External Libraries
from Bio import Motif

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def bitScoreMM(pwmFileName,genomeDict,mpbsDict,scoringMethod,tempLocation,pseudocounts=0.1,bitscore=12.0,fpr=0.01,precision=10**4,highCutoff=0.7,functionalDepth=0.9):
    """Performs basic motif matching algorithm and writes the results to a dictionary indexed by chromosome.

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
    pwm = createPwmDict(pwmFileName,pseudocounts)
    pwmName = pwmFileName.split("/")[-1].split(".")[0]
    pwmLen = len(pwm["A"])
    background = math.log(0.25,2)*pwmLen

    # Evaluating threshold
    pwmThreshold = 0.0
    if(scoringMethod == "bitscore"):
        pwmThreshold = bitscore
    elif(scoringMethod == "fpr"):
        bioPwm = biopythonMM.readPwmFile(pwmFileName,tempLocation,pseudocounts)
        sd = Motif.ScoreDistribution(bioPwm,precision=precision)
        pwmThreshold = sd.threshold_fpr(fpr)
    elif(scoringMethod == "boyle"):
        maxScore = 0.0
        minScore = 0.0 # TODO Boyle's rule is not suited for negative values.
        for i in range(0,pwmLen): maxScore += max(pwm["A"][i],pwm["C"][i],pwm["G"][i],pwm["T"][i])
        maxScore -= background
        pwmThreshold = min(highCutoff*maxScore,functionalDepth*(maxScore-minScore))
    else:
        sys.stderr.write("Choose a valid scoring method.\n")
        sys.exit(0)

    # Creating aditional parameters
    chrList = constants.getChromList(reference=[mpbsDict])
    tempMpbsDict = dict([(e,[]) for e in chrList])
    maxValue = -99.0
    revDict = dict([("A","T"),("T","A"),("C","G"),("G","C"),("N","N")])

    # Iterating on chromosomes
    for chrName in chrList:

        # Reading genome
        sequence = genomeDict[chrName].upper()

        # Performing motif matching
        for pos in xrange(0,len(sequence)-pwmLen+1):
            scoreF = -background; scoreR = -background;
            for i in range(0,pwmLen):
                scoreF += pwm[sequence[pos+i]][i]
                scoreR += pwm[revDict[sequence[pos+pwmLen-i-1]]][i]
            if(scoreF > pwmThreshold):
                if(scoreF > maxValue): maxValue = scoreF
                tempMpbsDict[chrName].append([pos,pos+pwmLen,pwmName,scoreF,"+"])
            if(scoreR > pwmThreshold):
                if(scoreR > maxValue): maxValue = scoreR
                tempMpbsDict[chrName].append([pos,pos+pwmLen,pwmName,scoreR,"-"])  

    # Update scores - new scores are within [0,1000]
    for chrName in chrList:
        for e in tempMpbsDict[chrName]:
            mpbsDict[chrName].append([e[0],e[1],e[2],int(1000*(e[3]-pwmThreshold)/(maxValue-pwmThreshold)),e[4]])
    
    return 0

def createPwmDict(pwmFileName, pseudocounts=0.0):
    """Reads a PWM file in Jaspar format and returns a PWM represented as a dictionary of nucleotides.

    Keyword arguments:
    pwmFileName -- The path + name of the PWM file.
    pseudocounts -- Amount of pseudocounts to add in each matrix cell. (default 0.0)

    Returns:
    pwm -- PWM as a dictionary of nucleotides.
    """
    # Creating PWM from input
    pwmFile = open(pwmFileName,"r")
    pwm = dict([("A",[]),("C",[]),("G",[]),("T",[]),("N",[])])
    nucList = ["A","C","G","T"]; counter = 0
    for line in pwmFile:
        pwm[nucList[counter]] = [float(e)+pseudocounts for e in line.strip().split(" ")]
        counter += 1
    motifLen = len(pwm["A"])
    for i in range(0,motifLen):
        total = pwm["A"][i]+pwm["C"][i]+pwm["G"][i]+pwm["T"][i]
        pwm["A"][i] = math.log(pwm["A"][i]/total,2); pwm["C"][i] = math.log(pwm["C"][i]/total,2); pwm["G"][i] = math.log(pwm["G"][i]/total,2); pwm["T"][i] = math.log(pwm["T"][i]/total,2)
        pwm["N"].append(0.0)
    return pwm


