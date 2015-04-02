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
import bitScoreMM

# Distal Libraries
from .. util import *

# External Libraries
from Bio import Motif

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def fimoMM(pwmFileName,genomeFile,mpbsDict,scoringMethod,tempLocation,pseudocounts=0.1,bitscore=12.0,fpr=0.01,precision=10**4,highCutoff=0.7,functionalDepth=0.9,threshold=0.0001):
    """Performs FIMO motif matching algorithm and writes the results to a dictionary indexed by chromosome.

    Keyword arguments:
    pwmFileName -- PWM file name.
    genomeFile -- Fasta file containing the regions to be analyzed
    mpbsDict -- Dictionary of MPBSs to insert the results.
    scoringMethod -- Method to evaluate which MPBSs are enriched.
    tempLocation -- Location to write temporary PWM files in order to help PWM creation and pseudocounting.
    pseudocounts -- Amount of pseudocounts to add in each PWM matrix's cell. (default 0.1)
    bitscore -- The cutoff bitscore value. (default 12.0)
    fpr -- False positive rate to determine the cutoff value. (default 0.01)
    precision -- Motif score distribution precision. (default 10**4)
    highCutoff -- High cutoff for Boyle's rule. (default 0.7)
    functionalDepth -- Functional depth for Boyle's rule. (default 0.9)
    threshold -- The cutoff threshold value. (default 0.0001)

    Returns:
    mpbsDict -- This method inserts entries on the mpbsDict.
    """

    # Converting jaspar to MEME
    memeFileName = jasparToMeme(pwmFileName, tempLocation, pseudocounts)
    tempPath = "/".join(memeFileName.split("/")[:-1])+"/"
    fimoFileName = tempPath+"results.txt"
    errorOutputName = tempPath+"error.txt"

    # Evaluating threshold
    pwmThreshold = 0.0
    if(scoringMethod == "bitscore"):
        pwmThreshold = bitscore
        threshold = 0.1
    elif(scoringMethod == "fpr"):
        bioPwm = biopythonMM.readPwmFile(pwmFileName,tempLocation,pseudocounts)
        sd = Motif.ScoreDistribution(bioPwm,precision=precision)
        pwmThreshold = sd.threshold_fpr(fpr)
        threshold = 0.1
        print bioPwm.max_score()
    elif(scoringMethod == "boyle"):
        maxScore = 0.0
        minScore = 0.0 # TODO Boyle's rule is not suited for negative values.
        pwmBoyle = bitScoreMM.createPwmDict(pwmFileName,pseudocounts)
        pwmLen = len(pwmBoyle["A"])
        for i in range(0,pwmLen): maxScore += max(pwmBoyle["A"][i],pwmBoyle["C"][i],pwmBoyle["G"][i],pwmBoyle["T"][i])
        background = math.log(0.25,2)*pwmLen
        maxScore -= background
        pwmThreshold = min(highCutoff*maxScore,functionalDepth*(maxScore-minScore))
        threshold = 0.1
    elif(scoringMethod == "fimo"):
        pass
    else:
        sys.stderr.write("Choose a valid scoring method.\n")
        sys.exit(0)

    # Performing FIMO
    os.system("fimo --text --verbosity 1 --max-stored-scores 1000000 --output-pthresh "+str(threshold)+" "+memeFileName+" "+genomeFile+" > "+fimoFileName+" 2> "+errorOutputName)

    # Reading FIMO output
    tempMpbsDict = dict()
    fimoFile = open(fimoFileName,"r")
    fimoFile.readline()
    maxValue = -999
    for line in fimoFile:
        ll = line.strip().split("\t")
        ll = [ll[0][0],ll[0][1:]]+ll[1:]
        if(scoringMethod != "fimo" and float(ll[5]) < pwmThreshold): continue
        if(float(ll[5]) > maxValue): maxValue = float(ll[5])
        if(ll[2] in tempMpbsDict.keys()):
            if(ll[0] == "+"): tempMpbsDict[ll[2]].append([int(ll[3])-1,int(ll[4]),ll[1],float(ll[5]),ll[0]])
            else: tempMpbsDict[ll[2]].append([int(ll[4])-1,int(ll[3]),ll[1],float(ll[5]),ll[0]])
        else: 
            if(ll[0] == "+"): tempMpbsDict[ll[2]] = [[int(ll[3])-1,int(ll[4]),ll[1],float(ll[5]),ll[0]]]
            else: tempMpbsDict[ll[2]] = [[int(ll[4])-1,int(ll[3]),ll[1],float(ll[5]),ll[0]]]
    fimoFile.close()

    # Update scores and remove MPBSs with score below pwmThreshold (if it is being used)
    for chrName in tempMpbsDict.keys():
        for e in tempMpbsDict[chrName]:
            if(chrName in mpbsDict.keys()): mpbsDict[chrName].append([e[0],e[1],e[2],int(1000*(e[3]-pwmThreshold)/(maxValue-pwmThreshold)),e[4]])
            else: mpbsDict[chrName] = [[e[0],e[1],e[2],int(1000*(e[3]-pwmThreshold)/(maxValue-pwmThreshold)),e[4]]]

    # Removing temporary PWM folder
    os.system("rm -rf "+"/".join(memeFileName.split("/")[:-1]))

    return 0

def jasparToMeme(pwmFileName, outputLocation, pseudocounts=0.0):
    """Converts a jaspar PWM to the meme format and returns the name of the file.

    Keyword arguments:
    pwmFileName -- The path + name of the PWM file.
    outputLocation -- Path to write temporary pseudocount-pwm PWM.
    pseudocounts -- Amount of pseudocounts to add in each matrix cell. (default 0.0)

    Returns:
    pwmFileName -- 
    <pwmName>.meme -- Meme version of the PWM.
    """

    # Adding pseudocounts
    pwmFile = open(pwmFileName,"r");
    tempFolder = outputLocation+pwmFileName.split("/")[-1].split(".")[0]+"/"
    os.system("mkdir -p "+tempFolder)
    tempFileName = tempFolder+pwmFileName.split("/")[-1].split(".")[0]+".pfm"
    pwmFileT = open(tempFileName,"w")
    for line in pwmFile: pwmFileT.write(" ".join([str(float(e)+pseudocounts) for e in line.strip().split(" ")])+"\n")
    pwmFile.close()
    pwmFileT.close()

    # Creating meme output
    memeFileName = tempFolder+pwmFileName.split("/")[-1].split(".")[0]+".meme"
    os.system("jaspar2meme -pfm "+tempFolder+" > "+memeFileName)
    
    return memeFileName


