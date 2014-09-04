#################################################################################################
# Perform operations on motifs.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Python Libraries
import os
import itertools

# Distal Libraries
from .. util import *

# External Libraries
from Bio import Motif

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def motifMatchingBiopython(combinationList,pwmList,coordDict,pwmLocation,genomeList,tempLocation,fpr=0.01,pseudocounts=0.0,precision=10**4,color="black"):
    """Performs Biopython based motif matching and returns a list containing the matches and
       writes the results on bed files.

    Keyword arguments:
    combinationList -- List of the number of cobinding combinations.
    pwmList -- List of PWMs where each entry represents the name of a PWM file.
    coordDict -- Dictionary of coordinates where the motif matching will be applied.
    pwmLocation -- Path containing the motif pwm files.
    genomeList -- List of fasta files containing the sequences to perform the motif matching, where the headers are the chromosomes.
    tempLocation -- Location to write temporary PWM files in order to help PWM creation and pseudocounting.
    fpr -- False positive rate to determine the cutoff value. (default 0.01)
    pseudocounts -- Amount of pseudocounts to add in each PWM matrix's cell. (default 0.0)
    precision -- Motif score distribution precision. (default 10**4)
    color -- Color of the bed entries. Can be 'green', 'red' or 'black'. (default 'black')

    Returns:
    mpbsDict -- Dictionary (for each PWM) of dictionaries (for each chromosome) of motif predicted binding sites.
    statDict -- Dictionary of statistics for Fisher test concerning the number of motifs inside enriched regions.
    geneDict -- Dictionary of genes (position NAME in bed file) that contains each motif.
    """
    
    # Reading PWM
    pwmDict = dict()
    for pwmName in pwmList: pwmDict[pwmName] = readPwmFile(pwmLocation+pwmName+".pwm","/".join(tempLocation.split("/")[:-1])+"/",pseudocounts)

    # Evaluating thresholds
    pwmThresholdDict = dict()
    for pwmName in pwmList:
        sd = Motif.ScoreDistribution(pwmDict[pwmName],precision=precision)
        pwmThresholdDict[pwmName] = sd.threshold_fpr(fpr)

    # Reading genome
    genomeDict = genome.readFastaFiles(genomeList)

    # Creating chromosome list
    chrList = constants.getChromList(reference=[coordDict])
    # Removing chrX, chrY and chrM
    # TODO Stop removing these chromosomes
    #chrListT = []
    #for e in chrList:
    #    if(e not in ["chrX", "chrY", "chrM"]): chrListT.append(e)
    #chrList = chrListT

    # Evaluating bed additionals
    if(color == "green"): color = "0,130,0"
    elif(color == "red"): color = "130,0,0"
    elif(color == "black"): color = "0,0,0"

    # Create combinations dictionary keys
    combKeys = []
    for c in combinationList:
        for b in [",".join(e) for e in itertools.combinations(pwmList,c)]: combKeys.append(b)

    # Iterating on chromosomes
    mpbsDict = dict([(e,dict()) for e in pwmDict.keys()])
    statDict = dict([(e,[0,0]) for e in combKeys]) # Left is evidence / Right is not evidence
    geneDict = dict([(e,[]) for e in combKeys])
    maxDict = dict([(e,-99.0) for e in pwmDict.keys()])
    ct=0
    for chrName in chrList:

        # Reading genome
        if(chrName not in genomeDict.keys()): continue
        sequence = genomeDict[chrName]

        # Iterating on coordinate dictionary
        for e in mpbsDict.keys(): mpbsDict[e][chrName] = []
        for coord in coordDict[chrName]:
            ct=ct+1
            #print "region", ct
            # Getting current sequence based on coordinates
            currSeq = sequence[coord[0]:coord[1]]

            # Keeping track of the factors found in this coordinate
            flagMotifs = dict([(e,False) for e in pwmDict.keys()])

            # Iterating on PWMs
            for pwmName in pwmDict.keys():
                pwmLen = len(pwmDict[pwmName])
                for pos, score in pwmDict[pwmName].search_pwm(currSeq,threshold=pwmThresholdDict[pwmName]):
                    if(score > maxDict[pwmName]): maxDict[pwmName] = score
                    if(pos >= 0): mpbsDict[pwmName][chrName].append([pos+coord[0],pos+coord[0]+pwmLen,pwmName,score,"+",pos+coord[0],pos+coord[0]+pwmLen,color])
                    else: mpbsDict[pwmName][chrName].append([-pos+coord[0],-pos+coord[0]+pwmLen,pwmName,score,"-",-pos+coord[0],-pos+coord[0]+pwmLen,color])
                    flagMotifs[pwmName] = True
            
            # Updating statistic counts and genes
            motifsFoundList = [k for k in pwmList if flagMotifs[k]]
            motifsFoundKeys = []
            motifsNotFoundKeys = [e for e in combKeys]
            for c in combinationList:
                for b in [",".join(e) for e in itertools.combinations(motifsFoundList,c)]:
                    motifsFoundKeys.append(b)
                    motifsNotFoundKeys.remove(b)
            for k in motifsFoundKeys:
                statDict[k][0] += 1
                for e in coord[2].split(":"): geneDict[k].append(e)
            for k in motifsNotFoundKeys:
                statDict[k][1] += 1

    # Update scores - new scores are within [0,1000]
    for pwmName in pwmDict.keys():
        for chrName in mpbsDict[pwmName].keys():
            for e in mpbsDict[pwmName][chrName]:
                e[3] = int(1000*(e[3]-pwmThresholdDict[pwmName])/(maxDict[pwmName]-pwmThresholdDict[pwmName]))

    # Remove repetitive genes from geneList
    for k in geneDict.keys(): geneDict[k] = list(set(geneDict[k]))
    
    return mpbsDict, statDict, geneDict

def readMotifMatching(combinationList,coordDict,pwmFileNameList,color="black",pwmReferenceList=None,coordFileName=None):
    """Reads motif predicted binding sites files and creates necessary structures for the statistical test.

    Keyword arguments:
    combinationList -- List of the number of cobinding combinations.
    coordDict -- Dictionary of coordinates where the motif matching was applied.
    pwmFileNameList -- List of PWMs files where each entry's name will represent the name of the motif.
                       Alternatively, it can be a single file containing all the MPBSs and their name on the NAME field.
    color -- Color of the bed entries. Can be 'green', 'red' or 'black'. (default 'black')
    pwmReferenceList -- Optional argument. In case pwmFileNameList is a single file (final motif matching file), this
                        parameter can be set to be a pwmList that will preserve the order of the pwmList. This is useful
                        in the case you want the same combinations of cobinding factors be created. (default None)
    coordFileName -- If the motif matching file entries are in bigBed format, it need to be converted to a bed in order
                     to be read. In this case where different executions are looking into the same set of bigBed files,
                     there will be a conflict during the deletions of the bed files. In order to remove this conflict, 
                     this additional argument is going to be passed ONLY to be used to name the created bed files. The
                     coordinates ARE NOT used in this function.

    Returns:
    mpbsDict -- Dictionary (for each PWM) of dictionaries (for each chromosome) of motif predicted binding sites.
    statDict -- Dictionary of statistics for Fisher test concerning the number of motifs inside enriched regions.
    geneDict -- Dictionary of genes (position NAME in bed file) that contains each motif.
    """

    # Reading all MPBSs
    pwmList = []; allMpbsDict = dict()
    if(isinstance(pwmFileNameList, list)):
        for pwmFileName in pwmFileNameList:
            pwmList.append(".".join(pwmFileName.split("/")[-1].split(".")[:-1]))
            pwmFileNameToRead = pwmFileName
            removeBed = False
            print pwmFileName, pwmFileName.split(".")[-1]
            if(pwmFileName.split(".")[-1] == "bb"): 
                coordName = ".".join(coordFileName.split("/")[-1].split(".")[:-1])
                bedFileName = pwmFileName[:-3]+"_"+coordName+".bed"
                bedFunctions.bigBedToBed(pwmFileName, bedFileName, removeBed=False)
                removeBed = True
                pwmFileNameToRead = bedFileName
            allMpbsDict[pwmList[-1]] = bedFunctions.createBedDictFromSingleFile(pwmFileNameToRead, features=[1,2,3,4,5])
            if(removeBed): os.system("rm "+bedFileName)
    else:
        if(pwmReferenceList): pwmList = pwmReferenceList
        pwmFile = open(pwmFileNameList,"r")
        for line in pwmFile:
            ll = line.strip().split("\t")
            if(ll[3] in allMpbsDict.keys()):
                if(ll[0] in allMpbsDict[ll[3]].keys()):
                    allMpbsDict[ll[3]][ll[0]].append([int(ll[1]),int(ll[2]),ll[3],int(ll[4]),ll[5]])
                else:
                    allMpbsDict[ll[3]][ll[0]] = [[int(ll[1]),int(ll[2]),ll[3],int(ll[4]),ll[5]]]
            else:
                if(not pwmReferenceList): pwmList.append(ll[3])
                allMpbsDict[ll[3]] = dict()
                allMpbsDict[ll[3]][ll[0]] = [[int(ll[1]),int(ll[2]),ll[3],int(ll[4]),ll[5]]]
        pwmFile.close()

    # Creating chromosome list
    chrList = constants.getChromList(reference=[coordDict])
    # Removing chrX, chrY and chrM
    #chrListT = []
    #for e in chrList:
    #    if(e not in ["chrX", "chrY", "chrM"]): chrListT.append(e)
    #chrList = chrListT

    # Evaluating bed additionals
    if(color == "green"): color = "0,130,0"
    elif(color == "red"): color = "130,0,0"
    elif(color == "black"): color = "0,0,0"

    # Create combinations dictionary keys
    combKeys = []
    for c in combinationList:
        for b in [",".join(e) for e in itertools.combinations(pwmList,c)]: combKeys.append(b)

    # Counting statistics
    mpbsDict = dict([(e,dict()) for e in pwmList])
    statDict = dict([(e,[0,0]) for e in combKeys]) # Left is evidence / Right is not evidence
    geneDict = dict([(e,[]) for e in combKeys])
    for chrName in coordDict.keys():

        for e in mpbsDict.keys(): mpbsDict[e][chrName] = [] # Creating chrName keys
        counter = dict([(e,0) for e in pwmList]) # Counters to iterate over all mpbs dict

        # Iterating on coordinates
        for coord in coordDict[chrName]:

            flagMotifs = dict([(e,False) for e in pwmList]) # Motifs found on this coordinate

            # Searching for MPBSs that overlapped this coordinate
            for factorName in pwmList:
                if(chrName not in allMpbsDict[factorName].keys()): continue
                while(counter[factorName] < len(allMpbsDict[factorName][chrName])):
                    currMpbs = allMpbsDict[factorName][chrName][counter[factorName]]
                    check = util.overlap(coord,currMpbs)
                    if(check == 0): # Contain overlap
                        flagMotifs[factorName] = True
                        mpbsDict[factorName][chrName].append(currMpbs+[currMpbs[0],currMpbs[1],color])
                    elif(check == -1): break # Motif is after coord
                    counter[factorName] += 1

            # Updating statistic counts and genes
            motifsFoundList = [k for k in pwmList if flagMotifs[k]]
            motifsFoundKeys = []
            motifsNotFoundKeys = [e for e in combKeys]
            for c in combinationList:
                for b in [",".join(e) for e in itertools.combinations(motifsFoundList,c)]:
                    motifsFoundKeys.append(b)
                    motifsNotFoundKeys.remove(b)
            for k in motifsFoundKeys:
                statDict[k][0] += 1
                for e in coord[2].split(":"): geneDict[k].append(e)
            for k in motifsNotFoundKeys:
                statDict[k][1] += 1

    # Remove repetitive genes from geneList
    for k in geneDict.keys(): geneDict[k] = list(set(geneDict[k]))
    
    return mpbsDict, statDict, geneDict

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


