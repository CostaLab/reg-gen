#################################################################################################
# Gene - coordinate associations.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Distal Libraries
from .. util import *

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def geneAssociationByPromoter(coordDict,geneList,geneAssocLocation,chromSizesLocation,promoterLength=1000,threshDist=50000):
    """Associates coordinates to genes given the following rules:
       1. If the peak is inside gene (promoter+coding) then this peak is associated with that gene.
       2. If a peak is inside overlapping genes, then the peak is annotated with both genes.
       3. If peak is between two genes (not overlapping neither), then both genes are annotated.
       4. If the distance between peak and gene is greater than a threshold distance, then it is not annotated.

    Keyword arguments:
    coordDict -- Dictionary of coordinates.
    geneList -- List of gene names. If None, then consider all genes to be enriched.
    geneAssocLocation -- Location + name of the file that contains the location of the genes.
    chromSizesLocation -- Location + name of the chromosome sizes file.
    promoterLength -- Length of the promoter region. (default 1000)
    threshDist -- Threshold maximum distance for a coordinate to be considered associated with a gene. (default 50000)

    Returns:
    aDictFinal -- Dictionary containing the associated bed entries.
    aDict -- Dictionary containing the associated bed entries with '_PROX' nad '_DIST' labels (for writing).
    """

    # Reading assocDict
    assocDictTemp = bedFunctions.createBedDictFromSingleFile(geneAssocLocation, features=[1,2,3,4,5])
    assocDict = dict()
    geneFlagDict = dict()
    for k in assocDictTemp.keys():
        assocDict[k] = []
        for e in assocDictTemp[k]:
            if(len(e[2]) < 2): continue
            geneFlagDict[e[2].upper()] = False
            assocDict[k].append([e[0],e[1],e[2].upper(),e[3],e[4]])
    assocDict = sort.sortBedDictionary(assocDict, field=0)

    # Updating geneFlagDict based on geneList
    if(geneList):
        geneList=[g.upper() for g in geneList]
        for e in geneList: geneFlagDict[e] = True
    else:
        for e in geneFlagDict.keys(): geneFlagDict[e] = True

    # Updating assocDict for keys in coordDict
    for k in coordDict.keys():
        if(k not in assocDict.keys()): assocDict[k] = []

    # Associating coordDict with assocDict
    aDict = dict() # Results dictionary
    for chrName in coordDict.keys():

        aDict[chrName] = [] # Adding new result list
        counter = 0 # Counter to run through gene list

        # Iterating on coordinate list (main list)
        for coord in coordDict[chrName]: 

            didBreak = False # Check wether it breaked or not.

            # Running linearly through gene list
            while(counter < len(assocDict[chrName])):

                # Extend the gene coordinates to encompass promoter region based on strandness
                if(assocDict[chrName][counter][4] == "+"): geneCoord = [assocDict[chrName][counter][0]-promoterLength,assocDict[chrName][counter][1]]
                else: geneCoord = [assocDict[chrName][counter][0],assocDict[chrName][counter][1]+promoterLength]

                check = util.overlap(coord,geneCoord) # Check overlap between coordinate and gene+promoter

                if(check == 0): # If contain overlap, then check if the coordinate also contains overlap with the next gene

                    # Verify if next gene+promoter also overlap
                    genesList = [assocDict[chrName][counter][2]+"_PROX"] # List of overlapping genes
                    if(counter < len(assocDict[chrName])-1): # If this is not the last gene (there is a 'next') then verify overlap
                        if(assocDict[chrName][counter+1][4] == "+"): geneCoordNext = [assocDict[chrName][counter+1][0]-promoterLength,assocDict[chrName][counter+1][1]]
                        else: geneCoordNext = [assocDict[chrName][counter+1][0],assocDict[chrName][counter+1][1]+promoterLength]
                        checkNext = util.overlap(coord,geneCoordNext) # Verify overlap between coordinate and next gene
                        if(checkNext == 0): genesList.append(assocDict[chrName][counter+1][2]+"_PROX") # If there is an overlap then add this gene to association list
                    #else: # If this is the last gene (there is no 'next') then dont verify overlap

                    # Verify if genes are enriched
                    for i in range(0,len(genesList)): # Check if genes in genesList are enriched 
                        if(not geneFlagDict[genesList[i][:-5]]): genesList[i] = "."+genesList[i] # If the gene is not enriched, put a dot (.) in front of it (will be used by motif matching)
                        #else: If gene is enriched than let the name without the dot (.) in front of it
                    didBreak = True
                    break

                elif(check == -1): # If gene is after coordinate, then check if they are within maximum distance for overlap. Also do the same check for the previous gene.

                    # Verify overlap again using maximum distance with current gene+promoter
                    genesList = [] # List of overlapping genes
                    maxGeneCoord = [geneCoord[0]-threshDist,geneCoord[1]+threshDist]
                    maxCheck = util.overlap(coord,maxGeneCoord)
                    if(maxCheck == 0): genesList.append(assocDict[chrName][counter][2]+"_DIST") # If it overlapped then put current gene in overlap list
                    
                    # Verify overlap again using maximum distance with previous gene+promoter
                    if(counter > 0): # Do this verification only if this is not the first gene
                        if(assocDict[chrName][counter-1][4] == "+"): geneCoordPrev = [assocDict[chrName][counter-1][0]-promoterLength,assocDict[chrName][counter-1][1]]
                        else: geneCoordPrev = [assocDict[chrName][counter-1][0],assocDict[chrName][counter-1][1]+promoterLength]
                        maxGeneCoordPrev = [geneCoordPrev[0]-threshDist,geneCoordPrev[1]+threshDist]
                        maxCheckPrev = util.overlap(coord,maxGeneCoordPrev)
                        if(maxCheckPrev == 0): genesList.append(assocDict[chrName][counter-1][2]+"_DIST") # If it overlapped then put previous gene in overlap list

                    # Verify if genes are enriched
                    if(len(genesList) == 0): genesList.append(".") # If genesList is empty then put a '.' to represent non-association
                    else: # If genesList is not empty then verify enriched genes
                        for i in range(0,len(genesList)): # Check if genes in genesList are enriched 
                            if(not geneFlagDict[genesList[i][:-5]]): genesList[i] = "."+genesList[i] # If the gene is not enriched, put a dot (.) in front of it (will be used by motif matching)
                            #else: If gene is enriched than let the name without the dot (.) in front of it
                    didBreak = True
                    break

                #elif(check == -1): If gene is before coordinate, dont do anything! Just move to the next gene.
                counter += 1

            # If not breaked, then the gene list is over and the coordinate can only overlap the last gene.
            if(not didBreak): 
                genesList = []
                if(len(assocDict[chrName]) == 0): # If gene list is empty then dont associate coordinate with any gene
                    genesList.append(".")
                else: # If gene list is not empty then try to verify overlap between coordinate and last gene
                    if(assocDict[chrName][counter-1][4] == "+"): geneCoordPrev = [assocDict[chrName][counter-1][0]-promoterLength,assocDict[chrName][counter-1][1]]
                    else: geneCoordPrev = [assocDict[chrName][counter-1][0],assocDict[chrName][counter-1][1]+promoterLength]
                    maxGeneCoordPrev = [geneCoordPrev[0]-threshDist,geneCoordPrev[1]+threshDist]
                    maxCheckPrev = util.overlap(coord,maxGeneCoordPrev)
                    if(maxCheckPrev == 0): # If it overlapped then put previous gene in overlap list and check for enrichment
                        genesList.append(assocDict[chrName][counter-1][2]+"_DIST")
                        if(not geneFlagDict[genesList[0][:-5]]): genesList[0] = "."+genesList[0]
                    else: genesList.append(".") # If does not overlap then put a '.' to represent non-association

            # Write the curent coordinate with its corresponding overlapping genes (enriched or not)
            aDict[chrName].append(coord[:5]) # Writing the raw coordinate until STRAND field only
            aDict[chrName][-1][2] = ":".join(genesList) # Write list of overlapping genes
            aDict[chrName][-1][3] = 0 # Force the SCORE field to be '0'
            aDict[chrName][-1][4] = "." # Force the STRAND field to be '.'
    # End for chrName in coordDict.keys()

    # Removing proximity information
    aDictFinal = removeProximityInformation(aDict)

    return aDictFinal, aDict

def removeProximityInformation(aDict):
    """Removes _PROX and _DIST from an association dictionary.

    Keyword arguments:
    aDict -- Association (bed) dictionary.

    Returns:
    resDict -- Resulting dictionary.
    """

    # Removing proximity information
    resDict = dict()
    for chrName in constants.getChromList(reference=[aDict]):
        resDict[chrName] = []
        for e in aDict[chrName]:
            newGeneNames = []
            for g in e[2].split(":"):
                if("_PROX" in g or "_DIST" in g): newGeneNames.append(g[:-5])
                else: newGeneNames.append(g)
            newGeneNames = ":".join(newGeneNames)
            resDict[chrName].append([e[0],e[1],newGeneNames,e[3],e[4]])
    return resDict


