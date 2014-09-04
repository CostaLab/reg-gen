#################################################################################################
# Auxiliary miscelaneous functions
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################


#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def instanceof(s,instance="int"):
    """Verifies if a string is an instance of other primary type.

    Keyword arguments:
    s -- String.
    instance -- the instance type. (default "int")

    Returns:
    True if string s is an instance of the given type and False if it is not.
    """
    try:
        if(instance == "int"): int(s)
        elif(instance == "float"): float(s)
        return True
    except (TypeError, ValueError):
        return False

def readList(fileName):
    """Reads a file of items split by line break and returns them in a list.

    Keyword arguments:
    fileName -- The path + name of the file.

    Returns:
    fileList -- List of items.
    """
    fileList = []
    inFile = open(fileName,"r")
    for line in inFile: fileList.append(line.strip())
    inFile.close()
    return fileList

def overlap(t1, t2):
    """Checks if one interval contains any overlap with another interval.

    Keyword arguments:
    t1 -- First tuple.
    t2 -- Second tuple.

    Returns:
    Returns -1 if i1 is before i2; 1 if i1 is after i2; and 0 if there is any overlap.
    """
    if(t1[1] <= t2[0]): return -1 # interval1 is before interval2
    if(t2[1] <= t1[0]): return 1 # interval1 is after interval2
    return 0 # interval1 overlaps interval2

def readInputParameters(inputList):
    """Reads input parameters and creates a dictionary of inputs.

    Keyword arguments:
    inputList -- sys.argv list of parameters.

    Returns:
    Returns dictionary of parameters.
    """
    paramDict = dict()
    for e in inputList:
        eSplit = e.split("=")
        paramDict[eSplit[0]] = eSplit[1]
    return paramDict


