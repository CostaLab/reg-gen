#################################################################################################
# Constants.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# Python Libraries
import os

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def getDataLocation():
    """Returns the location of the data directory on the original installation path.
    
    Returns:
    dataLocation -- Location of the data directory.
    """
    packagePathFile = open("/".join(os.path.dirname(__file__).split("/")[:-2])+"/packagePathFile.txt","r")
    dataLocation = packagePathFile.readline().strip()+"data/"
    packagePathFile.close()
    return dataLocation

def getGenome_MM9():
    """Returns the location of the mm9 genome on the original installation path.
    
    Returns:
    dataLocation -- Location of the mm9 genome.
    """
    return getDataLocation()+"mm9/genome.fa"

def getGenome_HG19():
    """Returns the location of the hg19 genome on the original installation path.
    
    Returns:
    dataLocation -- Location of the hg19 genome.
    """
    return getDataLocation()+"hg19/genome.fa"

def getAssociationFile_MM9():
    """Returns the location of the mm9 gene association file on the original installation path.
    
    Returns:
    dataLocation -- Location of the mm9 gene association file.
    """
    return getDataLocation()+"mm9/association_file.bed"

def getAssociationFile_HG19():
    """Returns the location of the hg19 gene association file on the original installation path.
    
    Returns:
    dataLocation -- Location of the hg19 gene association file.
    """
    return getDataLocation()+"hg19/association_file.bed"

def getChromSizes_MM9():
    """Returns the location of the mm9 chrom.sizes file on the original installation path.
    
    Returns:
    dataLocation -- Location of the mm9 chrom.sizes file.
    """
    return getDataLocation()+"mm9/chrom.sizes"

def getChromSizes_HG19():
    """Returns the location of the hg19 chrom.sizes file on the original installation path.
    
    Returns:
    dataLocation -- Location of the hg19 chrom.sizes file.
    """
    return getDataLocation()+"hg19/chrom.sizes"

def getPwmFolder():
    """Returns the location of PWM folder on the original installation path.
    
    Returns:
    dataLocation -- Location of the PWM folder.
    """
    return getDataLocation()+"motifs/jaspar_vertebrates/"

def getLogoFolder():
    """Returns the location of the logo folder on the original installation path.
    
    Returns:
    dataLocation -- Location of the logo folder.
    """
    return getDataLocation()+"logos/jaspar_vertebrates/"

def getPrecompFdr4Folder():
    """Returns the location of the pre-computed motif matches using an FDR = 10^(-4) on the original installation path.
    
    Returns:
    dataLocation -- Location of the pre-computed folder.
    """
    return getDataLocation()+"PrecomputedMotifs_fdr4/"

def getChromList(x=True, y=True, m=True, nChrom=22, reference=[]):
    """Returns a chromosome aliases list.

    Keyword arguments:
    x -- Wether the chrX will be present or not. (default True)
    y -- Wether the chrY will be present or not. (default True)
    m -- Wether the chrM will be present or not. (default True)
    nChrom -- Number of chromosomes. (default 22)
    reference -- List of dictionaries. The returned chromList will only contain entries that appear on any of these dictionaries.

    Returns:
    chromList -- List of chromosome aliases.
    """

    if(not reference): # Creating basic chromosome list
        chromList = ["chr"+str(e) for e in range(1,nChrom+1)]
        if(x): chromList.append("chrX")
        if(y): chromList.append("chrY")
        if(m): chromList.append("chrM")
    else: # Filtering by reference
        chromList = []
        for n in [str(e) for e in range(1,nChrom+1)] + ["X","Y","M"]:
            appears = False
            for d in reference:
                if("chr"+n in d.keys()):
                    appears = True
                    break
            if(appears): chromList.append("chr"+n)
    return chromList


