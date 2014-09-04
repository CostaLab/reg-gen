#################################################################################################
# Genome based functions.
#################################################################################################

#################################################################################################
##### LIBRARIES #################################################################################
#################################################################################################

# External Libraries
from Bio import SeqIO
from Bio.Seq import Seq

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################

def readFastaFiles(fastaFileList):
    """Reads a fasta file and outputs a dictionary of sequences where the keys are the fasta headers.

    Keyword arguments:
    fastaFileList -- List of fasta files

    Returns:
    fastaDict -- Dictionary of sequences.
    """
    fastaDict = dict()
    for fastaFileName in fastaFileList:
        fastaFile = open(fastaFileName,"r")
        for rawSeq in SeqIO.parse(fastaFile,"fasta"):
            fastaDict[rawSeq.name] = rawSeq.seq
        fastaFile.close()
    return fastaDict


