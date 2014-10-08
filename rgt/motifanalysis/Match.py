
###################################################################################################
# Libraries
###################################################################################################

# Internal
from .. Util import ErrorHandler
from .. GenomicRegion import GenomicRegion

# External
from MOODS import search

###################################################################################################
# Functions
###################################################################################################

def match((motif, sequence, genomic_region)):
    """ 
    Performs motif matching given sequence and the motif.pssm passed as parameter.
    The genomic_region is needed to evaluate the correct binding position.
    Please note that the arguments should be passed as a list, to allow for parallelization
    mapping function.

    Keyword arguments:
    motif -- A Motif.
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
        
    Return:
    genomic_region_list -- A list with MPBSs represented as GenomicRegions
    """

    # Initializing MPBS genomic region list
    genomic_region_list = []

    # Performing motif matching
    for search_result in search(sequence, [motif.pssm_list], motif.threshold, absolute_threshold=True, both_strands=True):
        for (position, score) in search_result:

            # If match forward strand
            if(position >= 0):
                p1 = genomic_region.initial + position
                strand = "+"
            # If match reverse strand
            elif(not motif.is_palindrome):
                p1 = genomic_region.initial - position
                strand = "-"
            else: continue

            # Evaluating p2 and normalized score (integer between 0 and 1000 -- needed for bigbed transformation)
            p2 = p1 + motif.len
            if(motif.max > motif.threshold):
                norm_score = int( ( (score - motif.threshold) * 1000.0) / (motif.max - motif.threshold) )
            else: norm_score = 1000
            genomic_region_list.append(GenomicRegion(genomic_region.chrom,p1,p2,name=motif.name,orientation=strand,data=norm_score))

    return genomic_region_list


