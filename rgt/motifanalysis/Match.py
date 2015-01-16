
###################################################################################################
# Libraries
###################################################################################################

# Python
from gc import collect

# Internal
from .. Util import ErrorHandler

# External
from MOODS import search

###################################################################################################
# Functions
###################################################################################################

def match_single(motif, sequence, genomic_region, output_file, unique_threshold=None, normalize_bitscore=True):
    """
    Performs motif matching given sequence and the motif.pssm passed as parameter.
    The genomic_region is needed to evaluate the correct binding position.
    Please note that the arguments should be passed as a list, to allow for parallelization
    mapping function.

    Keyword arguments:
    motif -- TODO.
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
    output_file -- TODO.  
    unique_threshold -- If this argument is provided, the motif search will be made using a threshold of 0 and
                        then accepting only the motif matches with bitscore/motif_length >= unique_threshold.
        
    Return:
    Print MPBSs to output_file.
    """

    # Establishing threshold
    if(unique_threshold):
        current_threshold = 0.0
        eval_threshold = unique_threshold
        motif_max = motif.max / motif.len
    else:
        current_threshold = motif.threshold
        eval_threshold = motif.threshold
        motif_max = motif.max
 
    # Performing motif matching
    for search_result in search(sequence, [motif.pssm_list], current_threshold, absolute_threshold=True, both_strands=True):
        for (position, score) in search_result:

            # Verifying unique threshold acceptance
            if(unique_threshold and score/motif.len < unique_threshold): continue

            # If match forward strand
            if(position >= 0):
                p1 = genomic_region.initial + position
                strand = "+"
            # If match reverse strand
            elif(not motif.is_palindrome):
                p1 = genomic_region.initial - position
                strand = "-"
            else: continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if(normalize_bitscore): # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if(motif_max > eval_threshold):
                    norm_score = int( ( (score - eval_threshold) * 1000.0) / (motif_max - eval_threshold) )
                else: norm_score = 1000
            else: # Keep the original bitscore
                if(unique_threshold): norm_score = score/motif.len
                else: norm_score = score

            # Writing motifs
            output_file.write("\t".join([genomic_region.chrom,str(p1),str(p2),motif.name,str(norm_score),strand])+"\n")

def match_multiple(motif_name_list, motif_pssm_list, motif_thresh_list, motif_ispalindrome_list, motif_max_list, motif_len_list,
          sequence, genomic_region, output_file):
    """
    Performs motif matching given sequence and the motif.pssm passed as parameter.
    The genomic_region is needed to evaluate the correct binding position.
    Please note that the arguments should be passed as a list, to allow for parallelization
    mapping function.

    Keyword arguments:
    motif_name_list -- TODO.
    motif_pssm_list -- TODO.
    motif_thresh_list -- TODO.
    motif_ispalindrome_list -- TODO.
    motif_max_list -- TODO.    
    motif_len_list -- TODO.    
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
    output_file -- TODO.    
        
    Return:
    Print MPBSs to output_file.
    """

    # Performing motif matching
    counter = 0
    for search_result in search(sequence, motif_pssm_list, motif_thresh_list, absolute_threshold=True, both_strands=True):

        # Evaluating structures
        curr_index = int(counter/2)
        motif_name = motif_name_list[curr_index]
        motif_thresh = motif_thresh_list[curr_index]
        motif_ispalindrome = motif_ispalindrome_list[curr_index]
        motif_max = motif_max_list[curr_index]
        motif_len = motif_len_list[curr_index]
        counter += 1
        print motif_name, motif_thresh, motif_ispalindrome, motif_max, motif_len

        for (position, score) in search_result:

            # If match forward strand
            if(position >= 0):
                p1 = genomic_region.initial + position
                strand = "+"
            # If match reverse strand
            elif(not motif_ispalindrome):
                p1 = genomic_region.initial - position
                strand = "-"
            else: continue

            # Evaluating p2 and normalized score (integer between 0 and 1000 -- needed for bigbed transformation)
            p2 = p1 + motif_len
            if(motif_max > motif_thresh):
                norm_score = int( ( (score - motif_thresh) * 1000.0) / (motif_max - motif_thresh) )
            else: norm_score = 1000
            output_file.write("\t".join([genomic_region.chrom,str(p1),str(p2),motif_name,str(norm_score),strand])+"\n")

