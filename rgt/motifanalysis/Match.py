
###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function

# Internal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion

# External
try:
    import MOODS.tools
    import MOODS.scan
except:
    import MOODS


###################################################################################################
# Functions
###################################################################################################

def match_single(motif, sequence, genomic_region, unique_threshold=None, normalize_bitscore=True, sort=False):
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
    if unique_threshold:
        current_threshold = 0.0
        eval_threshold = unique_threshold
        motif_max = motif.max / motif.len
    else:
        current_threshold = motif.threshold
        eval_threshold = motif.threshold
        motif_max = motif.max

    # Performing motif matching
    try:
        # old MOODS version
        results = MOODS.search(sequence, [motif.pssm_list], current_threshold,
                               absolute_threshold=True, both_strands=True)
    except:
        # TODO: we can expand this to use bg from sequence, for example,
        # or from organism.
        bg = MOODS.tools.flat_bg(4)
        results = MOODS.scan.scan_dna(sequence, [motif.pssm_list], bg, [current_threshold], 7)

    grs = GenomicRegionSet("mpbs")

    for search_result in results:
        for r in search_result:
            try:
                position = r.pos
                score = r.score
            except:
                (position, score) = r

            # Verifying unique threshold acceptance
            if unique_threshold and score/motif.len < unique_threshold:
                continue

            # If match forward strand
            if position >= 0:
                p1 = genomic_region.initial + position
                strand = "+"
            # If match reverse strand
            elif not motif.is_palindrome:
                p1 = genomic_region.initial - position
                strand = "-"
            else:
                continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if normalize_bitscore:
                # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if motif_max > eval_threshold:
                    norm_score = int(((score - eval_threshold) * 1000.0) / (motif_max - eval_threshold))
                else:
                    norm_score = 1000
            else:
                # Keep the original bitscore
                if unique_threshold:
                    norm_score = score/motif.len
                else:
                    norm_score = score

            grs.add(GenomicRegion(genomic_region.chrom, int(p1), int(p2),
                                  name=motif.name, orientation=strand, data=str(norm_score)))

    if sort:
        grs.sort()

    return grs
