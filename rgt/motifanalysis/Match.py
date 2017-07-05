
###################################################################################################
# Libraries
###################################################################################################

# Python
from __future__ import print_function

# Internal
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.GenomicRegion import GenomicRegion

# External
import MOODS.tools
import MOODS.scan


###################################################################################################
# Functions
###################################################################################################

def match_single(motif, sequence, genomic_region, unique_threshold=None, normalize_bitscore=False, output=None):
    """
    Performs motif matching given sequence and the motif.pssm passed as parameter.
    The genomic_region is needed to evaluate the correct binding position.

    Keyword arguments:
    motif -- TODO.
    sequence -- A DNA sequence (string).
    genomic_region -- A GenomicRegion.
    unique_threshold -- If this argument is provided, the motif search will be made using a threshold of 0 and
                        then accepting only the motif matches with bitscore/motif_length >= unique_threshold.
    normalize_bitscore -- If True, it normalises the scores between 0 and 1000. Necessary for bigbed conversion.
    output -- A GenomicRegionSet where all matching GenomicRegions will be appended.
        
    Return:
    Either the "output" GenomicRegionSet if provided, or a newly-instantiated one.
    """

    # Establishing threshold
    if unique_threshold:
        current_threshold = 0.0
        threshold = unique_threshold
        motif_max = motif.max / motif.len
    else:
        current_threshold = motif.threshold
        threshold = motif.threshold
        motif_max = motif.max

    # Performing motif matching
    # TODO: we can expand this to use bg from sequence, for example,
    # or from organism.
    bg = MOODS.tools.flat_bg(4)
    results = MOODS.scan.scan_dna(sequence, [motif.pssm_list], bg, [current_threshold], 7)

    if output is None:
        output = GenomicRegionSet("mpbs")

    pos_start = genomic_region.initial
    chrom = genomic_region.chrom

    for search_result in results:
        for r in search_result:
            position = r.pos
            score = r.score

            score_len = score / motif.len

            # Verifying unique threshold acceptance
            if unique_threshold and score_len < unique_threshold:
                continue

            # If match forward strand
            if position >= 0:
                p1 = pos_start + position
                strand = "+"
            # If match reverse strand
            elif not motif.is_palindrome:
                p1 = pos_start - position
                strand = "-"
            else:
                continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if normalize_bitscore:
                # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if motif_max > threshold:
                    score = int(((score - threshold) * 1000.0) / (motif_max - threshold))
                else:
                    score = 1000
            elif unique_threshold:
                score = score_len

            output.add(GenomicRegion(chrom, int(p1), int(p2), name=motif.name, orientation=strand, data=str(score)))

    return output


# TODO must add normalisation stuff (needed?)
# only small speed boost, no memory boost
def match_multiple(motifs, sequence, genomic_region, unique_threshold=None, normalize_bitscore=False, output=None):
    """
    More efficient than calling match_single on every motif.

    Keyword arguments:
    motif -- TODO.
    sequence -- A DNA sequence (string).
        genomic_region -- A GenomicRegion.
    unique_threshold -- If this argument is provided, the motif search will be made using a threshold of 0 and
                        then accepting only the motif matches with bitscore/motif_length >= unique_threshold.
    normalize_bitscore -- If True, it normalises the scores between 0 and 1000. Necessary for bigbed conversion.
    output -- A GenomicRegionSet where all matching GenomicRegions will be appended.

    Return:
    Either the "output" GenomicRegionSet if provided, or a newly-instantiated one.
    """

    pssm_lists = []
    thresholds = []
    for motif in motifs:
        if unique_threshold:
            thresholds.append(0.0)
        else:
            thresholds.append(motif.threshold)
        pssm_lists.append(motif.pssm_list)

    # Performing motif matching
    # TODO: we can expand this to use bg from sequence, for example,
    # or from organism.
    bg = MOODS.tools.flat_bg(4)
    results = MOODS.scan.scan_dna(sequence, pssm_lists, bg, thresholds, 7)

    if output is None:
        output = GenomicRegionSet("mpbs")

    pos_start = genomic_region.initial
    chrom = genomic_region.chrom

    for i, search_result in enumerate(results):
        motif = motifs[i]

        if unique_threshold:
            motif_max = motif.max / motif.len
            threshold = unique_threshold
        else:
            motif_max = motif.max
            threshold = motif.threshold

        for r in search_result:
            position = r.pos
            score = r.score

            score_len = score / motif.len

            # Verifying unique threshold acceptance
            if unique_threshold and score_len < unique_threshold:
                continue

            # If match forward strand
            if position >= 0:
                p1 = pos_start + position
                strand = "+"
            # If match reverse strand
            elif not motif.is_palindrome:
                p1 = pos_start - position
                strand = "-"
            else:
                continue

            # Evaluating p2
            p2 = p1 + motif.len

            # Evaluating score (integer between 0 and 1000 -- needed for bigbed transformation)
            if normalize_bitscore:
                # Normalized bitscore = standardize to integer between 0 and 1000 (needed for bigbed transformation)
                if motif_max > threshold:
                    score = int(((score - threshold) * 1000.0) / (motif_max - threshold))
                else:
                    score = 1000
            elif unique_threshold:
                score = score_len

            output.add(GenomicRegion(chrom, int(p1), int(p2), name=motif.name, orientation=strand, data=str(score)))

    return output
