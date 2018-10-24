###################################################################################################
# Libraries
###################################################################################################

from __future__ import division
# Python 3 compatibility
from __future__ import print_function

# Python
from os.path import basename

# External
from MOODS import tools, parsers
from numpy import argmax


# Internal


###################################################################################################
# Classes
###################################################################################################


class Motif:
    """
    Represent a DNA binding affinity motif.
    """

    def __init__(self, input_file_name, pseudocounts, threshold):
        """ 
        Initializes Motif.

        Fields:
        pfm -- Position Frequency Matrix.
        bg -- Background frequencies.
        pssm -- Position Specific Scoring Matrix.
        alphabet -- A list of letters, eg ["Aa", "Cc", "Gg", "Tt"]
        threshold -- Motif matching threshold.
        len -- Length of the motif.
        max -- Maximum PSSM score possible.
        is_palindrome -- True if consensus is biologically palindromic.
        """

        # Initializing name
        self.name = ".".join(basename(input_file_name).split(".")[:-1])
        #repository = input_file_name.split("/")[-2]

        # Creating PFM & PSSM
        self.pfm = parsers.pfm(str(input_file_name))
        self.bg = tools.flat_bg(len(self.pfm))  # total number of "points" to add, not per-row
        self.pssm = tools.log_odds(self.pfm, self.bg, pseudocounts, 2)
        self.pssm_rc = tools.reverse_complement(self.pssm)

        # how many bases this motif has
        self.len = len(self.pfm[0])

        # maximum value found in the whole PSSM
        self.max = max([max(e) for e in self.pssm])

        # we only support pure DNA or methylated DNA, for now.
        self.alphabet = ["Aa", "Cc", "Gg", "Tt"]
        if len(self.pfm) == 6:
            self.alphabet += ["m", "1"]

        self.threshold = threshold

        # Evaluating threshold
        # try:
        #     if pseudocounts != 1.0:
        #         raise ValueError()
        #     self.threshold = [fpr]
        # except Exception:
        #     # FIXME: this requires a modified version of MOODS. Not sure if we actually need it.
        #     # self.threshold = tools.threshold_from_p(self.pssm, self.bg, fpr, 2000.0)  # 10000.0 would take too long
        #     self.threshold = tools.threshold_from_p(self.pssm, self.bg, fpr)
        #     print(">>> recomputing threshold for %s: %f" % (self.name, self.threshold))

        #self.threshold_rc = tools.threshold_from_p(self.pssm_rc, self.bg, fpr)

        self.consensus = "".join([self.alphabet[i][0] for i in argmax(self.pssm, axis=0)])
        self.consensus_rc = "".join([self.alphabet[i][0] for i in argmax(self.pssm_rc, axis=0)])

        # Evaluating if motif is palindromic
        self.is_palindrome = self.consensus == self.consensus_rc
