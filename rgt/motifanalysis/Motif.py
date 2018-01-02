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


# Internal


###################################################################################################
# Classes
###################################################################################################


class Motif:
    """
    Represent a DNA binding affinity motif.
    """

    def __init__(self, input_file_name, pseudocounts, fpr, thresholds):
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
        repository = input_file_name.split("/")[-2]

        # Creating PFM & PSSM
        self.pfm = parsers.pfm(str(input_file_name))
        self.bg = tools.flat_bg(len(self.pfm))  # total number of "points" to add, not per-row
        self.pssm = tools.log_odds(self.pfm, self.bg, pseudocounts, 2)

        # how many bases this motif has
        self.len = len(self.pfm[0])

        # maximum value found in the whole PSSM
        self.max = max([max(e) for e in self.pssm])

        # we only support pure DNA or methylated DNA, for now.
        self.alphabet = ["Aa", "Cc", "Gg", "Tt"]
        if len(self.pfm) == 6:
            self.alphabet += ["m", "1"]

        # Evaluating threshold
        try:
            if pseudocounts != 1.0:
                raise ValueError()
            self.threshold = thresholds.dict[repository][self.name][fpr]
        except Exception:
            # FIXME: this requires a modified version of MOODS. Not sure if we actually need it.
            # self.threshold = tools.threshold_from_p(self.pssm, self.bg, fpr, 2000.0)  # 10000.0 would take too long
            self.threshold = tools.threshold_from_p(self.pssm, self.bg, fpr)
            print(">>> recomputing threshold for %s: %f" % (self.name, self.threshold))

        # Evaluating if motif is palindromic
        self.is_palindrome = [max(e) for e in self.pssm] == [max(e) for e in reversed(self.pssm)]


class ThresholdTable:
    """
    Container for all motif thresholds given default FPRs.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, motif_data):
        """ 
        Initializes Thresholds. Motif thresholds are stored in a dictionary:
        [repository] -> [motif name] -> [fpr] -> threshold float

        Parameters:
        motif_data -- MotifData object.

        Variables:
        
        """

        # Initializing dictionary level 0
        self.dict = dict()

        # Iterating over fpr files
        for fpr_file_name in motif_data.get_fpr_list():

            if not fpr_file_name:
                continue

            # Initializing dictionary level 1
            fpr_name = ".".join(fpr_file_name.split("/")[-1].split(".")[:-1])
            self.dict[fpr_name] = dict()

            # Iterating in fpr file
            fpr_file = open(fpr_file_name, "r")
            header = fpr_file.readline()
            fpr_values = [float(e) for e in header.strip().split("\t")[1:]]
            for line in fpr_file:
                ll = line.strip().split("\t")
                # Initializing dictionary level 2
                self.dict[fpr_name][ll[0]] = dict()
                for i in range(1, len(ll)):
                    # Updating dictionary
                    self.dict[fpr_name][ll[0]][fpr_values[i - 1]] = float(ll[i])
            fpr_file.close()
