
###################################################################################################
# Libraries
###################################################################################################

# Python
from os.path import basename

# Internal
from .. Util import ErrorHandler

# External
from Bio import motifs

###################################################################################################
# Classes
###################################################################################################

class Motif:
    """
    Represent a DNA binding affinity motif.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, (input_file_name, pseudocounts, precision, fpr)):
        """ 
        Initializes Motif.

        Variables:
        pfm -- Position Frequency Matrix.
        pwm -- Position Weight Matrix.
        pssm -- Position Specific Scoring Matrix.
        threshold -- Motif matching threshold.
        len -- Length of the motif.
        max -- Maximum PSSM score possible.
        is_palindrome -- True if consensus is biologically palindromic.
        """

        # Initializing error handler
        main_error_handler = ErrorHandler()       
 
        # Initializing name
        self.name = ".".join(basename(input_file_name).split(".")[:-1])

        # Creating PFM & PWM
        input_file = open(input_file_name,"r")
        self.pfm = motifs.read(input_file, "pfm")
        self.pwm = self.pfm.counts.normalize(pseudocounts)
        input_file.close()
        self.len = len(self.pfm)

        # Creating PSSM
        background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
        self.pssm = self.pwm.log_odds(background)
        self.max = self.pssm.max

        # Evaluating threshold
        try: distribution = self.pssm.distribution(background=background, precision=precision)
        except Exception: main_error_handler.throw_error("MM_PSEUDOCOUNT_0")
        self.threshold = distribution.threshold_fpr(fpr)

        # Evaluating if motf is palindromic
        if(str(self.pfm.consensus) == str(self.pfm.consensus.reverse_complement())): self.is_palindrome = True
        else: self.is_palindrome = False


