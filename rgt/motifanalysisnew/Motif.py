
###################################################################################################
# Libraries
###################################################################################################

# Internal
from .. Util import ErrorHandler

# External
from Bio import motifs
from Bio.Seq import Seq

###################################################################################################
# Classes
###################################################################################################

class Motif:
    """
    Represent a DNA binding affinity motif.

    Authors: Eduardo G. Gusmao.

    Methods:

    match(sequence):
    Performs motif matching.
    """

    def __init__(self, input_file_name, pseudocounts=0.1, precision=10**4, fpr=0.0001):
        """ 
        Initializes Motif.

        Variables:
        pfm -- Position Frequency Matrix.
        pwm -- Position Weight Matrix.
        pssm -- Position Specific Scoring Matrix.
        threshold -- Motif matching threshold.
        """

        # Creating PFM & PWM
        input_file = open(input_file_name,"r")
        self.pfm = motifs.read(input_file, "pfm")
        self.pwm = self.pfm.counts.normalize(pseudocounts)
        input_file.close()

        # Creating PSSM
        background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
        self.pssm = self.pwm.log_odds(background)

        # Evaluating threshold
        distribution = self.pssm.distribution(background=background, precision=precision)
        self.threshold = distribution.threshold_fpr(fpr)

    def match(self, sequence):
        """ 
        Performs motif matching.

        Keyword arguments:
        sequence -- Sequence in which to find the hits.
        
        Return:
        position -- Index of hits according to Biopython.
        score -- Score of hits.
        """

        # Perform Biopython matching
        seq = Seq(sequence, self.pssm.alphabet)
        for position, score in self.pssm.search(seq, threshold=self.threshold):
            yield (position, score)


