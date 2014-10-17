
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

    def __init__(self, input_file_name, pseudocounts, precision, fpr, thresholds):
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
        repository = input_file_name.split("/")[-2]

        # Creating PFM & PWM
        input_file = open(input_file_name,"r")
        self.pfm = motifs.read(input_file, "pfm")
        self.pwm = self.pfm.counts.normalize(pseudocounts)
        input_file.close()
        self.len = len(self.pfm)

        # Creating PSSM
        background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
        self.pssm = self.pwm.log_odds(background)
        self.pssm_list = [self.pssm[e] for e in ["A","C","G","T"]]
        self.max = self.pssm.max

        # Evaluating threshold
        try:
            if(pseudocounts != 0.1 or precision != 10000): 1/0 # Induce Exception
            self.threshold = thresholds.dict[repository][self.name][fpr]
        except Exception:
            try: distribution = self.pssm.distribution(background=background, precision=precision)
            except Exception: main_error_handler.throw_error("MM_PSEUDOCOUNT_0")
            self.threshold = distribution.threshold_fpr(fpr)

        # Evaluating if motf is palindromic
        if(str(self.pfm.consensus) == str(self.pfm.consensus.reverse_complement())): self.is_palindrome = True
        else: self.is_palindrome = False


class Thresholds:
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

            # Initializing dictionary level 1
            fpr_name = ".".join(fpr_file_name.split("/")[-1].split(".")[:-1])
            self.dict[fpr_name] = dict()

            # Iterating in fpr file
            fpr_file = open(fpr_file_name,"r")
            header = fpr_file.readline()
            fpr_values = [float(e) for e in header.strip().split("\t")[1:]]
            for line in fpr_file:
                ll = line.strip().split("\t")
                self.dict[fpr_name][ll[0]] = dict()# Initializing dictionary level 2
                for i in range(1,len(ll)): self.dict[fpr_name][ll[0]][fpr_values[i-1]] = float(ll[i]) # Updating dictionary
            fpr_file.close()


