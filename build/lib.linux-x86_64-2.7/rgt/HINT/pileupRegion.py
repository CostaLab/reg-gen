
###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings
warnings.filterwarnings("ignore")

# Internal
from .. Util import ErrorHandler

###################################################################################################
# Classes
###################################################################################################

class PileupRegion:
    """
    Represent an region in which a fragment pileup will be calculated.
    It is passed to pysam's functions (such as 'fetch') and then __call__
    gets called for every instance of the iterator that pysam's functions
    return.

    Authors: Eduardo G. Gusmao.

    Methods:

    __call__(alignment):
    It is called for every 'alignment' instance found by pysam's 'fetch' method.
    It loads self.vector with such alignment based on self.ext.
    """

    def __init__(self,start,end,ext,left_ext=0,right_ext=0,f_shift=0,r_shift=0):
        """ 
        Initializes PileupRegion.

        Variables:
        start -- Start point (in bp) of this PileupRegion (integer).
        end -- End point (in bp) of this PileupRegion (integer).
        length -- Length of this PileupRegion (integer).
        ext -- Number of bps to extend each fragment found by Pysam's 'fetch' method.
               The fragment is extended from its 5' position to the downstream
               direction (integer).
        vector -- Pileup vector. Each extended fragment will contribute with +1 to each
                  position of this vector, in which the length equals self.length (list).
        """
        self.start = start
        self.end = end
        self.length = end-start
        self.ext = ext
        self.left_ext = left_ext
        self.right_ext = right_ext
        self.f_shift = f_shift
        self.r_shift = r_shift
        self.vector = [0.0] * self.length

    def __call__(self, alignment):
        """ 
        Updates self.vector with the fragment obtained with the 'fetch' function from
        pysam package. This function will be called for all fragments fetched.

        Keyword arguments:
        alignment -- File location + name in HMM format. See manual for full description of such format.
        
        Return:
        None -- It updates self.vector only.
        """

        try:
            # Forward strand update
            if(not alignment.is_reverse):
                for i in range(max(alignment.pos,self.start),min(alignment.pos+self.ext,self.end-1)):
                    self.vector[i-self.start] += 1.0 
            # Reverse strand update
            else:
                for i in range(max(alignment.aend-self.ext,self.start),min(alignment.aend,self.end-1)):
                    self.vector[i-self.start] += 1.0 
        except Exception: pass

    def __call2__(self, alignment):
        try:
            if(not alignment.is_reverse):
                for i in range(max(alignment.pos-self.ext,self.start),min(alignment.pos+self.ext,self.end-1)):
                    self.vector[i-self.start] += 1.0 
            else:
                for i in range(max(alignment.aend-self.ext,self.start),min(alignment.aend+self.ext,self.end-1)):
                    self.vector[i-self.start] += 1.0
        except Exception: pass


