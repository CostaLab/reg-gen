
###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings
warnings.filterwarnings("ignore")

# Internal
from .. Util import ErrorHandler

# External
from pysam import __version__ as ps_version
from pysam import Samfile
from numpy import exp, array, abs, int, mat, linalg, convolve
from scipy.stats import scoreatpercentile

"""
Processes DNase-seq and histone modification signal for
HMM footprinting input.

Authors: Eduardo G. Gusmao.
"""

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

    def __init__(self,start,end,ext):
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

        # Forward strand update
        if(not alignment.is_reverse):
            for i in range(max(alignment.pos,self.start),min(alignment.pos+self.ext,self.end-1)):
                self.vector[i-self.start] += 1.0 
        # Reverse strand update
        else:
            for i in range(max(alignment.aend-self.ext,self.start),min(alignment.aend,self.end-1)):
                self.vector[i-self.start] += 1.0 

class BamFile:
    """
    Represents a bam file. It should be used to fetch normalized and slope
    signals from a bam file.
    Usage:
    1. Initialize class.
    2. Call load_sg_coefs once.
    3. Call get_signal as many times as needed.

    Authors: Eduardo G. Gusmao.

    Methods:

    load_sg_coefs(self, slope_window_size):
    Loads Savitzky-Golay coefficients into self.sg_coefs based on a slope_window_size.

    get_signal(self, ref, start, end, ext, initial_clip = 1000, per_norm = 98, per_slope = 98)
    Gets the signal associated with self.bam based on start, end and ext.
    initial_clip, per_norm and per_slope are used as normalization factors during the normalization
    and slope evaluation procedures.

    hon_norm(self, sequence, mean, std):
    Normalizes a sequence according to hon's criterion using mean and std.
    This represents a between-dataset normalization.

    boyle_norm(self, sequence):
    Normalizes a sequence according to Boyle's criterion.
    This represents a within-dataset normalization.

    savitzky_golay_coefficients(self, window_size, order, deriv):
    Evaluate the Savitzky-Golay coefficients in order to evaluate the slope of the signal.
    It uses a window_size (of the interpolation), order (of the polynomial), deriv (derivative needed).

    slope(self, sequence, sg_coefs):
    Evaluates the slope of sequence given the sg_coefs loaded.
    """

    def __init__(self, file_name):
        """ 
        Initializes BamFile.

        Variables:
        bam -- Pysam's bam representation.
        sg_coefs -- Savitzky-Golay coefficients (list). Should be loaded after class initialization.
        """
        self.file_name = file_name
        self.bam = Samfile(file_name,"rb")
        self.sg_coefs = None

    def load_sg_coefs(self, slope_window_size):
        """ 
        Loads Savitzky-Golay coefficients into self.sg_coefs based on a slope_window_size.

        Keyword arguments:
        slope_window_size -- Window size of Savitzky-Golay coefficients.
        
        Return:
        None -- It updates self.sg_coefs.
        """
        self.sg_coefs = self.savitzky_golay_coefficients(slope_window_size, 2, 1)

    def get_signal(self, ref, start, end, ext, initial_clip = 1000, per_norm = 98, per_slope = 98):
        """ 
        Gets the signal associated with self.bam based on start, end and ext.
        initial_clip, per_norm and per_slope are used as normalization factors during the normalization
        and slope evaluation procedures.

        Keyword arguments:
        ref -- Chromosome name.
        start -- Initial genomic coordinate of signal.
        end -- Final genomic coordinate of signal.
        ext -- Fragment extention. Eg. 1 for DNase and 200 for histone modifications.
        initial_clip -- Signal will be initially clipped at this level to avoid outliers.
        per_norm -- Percentile value for 'hon_norm' function of the normalized signal.
        per_slope -- Percentile value for 'hon_norm' function of the slope signal.
        
        Return:
        hon_signal -- Normalized signal.
        slopehon_signal -- Slope signal.
        """

        # Fetch raw signal
        pileup_region = PileupRegion(start,end,ext)
        if(ps_version == "0.7.5"):
            self.bam.fetch(reference=ref, start=start, end=end, callback = pileup_region)
        else:
            iter = self.bam.fetch(reference=ref, start=start, end=end)
            for alignment in iter: pileup_region.__call__(alignment)
        raw_signal = array([min(e,initial_clip) for e in pileup_region.vector])

        # Std-based clipping
        mean = raw_signal.mean()
        std = raw_signal.std()
        clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

        # Boyle normalization (within-dataset normalization)
        boyle_signal = array(self.boyle_norm(clip_signal))

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(boyle_signal, per_norm)
        std = boyle_signal.std()
        hon_signal = self.hon_norm(boyle_signal, perc, std)
        
        # Slope signal
        slope_signal = self.slope(hon_signal, self.sg_coefs)

        # Hon normalization on slope signal (between-dataset slope smoothing)
        abs_seq = array([abs(e) for e in slope_signal])
        perc = scoreatpercentile(abs_seq, per_slope)
        std = abs_seq.std()
        slopehon_signal = self.hon_norm(slope_signal, perc, std)

        # Returning normalized and slope sequences
        return hon_signal, slopehon_signal

    def hon_norm(self, sequence, mean, std):
        """ 
        Normalizes a sequence according to hon's criterion using mean and std.
        This represents a between-dataset normalization.

        Keyword arguments:
        sequence -- Input sequence.
        mean -- Global mean.
        std -- Global std.
        
        Return:
        norm_seq -- Normalized sequence.
        """

        norm_seq = []
        for e in sequence:
            if(e == 0.0): norm_seq.append(0.0)
            elif(e > 0.0): norm_seq.append(1.0/(1.0+(exp(-(e-mean)/std))))
            else: norm_seq.append(-1.0/(1.0+(exp(-(-e-mean)/std))))
        return norm_seq

    def boyle_norm(self, sequence):
        """ 
        Normalizes a sequence according to Boyle's criterion.
        This represents a within-dataset normalization.

        Keyword arguments:
        sequence -- Input sequence.
        
        Return:
        norm_seq -- Normalized sequence.
        """

        mean = array([e for e in sequence if e>0]).mean()
        norm_seq = [(float(e)/mean) for e in sequence]
        return norm_seq

    def savitzky_golay_coefficients(self, window_size, order, deriv):
        """ 
        Evaluate the Savitzky-Golay coefficients in order to evaluate the slope of the signal.
        It uses a window_size (of the interpolation), order (of the polynomial), deriv (derivative needed).

        Keyword arguments:
        window_size -- Size of the window for function interpolation.
        order -- Order of polynomial.
        deriv -- Derivative.
        
        Return:
        m[::-1] -- The Savitzky-Golay coefficients.
        """

        # Get statistics
        #try: # TODO Errors
        window_size = abs(int(window_size))
        order = abs(int(order))
        #except ValueError, msg:
        #    raise ValueError("windowSize and order have to be of type int")
        #if windowSize % 2 != 1 or windowSize < 1:
        #    raise TypeError("windowSize size must be a positive odd number")
        #if windowSize < order + 2:
        #    raise TypeError("windowSize is too small for the polynomials order")
        order_range = range(order+1)
        half_window = (window_size -1) // 2

        # Precompute Coefficients
        b = mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
        m = linalg.pinv(b).A[deriv]
        return m[::-1]

    def slope(self, sequence, sg_coefs):
        """ 
        Evaluates the slope of sequence given the sg_coefs loaded.

        Keyword arguments:
        sequence -- Input sequence.
        sg_coefs -- Savitzky-Golay coefficients.
        
        Return:
        slope_seq -- Slope sequence.
        """
        slope_seq = convolve(sequence, sg_coefs)
        slope_seq = [e for e in slope_seq[(len(sg_coefs)/2):(len(slope_seq)-(len(sg_coefs)/2))]]
        return slope_seq


