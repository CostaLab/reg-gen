
###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings
warnings.filterwarnings("ignore")
from math import log

# Internal
from .. Util import ErrorHandler
from .. Util import AuxiliaryFunctions
from pileupRegion import PileupRegion

# External
from pysam import __version__ as ps_version
from pysam import Samfile
from pysam import Fastafile
from numpy import exp, array, abs, int, mat, linalg, convolve
from scipy.stats import scoreatpercentile
from ngslib import BigWigFile

"""
Processes DNase-seq and histone modification signal for
HMM footprinting input.

Authors: Eduardo G. Gusmao.
"""

class GenomicSignal:
    """
    Represents a genomic signal. It should be used to fetch normalized and slope
    signals from a bam or bw file.
    Usage:
    1. Initialize class.
    2. Call load_sg_coefs once.
    3. Call get_signal as many times as needed.

    Authors: Eduardo G. Gusmao.

    Methods:

    load_sg_coefs(self, slope_window_size):
    Loads Savitzky-Golay coefficients into self.sg_coefs based on a slope_window_size.

    get_signal(self, ref, start, end, ext, initial_clip = 1000, per_norm = 98, per_slope = 98)
    Gets the signal associated with self.bam or self.bw based on start, end and ext.
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
        Initializes GenomicSignal.
        """
        self.file_name = file_name
        self.bam = None
        self.bw = None
        self.sg_coefs = None
        self.is_bam = False
        self.is_bw = False
        if(self.file_name.split(".")[-1].upper() == "BAM"):
            self.is_bam = True
            self.bam = Samfile(file_name,"rb")
        elif(self.file_name.split(".")[-1].upper() == "BW" or self.file_name.split(".")[-1].upper() == "BIGWIG"):
            self.is_bw = True
            self.bw = BigWigFile(file_name)
        else: pass # TODO ERROR

    def load_sg_coefs(self, slope_window_size):
        """ 
        Loads Savitzky-Golay coefficients into self.sg_coefs based on a slope_window_size.

        Keyword arguments:
        slope_window_size -- Window size of Savitzky-Golay coefficients.
        
        Return:
        None -- It updates self.sg_coefs.
        """
        self.sg_coefs = self.savitzky_golay_coefficients(slope_window_size, 2, 1)

    def get_tag_count(self, ref, start, end, ext, initial_clip = 1000, ext_both_directions=False):
        """ 
        Gets the tag count associated with self.bam based on start, end and ext.

        Keyword arguments:
        ref -- Chromosome name.
        start -- Initial genomic coordinate of signal.
        end -- Final genomic coordinate of signal.
        ext -- Fragment extention. Eg. 1 for DNase and 200 for histone modifications.
        initial_clip -- Signal will be initially clipped at this level to avoid outliers.
        
        Return:
        tag_count -- Total signal.
        """

        # Fetch raw signal
        pileup_region = PileupRegion(start,end,ext)
        if(self.is_bam):
            if(ps_version == "0.7.5"):
                self.bam.fetch(reference=ref, start=start, end=end, callback = pileup_region)
            else:
                iter = self.bam.fetch(reference=ref, start=start, end=end)
                if(not ext_both_directions):
                    for alignment in iter: pileup_region.__call__(alignment)
                else:
                    for alignment in iter: pileup_region.__call2__(alignment)
            raw_signal = array([min(e,initial_clip) for e in pileup_region.vector])
        elif(self.is_bw):
            signal = self.bw.pileup(ref, start, end)
            raw_signal = array([min(e,initial_clip) for e in signal])

        # Std-based clipping
        mean = raw_signal.mean()
        std = raw_signal.std()
        clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

        # Tag count
        try: tag_count = sum(clip_signal)
        except Exception: tag_count = 0

        return tag_count

    def get_signal(self, ref, start, end, ext, initial_clip = 1000, per_norm = 99.5, per_slope = 98, 
                   bias_table = None, genome_file_name = None, ext_both_directions=False, print_wig = None):
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
        bias_table -- Bias table to perform bias correction.
        
        Return:
        hon_signal -- Normalized signal.
        slopehon_signal -- Slope signal.
        """

        # Fetch raw signal
        pileup_region = PileupRegion(start,end,ext)
        if(self.is_bam):
            if(ps_version == "0.7.5"):
                self.bam.fetch(reference=ref, start=start, end=end, callback = pileup_region)
            else:
                iter = self.bam.fetch(reference=ref, start=start, end=end)
                if(not ext_both_directions):
                    for alignment in iter: pileup_region.__call__(alignment)
                else:
                    for alignment in iter: pileup_region.__call2__(alignment)
            raw_signal = array([min(e,initial_clip) for e in pileup_region.vector])
        elif(self.is_bw):
            signal = self.bw.pileup(ref, start, end)
            raw_signal = array([min(e,initial_clip) for e in signal])

        # Std-based clipping
        mean = raw_signal.mean()
        std = raw_signal.std()
        clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

        # Bias correction
        bias_corrected_signal = self.bias_correction(clip_signal, bias_table, genome_file_name, ref, start, end)

        # Boyle normalization (within-dataset normalization)
        boyle_signal = array(self.boyle_norm(bias_corrected_signal))

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

        # Writing signal
        if(print_wig):
            signal_file = open(print_wig+"signal.wig","a")
            norm_file = open(print_wig+"norm.wig","a")
            slope_file = open(print_wig+"slope.wig","a")
            signal_file.write("fixedStep chrom="+ref+" start="+str(start+1)+" step=1\n"+"\n".join([str(e) for e in clip_signal])+"\n")
            norm_file.write("fixedStep chrom="+ref+" start="+str(start+1)+" step=1\n"+"\n".join([str(e) for e in hon_signal])+"\n")
            slope_file.write("fixedStep chrom="+ref+" start="+str(start+1)+" step=1\n"+"\n".join([str(e) for e in slopehon_signal])+"\n")
            signal_file.close()
            norm_file.close()
            slope_file.close()

        # Returning normalized and slope sequences
        return hon_signal, slopehon_signal

    def bias_correction(self, signal, bias_table, genome_file_name, chrName, start, end):
        """ 
        Performs bias correction.

        Keyword arguments:
        signal -- Input signal.
        bias_table -- Bias table.
        
        Return:
        bias_corrected_signal -- Bias-corrected sequence.
        """

        if(not bias_table): return signal

        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fastaFile = Fastafile(genome_file_name)
        fBiasDict = bias_table.table[0]; rBiasDict = bias_table.table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start; p2 = end
        p1_w = p1 - (window/2); p2_w = p2 + (window/2)
        p1_wk = p1_w - (k_nb/2); p2_wk = p2_w + (k_nb/2)

        # Raw counts
        nf = [0.0] * (p2_w-p1_w); nr = [0.0] * (p2_w-p1_w)
        for r in self.bam.fetch(chrName, p1_w, p2_w):
            if((not r.is_reverse) and (r.pos > p1_w)): nf[r.pos-p1_w] += 1.0
            if((r.is_reverse) and ((r.aend-1) < p2_w)): nr[r.aend-1-p1_w] += 1.0

        # Smoothed counts
        Nf = []; Nr = [];
        fSum = sum(nf[:window]); rSum = sum(nr[:window]);
        fLast = nf[0]; rLast = nr[0]
        for i in range((window/2),len(nf)-(window/2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast; fSum += nf[i+(window/2)]; fLast = nf[i-(window/2)+1]
            rSum -= rLast; rSum += nr[i+(window/2)]; rLast = nr[i-(window/2)+1]

        # Fetching sequence
        currStr = str(fastaFile.fetch(chrName, p1_wk-1, p2_wk-2)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fastaFile.fetch(chrName,p1_wk+2, p2_wk+1)).upper())

        # Iterating on sequence to create signal
        af = []; ar = []
        for i in range((k_nb/2),len(currStr)-(k_nb/2)+1):
            fseq = currStr[i-(k_nb/2):i+(k_nb/2)]
            rseq = currRevComp[len(currStr)-(k_nb/2)-i:len(currStr)+(k_nb/2)-i]
            try: af.append(fBiasDict[fseq])
            except Exception: af.append(defaultKmerValue)
            try: ar.append(rBiasDict[rseq])
            except Exception: ar.append(defaultKmerValue)

        # Calculating bias and writing to wig file
        fSum = sum(af[:window]); rSum = sum(ar[:window]);
        fLast = af[0]; rLast = ar[0]
        bias_corrected_signal = []
        for i in range((window/2),len(af)-(window/2)):
            nhatf = Nf[i-(window/2)]*(af[i]/fSum)
            nhatr = Nr[i-(window/2)]*(ar[i]/rSum)
            zf = log(nf[i]+1)-log(nhatf+1)
            zr = log(nr[i]+1)-log(nhatr+1)
            bias_corrected_signal.append(zf+zr)
            fSum -= fLast; fSum += af[i+(window/2)]; fLast = af[i-(window/2)+1]
            rSum -= rLast; rSum += ar[i+(window/2)]; rLast = ar[i-(window/2)+1]

        # Termination
        fastaFile.close()
        return bias_corrected_signal

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


