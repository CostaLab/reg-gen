###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings

warnings.filterwarnings("ignore")
from math import log, ceil, floor, isnan

# Internal
from ..Util import ErrorHandler
from ..Util import AuxiliaryFunctions
from pileupRegion import PileupRegion

# External
from pysam import __version__ as ps_version
from pysam import Samfile
from pysam import Fastafile
from numpy import exp, array, abs, int, mat, linalg, convolve, nan, nan_to_num
from scipy.stats import scoreatpercentile

"""
Processes DNase-seq and histone modification signal for
HMM footprinting input.

Authors: Eduardo G. Gusmao.
"""


class GenomicSignal:
    """
    Represents a genomic signal. It should be used to fetch normalized and slope
    signals from a bam file.
    Usage:
    1. Initialize class.
    2. Call load_sg_coefs once.
    3. Call get_signal as many times as needed.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, file_name):
        """ 
        Initializes GenomicSignal.
        """
        self.file_name = file_name
        self.sg_coefs = None
        self.bam = Samfile(file_name, "rb")

    def load_sg_coefs(self, slope_window_size):
        """ 
        Loads Savitzky-Golay coefficients into self.sg_coefs based on a slope_window_size.

        Keyword arguments:
        slope_window_size -- Window size of Savitzky-Golay coefficients.

        Return:
        None -- It updates self.sg_coefs.
        """
        self.sg_coefs = self.savitzky_golay_coefficients(slope_window_size, 2, 1)

    def get_tag_count(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                      initial_clip=1000):
        """
        Gets the tag count associated with self.bam based on start, end and ext.

        Keyword arguments:
        ref -- Chromosome name.
        start -- Initial genomic coordinate of signal.
        end -- Final genomic coordinate of signal.
        downstream_ext -- Number of bps to extend towards the downstream region (right for forward strand and left for reverse strand).
        upstream_ext -- Number of bps to extend towards the upstream region (left for forward strand and right for reverse strand).
        forward_shift -- Number of bps to shift the reads aligned to the forward strand. Can be a positive number for a shift towards the downstream region (towards the inside of the aligned read) and a negative number for a shift towards the upstream region.
        reverse_shift -- Number of bps to shift the reads aligned to the reverse strand. Can be a positive number for a shift towards the upstream region and a negative number for a shift towards the downstream region (towards the inside of the aligned read).
        initial_clip -- Signal will be initially clipped at this level to avoid outliers.

        Return:
        tag_count -- Total signal.
        """

        # Fetch raw signal
        pileup_region = PileupRegion(start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift)
        if (ps_version == "0.7.5"):
            self.bam.fetch(reference=ref, start=start, end=end, callback=pileup_region)
        else:
            iter = self.bam.fetch(reference=ref, start=start, end=end)
            for alignment in iter: pileup_region.__call__(alignment)
        raw_signal = array([min(e, initial_clip) for e in pileup_region.vector])

        # Std-based clipping
        # mean = raw_signal.mean()
        # std = raw_signal.std()
        # clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

        # Tag count
        try:
            tag_count = sum(raw_signal)
        except Exception:
            tag_count = 0

        return tag_count

    def get_signal(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                   initial_clip=1000, per_norm=98, per_slope=98,
                   bias_table=None, genome_file_name=None, print_raw_signal=False,
                   print_bc_signal=False, print_norm_signal=False, print_slope_signal=False,
                   strands_specific=False):
        """
        Gets the signal associated with self.bam based on start, end and ext.
        initial_clip, per_norm and per_slope are used as normalization factors during the normalization
        and slope evaluation procedures.

        Keyword arguments:
        ref -- Chromosome name.
        start -- Initial genomic coordinate of signal.
        end -- Final genomic coordinate of signal.
        initial_clip -- Signal will be initially clipped at this level to avoid outliers.
        per_norm -- Percentile value for 'hon_norm' function of the normalized signal.
        per_slope -- Percentile value for 'hon_norm' function of the slope signal.
        bias_table -- Bias table to perform bias correction.
        genome_file_name -- Genome to perform bias correction.
        downstream_ext -- Number of bps to extend towards the downstream region
        (right for forward strand and left for reverse strand).
        upstream_ext -- Number of bps to extend towards the upstream region
        (left for forward strand and right for reverse strand).
        forward_shift -- Number of bps to shift the reads aligned to the forward strand.
        Can be a positive number for a shift towards the downstream region
        (towards the inside of the aligned read) and a negative number for a shift towards the upstream region.
        reverse_shift -- Number of bps to shift the reads aligned to the reverse strand.
        Can be a positive number for a shift towards the upstream region and a negative number
        for a shift towards the downstream region (towards the inside of the aligned read).

        Return:
        hon_signal -- Normalized signal.
        slopehon_signal -- Slope signal.
        """

        # Fetch raw signal
        pileup_region = PileupRegion(start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift)
        if (ps_version == "0.7.5"):
            self.bam.fetch(reference=ref, start=start, end=end, callback=pileup_region)
        else:
            iter = self.bam.fetch(reference=ref, start=start, end=end)
            for alignment in iter:
                pileup_region.__call__(alignment)
        raw_signal = array([min(e, initial_clip) for e in pileup_region.vector])

        # Std-based clipping
        mean = raw_signal.mean()
        std = raw_signal.std()
        clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

        # Cleavage bias correction
        bias_corrected_signal = self.bias_correction(clip_signal, bias_table, genome_file_name,
                                                     ref, start, end, forward_shift, reverse_shift, strands_specific)

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
        if (print_raw_signal):
            signal_file = open(print_raw_signal, "a")
            signal_file.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                [str(e) for e in nan_to_num(raw_signal)]) + "\n")
            signal_file.close()
        if (print_bc_signal):
            signal_file = open(print_bc_signal, "a")
            signal_file.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                [str(e) for e in nan_to_num(bias_corrected_signal)]) + "\n")
            signal_file.close()
        if (print_norm_signal):
            signal_file = open(print_norm_signal, "a")
            signal_file.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                [str(e) for e in nan_to_num(hon_signal)]) + "\n")
            signal_file.close()
        if (print_slope_signal):
            signal_file = open(print_slope_signal, "a")
            signal_file.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                [str(e) for e in nan_to_num(slope_signal)]) + "\n")
            signal_file.close()

        # Returning normalized and slope sequences
        return hon_signal, slopehon_signal

    def bias_correction(self, signal, bias_table, genome_file_name, chrName, start, end,
                        forward_shift, reverse_shift, strands_specific):
        """
        Performs bias correction.

        Keyword arguments:
        signal -- Input signal.
        bias_table -- Bias table.

        Return:
        bias_corrected_signal -- Bias-corrected sequence.
        """

        if (not bias_table): return signal

        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fastaFile = Fastafile(genome_file_name)
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(fBiasDict.keys()[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window / 2)
        p2_w = p2 + (window / 2)
        p1_wk = p1_w - int(floor(k_nb / 2.))
        p2_wk = p2_w + int(ceil(k_nb / 2.))
        if (p1 <= 0 or p1_w <= 0 or p1_wk <= 0): return signal

        # Raw counts
        nf = [0.0] * (p2_w - p1_w)
        nr = [0.0] * (p2_w - p1_w)
        for read in self.bam.fetch(chrName, p1_w, p2_w):
            if (not read.is_reverse):
                cut_site = read.pos + forward_shift
                if cut_site >= start and cut_site < end:
                    nf[cut_site - p1_w] += 1.0
                    # for i in range(max(read.pos + forward_shift, start), min(read.pos + forward_shift + 1, end - 1)):
                    #    nf[i - start] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if cut_site >= start and cut_site < end:
                    nr[cut_site - p1_w] += 1.0
                    # for i in range(max(read.aend + reverse_shift - 1, start), min(read.aend + reverse_shift, end - 1)):
                    #    nr[i - start] += 1.0

                    # if ((not read.is_reverse) and (read.pos > p1_w)): nf[read.pos - p1_w] += 1.0
                    # if ((read.is_reverse) and ((read.aend - 1) < p2_w)): nr[read.aend - 1 - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(nf[:window])
        rSum = sum(nr[:window])
        fLast = nf[0]
        rLast = nr[0]
        for i in range((window / 2), len(nf) - (window / 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += nf[i + (window / 2)]
            fLast = nf[i - (window / 2) + 1]
            rSum -= rLast
            rSum += nr[i + (window / 2)]
            rLast = nr[i - (window / 2) + 1]

        # Fetching sequence
        # currStr = str(fastaFile.fetch(chrName, p1_wk-1, p2_wk-2)).upper()
        # currRevComp = AuxiliaryFunctions.revcomp(str(fastaFile.fetch(chrName,p1_wk+2, p2_wk+1)).upper())
        currStr = str(fastaFile.fetch(chrName, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fastaFile.fetch(chrName, p1_wk + 1,
                                                                     p2_wk)).upper())

        # Iterating on sequence to create signal
        af = []
        ar = []
        for i in range(int(ceil(k_nb / 2.)), len(currStr) - int(floor(k_nb / 2)) + 1):
            fseq = currStr[i - int(floor(k_nb / 2.)):i + int(ceil(k_nb / 2.))]
            rseq = currRevComp[len(currStr) - int(ceil(k_nb / 2.)) - i:len(currStr) + int(floor(k_nb / 2.)) - i]
            try:
                af.append(fBiasDict[fseq])
            except Exception:
                af.append(defaultKmerValue)
            try:
                ar.append(rBiasDict[rseq])
            except Exception:
                ar.append(defaultKmerValue)

        # Calculating bias and writing to wig file
        fSum = sum(af[:window])
        rSum = sum(ar[:window])
        fLast = af[0]
        rLast = ar[0]
        bias_corrected_signal = []
        bias_corrected_signal_forward = []
        bias_corrected_signal_reverse = []
        for i in range((window / 2), len(af) - (window / 2)):
            nhatf = Nf[i - (window / 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window / 2)] * (ar[i] / rSum)
            zf = log(nf[i] + 1) - log(nhatf + 1)
            zr = log(nr[i] + 1) - log(nhatr + 1)
            bias_corrected_signal_forward.append(zf)
            bias_corrected_signal_reverse.append(zr)
            bias_corrected_signal.append(zf + zr)
            fSum -= fLast
            fSum += af[i + (window / 2)]
            fLast = af[i - (window / 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window / 2)]
            rLast = ar[i - (window / 2) + 1]

        # Fixing the negative number in bias corrected signal
        min_value = abs(min(bias_corrected_signal_forward))
        bias_fixed_signal_forward = [e + min_value for e in bias_corrected_signal_forward]

        min_value = abs(min(bias_corrected_signal_reverse))
        bias_fixed_signal_reverse = [e + min_value for e in bias_corrected_signal_reverse]

        min_value = abs(min(bias_corrected_signal))
        bias_fixed_signal = [e + min_value for e in bias_corrected_signal]

        # Termination
        fastaFile.close()
        if not strands_specific:
            return bias_fixed_signal
        else:
            return bias_fixed_signal_forward, bias_fixed_signal_reverse

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

        if std != 0:
            norm_seq = []
            for e in sequence:
                norm_seq.append(1.0 / (1.0 + (exp(-(e - mean) / std))))
            return norm_seq
        else:
            return sequence

    def boyle_norm(self, sequence):
        """
        Normalizes a sequence according to Boyle's criterion.
        This represents a within-dataset normalization.

        Keyword arguments:
        sequence -- Input sequence.

        Return:
        norm_seq -- Normalized sequence.
        """
        mean = array([e for e in sequence if e > 0]).mean()
        if isnan(mean):
            return sequence
        else:
            norm_seq = [(float(e) / mean) for e in sequence]
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
        # try: # TODO ERRORS
        window_size = abs(int(window_size))
        order = abs(int(order))
        # except ValueError, msg:
        #    raise ValueError("windowSize and order have to be of type int")
        # if windowSize % 2 != 1 or windowSize < 1:
        #    raise TypeError("windowSize size must be a positive odd number")
        # if windowSize < order + 2:
        #    raise TypeError("windowSize is too small for the polynomials order")
        order_range = range(order + 1)
        half_window = (window_size - 1) // 2

        # Precompute Coefficients
        b = mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
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
        slope_seq = [e for e in slope_seq[(len(sg_coefs) / 2):(len(slope_seq) - (len(sg_coefs) / 2))]]

        return slope_seq

    def get_signal_per_strand(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                              initial_clip=1000, per_norm=98, per_slope=98,
                              bias_table=None, genome_file_name=None, print_raw_signal=False,
                              print_bc_signal=False, print_norm_signal=False, print_slope_signal=False,
                              strands_specific=True):
        """

        :param ref: Chromosome name.
        :param start: Initial genomic coordinate of signal.
        :param end: Final genomic coordinate of signal.
        :param downstream_ext: Number of bps to extend towards the downstream region
        :param upstream_ext: Number of bps to extend towards the upstream region
        :param forward_shift: Number of bps to shift the reads aligned to the forward strand.
        :param reverse_shift: Number of bps to shift the reads aligned to the reverse strand.
        :param initial_clip: Signal will be initially clipped at this level to avoid outliers.
        :param per_norm: Percentile value for 'hon_norm' function of the normalized signal.
        :param per_slope: Percentile value for 'hon_norm' function of the slope signal.
        :param bias_table: Bias table to perform bias correction.
        :param genome_file_name: Genome to perform bias correction.
        :param print_raw_signal:
        :param print_bc_signal:
        :param print_norm_signal:
        :param print_slope_signal:
        :return: normalized and slope signal for each strand.
        """

        raw_signal_forward = [0.0] * (end - start)
        raw_signal_reverse = [0.0] * (end - start)

        reads = self.bam.fetch(reference=ref, start=start, end=end)
        for read in reads:
            if (not read.is_reverse):
                cut_site = read.pos + forward_shift
                if cut_site >= start and cut_site < end:
                    raw_signal_forward[cut_site - start] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if cut_site >= start and cut_site < end:
                    raw_signal_reverse[cut_site - start] += 1.0

        raw_signal_forward = array([min(e, initial_clip) for e in raw_signal_forward])
        raw_signal_reverse = array([min(e, initial_clip) for e in raw_signal_reverse])

        # Std-based clipping
        mean = raw_signal_forward.mean()
        std = raw_signal_forward.std()
        clip_signal_forward = [min(e, mean + (10 * std)) for e in raw_signal_forward]
        mean = raw_signal_reverse.mean()
        std = raw_signal_reverse.std()
        clip_signal_reverse = [min(e, mean + (10 * std)) for e in raw_signal_reverse]

        # Cleavage bias correction
        bc_signal_forward = None
        bc_signal_reverse = None
        if bias_table:
            bc_signal_forward, bc_signal_reverse = self.bias_correction(raw_signal_forward, bias_table,
                                                                        genome_file_name,
                                                                        ref, start, end, forward_shift, reverse_shift,
                                                                        strands_specific)
        else:
            bc_signal_forward = clip_signal_forward
            bc_signal_reverse = clip_signal_reverse

        # Boyle normalization (within-dataset normalization)
        boyle_signal_forward = array(self.boyle_norm(bc_signal_forward))
        boyle_signal_reverse = array(self.boyle_norm(bc_signal_reverse))

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(boyle_signal_forward, per_norm)
        std = boyle_signal_forward.std()
        hon_signal_forward = self.hon_norm(boyle_signal_forward, perc, std)

        perc = scoreatpercentile(boyle_signal_reverse, per_norm)
        std = boyle_signal_reverse.std()
        hon_signal_reverse = self.hon_norm(boyle_signal_reverse, perc, std)

        # Slope signal
        slope_signal_forward = self.slope(hon_signal_forward, self.sg_coefs)
        slope_signal_reverse = self.slope(hon_signal_reverse, self.sg_coefs)

        # Returning normalized and slope sequences
        return hon_signal_forward, slope_signal_forward, hon_signal_reverse, slope_signal_reverse