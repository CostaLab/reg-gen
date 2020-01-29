###################################################################################################
# Libraries
###################################################################################################
from math import log, ceil, floor, isnan
import numpy as np
from numpy import exp, array, abs, int, mat, linalg, convolve, nan_to_num
from pysam import Samfile, Fastafile
from pysam import __version__ as ps_version
from scipy.stats import scoreatpercentile

# Internal
from rgt.Util import AuxiliaryFunctions
from rgt.HINT.pileupRegion import PileupRegion

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

    def __init__(self, file_name=None):
        """ 
        Initializes GenomicSignal.
        """
        self.file_name = file_name
        self.sg_coefs = None
        if file_name is not None:
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
        tag_count = 0.0
        iter = self.bam.fetch(reference=ref, start=start, end=end)
        for alignment in iter:
            tag_count += 1

        return tag_count

    def get_signal(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                   initial_clip=1000, per_norm=98, per_slope=98,
                   bias_table=None, genome_file_name=None, print_raw_signal=False):
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
        if ps_version == "0.7.5":
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
        bc_signal = self.bias_correction_dnase(clip_signal, bias_table, genome_file_name, ref, start, end,
                                               forward_shift, reverse_shift)

        # Boyle normalization (within-dataset normalization)
        boyle_signal = array(self.boyle_norm(bc_signal))

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(boyle_signal, per_norm)
        std = boyle_signal.std()
        hon_signal = self.hon_norm_dnase(boyle_signal, perc, std)

        # Slope signal
        slope_signal = self.slope(hon_signal, self.sg_coefs)

        # Returning normalized and slope sequences
        return hon_signal, slope_signal

    def get_signal_atac(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                        initial_clip=50, per_norm=98, per_slope=98,
                        bias_table=None, genome_file_name=None):

        # Cleavage bias correction
        bc_signal_forward, bc_signal_reverse = self.bias_correction_atac(bias_table, genome_file_name,
                                                                         ref, start, end, forward_shift, reverse_shift)

        # Boyle normalization (within-dataset normalization)
        boyle_signal_forward = array(self.boyle_norm(bc_signal_forward))
        boyle_signal_reverse = array(self.boyle_norm(bc_signal_reverse))

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(boyle_signal_forward, per_norm)
        std = boyle_signal_forward.std()
        hon_signal_forward = self.hon_norm_atac(boyle_signal_forward, perc, std)

        perc = scoreatpercentile(boyle_signal_reverse, per_norm)
        std = boyle_signal_reverse.std()
        hon_signal_reverse = self.hon_norm_atac(boyle_signal_reverse, perc, std)

        # Slope signal
        slope_signal_forward = self.slope(hon_signal_forward, self.sg_coefs)
        slope_signal_reverse = self.slope(hon_signal_reverse, self.sg_coefs)

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(slope_signal_forward, per_norm)
        std = np.std(slope_signal_forward)
        slope_signal_forward = self.hon_norm_atac(slope_signal_forward, perc, std)

        perc = scoreatpercentile(slope_signal_reverse, per_norm)
        std = np.std(slope_signal_forward)
        slope_signal_reverse = self.hon_norm_atac(slope_signal_reverse, perc, std)

        # Returning normalized and slope sequences
        return hon_signal_forward, slope_signal_forward, hon_signal_reverse, slope_signal_reverse

    def get_signal_atac2(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                         initial_clip=50, per_norm=98, per_slope=98,
                         bias_table=None, genome_file_name=None):

        # Cleavage bias correction
        bc_signal = self.bias_correction_atac2(bias_table, genome_file_name,
                                               ref, start, end, forward_shift, reverse_shift)

        # Boyle normalization (within-dataset normalization)
        boyle_signal = array(self.boyle_norm(bc_signal))

        # Hon normalization (between-dataset normalization)
        perc = scoreatpercentile(boyle_signal, per_norm)
        std = boyle_signal.std()
        hon_signal = self.hon_norm_atac(boyle_signal, perc, std)

        # Slope signal
        slope_signal = self.slope(hon_signal, self.sg_coefs)

        # Hon normalization (between-dataset normalization)
        slope_signal = self.boyle_norm(slope_signal)

        perc = scoreatpercentile(slope_signal, per_slope)
        std = np.std(slope_signal)
        slope_signal = self.hon_norm_atac(slope_signal, perc, std)

        # Returning normalized and slope sequences
        return hon_signal, slope_signal

    def bias_correction_dnase(self, signal, bias_table, genome_file_name, chrName, start, end,
                              forward_shift, reverse_shift):

        if not bias_table: return signal
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fastaFile = Fastafile(genome_file_name)
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window // 2)
        p2_w = p2 + (window // 2)
        p1_wk = p1_w - int(floor(k_nb / 2.))
        p2_wk = p2_w + int(ceil(k_nb / 2.))
        if p1 <= 0 or p1_w <= 0 or p1_wk <= 0: return signal

        # Raw counts
        nf = [0.0] * (p2_w - p1_w)
        nr = [0.0] * (p2_w - p1_w)
        for read in self.bam.fetch(chrName, p1_w, p2_w):
            # check if the read is unmapped, according to issue #112
            if read.is_unmapped:
                continue

            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1_w <= cut_site < p2_w:
                    nf[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1_w <= cut_site < p2_w:
                    nr[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(nf[:window])
        rSum = sum(nr[:window])
        fLast = nf[0]
        rLast = nr[0]
        for i in range((window // 2), len(nf) - (window // 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += nf[i + (window // 2)]
            fLast = nf[i - (window // 2) + 1]
            rSum -= rLast
            rSum += nr[i + (window // 2)]
            rLast = nr[i - (window // 2) + 1]

        # Fetching sequence
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
        for i in range((window // 2), len(af) - (window // 2)):
            nhatf = Nf[i - (window // 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window // 2)] * (ar[i] / rSum)
            zf = log(nf[i] + 1) - log(nhatf + 1)
            zr = log(nr[i] + 1) - log(nhatr + 1)
            bias_corrected_signal.append(zf + zr)
            fSum -= fLast
            fSum += af[i + (window // 2)]
            fLast = af[i - (window // 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window // 2)]
            rLast = ar[i - (window // 2) + 1]

        # Termination
        fastaFile.close()
        return bias_corrected_signal

    def bias_correction_atac(self, bias_table, genome_file_name, chrName, start, end,
                             forward_shift, reverse_shift):

        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fastaFile = Fastafile(genome_file_name)
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window // 2)
        p2_w = p2 + (window // 2)
        p1_wk = p1_w - int(floor(k_nb / 2.))
        p2_wk = p2_w + int(ceil(k_nb / 2.))

        if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
            # Return raw counts
            nf = [0.0] * (p2 - p1)
            nr = [0.0] * (p2 - p1)
            for read in self.bam.fetch(chrName, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1 <= cut_site < p2:
                        nf[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1 <= cut_site < p2:
                        nr[cut_site - p1] += 1.0

            return nf, nr

        # Raw counts
        nf = [0.0] * (p2_w - p1_w)
        nr = [0.0] * (p2_w - p1_w)
        for read in self.bam.fetch(chrName, p1_w, p2_w):
            # check if the read is unmapped, according to issue #112
            if read.is_unmapped:
                continue

            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1_w <= cut_site < p2_w:
                    nf[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1_w <= cut_site < p2_w:
                    nr[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(nf[:window])
        rSum = sum(nr[:window])
        fLast = nf[0]
        rLast = nr[0]
        for i in range((window // 2), len(nf) - (window // 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += nf[i + (window // 2)]
            fLast = nf[i - (window // 2) + 1]
            rSum -= rLast
            rSum += nr[i + (window // 2)]
            rLast = nr[i - (window // 2) + 1]

        # Fetching sequence
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
        bias_corrected_signal_forward = []
        bias_corrected_signal_reverse = []
        for i in range((window // 2), len(af) - (window // 2)):
            nhatf = Nf[i - (window // 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window // 2)] * (ar[i] / rSum)
            bias_corrected_signal_forward.append(nhatf)
            bias_corrected_signal_reverse.append(nhatr)
            fSum -= fLast
            fSum += af[i + (window // 2)]
            fLast = af[i - (window // 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window // 2)]
            rLast = ar[i - (window // 2) + 1]

        # Termination
        fastaFile.close()
        return bias_corrected_signal_forward, bias_corrected_signal_reverse

    def bias_correction_atac2(self, bias_table, genome_file_name, chrName, start, end,
                              forward_shift, reverse_shift):

        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fastaFile = Fastafile(genome_file_name)
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window // 2)
        p2_w = p2 + (window // 2)
        p1_wk = p1_w - int(floor(k_nb / 2.))
        p2_wk = p2_w + int(ceil(k_nb / 2.))
        if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
            # Return raw counts
            signal = [0.0] * (p2 - p1)
            for read in self.bam.fetch(chrName, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0

            return signal

        # Raw counts
        nf = [0.0] * (p2_w - p1_w)
        nr = [0.0] * (p2_w - p1_w)
        for read in self.bam.fetch(chrName, p1_w, p2_w):
            # check if the read is unmapped, according to issue #112
            if read.is_unmapped:
                continue

            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1_w <= cut_site < p2_w:
                    nf[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1_w <= cut_site < p2_w:
                    nr[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(nf[:window])
        rSum = sum(nr[:window])
        fLast = nf[0]
        rLast = nr[0]
        for i in range((window // 2), len(nf) - (window // 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += nf[i + (window // 2)]
            fLast = nf[i - (window // 2) + 1]
            rSum -= rLast
            rSum += nr[i + (window // 2)]
            rLast = nr[i - (window // 2) + 1]

        # Fetching sequence
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
        bc_signal = []
        for i in range((window // 2), len(af) - (window // 2)):
            nhatf = Nf[i - (window // 2)] * (af[i] / fSum)
            nhatr = Nr[i - (window // 2)] * (ar[i] / rSum)
            bc_signal.append(nhatf + nhatr)
            fSum -= fLast
            fSum += af[i + (window // 2)]
            fLast = af[i - (window // 2) + 1]
            rSum -= rLast
            rSum += ar[i + (window // 2)]
            rLast = ar[i - (window // 2) + 1]

        # Termination
        fastaFile.close()
        return bc_signal

    def hon_norm_atac(self, sequence, mean, std):
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
                if e == 0:
                    norm_seq.append(e)
                else:
                    norm_seq.append(1.0 / (1.0 + (exp(-(e - mean) / std))))
            return norm_seq
        else:
            return sequence

    def hon_norm_dnase(self, sequence, mean, std):
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

        # if std != 0:
        #    norm_seq = []
        #    for e in sequence:
        #        norm_seq.append(1.0 / (1.0 + (exp(-(e - mean) / std))))
        #    return norm_seq
        # else:
        #    return sequence
        norm_seq = []
        for e in sequence:
            if e == 0.0:
                norm_seq.append(0.0)
            elif e > 0.0:
                norm_seq.append(1.0 / (1.0 + (exp(-(e - mean) / std))))
            else:
                norm_seq.append(-1.0 / (1.0 + (exp(-(-e - mean) / std))))
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
        order_range = list(range(order + 1))
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
        slope_seq = [e for e in slope_seq[(len(sg_coefs) // 2):(len(slope_seq) - (len(sg_coefs) // 2))]]

        return slope_seq

    def print_signal(self, ref, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift,
                     initial_clip=1000, per_norm=98, per_slope=98, bias_table=None, genome_file_name=None,
                     raw_signal_file=None, bc_signal_file=None, norm_signal_file=None, strand_specific=False):

        if raw_signal_file:
            pileup_region = PileupRegion(start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift)
            if ps_version == "0.7.5":
                self.bam.fetch(reference=ref, start=start, end=end, callback=pileup_region)
            else:
                iter = self.bam.fetch(reference=ref, start=start, end=end)
                for alignment in iter:
                    pileup_region.__call__(alignment)
            raw_signal = array([min(e, initial_clip) for e in pileup_region.vector])

            f = open(raw_signal_file, "a")
            f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                [str(e) for e in nan_to_num(raw_signal)]) + "\n")
            f.close()

        if bc_signal_file or norm_signal_file:
            # Parameters
            window = 50
            defaultKmerValue = 1.0

            # Initialization
            fasta = Fastafile(genome_file_name)
            fBiasDict = bias_table[0]
            rBiasDict = bias_table[1]
            k_nb = len(list(fBiasDict.keys())[0])
            p1 = start
            p2 = end
            p1_w = p1 - (window // 2)
            p2_w = p2 + (window // 2)
            p1_wk = p1_w - int(k_nb / 2.)
            p2_wk = p2_w + int(k_nb / 2.)

            currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
            currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

            # Iterating on sequence to create the bias signal
            signal_bias_f = []
            signal_bias_r = []
            for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
                fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
                rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
                try:
                    signal_bias_f.append(fBiasDict[fseq])
                except Exception:
                    signal_bias_f.append(defaultKmerValue)
                try:
                    signal_bias_r.append(rBiasDict[rseq])
                except Exception:
                    signal_bias_r.append(defaultKmerValue)

            # Raw counts
            signal_raw_f = [0.0] * (p2_w - p1_w)
            signal_raw_r = [0.0] * (p2_w - p1_w)
            for read in self.bam.fetch(ref, p1_w, p2_w):
                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1_w <= cut_site < p2_w:
                        signal_raw_f[cut_site - p1_w] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1_w <= cut_site < p2_w:
                        signal_raw_r[cut_site - p1_w] += 1.0

            # Smoothed counts
            Nf = []
            Nr = []
            fSum = sum(signal_raw_f[:window])
            rSum = sum(signal_raw_r[:window])
            fLast = signal_raw_f[0]
            rLast = signal_raw_r[0]
            for i in range((window // 2), len(signal_raw_f) - (window // 2)):
                Nf.append(fSum)
                Nr.append(rSum)
                fSum -= fLast
                fSum += signal_raw_f[i + (window // 2)]
                fLast = signal_raw_f[i - (window // 2) + 1]
                rSum -= rLast
                rSum += signal_raw_r[i + (window // 2)]
                rLast = signal_raw_r[i - (window // 2) + 1]

            # Calculating bias and writing to wig file
            fSum = sum(signal_bias_f[:window])
            rSum = sum(signal_bias_r[:window])
            fLast = signal_bias_f[0]
            rLast = signal_bias_r[0]
            signal_bc = []
            signal_bc_f = []
            signal_bc_r = []
            for i in range((window // 2), len(signal_bias_f) - (window // 2)):
                nhatf = Nf[i - (window // 2)] * (signal_bias_f[i] / fSum)
                nhatr = Nr[i - (window // 2)] * (signal_bias_r[i] / rSum)
                signal_bc.append(nhatf + nhatr)
                signal_bc_f.append(nhatf)
                signal_bc_r.append(nhatr)
                fSum -= fLast
                fSum += signal_bias_f[i + (window // 2)]
                fLast = signal_bias_f[i - (window // 2) + 1]
                rSum -= rLast
                rSum += signal_bias_r[i + (window // 2)]
                rLast = signal_bias_r[i - (window // 2) + 1]

            if bc_signal_file:
                f = open(bc_signal_file, "a")
                f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                    [str(e) for e in nan_to_num(signal_bc)]) + "\n")
                f.close()

                if strand_specific:
                    prefix = bc_signal_file.split(".")[0]
                    bc_signal_file_f = prefix + "_Forward" + ".bc.wig"
                    bc_signal_file_r = prefix + "_Reverse" + ".bc.wig"
                    f = open(bc_signal_file_f, "a")
                    f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                        [str(e) for e in nan_to_num(signal_bc_f)]) + "\n")
                    f.close()
                    f = open(bc_signal_file_r, "a")
                    f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                        [str(e) for e in nan_to_num(signal_bc_r)]) + "\n")
                    f.close()

            if norm_signal_file:
                norm_signal_bc = self.boyle_norm(signal_bc)
                perc = scoreatpercentile(norm_signal_bc, 98)
                std = np.std(norm_signal_bc)
                norm_signal_bc = self.hon_norm_atac(norm_signal_bc, perc, std)
                f = open(norm_signal_file, "a")
                f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                    [str(e) for e in nan_to_num(norm_signal_bc)]) + "\n")
                f.close()

                if strand_specific:
                    prefix = bc_signal_file.split(".")[0]
                    norm_signal_file_f = prefix + "_Forward" + ".norm.wig"
                    norm_signal_file_r = prefix + "_Reverse" + ".norm.wig"

                    signal_norm_f = self.boyle_norm(signal_bc_f)
                    perc = scoreatpercentile(signal_norm_f, 98)
                    std = np.std(signal_norm_f)
                    signal_norm_f = self.hon_norm_atac(signal_norm_f, perc, std)

                    signal_norm_r = self.boyle_norm(signal_bc_r)
                    perc = scoreatpercentile(signal_norm_r, 98)
                    std = np.std(signal_norm_r)
                    signal_norm_r = self.hon_norm_atac(signal_norm_r, perc, std)

                    f = open(norm_signal_file_f, "a")
                    f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                        [str(e) for e in nan_to_num(signal_norm_f)]) + "\n")
                    f.close()
                    f = open(norm_signal_file_r, "a")
                    f.write("fixedStep chrom=" + ref + " start=" + str(start + 1) + " step=1\n" + "\n".join(
                        [str(e) for e in nan_to_num(signal_norm_r)]) + "\n")
                    f.close()

    def get_raw_signal_by_fragment_length(self, ref, start, end, bam,
                                          forward_shift, reverse_shift, min_length=None, max_length=None,
                                          strand=True):

        p1 = start
        p2 = end
        raw_f = [0.0] * (p2 - p1)
        raw_r = [0.0] * (p2 - p1)

        if min_length is None and max_length is None:
            for read in bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1 <= cut_site < p2:
                        raw_f[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1 <= cut_site < p2:
                        raw_r[cut_site - p1] += 1.0
        elif min_length is None and max_length is not None:
            for read in bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if abs(read.template_length) <= max_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1 <= cut_site < p2:
                            raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1 <= cut_site < p2:
                            raw_r[cut_site - p1] += 1.0
        elif min_length is not None and max_length is None:
            for read in bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if abs(read.template_length) > min_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1 <= cut_site < p2:
                            raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1 <= cut_site < p2:
                            raw_r[cut_site - p1] += 1.0
        elif min_length is not None and max_length is not None:
            for read in bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if min_length <= abs(read.template_length) <= max_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1 <= cut_site < p2:
                            raw_f[cut_site - p1] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1 <= cut_site < p2:
                            raw_r[cut_site - p1] += 1.0
        if strand:
            return np.array(raw_f), np.array(raw_r)
        else:
            return np.add(np.array(raw_f), np.array(raw_r))

    def get_bc_signal_by_fragment_length(self, ref, start, end, bam, fasta, bias_table,
                                         forward_shift, reverse_shift, min_length=None, max_length=None,
                                         strand=True):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window // 2)
        p2_w = p2 + (window // 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        if (p1 <= 0 or p1_w <= 0 or p2_wk <= 0):
            # Return raw counts
            signal = [0.0] * (p2 - p1)
            for read in self.bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0

            return signal

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        raw_f = [0.0] * (p2_w - p1_w)
        raw_r = [0.0] * (p2_w - p1_w)

        if min_length is None and max_length is None:
            for read in bam.fetch(ref, p1_w, p2_w):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1_w <= cut_site < p2_w:
                        raw_f[cut_site - p1_w] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1_w <= cut_site < p2_w:
                        raw_r[cut_site - p1_w] += 1.0
        elif min_length is None and max_length is not None:
            for read in bam.fetch(ref, p1_w, p2_w):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if abs(read.template_length) <= max_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1_w <= cut_site < p2_w:
                            raw_f[cut_site - p1_w] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1_w <= cut_site < p2_w:
                            raw_r[cut_site - p1_w] += 1.0
        elif min_length is not None and max_length is None:
            for read in bam.fetch(ref, p1_w, p2_w):
                if abs(read.template_length) > min_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1_w <= cut_site < p2_w:
                            raw_f[cut_site - p1_w] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1_w <= cut_site < p2_w:
                            raw_r[cut_site - p1_w] += 1.0
        elif min_length is not None and max_length is not None:
            for read in bam.fetch(ref, p1_w, p2_w):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if min_length < abs(read.template_length) <= max_length:
                    if not read.is_reverse:
                        cut_site = read.pos + forward_shift
                        if p1_w <= cut_site < p2_w:
                            raw_f[cut_site - p1_w] += 1.0
                    else:
                        cut_site = read.aend + reverse_shift - 1
                        if p1_w <= cut_site < p2_w:
                            raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(raw_f[:window])
        rSum = sum(raw_r[:window])
        fLast = raw_f[0]
        rLast = raw_r[0]
        for i in range((window // 2), len(raw_f) - (window // 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += raw_f[i + (window // 2)]
            fLast = raw_f[i - (window // 2) + 1]
            rSum -= rLast
            rSum += raw_r[i + (window // 2)]
            rLast = raw_r[i - (window // 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        bc_f = []
        bc_r = []
        for i in range((window // 2), len(signal_bias_f) - (window // 2)):
            nhatf = Nf[i - (window // 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window // 2)] * (signal_bias_r[i] / rSum)
            bc_f.append(nhatf)
            bc_r.append(nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window // 2)]
            fLast = signal_bias_f[i - (window // 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window // 2)]
            rLast = signal_bias_r[i - (window // 2) + 1]

        if strand:
            return np.array(bc_f), np.array(bc_r)
        else:
            return np.add(np.array(bc_f), np.array(bc_r))

    def get_bias_raw_bc_signal(self, ref, start, end, bam, fasta, bias_table, forward_shift, reverse_shift,
                               strand=False):
        # Parameters
        window = 50
        defaultKmerValue = 1.0

        # Initialization
        fBiasDict = bias_table[0]
        rBiasDict = bias_table[1]
        k_nb = len(list(fBiasDict.keys())[0])
        p1 = start
        p2 = end
        p1_w = p1 - (window // 2)
        p2_w = p2 + (window // 2)
        p1_wk = p1_w - int(k_nb / 2.)
        p2_wk = p2_w + int(k_nb / 2.)

        if p1 <= 0 or p1_w <= 0 or p2_wk <= 0:
            # Return raw counts
            signal = [0.0] * (p2 - p1)
            for read in self.bam.fetch(ref, p1, p2):
                # check if the read is unmapped, according to issue #112
                if read.is_unmapped:
                    continue

                if not read.is_reverse:
                    cut_site = read.pos + forward_shift
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0
                else:
                    cut_site = read.aend + reverse_shift - 1
                    if p1 <= cut_site < p2:
                        signal[cut_site - p1] += 1.0

            return signal

        currStr = str(fasta.fetch(ref, p1_wk - 1 + forward_shift, p2_wk - 2 + forward_shift)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + reverse_shift + 2,
                                                                 p2_wk + reverse_shift + 1)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        # Raw counts
        signal_raw_f = [0.0] * (p2_w - p1_w)
        signal_raw_r = [0.0] * (p2_w - p1_w)
        for read in bam.fetch(ref, p1_w, p2_w):
            # check if the read is unmapped, according to issue #112
            if read.is_unmapped:
                continue

            if not read.is_reverse:
                cut_site = read.pos + forward_shift
                if p1_w <= cut_site < p2_w:
                    signal_raw_f[cut_site - p1_w] += 1.0
            else:
                cut_site = read.aend + reverse_shift - 1
                if p1_w <= cut_site < p2_w:
                    signal_raw_r[cut_site - p1_w] += 1.0

        # Smoothed counts
        Nf = []
        Nr = []
        fSum = sum(signal_raw_f[:window])
        rSum = sum(signal_raw_r[:window])
        fLast = signal_raw_f[0]
        rLast = signal_raw_r[0]
        for i in range((window // 2), len(signal_raw_f) - (window // 2)):
            Nf.append(fSum)
            Nr.append(rSum)
            fSum -= fLast
            fSum += signal_raw_f[i + (window // 2)]
            fLast = signal_raw_f[i - (window // 2) + 1]
            rSum -= rLast
            rSum += signal_raw_r[i + (window // 2)]
            rLast = signal_raw_r[i - (window // 2) + 1]

        # Calculating bias and writing to wig file
        fSum = sum(signal_bias_f[:window])
        rSum = sum(signal_bias_r[:window])
        fLast = signal_bias_f[0]
        rLast = signal_bias_r[0]
        bias_f = []
        bias_r = []
        raw = []
        raw_f = []
        raw_r = []
        bc = []
        bc_f = []
        bc_r = []
        for i in range((window // 2), len(signal_bias_f) - (window // 2)):
            nhatf = Nf[i - (window // 2)] * (signal_bias_f[i] / fSum)
            nhatr = Nr[i - (window // 2)] * (signal_bias_r[i] / rSum)
            bias_f.append(signal_bias_f[i])
            bias_r.append(signal_bias_r[i])
            raw.append(signal_raw_f[i] + signal_raw_r[i])
            raw_f.append(signal_raw_f[i])
            raw_r.append(signal_raw_r[i])
            # zf = (signal_raw_f[i]) / (signal_bias_f[i])
            # zr = (signal_raw_r[i]) / (signal_bias_r[i])
            bc.append(nhatf + nhatr)
            bc_f.append(nhatf)
            bc_r.append(nhatr)
            fSum -= fLast
            fSum += signal_bias_f[i + (window // 2)]
            fLast = signal_bias_f[i - (window // 2) + 1]
            rSum -= rLast
            rSum += signal_bias_r[i + (window // 2)]
            rLast = signal_bias_r[i - (window // 2) + 1]

        currStr = str(fasta.fetch(ref, p1_wk, p2_wk - 1)).upper()
        currRevComp = AuxiliaryFunctions.revcomp(str(fasta.fetch(ref, p1_wk + 1, p2_wk)).upper())

        # Iterating on sequence to create the bias signal
        signal_bias_f = []
        signal_bias_r = []
        for i in range(int(k_nb / 2.), len(currStr) - int(k_nb / 2) + 1):
            fseq = currStr[i - int(k_nb / 2.):i + int(k_nb / 2.)]
            rseq = currRevComp[len(currStr) - int(k_nb / 2.) - i:len(currStr) + int(k_nb / 2.) - i]
            try:
                signal_bias_f.append(fBiasDict[fseq])
            except Exception:
                signal_bias_f.append(defaultKmerValue)
            try:
                signal_bias_r.append(rBiasDict[rseq])
            except Exception:
                signal_bias_r.append(defaultKmerValue)

        bias_f = []
        bias_r = []
        for i in range((window // 2), len(signal_bias_f) - (window // 2)):
            bias_f.append(signal_bias_f[i])
            bias_r.append(signal_bias_r[i])

        if strand:
            return bias_f, bias_r, raw, raw_f, raw_r, bc, bc_f, bc_r
        else:
            return bias_f, bias_r, raw, bc
