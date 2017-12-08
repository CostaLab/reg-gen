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

    def __init__(self, start, end, downstream_ext, upstream_ext, forward_shift, reverse_shift):
        """ 
        Initializes PileupRegion.

        Variables:
        start -- Start point (in bp) of this PileupRegion (integer).
        end -- End point (in bp) of this PileupRegion (integer).
        length -- Length of this PileupRegion (integer).
        downstream_ext -- Number of bps to extend towards the downstream region (right for forward strand and left for reverse strand).
        upstream_ext -- Number of bps to extend towards the upstream region (left for forward strand and right for reverse strand).
        forward_shift -- Number of bps to shift the reads aligned to the forward strand.
        Can be a positive number for a shift towards the downstream region
        (towards the inside of the aligned read) and a negative number for a shift towards the upstream region.
        reverse_shift -- Number of bps to shift the reads aligned to the reverse strand.
        Can be a positive number for a shift towards the upstream region and a negative number for a shift towards the downstream region
        (towards the inside of the aligned read).
        vector -- Pileup vector. Each extended fragment will contribute with +1 to each
                  position of this vector, in which the length equals self.length (list).
        """
        self.start = start
        self.end = end
        self.length = end - start
        self.downstream_ext = downstream_ext
        self.upstream_ext = upstream_ext
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift
        self.vector = [0.0] * self.length

    def __call__(self, alignment):
        try:
            if not alignment.is_reverse:
                cut_site = alignment.pos + self.forward_shift
                if self.start <= cut_site < self.end:
                    self.vector[cut_site - self.start] += 1.0
                #for i in range(max(alignment.pos + self.forward_shift - self.upstream_ext, self.start),
                #               min(alignment.pos + self.forward_shift + self.downstream_ext, self.end - 1)):
            else:
                cut_site = alignment.aend + self.reverse_shift - 1
                if self.start <= cut_site < self.end:
                    self.vector[cut_site - self.start] += 1.0
                #for i in range(max(alignment.aend + self.reverse_shift - self.downstream_ext, self.start),
                #               min(alignment.aend + self.reverse_shift + self.upstream_ext, self.end - 1)):
                #    self.vector[i - self.start] += 1.0
        except Exception:
            pass

    """
    def __call__(self, alignment):
        ""
        Updates self.vector with the fragment obtained with the 'fetch' function from
        pysam package. This function will be called for all fragments fetched.

        Keyword arguments:
        alignment -- File location + name in HMM format. See manual for full description of such format.
        
        Return:
        None -- It updates self.vector only.
        ""

        try:
            # Forward strand update
            if not alignment.is_reverse:
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
    """
