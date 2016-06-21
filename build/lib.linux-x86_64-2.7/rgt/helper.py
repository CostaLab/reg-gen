"""
helper codes
===================
some other extra codes
"""

import os
import sys
import numpy
from GenomicRegionSet import GenomicRegionSet
from GenomicRegion import GenomicRegion
# if sys.platform.startswith("linux"):
#     import wWigIO

def get_chrom_sizes_as_genomicregionset(chrom_size_path):
    regionset = GenomicRegionSet('')
    with open(chrom_size_path) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            chrom, end = line[0], int(line[1])
            regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
    
    return regionset

def pretty(d, indent=0):
    for key, value in d.iteritems():
        print('\t' * indent + str(key))
        if isinstance(value, dict):
            pretty(value, indent+1)
        else:
            for t in value:
                print('\t' * (indent+1) + str(t))
# class BigWigFile():
#     """Fast reader of BigWig format file.

#     *Usage:*
    
#         Open file:
#             fh=BigWigFile("test.bw")
#         Get chromosome sizes:
#             chroms=fh.chromSizes() # (chroms,sizes)
#         Fetch regions:
#             wigs=fh.fetch(chrom="chr1",start=100,stop=200)
#             for wig in wigs:
#                 #do some thing with wig
#                 print wig.chrom,wig.start,wig.stop,wig.score
#         Close file:
#             fh.close()
    
#     *Parameters:*

#         chrom=None: return empty list.
#         start=None: start at first position.
#         stop=None:  stop at the end of chromosome.

#     (Source: ngslib source code https://pypi.python.org/pypi/ngslib/1.1.14 retrieved on, Jan 12, 2016.)
#     """
#     def __init__(self,infile,chrom_size=None):
#         """Open BigWig file. """
#         # Check if file exists
#         if sys.platform.startswith("linux"):
#             self.closed = True
#             if not os.path.exists(infile) and infile != 'stdin':
#                 raise IOError("ERROR: input file '{0}' dose not exist.".format(infile))
#             if infile.endswith('.wig'):
#                 bwfile = os.path.splitext(infile)[0]+".bw"
#                 if os.path.exists(bwfile):
#                     self.infile = bwfile
#                 else:
#                     if chrom_size is None:
#                         raise IOError("ERROR: 'chrom_size' file is required for wig -> bigwig conversion.")
#                     BigWigFile.wigToBigWig(infile,chrom_size,bwfile)
#                 self.infile=bwfile
#             else:
#                 self.infile=infile
#             wWigIO.open(self.infile)
#             self.closed = False
#             self.chroms, self.sizes = self.chromSizes()
#             self.closed = False
#         else:
#             print("Error: RGT doesn't support reading BigWig file on Mac or windows.")
#             sys.exit(0)

#     def chromSizes(self):
#         """Get chromosome sizes."""
#         chroms,sizes=wWigIO.getChromSize(self.infile)
#         return chroms,sizes

#     def fetch(self, chrom,start=None,stop=None, strand="+", zerobased=True,**kwargs):
#         """Fetch intervals in a given region. 
#         Note: strand is useless here.
#         Parameters:
#             start: int or None
#                 start position, None for start=0
#             end: int or None
#                 end position, None for stop = end of chrom
#             strand: string
#                 choice from '+', '-' and '.'
#             zerosbased: bool
#                 indices are zerobased or not. Useless here.
#         Dictionary parameters:
#             chunksize: int
#                 chunk size when fetch items from a large region.
#         Generates:
#             wig: tuple
#                 (start, stop, value)
#         """
#         if chrom is None: raise ValueError("ERROR: chrom name is required.")
#         if start is None: start = 0
#         if not zerobased: start -= 1
#         if start <0:
#             raise ValueError('ERROR: start should be >=0 (zerobased=True) or >=1 (zerobased=False).')
#         if stop is None: stop = self.sizes[self.chroms.index['chrom']]
#         chunksize = kwargs.get('chunksize',10000)
#         try:
#             for cstart in xrange(start,stop,chunksize):
#                 cend = cstart + chunksize
#                 if cend > stop: cend = stop
#                 for wig in wWigIO.getIntervals(self.infile,chrom,cstart,cend):
#                     if wig[0] < cstart: wig[0] = cstart
#                     if wig[1] > cend: wig[1] = cend
#                     yield wig
#         except:
#             # check result
#             if wigs == 1:
#                 raise ValueError("ERROR: wWigIO.getIntervals doesn't have correct parameters.")
#             if wigs == 2:
#                 raise ValueError("ERROR: BigWig file cannot be opened.")

#     def pileup(self,chrom,start=None,stop=None,strand="+",zerobased=True,**kwargs):
#         """Fetch depth for a genomic region."""
#         if chrom is None: raise ValueError("ERROR: chrom name is required.")
#         if start is None: start = 0
#         if not zerobased: start -= 1
#         if stop is None: stop = self.sizes[self.chroms.index(chrom)]
#         vals = numpy.zeros(stop-start)
#         for wstart,wstop,depth in self.fetch(chrom,start,stop,strand,zerobased):
#             vals[(wstart-start):(wstop-start)]+=depth
#         if strand == "-":
#             vals = vals[::-1]
#         return vals

#     def fetchBed(self,tbed,byPos=False,forcestrand=True):
#         """Fetch intervals for Bed."""
#         wigs=wWigIO.getIntervals(self.infile,tbed.chrom,tbed.start,tbed.stop)
#         if not byPos:
#             return wigs
#         # get depth by position, return a numpy array.
#         vals = numpy.zeros(tbed.length())
#         for start,stop,depth in wigs:
#             start = max(start,tbed.start)
#             stop  = min(stop,tbed.stop)
#             vals[(start-tbed.start):(stop-tbed.start)]+=depth
#         if forcestrand and tbed.strand=="-":
#             vals=vals[::-1]
#         return vals

#     def __enter__(self):
#         """Enter instance."""
#         return self

#     def __exit__(self, type, value, traceback):
#         """Exit instance."""
#         self.close()

#     def close(self):
#         """Close BigWig file."""
#         if not self.closed:
#             wWigIO.close(self.infile)
#             self.closed = True

#     def __del__(self): 
#         """Close BigWig file. Avoid memory leaks."""
#         self.close()

#     def wigToBigWig(wigfile, sizefile, bwfile):
#         """Convert wiggle file to BigWig file."""
#         wWigIO.wigToBigWig(wigfile, sizefile, bwfile)
#     wigToBigWig=staticmethod(wigToBigWig)
    
#     def bigWigToWig(bwfile, wigfile):
#         """Convert BigWig file to Wiggle file."""
#         wWigIO.bigWigToWig(bwfile, wigfile)
#     bigWigToWig=staticmethod(bigWigToWig)
