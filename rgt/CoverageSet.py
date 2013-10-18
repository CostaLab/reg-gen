from __future__ import print_function
from rgt.GenomicRegionSet import *
import pysam, sys  # @UnusedImport @UnresolvedImport
import numpy as np
import numpy.ma
import matplotlib  # @UnresolvedImport
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # @UnresolvedImport

"""
Represent coverage data.

Authors: Ivan G. Costa, Manuel Allhoff

CoverageSet has an instance of GenomicRegionSet. 
It has the coverage data for its GenomicRegionSet.
The coverage depends on the binsize.

Methods:

coverageFromBam(filename):
Compute coverage of GenomicRegionSet based on BAM file.

writeBed:
Output coverage in BED format.
"""

class CoverageSet:
    def __init__(self, name, GenomicRegionSet):
        """Initialize CoverageSet <name>."""
        self.name = name
        self.genomicRegions = GenomicRegionSet
        self.values = [] #coverage data for genomicRegions
        self.binsize = 1
        self.mapped_reads = 0 #number of mapped read
        self.reads = 0 #number of reads

    def coverage_from_bam(self, bamFile, readSize = 200, binsize = 50):
        """Return list of arrays describing the coverage of each genomicRegions from <bamFile>. 
        Consider reads in <bamFile> with a length of <readSize>.
        Divide the genomic regions in bins with a width of <binsize>."""
        self.binsize = binsize
        bam = pysam.Samfile(bamFile, "rb" )
        self.mapped_reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:3]) ) for l in pysam.idxstats(bamFile) ])
        self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile) ])
        
        for region in self.genomicRegions:
            cov = [0] * (len(region) / binsize)
            for read in bam.fetch(region.chrom, region.initial - readSize, region.final + readSize):
                if read.is_reverse is False:
                    pos = read.pos 
                    for i in range( max(0, pos - region.initial) / binsize, min(len(region), pos + readSize - region.initial) / binsize ):
                        cov[i] += 1
                else:
                    pos = read.pos + read.rlen
                    for i in range( max(0, pos - readSize - region.initial) / binsize, min(len(region), pos - region.initial) / binsize ):
                        cov[i] += 1

            self.values.append(np.array(cov))

        self.valuesorig = self.values[:]

    def write_bed(self, filename):
        """Output coverage in BED format"""
        with open(filename, 'w') as f:
            i = 0
            for region in self.genomicRegions:
                coverage = self.values[i]
                i += 1
                assert len(coverage) == (region.final - region.initial) / self.binsize
                for j in range(1, len(coverage) + 1):
                    if coverage[j-1] == 0:
                        continue
                    print(region.chrom, region.initial+ (j-1)*self.binsize, min(region.initial+j*self.binsize, region.final), \
                          coverage[j-1], sep='\t', file=f)

#    def normRPM(self):
#        self.values = self.values * (1000000.0 / self.mapped_reads)
#    
#    def normFPKM(self):
#        self.values = self.values * (1000000000.0 / self.mapped_reads)
#    
#    def normFactor(self, factor):
#        self.values = self.values * (1000000.0 / (self.step * factor))
#
#    def log(self):
#        print(self.name, numpy.min(self.values), file=sys.stderr)
#        self.values=numpy.log2(self.values + 0.01)
#
#    def mean(self):
#        return numpy.mean(self.values, axis=0)
#    
#    def diff(self, coverage):
#        self.valuesorig = self.valuesorig + coverage.valuesorig
#        self.values = self.values - coverage.values
#    
#    def abs(self):
#        self.values = abs(self.values)
#
#    def plot(self, log = False, name = None, outTxt = True):
#        if name is None:
#            name = self.name
#        mean = self.mean()
#        plt.plot(range(0, self.step * len(mean), self.step), mean)
#        plt.title("Mean Expression "+name)
#        plt.axis([0, len(mean) * self.step, mean.min(), mean.max()])
#        plt.savefig(name + ".pdf")
#        plt.close("all")
#        if outTxt:
#            f = open(name + ".txt", "w")
#            f.write("Pos\tValue\n")
#        for i, m in enumerate(mean):
#            f.write(str(i * self.step) + '\t' + str(m) + '\n')
#        f.close()
    

#def statisticsFromDiffProfiles(self,coverage,log=False,norm="RPM"):
#  compvalues=numpy.array(coverage.values,numpy.float)
#  values=numpy.array(self.values,numpy.float) 
#  if norm=="RPM":
#    values= values*(1000000.0/self.sizeBam)
#    compvalues= compvalues*(1000000.0/coverage.sizeBam)      
#  if norm=="FPKM":
#    values= values*(1000000000.0/self.sizeBam)
#    compvalues= compvalues*(1000000000.0/coverage.sizeBam)
#  values=abs((values)-(compvalues))
#  if log:     
#    values= numpy.log2(values+1)
#  return numpy.mean(values,axis=0), numpy.std(values,axis=0)      

        #for pileupcolumn in bam.pileup(region.chr, v-window, v+step+window,stepper="all"):  
        #  if (pileupcolumn.pos-region.initial-window) <= len(region):
        #    if (pileupcolumn.pos-region.initial-window) >= 0:           
        #      cov[j]+=pileupcolumn.n

        #for pileupcolumn in bam.pileup(region.chr, region.initial, region.final,stepper="all"):
        #  try:
        #    if (pileupcolumn.pos-region.initial) <= len(region):
        #      if (pileupcolumn.pos-region.initial) >= 0:           
        #        cov[(pileupcolumn.pos-region.initial)/step]+=pileupcolumn.n
        #  except:
        #    #print "Warning: pile up returned outside read",pileupcolumn.pos,region
        #    pass
#         self.window=window