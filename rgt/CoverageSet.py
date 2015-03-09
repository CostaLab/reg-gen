from __future__ import print_function
from rgt.GenomicRegionSet import *
import pysam, sys  # @UnresolvedImport
import numpy as np
import numpy.ma
import os
import tempfile

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

subtract:
Subtract coverage from CoverageSet object

scale:
Scale coverage by factor

"""

class CoverageSet:
    def __init__(self, name, GenomicRegionSet):
        """Initialize CoverageSet <name>."""
        self.name = name
        self.genomicRegions = GenomicRegionSet
        self.coverage = [] #coverage data for genomicRegions
        self.binsize = 100
        self.mapped_reads = None #number of mapped read
        self.reads = None #number of reads
        self.stepsize = 50
        

    
    def subtract(self, cs):
        """Substract CoverageSet <cs>, set negative values to 0."""
        cs_chroms = cs.genomicRegions.get_chrom()
        assert len(cs_chroms) == len(set(cs_chroms)) #no double entries
        assert len(self.genomicRegions.get_chrom()) == len(set(self.genomicRegions.get_chrom()))
        
        i = 0
        for c in self.genomicRegions.get_chrom(): #c corresponds to self.coverage[i]
            try:
                j = cs_chroms.index(c)
                assert len(self.coverage[i]) == len(cs.coverage[j])
                self.coverage[i] -= cs.coverage[j]
                self.coverage[i] = self.coverage[i].clip(0, max(max(self.coverage[i]), 0)) #neg. values to 0
            except ValueError:
                pass
            i += 1
    
    def add(self, cs):
        """Substract CoverageSet <cs>, set negative values to 0."""
        cs_chroms = cs.genomicRegions.get_chrom()
        assert len(cs_chroms) == len(set(cs_chroms)) #no double entries
        assert len(self.genomicRegions.get_chrom()) == len(set(self.genomicRegions.get_chrom()))
        
        i = 0
        for c in self.genomicRegions.get_chrom(): #c corresponds to self.coverage[i]
            try:
                j = cs_chroms.index(c)
                assert len(self.coverage[i]) == len(cs.coverage[j])
                self.coverage[i] += cs.coverage[j]
            except ValueError:
                pass
            i += 1
    
    def scale(self, factor):
        """Scale coverage with <factor>"""
#        print(factor)
        for i in range(len(self.coverage)):
            self.coverage[i] = np.rint(self.coverage[i] * float(factor)).astype(int)

    def normRPM(self):
        factor=1000000/float(self.reads)
        self.coverage=np.array(self.coverage)*factor

    def write_bed(self, filename, zero=False):
        """Output coverage in BED format"""
        with open(filename, 'w') as f:
            i = 0
            for region in self.genomicRegions:
                c = self.coverage[i]
                i += 1
                for j in range(len(c)):
                    if zero:
                        print(region.chrom, j * self.stepsize + ((self.binsize-self.stepsize)/2) + region.initial, \
                              j * self.stepsize + ((self.binsize+self.stepsize)/2) + region.initial, c[j], sep='\t', file=f)
                    else:
                        if c[j] != 0:
                            print(region.chrom, j * self.stepsize + ((self.binsize-self.stepsize)/2) + region.initial, \
                                  j * self.stepsize + ((self.binsize+self.stepsize)/2) + region.initial, c[j], sep='\t', file=f)
                        
    def write_wig(self, filename):
        with open(filename, 'w') as f:
            i = 0
            for region in self.genomicRegions:
                print('variableStep chrom=' + str(region.chrom) + ' span=' +str(self.stepsize), file=f)
                c = self.coverage[i]
                i += 1
                for j in range(len(c)):
                    if c[j] != 0:
                        print(j * self.stepsize + ((self.binsize-self.stepsize)/2), c[j], file=f)
    
    def write_bigwig(self, filename, chrom_file, save_wig=False):
        if save_wig:
            tmp_path = filename + '.wig'
            self.write_wig(tmp_path)
            t = ['wigToBigWig', "-clip", tmp_path, chrom_file, filename] #TODO: something is wrong here, call only wigToBigWig
            c = " ".join(t)
            os.system(c)
        else:
            _, tmp_path = tempfile.mkstemp()
            self.write_wig(tmp_path)
            t = ['wigToBigWig', "-clip", tmp_path, chrom_file, filename] #TODO: something is wrong here, call only wigToBigWig
            c = " ".join(t)
            #print(c)
            os.system(c)
            os.remove(tmp_path)

    def coverage_from_genomicset(self,bamFile,readSize=200):
        bam = pysam.Samfile(bamFile, "rb" )
        self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile)])
        cov=[0]*len(self.genomicRegions)
        for i,region in enumerate(self.genomicRegions):
            try:
                for r in bam.fetch(region.chrom,region.initial-readSize,region.final+readSize):
                    cov[i] += 1
            except:
                print("\tError: "+str(region))
        self.coverage=cov 
        self.coverageOrig=cov

    def _get_bedinfo(self, l):
        if len(l) > 1:
            l.strip()
            l = l.split('\t')
            return l[0], int(l[1]), int(l[2]), True
        else:
            return -1, -1, -1, False

    def coverage_from_bam(self, bam_file, read_size = 200, binsize = 100, stepsize = 50, rmdup = True, mask_file = None):
        """Return list of arrays describing the coverage of each genomicRegions from <bam_file>. 
        Consider reads in <bam_file> with a extension size of <read_size>.
        Remove duplicates (read with same position) with rmdup=True (default).
        Divide the genomic regions in bins with a width of <binsize> and use <stepsize> to smooth the signal."""
        self.binsize = binsize
        self.stepsize = stepsize
        
        bam = pysam.Samfile(bam_file, "rb" )
        for read in bam.fetch():
            read_size += read.rlen
            break
        self.mapped_reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:3]) ) for l in pysam.idxstats(bam_file) ])
        self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bam_file) ])
        #print("Loading reads of %s..." %self.name, file=sys.stderr)
        
        #check whether one should mask
        next_it = True
        if mask_file is not None and os.path.exists(mask_file):
            mask = True
            f = open(mask_file, 'r')
            c_help, s_help, e_help = self.genomicRegions.sequences[0].chrom, -1, -1
        else:
            mask = False
        
        chrom_regions = [r.chrom for r in self.genomicRegions.sequences] #chroms by regions
        
        for region in self.genomicRegions:
            cov = [0] * (len(region) / stepsize)

            positions = []
            j = 0
            read_length = -1
            try:
                for read in bam.fetch(region.chrom, max(0, region.initial-read_size), region.final+read_size):
                    
                    j += 1
                    read_length = read.rlen 
                    if not read.is_unmapped:
                        pos = read.pos - read_size if read.is_reverse else read.pos
                        pos_help = read.pos - read.qlen if read.is_reverse else read.pos
                        
                        #if position in mask region, then ignore
                        if mask:
                            while next_it and c_help not in chrom_regions: #do not consider this deadzone
                                c_help, s_help, e_help, next_it = self._get_bedinfo(f.readline())
                            if c_help != -1 and chrom_regions.index(region.chrom) >= chrom_regions.index(c_help): #deadzones behind, go further
                                while next_it and c_help != region.chrom: #get right chromosome
                                    c_help, s_help, e_help, next_it = self._get_bedinfo(f.readline())
                            while next_it and e_help <= pos_help and c_help == region.chrom: #check right position
                                c_help, s_help, e_help, next_it = self._get_bedinfo(f.readline())
                            if next_it and s_help <= pos_help and c_help == region.chrom:
                                continue #pos in mask region
                        
                        positions.append(pos)
            except ValueError:
                pass
            if rmdup:
                positions = list(set(positions))
            positions.sort()
            positions.reverse()
            
            i = 0
            while positions:
                win_s = max(0, i * stepsize - binsize*0.5) + region.initial
                win_e = i * stepsize + binsize*0.5 + region.initial 
                c = 0
                taken = []
                while True:
                    s = positions.pop()
                    
                    taken.append(s)
                    if s < win_e: #read within window
                        c += 1
                    if s >= win_e or not positions:
                        taken.reverse()
                        for s in taken:
                            if s + read_size + read_length >= win_s: #consider read in next iteration
                                positions.append(s)
                            else:
                                break #as taken decreases monotonously
                        taken = []
                        break
                
                if i < len(cov):
                    cov[i] = c

                i += 1

            self.coverage.append(np.array(cov))

        self.coverageorig = self.coverage[:]
        
        
#     def coverage_from_bam(self, bamFile, readSize = 200, binsize = 100, stepsize = 50, rmdup = False, mask_file=None, get_pos=False):
#         """Return list of arrays describing the coverage of each genomicRegions from <bamFile>. 
#         Consider reads in <bamFile> with a length of <readSize>.
#         Remove duplicates (read with same position) with rmdup=True, else rmdup=False (default).
#         Do not consider reads, that originate from positions described by ,mask_file>.
#         If <get_pos> computes list of forward and backward reads for each bin.
#         Divide the genomic regions in bins with a width of <binsize>."""
#         self.binsize = binsize
#         bam = pysam.Samfile(bamFile, "rb" )
#         self.mapped_reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:3]) ) for l in pysam.idxstats(bamFile) ])
#         self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile) ])
#         
#         last_pos, last_chrom, next_it = -1, -1, True
#         
#         #check whether one should mask
#         if mask_file is not None and os.path.exists(mask_file):
#             mask = True
#             f = open(mask_file, 'r')
#             c, s, e = self.genomicRegions.sequences[0].chrom, -1, -1
#         else:
#             mask = False
# 
#         step = 100
#         self.step = step
#         for region in self.genomicRegions:
#             cov = [0] * (len(region) / step)
#                 
#             print("load reads of %s" % region.chrom, file=sys.stderr)
#             
#             for i in range(int(self.binsize/self.step / 2), len(cov)):
#                 f_pos, b_pos = set(), set()
#                 if i % 10000 ==0 :
#                     print(i, len(cov), file=sys.stderr)
#                 for read in bam.fetch(region.chrom, i*step - self.binsize/2, i*step + self.binsize/2):
#                     if read.is_reverse:
#                         b_pos.add(read.pos)
#                     else:
#                         f_pos.add(read.pos)
#                     
#                     cov[i] = len(f_pos) + len(b_pos)
#             
#             self.coverage.append(np.array(cov))
            
#             for read in bam.fetch(region.chrom, max(0, region.initial - readSize), region.final + readSize):
#                 pos = read.pos if not read.is_reverse else read.pos + read.rlen #get right or left position of read
#                 
#                 if rmdup and last_chrom == region.chrom and last_pos == pos: 
#                     continue #rmdup
#                 
#                 if mask:
#                     while next_it and c != region.chrom: #get right chromosome
#                         c, s, e, next_it = self._get_bedinfo(f.readline())
#                     while next_it and e <= pos: #check right position
#                         c, s, e, next_it = self._get_bedinfo(f.readline())
#                     if next_it and s <= pos:
#                         continue #pos in mask region
# 
#                 if read.is_reverse is False:
#                     for i in range( max(0, pos - region.initial) / binsize, min(len(region), pos + readSize - region.initial) / binsize + 1 ):
#                         cov[i] += 1
#                         if get_pos:
#                             positions[i][0].append(pos)
#                 else:
#                     for i in range( max(0, pos - readSize - region.initial) / binsize, min(len(region), pos - region.initial) / binsize + 1 ):
#                         cov[i] += 1
#                         if get_pos:
#                             positions[i][1].append(pos)
# 
#                 last_pos, last_chrom = pos, region.chrom
# 
#             self.coverage.append(np.array(cov))
#             
#             if get_pos:
#                 self.read_pos.append( np.array(positions) )
#         
#         self.coverageorig = self.coverage[:]
       
#     def _get_bedinfo(self, l):
#         if l != "":
#             l.strip()
#             l = l.split('\t')
#             return l[0], int(l[1]), int(l[2]),True
#         else:
#             return -1, -1, -1, False
              
#     def coverage_from_bam(self, bamFile, readSize = 200, binsize = 50, rmdup = False, mask_file=None, get_pos=False):
#         """Return list of arrays describing the coverage of each genomicRegions from <bamFile>. 
#         Consider reads in <bamFile> with a length of <readSize>.
#         Remove duplicates (read with same position) with rmdup=True, else rmdup=False (default).
#         Do not consider reads, that originate from positions described by ,mask_file>.
#         If <get_pos> computes list of forward and backward reads for each bin.
#         Divide the genomic regions in bins with a width of <binsize>."""
#         self.binsize = binsize
#         bam = pysam.Samfile(bamFile, "rb" )
#         self.mapped_reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:3]) ) for l in pysam.idxstats(bamFile) ])
#         self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile) ])
#         
#         last_pos, last_chrom, next_it = -1, -1, True
#         
#         #check whether one should mask
#         if mask_file is not None and os.path.exists(mask_file):
#             mask = True
#             f = open(mask_file, 'r')
#             c, s, e = self.genomicRegions.sequences[0].chrom, -1, -1
#         else:
#             mask = False
#         
#         for region in self.genomicRegions:
#             cov = [0] * (len(region) / binsize)
#             if get_pos:
#                 positions = [([],[]) for _ in range(len(region) / binsize) ] #positions of forward/backward reads within bin
#             
#             if len(region) % binsize != 0: #end of chromosome
#                 cov += [0]
#                 if get_pos:
#                     positions += [([],[])]
#                 
#             print("load reads of %s" % region.chrom, file=sys.stderr)
#             
#             for read in bam.fetch(region.chrom, max(0, region.initial - readSize), region.final + readSize):
#                 pos = read.pos if not read.is_reverse else read.pos + read.rlen #get right or left position of read
#                 
#                 if rmdup and last_chrom == region.chrom and last_pos == pos: 
#                     continue #rmdup
#                 
#                 if mask:
#                     while next_it and c != region.chrom: #get right chromosome
#                         c, s, e, next_it = self._get_bedinfo(f.readline())
#                     while next_it and e <= pos: #check right position
#                         c, s, e, next_it = self._get_bedinfo(f.readline())
#                     if next_it and s <= pos:
#                         continue #pos in mask region
# 
#                 if read.is_reverse is False:
#                     for i in range( max(0, pos - region.initial) / binsize, min(len(region), pos + readSize - region.initial) / binsize + 1 ):
#                         cov[i] += 1
#                         if get_pos:
#                             positions[i][0].append(pos)
#                 else:
#                     for i in range( max(0, pos - readSize - region.initial) / binsize, min(len(region), pos - region.initial) / binsize + 1 ):
#                         cov[i] += 1
#                         if get_pos:
#                             positions[i][1].append(pos)
# 
#                 last_pos, last_chrom = pos, region.chrom
# 
#             self.coverage.append(np.array(cov))
#             
#             if get_pos:
#                 self.read_pos.append( np.array(positions) )
#         
#         self.coverageorig = self.coverage[:]

#     def write_bed(self, filename):
#         """Output coverage in BED format"""
#         with open(filename, 'w') as f:
#             i = 0
#             for region in self.genomicRegions:
#                 c = self.coverage[i]
#                 i += 1
# #                 assert len(c) == (region.final - region.initial) / self.binsize + 1
#                 for j in range(1, len(c) + 1):
#                     if c[j-1] == 0:
#                         continue
#                     print(region.chrom, region.initial+ (j-1)*self.step, min(region.initial+j*self.step, region.final), \
#                           c[j-1], sep='\t', file=f)

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


#     def mask(self, mask_region):
#         """mask coverage with 0 in GenomicRegion mask_region"""
#         chrom = 0
#         for region in self.genomicRegions:
#             if region.chrom != mask_region.chrom: #jump over chromosomes
#                 chrom += 1 
#             else:
#                 break
#             
#         start = mask_region.initial / self.binsize
#         end = mask_region.final / self.binsize
#         self.coverage[chrom][max(0,start) : end] = 0






