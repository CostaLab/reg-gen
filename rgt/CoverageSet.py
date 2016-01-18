"""
CoverageSet
===================
CoverageSet represents coverage data. CoverageSet has an instance of GenomicRegionSet. It has the coverage data for its GenomicRegionSet.
The coverage depends on the binsize.

Author: Ivan G. Costa, Manuel Allhoff
"""

from __future__ import print_function
from rgt.GenomicRegionSet import *
import pysam, sys  # @UnresolvedImport
import numpy as np
import numpy.ma
import os
import tempfile
import wWigIO
#from ngslib import BigWigFile

class BigWigFile():
    '''
    Fast reader of BigWig format file.
    Usage:
        Open file:
            fh=BigWigFile("test.bw")
        Get chromosome sizes:
            chroms=fh.chromSizes() # (chroms,sizes)
        Fetch regions:
            wigs=fh.fetch(chrom="chr1",start=100,stop=200)
            for wig in wigs:
                #do some thing with wig
                print wig.chrom,wig.start,wig.stop,wig.score
        Close file:
            fh.close()
    
    Parameters:
        chrom=None: return empty list.
        start=None: start at first position.
        stop=None:  stop at the end of chromosome.
    '''
    def __init__(self,infile,chrom_size=None):
        ''' Open BigWig file. '''
        # Check if file exists
        self.closed = True
        if not os.path.exists(infile) and infile != 'stdin':
            raise IOError("ERROR: input file '{0}' dose not exist.".format(infile))
        if infile.endswith('.wig'):
            bwfile = os.path.splitext(infile)[0]+".bw"
            if os.path.exists(bwfile):
                self.infile = bwfile
            else:
                if chrom_size is None:
                    raise IOError("ERROR: 'chrom_size' file is required for wig -> bigwig conversion.")
                BigWigFile.wigToBigWig(infile,chrom_size,bwfile)
            self.infile=bwfile
        else:
            self.infile=infile
        wWigIO.open(self.infile)
        self.closed = False
        self.chroms, self.sizes = self.chromSizes()
        self.closed = False

    def chromSizes(self):
        ''' Get chromosome sizes.'''
        chroms,sizes=wWigIO.getChromSize(self.infile)
        return chroms,sizes

    def fetch(self, chrom,start=None,stop=None, strand="+", zerobased=True,**kwargs):
        '''
        Fetch intervals in a given region. 
        Note: strand is useless here.
        Parameters:
            start: int or None
                start position, None for start=0
            end: int or None
                end position, None for stop = end of chrom
            strand: string
                choice from '+', '-' and '.'
            zerosbased: bool
                indices are zerobased or not. Useless here.
        Dictionary parameters:
            chunksize: int
                chunk size when fetch items from a large region.
        Generates:
            wig: tuple
                (start, stop, value)
        '''
        if chrom is None: raise ValueError("ERROR: chrom name is required.")
        if start is None: start = 0
        if not zerobased: start -= 1
        if start <0:
            raise ValueError('ERROR: start should be >=0 (zerobased=True) or >=1 (zerobased=False).')
        if stop is None: stop = self.sizes[self.chroms.index['chrom']]
        chunksize = kwargs.get('chunksize',10000)
        try:
            for cstart in xrange(start,stop,chunksize):
                cend = cstart + chunksize
                if cend > stop: cend = stop
                for wig in wWigIO.getIntervals(self.infile,chrom,cstart,cend):
                    if wig[0] < cstart: wig[0] = cstart
                    if wig[1] > cend: wig[1] = cend
                    yield wig
        except:
            # check result
            if wigs == 1:
                raise ValueError("ERROR: wWigIO.getIntervals doesn't have correct parameters.")
            if wigs == 2:
                raise ValueError("ERROR: BigWig file cannot be opened.")

    def pileup(self,chrom,start=None,stop=None,strand="+",zerobased=True,**kwargs):
        '''
        Fetch depth for a genomic region.        
        '''
        if chrom is None: raise ValueError("ERROR: chrom name is required.")
        if start is None: start = 0
        if not zerobased: start -= 1
        if stop is None: stop = self.sizes[self.chroms.index(chrom)]
        vals = numpy.zeros(stop-start)
        for wstart,wstop,depth in self.fetch(chrom,start,stop,strand,zerobased):
            vals[(wstart-start):(wstop-start)]+=depth
        if strand == "-":
            vals = vals[::-1]
        return vals

    def fetchBed(self,tbed,byPos=False,forcestrand=True):
        '''Fetch intervals for Bed.'''
        wigs=wWigIO.getIntervals(self.infile,tbed.chrom,tbed.start,tbed.stop)
        if not byPos:
            return wigs
        # get depth by position, return a numpy array.
        vals = numpy.zeros(tbed.length())
        for start,stop,depth in wigs:
            start = max(start,tbed.start)
            stop  = min(stop,tbed.stop)
            vals[(start-tbed.start):(stop-tbed.start)]+=depth
        if forcestrand and tbed.strand=="-":
            vals=vals[::-1]
        return vals

    def __enter__(self):
        ''' Enter instance. '''
        return self

    def __exit__(self, type, value, traceback):
        ''' Exit instance. '''
        self.close()

    def close(self):
        ''' Close BigWig file. '''
        if not self.closed:
            wWigIO.close(self.infile)
            self.closed = True

    def __del__(self): 
        ''' Close BigWig file.  Avoid memory leaks.'''
        self.close()

    def wigToBigWig(wigfile, sizefile, bwfile):
        ''' Convert wiggle file to BigWig file. '''
        wWigIO.wigToBigWig(wigfile, sizefile, bwfile)
    wigToBigWig=staticmethod(wigToBigWig)
    
    def bigWigToWig(bwfile, wigfile):
        ''' Convert BigWig file to Wiggle file. '''
        wWigIO.bigWigToWig(bwfile, wigfile)
    bigWigToWig=staticmethod(bigWigToWig)


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
                        
    def write_wig(self, filename, end):
        f = open(filename, 'w')
        i = 0
        for region in self.genomicRegions:
            print('variableStep chrom=' + str(region.chrom) + ' span=' +str(self.stepsize), file=f)
            c = self.coverage[i]
            i += 1
            for j in range(len(c)):
                if c[j] != 0:
                    print(j * self.stepsize + ((self.binsize-self.stepsize)/2), c[j], file=f)
        f.close()
    
    def write_bigwig(self, filename, chrom_file, end=True, save_wig=False):
        if save_wig:
            tmp_path = filename + '.wig'
            self.write_wig(tmp_path)
            t = ['wigToBigWig', "-clip", tmp_path, chrom_file, filename] #TODO: something is wrong here, call only wigToBigWig
            c = " ".join(t)
            os.system(c)
        else:
            _, tmp_path = tempfile.mkstemp()
            self.write_wig(tmp_path, end)
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

    def coverage_from_bam(self, bam_file, read_size = 200, binsize = 100, stepsize = 50, rmdup = True, mask_file = None, get_strand_info = False):
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
        
        if get_strand_info:
            self.cov_strand_all = []
        
        for region in self.genomicRegions:
            cov = [0] * (len(region) / stepsize)
            
            if get_strand_info:
                cov_strand = [[0,0]] * (len(region) / stepsize)
                strand_info = {}
            
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
                        
                        if get_strand_info:
                            if pos not in strand_info:
                                strand_info[pos] = (1,0) if not read.is_reverse else (0,1)
                        
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
                if get_strand_info:
                    sum_strand_info = [0,0]
                
                taken = []
                while True:
                    s = positions.pop()
                    
                    taken.append(s)
                    if s < win_e: #read within window
                        c += 1
                        if get_strand_info:
                            sum_strand_info[0] += strand_info[s][0]
                            sum_strand_info[1] += strand_info[s][1]
                        
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
                    if get_strand_info:
                        cov_strand[i] = sum_strand_info
                i += 1

            self.coverage.append(np.array(cov))
            if get_strand_info:
                self.cov_strand_all.append(np.array(cov_strand))
            
        self.coverageorig = self.coverage[:]
        self.overall_cov = reduce(lambda x,y: np.concatenate((x,y)), [self.coverage[i] for i in range(len(self.genomicRegions))])
    
    def index2coordinates(self, index, r):
        """Translate index within coverage array to genomic coordinates."""
        iter = r.__iter__()
        r = iter.next()
        sum = r.final
        last = 0
        i = 0
        while sum <= index * self.stepsize:
            last += len(self.coverage[i])
            try:
                r = iter.next()
            except StopIteration:
                sum += r.final
                i += 1
                break
            sum += r.final
            i += 1
        
        return r.chrom, (index-last) * self.stepsize, \
            min((index-last) * self.stepsize + self.stepsize, r.final)
    
    def coverage_from_bigwig(self, bigwig_file, stepsize=100):
        """Return list of arrays describing the coverage of each genomicRegions from <bigwig_file>."""
        self.coverage = []
        bwf = BigWigFile(bigwig_file)
        #ds = []
        for gr in self.genomicRegions:
            #print(".", end="")
            depth = bwf.pileup(gr.chrom, gr.initial-stepsize/2, gr.final+stepsize/2)
            #ds = []
            ds = [depth[d] for d in range(0, gr.final-gr.initial, stepsize)]
            #for i in range(0, gr.final-gr.initial):
            #    d = [ depth[j] for j in range(i,i+stepsize) ]
            #    ds.append(sum(d)/len(d))

            #if gr.orientation == "-":
            #    self.coverage.append( np.array(list(reversed(ds))) )
            #else:
            self.coverage.append( np.array(ds) )
        #print(len(ds))
        bwf.close()
        
    def phastCons46way_score(self, stepsize=100):
        """Load the phastCons46way bigwig files to fetch the scores as coverage"""
        self.coverage = []
        phastCons46way_dir = "/data/phastCons46way/"
        for gr in self.genomicRegions:
            bwf = BigWigFile(os.path.join(phastCons46way_dir, gr.chrom+".phastCons46way.bw"))
            depth = bwf.pileup(gr.chrom, gr.initial-stepsize/2, gr.final+stepsize/2)
            ds = []
            for i in range(0, gr.final-gr.initial):
                d = [ depth[j] for j in range(i,i+stepsize) ]
                ds.append(sum(d)/len(d))
                
            if gr.orientation == "-":
                self.coverage.append( np.array(list(reversed(ds))) )
            else:
                self.coverage.append( np.array(ds) )

            bwf.close()

    def count_unique_reads(self, bamFile):
        """Count the number of unique reads on the given GenomicRegionSet"""
        bam = pysam.Samfile(bamFile, "rb" )
        #self.reads = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile)])
        reads = []
        for i,region in enumerate(self.genomicRegions):
            #print(bam.fetch(region.chrom,region.initial,region.final))
            #try:
            for r in bam.fetch(region.chrom,region.initial,region.final):
                reads.append(r.qname)
            #except:
            #    print("\tError: "+str(region))

        reads = list(set(reads))
        return len(reads)
        


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






