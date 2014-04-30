from __future__ import print_function
from __future__ import division
from GenomicRegion import *
import random
import copy
import os
import time
from scipy import stats
from Util import GenomeData
from Util import OverlapType

"""
Represent list of GenomicRegions.

Authors: Ivan G. Costa, Manuel Allhoff, Joseph Kuo

Define a list of GenomicRegions. Perform operations of this set.

Methods:

readBed(filename):
Add GenomicRegions defined by row in bed file.

intersect(region):
Return intersection with region.

"""

class GenomicRegionSet:

    def __init__(self, name):
        """Initialize a set of GenomicRegions"""
        self.name = name
        self.sequences = []
        self.sorted = False
        self.fileName=""
        self.genome_path = ""
    
    def get_chrom(self):
        """Return all chromosomes"""
        return [ r.chrom for r in self.sequences ]
    
    def add(self, region):
        """Add GenomicRegion"""
        self.sequences.append(region)
        self.sorted = False
    
    def __len__(self):
        return len(self.sequences)
    
    def __iter__(self):
        return iter(self.sequences)

    def extend(self,left,right):
        """Perform extend step for every element"""
        for s in self.sequences:
            s.extend(left,right)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicRegion.__cmp__)

    def read_bed(self, filename):
        """Read BED file and add every row as a GenomicRegion. 
        Chrom (1), start (2), end (2), name (4) and orientation (6) is used for GenomicRegion. 
        All other columns (5, 7, 8, ...) are put to the data variable of the GenomicRegion.
        The numbers in parentheses are the columns of the BED format.
        See BED format at: http://genome.ucsc.edu/FAQ/FAQformat.html#format1 """
        self.fileName=filename
        with open(filename) as f:
            for line in f:
                try:
                    name, orientation, data = None, None, None
                    line = line.strip("\n")
                    line = line.split("\t")
                    size = len(line)
                    chrom = line[0]
                    start, end = int(line[1]), int(line[2])
                    if start > end:
                        start, end =  end, start
                    if size > 3:
                        name = line[3]
                    if size > 5:
                        orientation = line[5]
                    if size > 5:
                        data = "\t".join( [line[4]] + line[6:] )
                    self.add( GenomicRegion(chrom, start, end, name, orientation, data) )
                except:
                    print("Error at line",line,self.fileName)
            self.sort()
  
    def read_bedgraph(self, filename):
        """Read BEDGRAPH file and add every row as a GenomicRegion. 
        See BEDGRAPH format at: http://genome.ucsc.edu/goldenPath/help/bedgraph.html """
        self.fileName=filename
        with open(filename) as f:
            for line in f:
                try:
                    line = line.strip("\n")
                    line = line.split("\t")
                    assert len(line) == 4

                    chrom, start, end, data = line[0], int(line[1]), int(line[2]), float(line[3])

                    self.add( GenomicRegion(chrom=chrom, initial=start, final=end, data=data) )
                except:
                    print("Error at line", line, self.fileName)
                    
            self.sort()
            
    def random_subregions(self,size):
        """Return a subsampling of the genomic region set with a specific number of regions"""
        z = GenomicRegionSet(self.name + '_random')
        samp = random.sample(range(len(self)),size)
        for i in samp:
            z.add(self.sequences[i])
        z.sort()
        return z                

    def write_bed(self,filename):
        """Write GenomicRegions to BED file"""
        with open(filename, 'w') as f:
            for s in self:
                print(s, file=f)

    def filter_by_gene_association(self,fileName,geneSet,geneAnnotation,genomeSize,promoterLength=1000,threshDist=50000):
        """code based on eduardos functions. This should be  integrated in the core framework soon
        TODO: Eduardo should check this!"""
        from rgt.motifanalysis.util import bedFunctions, sort
        from rgt.motifanalysis.enrichment.geneAssociation import *
        self.fileName=fileName
        regionsToGenes={}
        coordDict = bedFunctions.createBedDictFromSingleFile(fileName, features=[1,2,3,4,5]) 
        coordDict = sort.sortBedDictionary(coordDict, field=0)
        [dictBed,allBed] = geneAssociationByPromoter(coordDict,geneSet,geneAnnotation,genomeSize,promoterLength,threshDist)  
        #print dictBed
        genes=[]
        totalPeaks=0
        allgenes=[]
        for chr in dictBed.keys():
            for (v1,v2,name,orientation,data) in dictBed[chr]:
                totalPeaks+=1
                names=name.split(":")
                keep=[n for n in names if "." not in n]
                if len(keep) > 0:
                    self.add(GenomicRegion(chr,v1,v2,":".join(keep)))
                genes = genes+keep
                allgenes=allgenes+[n.strip(".") for n in names]
                regionsToGenes[chr+":"+str(v1)+"-"+str(v2)]=[n.strip(".") for n in names]
        #print "Total Peaks", total
        mappedGenes=len(list(set(allgenes)))
        self.sort()
        self.genes=list(set(genes))
        if geneSet == None:
          le=0
        else:
          le=len(geneSet)
        return le, len(self.genes), mappedGenes, totalPeaks,regionsToGenes

    def intersect(self,y,mode=OverlapType.OVERLAP):
        """Return the overlapping regions with three different modes.
        
        (mode = OverlapType.OVERLAP) 
        Return new GenomicRegionSet including only the overlapping regions with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
        
            Return:
            z -- a new GenomicRegionSet including only the overlapping regions
        
            Graphical explanation:
            self       ----------              ------
            y                 ----------                    ----
            Result            ---
            
        (mode = OverlapType.ORIGINAL)
        Return the regions of original GenomicRegionSet which have any intersections with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- the regions of original GenomicRegionSet which have any intersections with y
            
            Graphical explanation:
            self       ----------              ------
            y              ----------                    ----
            Result     ----------
            
        (mode = OverlapType.COMP_INCL)
        Return region(s) of the GenomicRegionSet which are 'completely' included by y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- region(s) of self which are 'completely' included by y
            
            Graphical explanation:
            self        -------------             ------
            y              ----------      ---------------              ----
            Result                                ------
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet('None region')
        z = GenomicRegionSet(self.name + '_' + y.name)
        
        # If there is overlap within self or y, they should be merged first. 
        self.sort()
        y.sort()
        # Iteration
        con_self = self.__iter__()
        con_y = y.__iter__()
        s = con_self.next()
        ss = con_y.next()
        cont_loop = True
        pre_z = GenomicRegion(chrom="chr1",initial=0,final=0) # Store the previous region which add in z
        while cont_loop:
            # When the regions overlap
            if s.overlap(ss):
                if mode == OverlapType.OVERLAP:
                    sss = GenomicRegion(chrom=s.chrom, 
                                        initial=max(s.initial, ss.initial), 
                                        final=min(s.final, ss.final))
                    z.add(sss)
                elif mode == OverlapType.ORIGINAL:
                    if s == pre_z:
                        pass
                    else:
                        z.add(s)
                        pre_z = s
                elif mode == OverlapType.COMP_INCL:
                    if s.initial >= ss.initial and s.final <= ss.final:
                        z.add(s)
                    else:
                        pass    
                if s.final >= ss.final:
                    try:
                        ss = con_y.next()
                    except:
                        cont_loop = False
                else:
                    try:
                        s = con_self.next()
                    except:
                        cont_loop = False   
            # When the region have no overlap
            elif s < ss:
                try:
                    s = con_self.next()
                except:
                    cont_loop = False
            elif s > ss:
                try:      
                    ss = con_y.next()
                except:
                    cont_loop = False
        z.sort()
        return z
    
    def closest(self,y):
        """Return a new GenomicRegionSet including the region(s) of y which is closest to any self region. 
           If there are intersection, return False.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        
        Return:
        z -- the region(s) which is nearest to self
        
        """
        if self.__len__() == 0 or y.__len__() == 0:
            return GenomicRegionSet('Empty set') 
        elif self.intersect(y).__len__() != 0:
            return False
        else:
            z_dict = {}  # For storing the distance and the regions
            self.sort()
            y.sort()
            con_self = self.__iter__()
            con_y = y.__iter__()
            s = con_self.next()
            ss = con_y.next()
            cont_loop = True
            while cont_loop:
                if s.chrom == ss.chrom:
                    # ----
                    #        ------
                    if s < ss:
                        z_dict[ss.initial - s.final] = ss
                        try: s = con_self.next()
                        except: cont_loop = False
                    #          ------
                    #  -----
                    elif s > ss:
                        z_dict[s.initial - ss.final] = ss
                        try: ss = con_y.next()
                        except: cont_loop = False
                elif s.chrom != ss.chrom:
                    if s < ss:
                        try: s = con_self.next()
                        except: cont_loop = False
                    elif s > ss:
                        try: ss = con_y.next()
                        except: cont_loop = False
            
            if len(z_dict.keys()) == 0:
                return GenomicRegionSet('Empty set')
            else:
                minimum_distance = min(z_dict.keys())
                z = GenomicRegionSet('Closest region')
                z.add(z_dict[minimum_distance])
                return z
            
    def remove_duplicates(self):
        """Remove the duplicate regions and remain the unique regions. (No return)"""
        self.sort()
        for i in range(len(self.sequences) - 1):
            try:
                if self.sequences[i] == self.sequences[i+1]:
                    del self.sequences[i+1]
            except:
                continue
            
    def window(self,y,adding_length = 1000):
        """Return the overlapping regions of self and y with adding a specified number 
        (1000, by default) of base pairs upstream and downstream of each region in self. 
        In effect, this allows regions in y that are near regions in self to be detected.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        adding_length -- the length of base pairs added to upstream and downstream of self (default 1000)
        
        Return:
        z -- a GenomicRegionSet including the regions of overlapping between extended self and original y.
        
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet('None region')
        # Establish an extended GenomicRegionSet
        extended_self = copy.deepcopy(self)
        for i in extended_self:
            i.extend(adding_length,adding_length)
        # Find their intersections
        return extended_self.intersect(y)
    
    def subtract(self,y):
        """Return a GenomicRegionSet excluded the overlapping regions with y.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to subtract by
        
        Return:
        z -- the remaining regions of self after subtraction
        
        Graphical explanation:
        self     ----------              ------
        y               ----------                    ----
        Result   -------                 ------
        
        """
        if len(self) == 0:
            return GenomicRegionSet('None region')
        elif len(y) == 0:
            return self
        else:
            z = GenomicRegionSet('Subtracted RegionSet')
            self.sort()
            y.sort()
            con_self = self.__iter__()
            con_y = y.__iter__()
            s = con_self.next()
            ss = con_y.next()
            cont_loop = True
            while cont_loop:
                # When the regions overlap
                if s.overlap(ss):
                    if s.initial < ss.initial:
                        z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=ss.initial))
                    if s.final > ss.final:
                        s.initial = ss.final
                        try:
                            ss = con_y.next()
                            pass
                        except:
                            z.add(GenomicRegion(chrom=s.chrom, initial=ss.final, final=s.final))
                            try:
                                s = con_self.next()
                                pass
                            except:
                                cont_loop = False
                    else:
                        try:
                            s = con_self.next()
                            pass
                        except:
                            cont_loop = False 
                # When the region have no overlap
                elif s < ss:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final))
                    try:
                        s = con_self.next()
                        pass
                    except:
                        try:
                            ss = con_y.next()
                            pass
                        except:
                            cont_loop = False
                elif s > ss:
                    try:
                        ss = con_y.next()
                        pass
                    except:
                        z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final))
                        try:
                            s = con_self.next()
                            pass
                        except:
                            cont_loop = False
            z.sort()
            return z
    
    def merge(self):
        """Merge the regions within the GenomicRegionSet"""
        self.sort()
        
        if len(self.sequences) in [0, 1]:
            return self
        else:
            z = GenomicRegionSet(name = self.name)
            prev_region = self.sequences[0]
            for cur_region in self.sequences[1:]:
                if prev_region.overlap(cur_region):
                    prev_region.initial = min(prev_region.initial, cur_region.initial)
                    prev_region.final = max(prev_region.final, cur_region.final)
                    prev_region.data += '_$_' + cur_region.data #use extra character to distinguish data sets
                else:
                    z.add(prev_region)
                    prev_region = cur_region
            z.add(prev_region)
            self.sequences = z.sequences
    
    def combine(self,region_set):
        """ Adding another GenomicRegionSet without merging the overlapping regions. """
        self.sequences.extend(region_set.sequences)
        self.sorted = False
        
    def cluster(self,max_distance):
        """Cluster the regions with a certain distance and return the result as a new GenomicRegionSet.
        
        Keyword arguments:
        max_distance -- the maximum distance between regions within the same cluster
        
        Return:
        z -- a GenomicRegionSet including clusters
        
        Graphical explanation:
        self           ----           ----            ----
                          ----             ----                 ----
        Result(d=1)    -------        ---------       ----      ----
        Result(d=10)   ---------------------------------------------        
        
        """
        self.sort()
        if len(self) == 0:
            return GenomicRegionSet('None region')
        elif len(self) == 1:
            return self
        else:
            z = GenomicRegionSet('Clustered region set')
            previous = self.sequences[0]
            for s in self.sequences[1:]:
                s_ext = copy.deepcopy(s)
                s_ext.extend(max_distance,max_distance)
                if s_ext.overlap(previous):
                    previous.initial = min(previous.initial, s.initial)
                    previous.final = max(previous.final, s.final)
                else:
                    z.add(previous)
                    previous = s
            z.add(previous)
            return z
        
    def flank(self,size):
        """Return two flanking intervals with given size from both ends of each region.
        
        Keyword arguments:
        size -- the length of flanking intervals (default = SAME length as the region)
        
        Return:
        z -- a GenomicRegionSet including all flanking intervals
        
        Graphical explanation:
        self        -----           --            ---
        Result -----     -----    --  --       ---   ---
        
        """
        if self.__len__() == 0:
            return GenomicRegionSet("Empty")
        else:
            z = GenomicRegionSet("Flanking intervals")
            for s in self:
                s1 = GenomicRegion(name='upstream',chrom=s.chrom,
                                   initial=max(0, s.initial - size),
                                   final=s.initial)
                s2 = GenomicRegion(name='downstream',chrom=s.chrom,
                                   initial=max(0, s.final),
                                   final=s.final + size)  # Adding the limit of chromosome length
                z.add(s1)
                z.add(s2)
            return z
    
    def jaccard(self,query):
        """Return jaccard index, a value of similarity of these two GenomicRegionSet
        
        Keyword arguments:
        query -- the GenomicRegionSet which to compare with
        
        Return:
        similarity -- (Total length of overlapping regions)/(Total length of original regions)
        
        Graphical explanation:
        self           --8--      ---10---    -4-
        query     ---10---             ---10---
        intersect      -5-             -4-    2
        similarity:   ( 5 + 4 + 2)/[(8 + 10 + 4) + (10 +10) - (5 + 4 + 2)]
                      = 11/31
        
        """
        intersects = self.intersect(query)
        #print(intersects.total_coverage(),self.total_coverage(), query.total_coverage(),sep="\t")
        similarity = intersects.total_coverage()/(self.total_coverage() + query.total_coverage() - intersects.total_coverage())
        return similarity
    
    def within_overlap(self):
        """ Check whether there is overlapping within or not. """
        refer_posi = GenomicRegion(name="reference",chrom="chr1",initial=0,final=0)
        self.sort()
        for s in self:
            if s.overlap(refer_posi):
                return True
            refer_posi = s
        return False

    def total_coverage(self):
        """ Return the sum of all lengths of regions. """
        length = 0
        for s in self:
            length = length + s.__len__()
        return length
    
    def get_genome_data(self,organism, chrom_X=False, chrom_M=False):
        """ Add genome data from database into the GenomicRegionSet. """
        genome = Genomedata(organism)
        chromosome_file = open(genome.get_chromosome_sizes(),'r')
        for line in chromosome_file:
            if chrom_X and "chrX" in line and "random" not in line:
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            elif chrom_M and "chrM" in line and "random" not in line:
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            elif "random" not in line and "_" not in line and "chrX" not in line and "chrM" not in line: 
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
            else:
                continue
        chromosome_file.close()
        
    def random_regions(self, organism, total_size=None, multiply_factor=0, 
                       overlap_result=False, overlap_input=False, chrom_X=False, chrom_M=False, 
                       filter_path=None):
        """Return a GenomicRegionSet which contains the random regions generated by given entries and given number on the given organism.
        
        Keyword arguments:
        organism -- which organism's genome to use. (hg19, mm9)
        total_size -- Given the number of result random regions.
        multiply_factor -- this factor multiplies to the number of entries is the number of exporting random regions.
                           ** total_size has higher priority than multiply_factor. **
        overlap_result -- the results whether overlap with each other or not. (True/False)
        overlap_input -- the results whether overlap with input entries or not. (True/False)
        chrom_M -- the result should cover mitochondria chromosome or not. (True/False)
        
        Return:
        z -- a GenomicRegionSet which contains the random regions
        
        """

        # Fetching the chromosome length from data
        chrom_map = GenomicRegionSet("chrom_map")
        chrom_map.get_genome_data(organism, chrom_X=False, chrom_M=False)
        chrom_list = chrom_map.get_chrom()
        if filter_path:
            filter_map = GenomicRegionSet('filter')
            filter_map.read_bed(filter_path)
            chrom_map = chrom_map.subtract(filter_map)
            print("chrom_map length: ", len(chrom_map.sequences))
        # Defining input_map, result_list (containing lengths of all result regions)
        result_list = []
        input_map = self
        input_num = self.__len__()
        input_list = [x.__len__() for x in self.sequences] # Contain a list of length of inputs
            
        if type(total_size) == int:
            for i in range(total_size):
                result_list.append(input_list[i % input_num])
        elif multiply_factor > 0:
            for i in range(int(multiply_factor * input_num)):
                result_list.append(input_list[i % input_num])
        
        # Dealing with positions of random regions
        z = GenomicRegionSet(name="random regions")
        if overlap_input == False and overlap_result == False:
            result_map = chrom_map.subtract(input_map)
            for length in result_list:
                cont_loop = True
                while cont_loop:
                    sample = random.choice(result_map.any_chrom(random.choice(chrom_list)))
                    try:
                        random_posi = random.randint(sample.initial, sample.final - length)
                        cont_loop = False
                    except:
                        pass
                new_region = GenomicRegion(chrom=sample.chrom,
                                           initial=random_posi,
                                           final=random_posi + length)
                z.add(new_region)
                x = GenomicRegionSet("temp")
                x.add(new_region)
                result_map = result_map.subtract(x)
            
        elif overlap_input == False and overlap_result == True:
            result_map = chrom_map.subtract(input_map)
            for length in result_list:
                cont_loop = True
                while cont_loop:
                    sample = random.choice(result_map.any_chrom(random.choice(chrom_list)))
                    try:
                        random_posi = random.randint(sample.initial, sample.final - length)
                        cont_loop = False
                    except:
                        pass
                z.add(GenomicRegion(chrom=sample.chrom,
                                    initial=random_posi,
                                    final=random_posi + length))
        
        elif overlap_input == True and overlap_result == False:
            result_map = chrom_map
            for length in result_list:
                cont_loop = True
                while cont_loop:
                    sample = random.choice(result_map.any_chrom(random.choice(chrom_list)))
                    try:
                        random_posi = random.randint(sample.initial, sample.final - length)
                        cont_loop = False
                    except:
                        pass
                new_region = GenomicRegion(chrom=sample.chrom,
                                           initial=random_posi,
                                           final=random_posi + length)
                z.add(new_region)
                x = GenomicRegionSet("temp")
                x.add(new_region)
                result_map = result_map.subtract(x)
                
        elif overlap_result == True and overlap_result == True:
            result_map = chrom_map
            for length in result_list:
                cont_loop = True
                while cont_loop:
                    sample = random.choice(result_map.any_chrom(random.choice(chrom_list)))
                    try:
                        random_posi = random.randint(sample.initial, sample.final - length)
                        cont_loop = False
                    except:
                        pass
                z.add(GenomicRegion(chrom=sample.chrom,
                                    initial=random_posi,
                                    final=random_posi + length))
        return z
    
    def projection_test(self, query, organism):
        """" Return the p value of binomial test. """
        chrom_map = GenomicRegionSet("Genome")
        chrom_map.get_genome_data(organism=organism)
        #print("coverage of reference: ",self.total_coverage(),"\tcoverage of genome: ",chrom_map.total_coverage())
        possibility = self.total_coverage() / chrom_map.total_coverage() # The average likelihood
        #print("The average likelihood: ", possibility)
        
        intersect_regions = self.intersect(query,mode=OverlapType.OVERLAP)
        n = query.__len__()
        k = intersect_regions.__len__()
        #print("intersections: ",k,"\tnumber of query",n,"\tgenetic coverage: ",possibility)
        p = float(stats.binom_test(k, n, possibility))
        return p

    def any_chrom(self,chrom):
        """ Return a list of regions which belongs to given chromosome. (used in random_regions) """
        return [s for s in self if s.chrom == chrom] 
    
    def relocate_regions(self,center='midpoint', left_length=2000, right_length=2000):
        """ Return a new GenomicRegionSet which relocates the regions by given center and extend length.
        
        Parameters:
        center:    
            midpoint  -- locate the new region's center as original region's midpoint
            leftend      -- locate the new region's center as original region's 5' end (if no orientation information, default is left end)
            rightend      -- locate the new region's center as original region's 3' end (if no orientation information, default is right end)
            bothends  -- locate the new region's center as original region's both ends
        left_length:
        right_length:
            The length to extend from the center
        """
        new_regions = GenomicRegionSet('new'+self.name)
        for r in self.sequences:
            # Define the position
            if center == 'midpoint':
                mp = int(0.5*(r.initial + r.final))
                nr = GenomicRegion(chrom=r.chrom, initial=mp, final=mp)
            elif center == 'leftend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial)
            elif center == 'rightend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final)
            elif center == 'bothends':
                nrL = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial) # newbed on left end
                nrR = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final) # newbed on right end
            try:new_regions.add(nr)
            except: 
                new_regions.add(nrL)
                new_regions.add(nrR)
                
        # Extend the region
        new_regions.extend(left_length, right_length)
        return new_regions
