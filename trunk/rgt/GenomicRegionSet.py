from __future__ import print_function
from rgt.GenomicRegion import *
import random
import copy

"""
Represent list of GenomicRegions.

Authors: Ivan G. Costa, Manuel Allhoff

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
    
    def get_chrom(self):
        """Return all chromosomes"""
        return [ r.chrom for r in self.sequences ]
    
    def add(self, region):
        """Add GenomicRegion"""
        self.sequences.append(region)
        self.sorted = False
    
    def __len__(self):
        return len(self.sequences)
    
    def  __iter__(self):
        return iter(self.sequences)

    def extend(self,left,right):
        """Perform extend step for every element"""
        for s in self.sequences:
            s.extend(left,right)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicRegion.__cmp__)
        self.sorted = True

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
  
 
    def randomRegions(self,size):
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
        from rgt.motifanalysis.util import bedFunctions,sort
        from rgt.motifanalysis.enrichment.geneAssociation import *
        self.fileName=fileName
        de_genes=geneSet.genes
        coordDict = bedFunctions.createBedDictFromSingleFile(fileName, features=[1,2,3,4,5]) 
        coordDict = sort.sortBedDictionary(coordDict, field=0)
        [dictBed,allBed]=geneAssociationByPromoter(coordDict,de_genes,geneAnnotation,genomeSize,promoterLength,threshDist)  
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
                genes=genes+keep
                allgenes=allgenes+[n.strip(".") for n in names]
        #print "Total Peaks", total
        mappedGenes=len(list(set(allgenes)))
        self.sort()
        self.genes=list(set(genes))
        return len(de_genes), len(self.genes), mappedGenes, totalPeaks 

    def intersect(self,y,mode='overlap'):
        """Return the overlapping regions with three different modes.
        
        (mode='overlap') 
        Return new GenomicRegionSet including only the overlapping regions with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
        
            Return:
            z -- a new GenomicRegionSet including only the overlapping regions
        
            Graphical explanation:
            self       ----------              ------
            y                 ----------                    ----
            Result            ---
            
        (mode='original')
        Return the regions of original GenomicRegionSet which have any intersections with y.
        
            Keyword arguments:
            y -- the GenomicRegionSet which to compare with
            
            Return:
            z -- the regions of original GenomicRegionSet which have any intersections with y
            
            Graphical explanation:
            self       ----------              ------
            y              ----------                    ----
            Result     ----------
            
        (mode='comp_incl')
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
        cont_loop = True
        s = con_self.next()
        ss = con_y.next()
        while cont_loop:
            # When the regions overlap
            if s.overlap(ss):
                if mode == 'overlap':
                    sss = GenomicRegion(chrom=s.chrom, 
                                        initial=max(s.initial, ss.initial), 
                                        final=min(s.final, ss.final))
                    z.add(sss)
                elif mode == 'original':
                    z.add(s)
                elif mode == 'comp_incl':
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
        if mode == 'original':
            z.remove_duplicates()
        return z
    
    def closest(self,y):
        """Return a new GenomicRegionSet including the region(s) of y which is closest to any self region. 
           If there are intersection, return False.
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        
        Return:
        z -- the region(s) which is nearest to self
        
        """
        if len(self.sequences) == 0 or len(y.sequences) == 0:
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
        for i in self.sequences:
            if self.sequences.count(i) > 1:
                self.sequences.remove(i)  # remove the first item with value i
    
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
        """Merge the regions within the GenomicRegionSet (no return)"""
        #self.sort()
        if len(self) == 0:
            return GenomicRegionSet('None region')
        elif len(self) == 1:
            return self
        else:
            z = GenomicRegionSet(name=self.name)
            previous = self.sequences[0]
            for s in self.sequences[1:]:
                if previous.overlap(s):
                    previous.initial = min(previous.initial, s.initial)
                    previous.final = max(previous.final, s.final)
                else:
                    z.add(previous)
                    previous = s
            z.add(previous)
            self.sequences = z.sequences
            return
    
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
        if len(self.sequences) == 0:
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
    
    def jaccard(self,y):
        """Return a value of similarity of these two GenomicRegionSet
        
        Keyword arguments:
        y -- the GenomicRegionSet which to compare with
        
        Return:
        similarity -- (Total length of overlapping regions)/(Total length of original regions)
        
        Graphical explanation:
        self           --8--      ---10---    -4-
        y         ---10---             ---10---
        intersect      -5-             -4-    2
        similarity:   ( 5 + 4 + 2)/[(8 + 10 + 4) + (10 +10) - (5 + 4 + 2)]
                      = 11/31
        
        """
        sum_self = sum(s.__len__() for s in self)
        sum_y = sum(s.__len__() for s in y)
        intersects = self.intersect(y)
        sum_inter = sum(s.__len__() for s in intersects)
        similarity = sum_inter/(sum_self + sum_y - sum_inter)
        return similarity
    
    def within_overlap(self):
        refer_posi = GenomicRegion(name="reference",chrom="chr1",initial=0,final=0)
        self.sort()
        for s in self:
            if s.overlap(refer_posi):
                return True
            refer_posi = s
        return False
