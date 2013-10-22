from __future__ import print_function
from rgt.GenomicRegion import *

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
        """Initialize a set of GenmocRegions"""
        self.name = name
        self.sequences = []
        self.sorted = False
        self.fileName=""
    
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
            self.sort()
  
    def intersect(self,y):
        """Return new GenomicRegionSet as the intersection with y"""
        z = GenomicRegionSet(self.name + '&' + y.name)
        x = self.__iter__()
        y = y.__iter__()
        cont_loop = True
        s = x.next()
        ss = y.next()
        while cont_loop:
            if s.overlapp(ss):
                z.add(s)
                try:
                    x.next()
                except:
                    cont_loop = False
            elif s >= ss:
                try:
                    ss = y.next()
                except:
                    cont_loop = False
            else:
                try:
                    s = x.next()
                except:
                    cont_loop = False
        z.sort()
        return z
        

    def write_bed(self,filename):
        """Write GenomicRegions to BED file"""
        with open(filename, 'w') as f:
            for s in self:
                print(s, file=f)

    def filter_by_gene_association(self,fileName,geneSet,geneAnnotation,genomeSize):
        """code based on eduardos functions. This should be  integrated in the core framework soon
        TODO: Eduardo should check this!"""
        from rgt.motifanalysis.util import bedFunctions,sort
        from rgt.motifanalysis.enrichment.geneAssociation import *
        self.fileName=fileName
        de_genes=geneSet.genes
        coordDict = bedFunctions.createBedDictFromSingleFile(fileName, features=[1,2,3,4,5]) 
        coordDict = sort.sortBedDictionary(coordDict, field=0)
        [dictBed,allBed]=geneAssociationByPromoter(coordDict,de_genes,geneAnnotation,genomeSize)  
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
  
