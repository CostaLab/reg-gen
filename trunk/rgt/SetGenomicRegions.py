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

class SetGenomicRegions:
  
    def __init__(self, name):
        """Initialize a set of GenmocRegions"""
        self.name = name
        self.sequences = []
        self.sorted = False
    
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

    def readBed(self, filename):
        """Read BED file (see format: http://genome.ucsc.edu/FAQ/FAQformat.html#format1)
        and add every row as a GenomicRegion."""
        with open(filename) as line:
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
                data = "\t".join( [line[4], line[6:]] )
            self.add( GenomicRegion(chrom, start, end, name, orientation, data) )
        self.sort()
  
    def intersect(self,y):
        """Return new SetGenomicRegions as the intersection with y"""
        z = SetGenomicRegions(self.name + '&' + y.name)
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
        

    def writeBed(self,filename):
        """Write GenomicRegions to BED file"""
        with open(filename, 'w') as f:
            for s in self:
                print(s, file=f)

    def filterByGeneAssociation(self,fileName,geneListFile,geneAnnotation,genomeSize):
        """code based on eduardos functions. This should be  integrated in the core framework soon
        TODO: Eduardo should check this!"""
        from rgt.motifanalysis.util import bedFunctions,sort
        from rgt.motifanalysis.enrichment.geneAssociation import *
        self.fileName=fileName
        de_genes=list(set([l.strip("\n") for l in open(geneListFile)]))
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
  
