from rgt.GenomicRegion import *

class SetGenomicRegions:
  
  def __init__(self,name):
    self.name=name
    self.sequences=[]
    self.sorted=False
    
  def add(self,sequence):
    self.sequences.append(sequence)
    self.sorted=False
    
  def __len__(self):
    return len(self.sequences)
    
  def  __iter__(self):
    return iter(self.sequences)

  def extend(self,left,right):
    for s in self.sequences:
      s.extend(left,right)

  def sort(self):
    self.sequences.sort(cmp=GenomicRegion.__cmp__)
    self.sorted=True

  def readBed(self,file):
    self.fileName=file
    chrs={}
    ct=0
    for l in open(file):
      #try:
        name=None
        orientation=None
        data=None
        l = l.strip("\n")
        l = l.split("\t")
        size=len(l)
        chr=l[0]
        v1= int(l[1])
        v2 = int(l[2])
        if v1>v2:
            aux=v1
            v1=v2
            v2=aux
        if size > 3:
          name = l[3]
        if size > 4:
          orientation= l[4]
        if size > 5:
          data="\t".join(l[5:])
        self.add(GenomicRegion(chr,v1,v2,name,orientation,data))
      #except:
      #  print file,l
    self.sort()

  
  def intersect(self,y):
    chr=None
    z=SetGenomicRegions(self.name+'_'+y.name)
    x=self.__iter__()
    y=y.__iter__()
    iter=True
    s=x.next()
    ss=y.next()
    while iter:
      if s.overlapp(ss):
        z.add(s)
        try:
          s=x.next()
        except:
          iter=False
      elif cmp(s,ss)>=0:
        try:
          ss=y.next()
        except:
          iter=False
      else:
        try:
          s=x.next()
        except:
          iter=False
    z.sort()
    return z
        

  def writeBed(self,fileName):
    file=open(fileName,'w')
    for s in self:
       file.write(str(s)+'\n')
    file.close()

  def filterByGeneAssociation(self,fileName,genes,geneAnnotation,genomeSize):
    """code based on eduardos functions. This should be  integrated in the core framework soon"""
    "XXX"
    from rgt.motifanalysis.util import bedFunctions,sort
    from rgt.motifanalysis.enrichment.geneAssociation import *
    self.fileName=fileName
    #de_genes=list(set([l.strip("\n") for l in open(geneListFile)]))
    
    coordDict = bedFunctions.createBedDictFromSingleFile(fileName, features=[1,2,3,4,5]) 
    coordDict = sort.sortBedDictionary(coordDict, field=0)
    [dictBed,allBed]=geneAssociationByPromoter(coordDict,genes.genes,geneAnnotation,genomeSize)
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
  



   
   
      
      




    

       




    
