import sys
import os.path
from rgt.SetGenomicRegions import *
from fisher import pvalue



#geneList = sys.argv[2]
anotationPath = sys.argv[2]
genomeFile=anotationPath+"chrom.sizes"
geneFile=anotationPath+"association_file.bed"

print os.path.realpath('.')
#get experiments
beds=[]
geneLists=[]

bedGenes = SetGenomicRegions(geneFile)
bedGenes.readBed(geneFile)
allgenes=[]
for r in bedGenes:
 allgenes.append(r.name)
allgenes=list(set(allgenes))

for l in open(sys.argv[1]):
      l=l.strip('\n')
      l=l.split('\t')
      name1=l[0]
      name2=l[1]
      geneList=l[3]
      bedFile=l[2]


      bed = SetGenomicRegions(bedFile)
      [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filterByGeneAssociation(bedFile,geneList,geneFile,genomeFile)
      #print degenes,de_peak_genes, mappedGenes, totalPeaks
      prop_de=de_peak_genes/float(degenes)
      prop_back=mappedGenes/float(len(allgenes))
      a=de_peak_genes
      b=degenes-de_peak_genes
      c=mappedGenes-de_peak_genes
      d=len(allgenes)-b-c-a
      p= pvalue(a,b,c,len(allgenes)-b-c-a)
      print name1,name2,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail

#print "No Peaks", len(bed)
#print "No Genes", len(bed.genes)
#print bed.genes
#bed.writeBed(bedFile+"_gene.bed")

