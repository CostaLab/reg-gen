"""
A intersectGenomicSets perform intersection statistics for a set of bed files.

Authors: Ivan G. Costa, Manuel Allhoff

It recieves as input a experimental matrix with a list of bed files and outputs simple overlap statistics.

"""


import sys
from rgt.ExperimentalMatrix import *
from rgt.GenomicRegionSet import *
from rgt.CoverageSet import *
import numpy

def averageExpression(bed,expression,regionsToGenes):
    values=numpy.zeros((len(bed),len(expression.cond)))
    labels=[]
    for i,b in enumerate(bed):
        region=b.toString()
        labels.append(region)
        try:
          genes=regionsToGenes[region]
          #print genes
          for g in genes:
            values[i,]+=expression.values[g.upper()]
            #print expression.values[g.upper()]
          values[i,]=values[i,]/float(len(genes))
        except:
          values[i,]='NaN'
    return values,labels

def printTable(namesCol,namesLines,table,fileName):
    f=open(fileName,"w")
    f.write("\t"+("\t".join(namesCol))+"\n")
    for i,line in enumerate(table):
      f.write(namesLines[i]+"\t"+("\t".join([str(j) for j in line]))+"\n")


out=""
experimentalFile = sys.argv[1]
exps=ExperimentalMatrix()
exps.read(experimentalFile)
beds = exps.get_regionsets()
geneExp= sys.argv[2]
outputFile = sys.argv[3]
anotationPath = sys.argv[4]
genomeFile=anotationPath+"chrom.sizes"
geneFile=anotationPath+"association_file.bed"


genes= GeneSet("Expression")
genes.readExpression(geneExp)

for bed in beds:
  bedNew = GenomicRegionSet("")
  [degenes,de_peak_genes, mappedGenes, totalPeaks,regionsToGenes] = bedNew.filter_by_gene_association(bed.fileName,genes,geneFile,genomeFile,threshDist=5000)
  [ct,labels]=averageExpression(bed,genes,regionsToGenes)
  aux=bed.fileName.split("/")
  fileName=aux[-1]
  printTable(genes.cond,labels,ct,outputFile+"/"+fileName[:-4]+".txt")


