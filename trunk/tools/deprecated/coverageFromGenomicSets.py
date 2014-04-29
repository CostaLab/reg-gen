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

def bedCoverage(bed,reads):
  c=[]
  for r in reads:
      cov=CoverageSet(r,bed)
      cov.coverage_from_genomicset(r)
      #cov.normRPM()
      c.append(cov.coverage)
  return numpy.transpose(c)   


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
reads = exps.get_readsfiles()
readsnames = exps.get_readsnames()
outputDir = sys.argv[2]
if len(sys.argv) > 3:
  experimentalFile2 = sys.argv[3]
  exps2=ExperimentalMatrix()
  exps2.read(experimentalFile2)
  reads = exps2.get_readsfiles()
  readsnames = exps2.get_readsnames()
  out=outputDir

for bed in beds:
  bednames=[r.chrom+":"+str(r.initial)+"-"+str(r.final) for r in bed]
  c=bedCoverage(bed,reads)
  printTable(readsnames,bednames,c,outputDir+"/"+bed.name+".txt")


