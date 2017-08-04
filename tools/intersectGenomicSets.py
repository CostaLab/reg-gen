"""
A intersectGenomicSets perform intersection statistics for a set of bed files.

Authors: Ivan G. Costa, Manuel Allhoff

It recieves as input a experimental matrix with a list of bed files and outputs simple overlap statistics.

"""


import sys
from rgt.ExperimentalMatrix import *
from rgt.GenomicRegionSet import *
import numpy

def bedOverllap(beds,beds2,outPath):
    namesCol=[b.name for b in beds2]
    namesLine=[b.name for b in beds]
    res=[]
    totalCol=[len(b) for b in beds2]
    totalLine=[len(b) for b in beds]
    #print len(beds), len(beds2)
    for b in beds:
        resLine=[]
        for b2 in beds2:
            inter=b.intersect(b2,OverlapType.ORIGINAL)
            if len(outPath) > 0:
              inter.write(outPath+"/"+inter.name+".bed")
            resLine.append(len(inter))
        res.append(resLine)
    res=numpy.array(res,float)
    #if len(totalCol)>1:
    resprop1=res/numpy.array(totalCol)
    #else:
    #resprop1=numpy.divide(res,numpy.array([totalLine]).transpose())
    resprop2=res.transpose()/numpy.array(totalLine)
    resprop2=resprop2.transpose()
    return namesCol, namesLine,res,resprop1,resprop2


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
beds2=beds
outputDir = sys.argv[2]
if len(sys.argv) > 3:
  experimentalFile2 = sys.argv[3]
  exps2=ExperimentalMatrix()
  exps2.read(experimentalFile2)
  beds2 = exps2.get_regionsets()
  out=outputDir


[namesCol,namesLine,res,resprop1,resprop2]=bedOverllap(beds,beds2,outPath=out)
printTable(namesCol,namesLine,res,outputDir+"/count.table")
printTable(namesCol,namesLine,resprop1,outputDir+"/propline.table")
printTable(namesCol,namesLine,resprop2,outputDir+"/propcollumn.table")



