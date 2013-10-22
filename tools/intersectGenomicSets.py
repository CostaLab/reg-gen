"""
A intersectGenomicSets perform intersection statistics for a set of bed files.

Authors: Ivan G. Costa, Manuel Allhoff

It recieves as input a experimental matrix with a list of bed files and outputs simple overlap statistics.

"""


import sys
from rgt.ExperimentalMatrix import *
from rgt.GenomicRegionSet import *

def bedOverllap(beds):
  names=[b.name for b in beds]
  res=[]
  for b in beds:
    resLine=[]
    for b2 in beds:
        inter=b.intersect(b2)
        res.append(len(inter))
    res.append(resLine)
  




designFile = sys.argv[1]
exps=ExperimentalMatrix()
exps.read(designFile)
beds = exps.get_regionsets()

print "\t"+("\t".join())

