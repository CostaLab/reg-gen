import sys
from rgt.ExperimentalMatrix import *
from rgt.GenomicRegionSet import *

designFile = sys.argv[1]
exps=ExperimentalMatrix()
exps.read(designFile)
beds = exps.get_regionsets()

print "\t"+("\t".join([b.name for b in beds]))
for b in beds:
    res=[b.name]
    for b2 in beds:
        inter=b.intersect(b2)
        res.append(str(len(inter)))
    print "\t".join(res)

