import sys
from SetGenomicRegions import *

tfs={}
for l in open(sys.argv[1]):
      l=l.strip('\n')
      l=l.split('\t')
      tfs[l[0]]=l[1]

histones={}
for l in open(sys.argv[2]):
    l=l.strip('\n')
    l=l.split('\t')
    histones[l[0]]=l[1]

bedtfs={}
for n in tfs.keys():
  bedtfs[n] = SetGenomicRegions(n)
  bedtfs[n].readBed(tfs[n])
nametfs = tfs.keys()

bedhistones={}
for n in histones.keys():
  bedhistones[n] = SetGenomicRegions(n)
  bedhistones[n].readBed(histones[n])
namehistones =  bedhistones.keys()
                
print "tf"+"\t"+"\t".join(namehistones)+"\t"+"cts"+"\t"+"\t".join(namehistones)
ct=0
for tf in bedtfs.keys():
    btf=bedtfs[tf]
    cttf=len(btf)
    res=[]
    for hist in bedhistones.keys():
      bhis=bedhistones[hist]
      chist=len(bhis)
      inter=btf.intersect(bhis)
      res.append(len(inter))
    print tf+"\t"+"\t".join([str(o) for o in res])+"\t"+str(cttf)+"\t"+"\t".join([str(float(o)/cttf) for o in res])
    ct+=1
#print bedtfs
