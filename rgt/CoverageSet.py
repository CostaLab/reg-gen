from rgt.SetGenomicRegions import *
import pysam
import numpy
import numpy.ma
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class CoverageSet:
  
  def __init__(self,name,genomicRegions):
    self.name=name
    self.genomicRegions=genomicRegions
    self.values=[]
    self.sizeBam=0
    self.step=1
    self.window=1

  def coverageFromBam(self,bamFile,readSize=200,step=50,window=50):
    self.step=step
    bam = pysam.Samfile(bamFile, "rb" )
    #size of bamfile
    self.sizeBam = reduce(lambda x, y: x + y, [ eval('+'.join(l.rstrip('\n').split('\t')[2:]) ) for l in pysam.idxstats(bamFile)])
      #print self.sizeBam
    for region in self.genomicRegions:
      cov=[0]*(len(region)/step)
      #cov=[0]*len(region)
      for r in bam.fetch(region.chr,region.initial-readSize,region.final+readSize):
      #for j,v in enumerate(range(region.initial-readExt,region.final+ReadExt,step)):
        if r.is_reverse:
          for i in range(max(0,r.pos-readSize-region.initial)/step,min(len(region),r.pos-region.initial)/step):
#         #   print i, len(cov)
            cov[i]=cov[i]+1
          pass
        else:
          for i in range(max(0,r.pos-region.initial)/step,min(len(region),r.pos+readSize-region.initial)/step):
#            print i, len(cov)
            cov[i]=cov[i]+1
          
        #for pileupcolumn in bam.pileup(region.chr, v-window, v+step+window,stepper="all"):  
        #  if (pileupcolumn.pos-region.initial-window) <= len(region):
        #    if (pileupcolumn.pos-region.initial-window) >= 0:           
        #      cov[j]+=pileupcolumn.n

      #for pileupcolumn in bam.pileup(region.chr, region.initial, region.final,stepper="all"):
      #  try:
      #    if (pileupcolumn.pos-region.initial) <= len(region):
      #      if (pileupcolumn.pos-region.initial) >= 0:           
      #        cov[(pileupcolumn.pos-region.initial)/step]+=pileupcolumn.n
      #  except:
      #    #print "Warning: pile up returned outside read",pileupcolumn.pos,region
      #    pass

      self.values.append(cov)

    self.values=numpy.array(self.values,numpy.float)
    self.valuesorig=self.values.copy()
    self.window=window


  def normRPM(self):
    self.values= self.values*(1000000.0/self.sizeBam)

  def normFPKM(self):
    self.values= self.values*(1000000000.0/(self.sizeBam))

  def normFactor(self,factor):
    self.values= self.values*(1000000.0/(self.step*factor))    

  def log(self):
    print self.name, numpy.min(self.values)
    self.values=numpy.log2(self.values+0.01)

  def mean(self):
    #values= numpy.ma.array(self.values, mask=self.valuesorig==0)
    #sum=numpy.sum(values,axis=0)
    #total=numpy.sum(self.valuesorig>0,axis=0)
    #print "sum",sum
    #print "total", total
    #print sum/(total+0000000000000000.1)
    #return sum/(total+0000000000000000.1)
    return numpy.mean(self.values,axis=0)
    
  def diff(self,coverage):
    self.valuesorig=self.valuesorig+coverage.valuesorig
    self.values=self.values-coverage.values

  def abs(self):
    self.values=abs(self.values)

  def plot(self,log=False,name=None,outTxt=True):
    if name==None:
      name=self.name
    mean=self.mean()
    plt.plot(range(0,self.step*len(mean),self.step),mean)
    plt.title("Mean Expression "+name)
    plt.axis([0, len(mean)*self.step, mean.min(), mean.max()])
    plt.savefig(name+".pdf")
    plt.close("all")
    if outTxt:
      f=open(name+".txt","w")
      f.write("Pos\tValue\n")
      for i,m in enumerate(mean):
        f.write(str(i*self.step)+'\t'+str(m)+'\n')
      f.close()
    

  #def statisticsFromDiffProfiles(self,coverage,log=False,norm="RPM"):
  #  compvalues=numpy.array(coverage.values,numpy.float)
  #  values=numpy.array(self.values,numpy.float) 
  #  if norm=="RPM":
  #    values= values*(1000000.0/self.sizeBam)
  #    compvalues= compvalues*(1000000.0/coverage.sizeBam)      
  #  if norm=="FPKM":
  #    values= values*(1000000000.0/self.sizeBam)
  #    compvalues= compvalues*(1000000000.0/coverage.sizeBam)
  #  values=abs((values)-(compvalues))
  #  if log:     
  #    values= numpy.log2(values+1)
  #  return numpy.mean(values,axis=0), numpy.std(values,axis=0)      

