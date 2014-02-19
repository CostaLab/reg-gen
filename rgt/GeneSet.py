class GeneSet:
    
    def __init__(self,name):
        self.name=name
        self.genes=[]
        self.values={}
        self.cond=[]

    def read(self,geneListFile):
        self.genes=list(set([l.strip("\n") for l in open(geneListFile)]))

    def readExpression(self,geneListFile,header=True):
        f=open(geneListFile)
        if header:
          l=f.readline()
          l=l.strip("\n")
          l=l.split("\t")
          self.cond=l[1:len(l)]
        for l in f.readlines():
          l=l.strip("\n")
          l=l.split("\t")
          self.genes=l[0].upper()
          self.values[l[0].upper()]= [float(v) for v in l[1:len(l)]]

			
