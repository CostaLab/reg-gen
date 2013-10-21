class GeneSet:
    
    def __init__(self,name):
        self.name=name
        self.genes=[]

    def read(self,geneListFile):
        self.genes=list(set([l.strip("\n") for l in open(geneListFile)]))
