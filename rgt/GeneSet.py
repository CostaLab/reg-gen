class SetGenes:
    
    def __init__(self,name):
        self.name=name
        self.genes=[]

    def read(self,file):
        self.genes=list(set([l.strip("\n") for l in open(geneListFile)]))