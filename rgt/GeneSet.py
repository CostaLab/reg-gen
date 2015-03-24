from Util import GenomeData
"""
Describe genes and their expression.

Authors: Ivan G. Costa, Manuel Allhoff

"""

class GeneSet:
    
    def __init__(self,name):
        self.name = name
        self.genes = [] #list of genes to consider
        self.values = {} #keys: gene, value: expression data as a list
        self.cond = []

    def __len__(self):
        return len(self.genes)
    
    def __iter__(self):
        return iter(self.genes)

    def read(self, geneListFile, extra_column=False):
        """Read genes"""
        with open(geneListFile) as f:
            if extra_column:
                self.extra_column = []
                for line in f:            
                    l = line.strip()
                    self.genes.append(l[0])
                    self.extra_column.append(l[1])
            else:
                lines = (line.rstrip() for line in f)
                self.genes = list(set((line for line in lines if line))) # Non-blank lines in a list
            
    def read_expression(self, geneListFile, header = True):
        """Read gene expression data"""
        with open(geneListFile) as f:
            if header:
                l = f.readline()
                l = l.strip("\n")
                l = l.split("\t")
                self.cond = l[1:len(l)]
            for l in f.readlines():
                l = l.strip("\n")
                l = l.split("\t")
                self.genes.append(l[0].upper())
                self.values[l[0].upper()] = [float(v) for v in l[1:len(l)]]

    def get_all_genes(self, organism):
        """Get all gene names for a given organism"""
        genome = GenomeData(organism=organism)
        self.genes = []
        f = open(genome.get_association_file())
        for l in f.readlines():
            l = l.strip("\n")
            l = l.split("\t")
            self.genes.append(l[3].upper())
        f.close()

        self.genes = list(set(self.genes))

    def subtract(self, gene_set):
        """Subtract another GeneSet"""
        self.genes = [gene for gene in self.genes if gene not in gene_set.genes]

