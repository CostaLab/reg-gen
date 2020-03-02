"""
GeneSet
===================
GeneSet describes genes and their expression.

"""

###############################################################################
# Libraries
###############################################################################

# Python


# Internal
from .Util import GenomeData

# External

###############################################################################
# Class
###############################################################################


class GeneSet:
    """*Keyword arguments:*

        - name -- Name of the GeneSet
    """

    def __init__(self, name):
        self.name = name
        self.genes = []  # list of genes to consider
        self.values = {}  # keys: gene, value: expression data as a list
        self.cond = []
        self.symbol_dict = {}

    def __len__(self):
        """Return the number of genes."""
        return len(self.genes)

    def __iter__(self):
        """Iterate this GeneSet."""
        return iter(self.genes)

    def add(self, gene_name, value=None):
        self.genes.append(gene_name)
        if value:
            self.values[gene_name] = value

    def read(self, gene_list_file, score=False):
        """Read genes from the file.

        *Keyword arguments:*

            - geneListFile -- Path to the file which contains a list of genes.
        """
        with open(gene_list_file) as f:
            for line in f:
                line = line.strip()
                if line:
                    l = line.split()
                    if l[0] != "":
                        if "." in l[0]:
                            # gene_name = l[0].upper()
                            gene_name = l[0].partition(".")[0].upper()
                            # self.genes.append(l[0].partition(".")[0].upper())
                        elif "," in l[0]:
                            gene_name = l[0].partition(",")[0].upper()
                        else:
                            gene_name = l[0].upper()
                        self.genes.append(gene_name)
                        if score:
                            self.values[gene_name] = l[1]

    def read_expression(self, geneListFile, header=False, valuestr=False):
        """Read gene expression data.

        *Keyword arguments:*
        
            - geneListFile -- Path to the file which contains genes and expression value.
            - header -- Read first line as header.
            - valuestr -- Keep the value as a string, otherwise convert to float number.
        """
        with open(geneListFile) as f:
            if header:
                l = f.readline()
                l = l.strip("\n")
                l = l.split("\t")
                self.cond = l[1:len(l)]
            else:
                l = f.readline()
                l = l.strip("\n")
                l = l.split("\t")
                self.cond = [str(e) for e in range(len(l) - 1)]
            for line in f.readlines():
                line = line.strip("\n")
                l = line.split("\t")
                if l[0] != "":
                    try:
                        # print(l[0])
                        if "." in l[0] and l[0][0:2] == "EN":
                            na = l[0].partition(".")[0].upper()
                        else:
                            na = l[0].upper()
                        self.genes.append(na)
                        # self.values[l[0].upper()] = [float(v) for v in l[1:len(l)]]
                        if not valuestr:
                            self.values[na] = float(l[1])
                        else:
                            self.values[na] = l[1]
                    except:
                        print("*** error in loading gene: " + line)

    def get_all_genes(self, organism):
        """Get all gene names for a given organism.
        
        *Keyword arguments:*
        
            - organism -- Define the organism.
        """
        genome = GenomeData(organism=organism)
        self.genes = []
        f = open(genome.get_gene_regions())
        for l in f.readlines():
            l = l.strip("\n")
            l = l.split("\t")
            self.genes.append(l[3].upper())
        f.close()

        self.genes = list(set(self.genes))

    def subtract(self, gene_set):
        """Subtract another GeneSet.
        
        *Keyword arguments:*

            - gene_set -- Another GeneSet for subtracting with.
        """
        self.genes = [gene for gene in self.genes if gene not in gene_set.genes]

    def check(self, a_gene):
        """Check a gene is in the list or not

        *Keyword arguments:*

            - a_gene -- A gene symbol.
        """
        # g = a_gene.partition(".")[0].upper()
        g = a_gene.upper()

        return g in self.genes

    def save(self, filename):
        """Save gene list into the given filename."""
        if self.values:
            with open(filename, "w") as f:
                for g in self.genes:
                    print("\t".join([g, str(self.values[g])]), file=f)
        else:
            with open(filename, "w") as f:
                for g in self.genes:
                    print(g, file=f)

    # def load_alias(self, organism):
    #     """Load alias file for gene names"""

    #     genome_data = GenomeData(organism)
    #     # Opening alias file
    #     alias_file = open(genome_data.get_gene_alias(),"r")
    #     # Iterating over alias file entries
    #     for line in alias_file:
    #         ll = line.strip().split("\t")
    #         ensembl_id = ll[0]
    #         official_name = ll[1]
    #         alias_vec = ll[2].split("&")
    #         self.symbol_dict[ensembl_id] = official_name
    #         for e in alias_vec:
    #             try: self.alias_dict[e].append(ensembl_id)
    #             except Exception: self.alias_dict[e] = [ensembl_id]
    #     # Termination
    #     alias_file.close()

    # def ensembl2symbol(self, organism):
    #     """Convert ensembl ID to gene symbol"""
    #     if not self.symbol_dict:
    #         self.load_alias()
    #     curr_list = self.genes
    #     mapped_list = []
    #     unmapped_list = []
    #     for e in curr_list:
    #         try: mapped_list.append(self.symbol_dict[e])
    #         except Exception: unmapped_list.append(self.symbol_dict[e])

    #     self.genes = mapped_list
