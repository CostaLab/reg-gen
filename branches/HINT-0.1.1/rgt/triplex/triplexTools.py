import os
from rgt.GeneSet import *
from rgt.GenomicRegionSet import *
from rgt.RNADNABindingSet import *

class FischerTest:
    def __init__(self, gene_list_file, organism, promoterLength):
        """Initiation"""
        self.organism = organism

        # DE gene regions
        self.df_gene = GeneSet("df genes")
        self.df_gene.read(gene_list_file)
        df_regions = GenomicRegionSet("df genes")
        self.df_regions = df_regions.get_from_genes(gene_list=self.df_gene, organism=organism, 
                                                    promoterLength=promoterLength)

        # All gene regions
        self.all_gene = GeneSet("all genes")
        self.all_gene.get_all_genes(organism=organism)
        all_regions = GenomicRegionSet("all genes")
        self.all_regions = all_regions.get_from_genes(gene_list=self.all_gene, organism=organism, 
                                                      promoterLength=promoterLength)

    def search_triplex(self, rna, temp):
        """Perform triplexator on DE genes and all genes to find triplexes"""
        # DE gene regions
        self.df_regions.write_bed(os.path.join(temp,"de_gene.bed"))
        bed = os.path.join(temp, "de_gene.bed")
        fasta = os.path.join(temp, "de_gene.fasta")
        txp = os.path.join(temp, "de.txp")
        os.system("bedtools getfasta -fi /data/genome/hg19genome.fa -bed "+bed+" -fo "+fasta)
        os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -ss "+rna+" -ds "+fasta+" > "+txp)
        
        # All gene regions
        self.all_regions.write_bed(os.path.join(temp,"all_gene.bed"))
        bed = os.path.join(temp, "all_gene.bed")
        fasta = os.path.join(temp, "all_gene.fasta")
        txp = os.path.join(temp, "all.txp")
        os.system("bedtools getfasta -fi /data/genome/hg19genome.fa -bed "+bed+" -fo "+fasta)
        os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -ss "+rna+" -ds "+fasta+" > "+txp)

    def load_txp(self,temp):
        """Loading the txp files from temp directory"""
        de_binding = RNADNABindingSet(organism=self.organism, filename=os.path.join(temp,"de.txp"))
        
        all_binding = RNADNABindingSet(organism=self.organism, filename=os.path.join(temp,"all.txp"))
        





class RandomTest:
    def __init__(self, txp_path, organism):
    	self.txp = RNADNABindingSet(organism=organism, filename=txp_path)
    	
    def target_count(self):
    	"""Count the number of TFFs for each TFO

    	The count number is stored in a dictionary as self.merged_TFO
    	"""
    	self.txp.merged_TFO()
    	self.merged_TFO = OderedDict()
    	for tfo, tffs in iteritems(self.txp.merged):
    		self.merged_TFO[tfo] = len(tffs)
        
    #def randomization(repeat):

