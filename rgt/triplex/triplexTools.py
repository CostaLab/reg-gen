# Python Libraries
from __future__ import print_function
#from collections import *
import os
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from SequenceSet import Sequence, SequenceSet
from Util import SequenceType
#import multiprocessing

####################################################################################
####################################################################################

class TriplexSearch:
    """Contains functions for potential triplex forming sites on DNA or RNA.

    Methods:
        1. For RNA as input: define the TFO (triplex-forming oligonucleotides) on RNA
        2. For DNA as input: define the TTS (triplex target sites) on DNA
        3. For RNA and DNA as input: define the possible triplex region between RNA and DNA

    Some parameters:

        motifs:
            R - the purine motif that permit guanines (G) and adenines (A)
            Y - the pyrimidine motif that permit cytosines (C) and thymines (T)
            M - the mixed motif, purine-pyrimdine, that permit guanines (G) and thymines (T)
            P - parallel binding, i.e. motifs facilitating Hoogsten bonds; 
                this covers the pyrimidine motif and the purine-pyrimidine motif in parallel configuration
            A - anti-parallel binding, i.e. motifs facilitating reverse Hoogsten bonds; 
                this covers the purine motif and the purine-pyrimidine motif in anti-parallel configuration
    """
    
    def search_bindingsites(self, sequence_set, seq_type, motif, min_len, max_len):
        """Find binding sites on RNA or DNA from given SequenceSet
        
        Parameters:
            sequence_set:  A SequenceSet object
            seq_type:      RNA or DNA
            motif:         R, Y, M, P, or A
            min_len:       Define the minimum length of RBS (default is 10 bp)
            max_len:       Define the maximum length of RBS (default is infinite)
        """
        
        all_binding_sites = BindingSiteSet(name=sequence_set.name)
        for s in sequence_set:
            if seq_type == SequenceType.RNA:
                bs = self.find_rbs(s, motif, min_len, max_len)
            elif seq_type == SequenceType.DNA:
                bs = self.find_dbs(s, motif, min_len, max_len)
            all_binding_sites.concatenate(bs)

        print(len(all_binding_sites))
        return all_binding_sites

    def find_rbs(self, a_sequence, motif="RYMPA", min_len=10, max_len=None):
        """Return a RNABindingSet with the potential RNA binding sites

        Parameters:
            a_sequence:  A Sequence object
            min_len:     Define the minimum length of RBS (default is 10 bp)
            max_len:     Define the maximum length of RBS (default is infinite)
            motif:       R: G, A
                         Y: C, T
                         M: G, T
                         P: C, T, G
                         A: A, G, T
        """
        rbss = BindingSiteSet(name=a_sequence.name)

        for m in motif:
            # Motif choice for tfo on RNA
            if m == "R":
                targets = ["G", "A"]   
            elif m == "Y":
                targets = ["C", "T"]
            elif m == "M":
                targets = ["G", "T"] 
            elif m == "P":
                targets = ["C", "T", "G"] 
            elif m == "A":
                targets = ["A", "G", "T"] 

            # Parsing the sequence

            count = 0
            sample = ""
            con_rbs = False

            for i, a in enumerate(a_sequence.seq):

                if a in targets:
                    if count < min_len:
                        sample += a
                        count += 1
                    elif min_len <= count:
                        sample += a
                        count += 1
                        rbs = BindingSite(chrm=a_sequence.name, initial=i-count+1, final=i+1, 
                                          name=a_sequence.name, score=count,  errors_bp=0, motif=m)
                        con_rbs = True

                    elif max_len and count > max_len:
                        rbss.add(rbs)
                        con_rbs = False
                else:  
                    sample = ""
                    count = 0
                    if con_rbs: 
                        rbss.add(rbs)
                        con_rbs = False
        return rbss


    def find_dbs(self,a_sequence, motif="RYMPA", min_len=10, max_len=None, organism="hg19"):
        """Calculate the potential DBS (DNA binding sites) on DNA in the given sequence

        Parameters:
            a_sequence:  A Sequence object
            min_len:     Define the minimum length of RBS (default is 10 bp)
            max_len:     Define the maximum length of RBS (default is infinite)
            motif:
                         R: G, A
                         Y: C, T
                         M: G, T
                         P: C, T, G
                         A: A, G, T

            organism:   Define the organism (hg19 as default)
        """
        
        all_dbs = RNABindingSet(name=a_sequence.name)

        
        targets = ["A", "G"] 

        # Parsing the sequence
        count = 0
        sample = ""
        con_dbs = False

        for i, a in enumerate(a_sequence.seq):
            if a in targets:
                if count < min_len:
                    sample += a
                    count += 1
                elif min_len <= count:
                    sample += a
                    count += 1
                    dbs = RNABinding(name=a_sequence.name, initial=i-count+1, final=i+1, 
                                     score=count, errors=0, motif=m, 
                                     orientation=None, seq=sample)
                    con_dbs = True

                elif max_len and count > max_len:
                    all_dbs.add(rbs)
                    con_dbs = False
            else:  
                sample = ""
                count = 0
                if con_rbs: 
                    all_dbs.add(dbs)
                    con_dbs = False
        return all_dbs
        

    #def find_triplex(self):

####################################################################################
####################################################################################

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
        



####################################################################################
####################################################################################

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


if __name__ == '__main__':
    a = SequenceSet(name="test", seq_type=SequenceType.RNA)
    a.read_fasta("/projects/lncRNA/data/hotair.fasta")
    print(len(a))


    triplex = TriplexSearch()
    b = triplex.find_rbs(a[0])
    print(b)
    print(len(b))
    for s in b:
        print(s.seq + "\t\t" + s.motif)
    b.write_rbs("test")
