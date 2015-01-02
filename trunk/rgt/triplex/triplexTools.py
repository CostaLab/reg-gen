# Python Libraries
from __future__ import print_function
#from collections import *
import os
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from RNABindingSet import RNABinding, RNABindingSet
from SequenceSet import Sequence, SequenceSet
from Util import SequenceType
#import multiprocessing

####################################################################################
####################################################################################

class RNADNABindingSet:
    """A framework to map DNA binding sites to corresponding RNA binding sites """

    def __init__(self, organism, name=None, filename=None):
        """Initialize"""
        
        self.name = name  # RNA name
        self.dict = OrderedDict()    # All RNA/DNA interactions and extra information 
        self.sorted_dna = False
        self.sorted_rna = False
        self.rna = RNABindingSet(name)
        self.dna = GenomicRegionSet("DNA binding sites")
        self.organism = organism
        if filename:
            self.fileName = filename
            self.read_txp(filename)
    
    def read_txp(self, filename):
        """Read txp file to load all interactions. """
        
        with open(filename) as f:
            for line in f:
                if line[0] == "#": continue # skip the comment line
                
                line = line.strip("\n")
                line = line.split()
                
                if len(line) < 12: continue # skip the unimportant lines in txp
                if len(line) == 12:
                    line.insert(8,"_")
                # Format of each line in txp
                #     [0] Sequence-ID   
                #     [1] TFO start   
                #     [2] TFO end 
                #     [3] Duplex-ID   
                #     [4] TTS start   
                #     [5] TTS end 
                #     [6] Score   
                #     [7] Error-rate  
                #     [8] Errors  
                #     [9] Motif   
                #     [10] Strand  
                #     [11] Orientation 
                #     [12] Guanine-rate

                # RNA binding site
                if not self.name: self.name = line[0]
                
                rna_start, rna_end = int(line[1]), int(line[2])
                if rna_start > rna_end: rna_start, rna_end =  rna_end, rna_start
                rna = RNABinding(name=line[0], initial=rna_start, final=rna_end, score=line[6], errors=line[8], 
                                 motif=line[9], orientation=line[11], guanine_rate=line[12])
                
                # DNA binding site
                dna_start = int(line[3].split(":")[1].split("-")[0]) + int(line[4])
                dna_end = int(line[3].split(":")[1].split("-")[0]) + int(line[5])
                dna = GenomicRegion(chrom=line[3].split(":")[0], 
                                   initial=dna_start, 
                                   final=dna_end, 
                                   name=line[3],
                                   orientation=line[10])
                
                # Map RNA binding site to DNA binding site
                self.dict[rna] = dna

    def __len__(self):
        """Return the number of interactions between DNA and RNA """
        return len(self.dict.keys())

    def __iter__(self):
        return self.dict.iteritems()

    def get_rna(self, sort=False):
        """Return RNABindingSet which contains all RNA binding sites"""
        rna_set = RNABindingSet(name=self.name)
        for rna in self.dict.keys():
            rna_set.add(rna)
        if sort: rna_set.sort()
        return rna_set

    def get_dna(self, sort=False):
        """Return GenomicRegionSet which contains all DNA binding sites"""
        dna_set = GenomicRegionSet(name="DNA_binding_sites")
        for dna in self.dict.values():
            dna_set.add(dna)
        if sort: dna_set.sort()
        return dna_set

    def sort_rna(self):
        """Sort the dictionary by RNA"""
        self.dict = OrderedDict(sorted(self.dict.iteritems(), key=lambda x: x[0], cmp=GenomicRegion.__cmp__))
        
    def sort_dna(self):
        """Sort the dictionary by DNA"""
        self.dict = OrderedDict(sorted(self.dict.iteritems(), key=lambda x: x[1], cmp=GenomicRegion.__cmp__))
    
    def merge_rna(self):
        """Merge the RNA binding regions which have overlap to each other and 
           combine their corresponding DNA binding regions.
        
        extend -> Define the extending length in basepair of each RNA binding regions
        perfect_match -> Merge only the exactly same RNA binding regions
        """
        rna_merged = self.get_rna()
        rna_merged.merge()
        new_dict = OrderedDict()
        
        con_rna, con_dna = iter(self)
        rna, dna = con_rna.next(), con_dna.next()
        
        con_rna_m = iter(rna_merged)
        rna_m = con_rna_m.next()

        con_loop = True
        while con_loop:
            if rna_m.overlap(rna):
                # Add to new dictionary
                try:
                    new_dict[rna_m].add(dna)
                except:
                    z = GenomicRegionSet(str(rna_m))
                    new_dict[rna_m] = z.add(dna)
                # next iteration
                try: rna, dna = con_rna.next(), con_dna.next()
                except: con_loop = False
            elif rna < rna_m:
                # next RNA and DNA
                try: rna, dna = con_rna.next(), con_dna.next()
                except: con_loop = False
            else:
                # next RNA_m
                try: rna_m = con_rna_m.next()
                except: con_loop = False
        self.merged_dict = new_dict

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
    
    def find_rbs(self, a_sequence, motif="RYMPA", min_len=10, max_len=None, organism="hg19"):
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
            organism:   Define the organism (hg19 as default)
        """
        tfos = RNABindingSet(name=a_sequence.name)

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
                        rbs = RNABinding(name=a_sequence.name, initial=i-count+1, final=i+1, 
                                         score=count, errors=0, motif=m, 
                                         orientation=None, seq=sample)
                        con_rbs = True

                    elif max_len and count > max_len:
                        tfos.add(rbs)
                        con_rbs = False
                else:  
                    sample = ""
                    count = 0
                    if con_rbs: 
                        tfos.add(rbs)
                        con_rbs = False
        return tfos


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
