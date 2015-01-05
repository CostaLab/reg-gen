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
                bs = self.find_dbs(s, min_len, max_len)
            all_binding_sites.concatenate(bs)

        return all_binding_sites

    def find_rbs(self, a_sequence, motif, min_len, max_len):
        """Return a BindingSiteSet with the potential RNA binding sites

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
                                          name=a_sequence.name, score=count,  errors_bp=0, 
                                          motif=m, seq=Sequence(seq=sample, strand="RNA"))
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


    def find_dbs(self, a_sequence, min_len, max_len):
        """Return a BindingSiteSet with the potential DNA triplex binding sites
        (including both strands)

        Parameters:
            a_sequence:  A Sequence object
            min_len:     Define the minimum length of RBS (default is 10 bp)
            max_len:     Define the maximum length of RBS (default is infinite)
        """
        
        all_dbs = BindingSiteSet(name=a_sequence.name)
        targets = ["A", "G"]
        if a_sequence.strand == "+": c_strand = "-"
        else: c_strand = "+"
        # Parsing the sequence
        count_1 = 0
        count_2 = 0
        sample_1 = ""
        sample_2 = ""
        con_1 = False
        con_2 = False
        
        for i, a in enumerate(a_sequence.seq):
            if a in targets:
                # Strand 2
                sample_2 = ""
                count_2 = 0
                if con_2: 
                    all_dbs.add(dbs)
                    con_2 = False
                    
                # Strand 1 
                if count_1 < min_len:
                    sample_1 += a
                    count_1 += 1
                elif min_len <= count_1:
                    sample_1 += a
                    count_1 += 1
                    dbs_1 = BindingSite(chrm=a_sequence.name, initial=i-count_1+1, final=i+1, 
                                        score=count_1, errors=0, strand=a_sequence.strand, 
                                        seq=Sequence(seq=sample, strand=a_sequence.strand))
                    con_1 = True
                elif max_len and count_1 > max_len:
                    all_dbs.add(dbs_1)
                    con_1 = False
            else:
                # Strand 1
                sample_1 = ""
                count_1 = 0
                if con_1: 
                    all_dbs.add(dbs_1)
                    con_1 = False
                    
                # Strand 2 
                if count_2 < min_len:
                    sample_2 += a
                    count_2 += 1
                elif min_len <= count_2:
                    sample_2 += a
                    count_2 += 1
                    dbs_2 = BindingSite(chrm=a_sequence.name, initial=i-count_2+1, final=i+1, 
                                        score=count_2, errors=0, strand=c_strand, 
                                        seq=Sequence(seq=sample, strand=c_strand))
                    con_2 = True
                elif max_len and count_2 > max_len:
                    all_dbs.add(dbs_2)
                    con_2 = False
        return all_dbs
        

    #def find_triplex(self):

    def compare_rb_db(self, rna_seq, dna_seq, match_dict):
        """Find the highest score between two sequences
        Parameters:
            rna_seq:     A RNA sequence
            dan_seq:     A DNA sequence
            match_dict:  A dictionary with defined matching cases
        """
        if not rna_seq.seq or not dna_seq.seq:
            print("*** Error: Compare empity sequence for local alignment score.")
            print("*** "+str(rna_seq))
            print("*** "+str(dna_seq))
            return
        else:
            m = len(rna_seq.seq)
            n = len(dna_seq.seq)
            
            if m == n:
                score = 0
                match = ""
                for i, r in enumerate(rna_seq.seq):
                    if dna_seq.seq[i] in match_dict[r]:
                        score += 1
                        match += "|"
                    else:
                        match += "*"
                        
                return score, rna_seq.seq, match, dna_seq.seq
            
            elif m > n:
                scores = []
                matches = []
                for i in range(m-n):
                    score = 0
                    match = ""
                    for j in range(n):
                        if dna_seq.seq[i] in match_dict[rna_seq.seq[i+j]]:
                            score += 1
                            match += "|"
                        else:
                            match += "*"
                    scores.append(score)
                    matches.append(match)
                hp = [ i for i, j in enumerate(scores) if j == max(scores) ]
                #######################################
                ####
                ####  ( O O )==<   ppbbppbbppbbpbpbpb
                ####    ---
                #######################################
                return 

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
