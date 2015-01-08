# Python Libraries
from __future__ import print_function
from collections import *
import os
import multiprocessing
# Local Libraries
from scipy import stats
# Distal Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from SequenceSet import Sequence, SequenceSet
from RNADNABindingSet import RNADNABinding, RNADNABindingSet
from Util import SequenceType, Html
#import multiprocessing

####################################################################################
####################################################################################
def mp_find_rbs(s, motif, min_len, max_len):
    triplex = TriplexSearch()
    return triplex.find_rbs(s, motif, min_len, max_len)
def mp_find_dbs(s, min_len, max_len):
    triplex = TriplexSearch()
    return triplex.find_dbs(s, min_len, max_len)

def value2str(value):
    if (isinstance(value,str)): return value
    if value == 0: return "0"
    if(isinstance(value,int)): return str(value)
    elif(isinstance(value,float)):
        if value >= 1000: 
            try: r = "{}".format(int(value))
            except: r = "Inf"
        elif 1000 > value > 10: r = "{:.1f}".format(value)
        elif 10 > value >= 1: r = "{:.2f}".format(value)
        elif 1 > value > 0.0001: r = "{:.4f}".format(value)
        else: r = "{:.1e}".format(value)
        return r

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
    
    def search_bindingsites(self, sequence_set, seq_type, motif, min_len, max_len, multiprocess=False):
        """Find binding sites on RNA or DNA from given SequenceSet
        
        Parameters:
            sequence_set:  A SequenceSet object
            seq_type:      RNA or DNA
            motif:         R, Y, M, P, or A
            min_len:       Define the minimum length of RBS (default is 10 bp)
            max_len:       Define the maximum length of RBS (default is infinite)
        """
        all_binding_sites = BindingSiteSet(name=sequence_set.name)
        if not multiprocess or len(sequence_set) == 1:
            for i,s in enumerate(sequence_set):
                if seq_type == SequenceType.RNA:
                    bs = self.find_rbs(s, motif, min_len, max_len)
                elif seq_type == SequenceType.DNA:
                    bs = self.find_dbs(s, min_len, max_len)
                all_binding_sites.concatenate(bs)
        else:
            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
            mp_input = [ [s, motif, min_len, max_len] for s in sequence_set ]
            print(mp_input)
            if seq_type == SequenceType.RNA:
                bs = pool.map(mp_find_rbs, mp_input)
            elif seq_type == SequenceType.DNA:
                bs = pool.map(mp_find_dbs, mp_input)
            pool.close()
            pool.join()
            for b in bs:
                all_binding_sites.concatenate(b)

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
            region:      A BindingSite or GenomicRegion
        """
        
        all_dbs = BindingSiteSet(name=a_sequence.name)
        if len(a_sequence) < min_len: return
        
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
        #l = len(a_sequence)

        for i, a in enumerate(a_sequence.seq):
        #    print(str(l)+"\t"+str(i))
            if a in targets:
                # Strand 2
                if con_2:
                    sample_2 = sample_2.replace("C","G")
                    sample_2 = sample_2.replace("T","A")
                    dbs_2 = BindingSite(chrm=a_sequence.name, initial=i-count_2+1, final=i+1, 
                                        score=count_2, errors_bp=0, strand=c_strand, 
                                        seq=Sequence(seq=sample_2, strand=c_strand))
                    all_dbs.add(dbs_2)
                    con_2 = False
                sample_2 = ""
                count_2 = 0

                # Strand 1 
                if count_1 < min_len:
                    sample_1 += a
                    count_1 += 1
                elif min_len <= count_1:
                    sample_1 += a
                    count_1 += 1
                    
                    con_1 = True
                elif max_len and count_1 > max_len:
                    dbs_1 = BindingSite(chrm=a_sequence.name, initial=i-count_1+1, final=i+1, 
                                        score=count_1, errors_bp=0, strand=a_sequence.strand, 
                                        seq=Sequence(seq=sample_1, strand=a_sequence.strand))
                    all_dbs.add(dbs_1)
                    con_1 = False
            else:
                # Strand 1
                if con_1: 
                    dbs_1 = BindingSite(chrm=a_sequence.name, initial=i-count_1+1, final=i+1, 
                                        score=count_1, errors_bp=0, strand=a_sequence.strand, 
                                        seq=Sequence(seq=sample_1, strand=a_sequence.strand))
                    all_dbs.add(dbs_1)
                    con_1 = False
                sample_1 = ""
                count_1 = 0
                # Strand 2
                if count_2 < min_len:
                    sample_2 += a
                    count_2 += 1
                elif min_len <= count_2:
                    sample_2 += a
                    count_2 += 1
                    con_2 = True
                elif max_len and count_2 > max_len:
                    sample_2 = sample_2.replace("C","G")
                    sample_2 = sample_2.replace("T","A")
                    dbs_2 = BindingSite(chrm=a_sequence.name, initial=i-count_2+1, final=i+1, 
                                        score=count_2, errors_bp=0, strand=c_strand, 
                                        seq=Sequence(seq=sample_2, strand=c_strand))
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

class FisherTest:
    """Test the association between given triplex and differential expression genes"""
    def __init__(self, gene_list_file, organism, promoterLength):
        """Initiation"""
        self.organism = organism

        # DE gene regions
        self.de_gene = GeneSet("de genes")
        self.de_gene.read(gene_list_file)
        self.de_regions = GenomicRegionSet("de genes")
        self.de_regions.get_promotors(gene_set=self.de_gene, organism=organism, 
                                      promoterLength=promoterLength)
        
        # nonDE gene regions
        self.nde_gene = GeneSet("nde genes")
        self.nde_gene.get_all_genes(organism=organism)
        self.nde_gene.subtract(self.de_gene)
        self.nde_regions = GenomicRegionSet("nde genes")
        self.nde_regions.get_promotors(gene_set=self.nde_gene, organism=organism, 
                                       promoterLength=promoterLength)

    def search_triplex(self, rna, temp, remove_temp=False):
        # DE
        self.de_regions.write_bed(os.path.join(temp,"de_regions.bed"))
        os.system("bedtools getfasta -fi /data/genome/hg19genome.fa -bed "+\
                  os.path.join(temp,"de_regions.bed")+" -fo "+os.path.join(temp,"de.fa"))
        os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -l 15 -e 15 -c 2 -fr off -fm 0 -of 1 -mf -ss "+\
                  rna+" -ds "+os.path.join(temp,"de.fa")+" > "+os.path.join(temp, "de.txp"))


        # non-DE
        self.nde_regions.write_bed(os.path.join(temp,"nde_regions.bed"))
        os.system("bedtools getfasta -fi /data/genome/hg19genome.fa -bed "+\
                  os.path.join(temp,"nde_regions.bed")+" -fo "+os.path.join(temp,"nde.fa"))
        os.system("/projects/lncRNA/bin/triplexator/bin/triplexator -ss "+rna+" -ds "+os.path.join(temp,"nde.fa")+\
                  " > "+os.path.join(temp, "nde.txp"))
        if remove_temp:
            os.remove(os.path.join(temp,"de_regions.bed"))
            os.remove(os.path.join(temp,"de.fa"))
            os.remove(os.path.join(temp,"nde_regions.bed"))
            os.remove(os.path.join(temp,"nde.fa"))
        
    def count_frequency(self, temp):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""
        
        # Read txp and merge RBS
        txp = RNADNABindingSet("all")
        txp.read_txp(os.path.join(temp, "de.txp"))
        txp.read_txp(os.path.join(temp, "nde.txp"))
        txp.merge_rbs()
        
        self.de_frequency = OrderedDict()
        self.nde_frequency = OrderedDict()
        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)
        for rbs, regions in txp.merged_dict.iteritems():
            # DE
            inter = len(self.de_regions.intersect(regions))
            self.de_frequency[rbs] = [inter, len_de - inter]
            # non-DE
            inter = len(self.nde_regions.intersect(regions))
            self.nde_frequency[rbs] = [inter, len_nde - inter]
        self.txp = txp

    def fisher_exact(self):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        for rbs in self.txp.merged_dict.keys():
            self.oddsratio[rbs], self.pvalue[rbs] = stats.fisher_exact([self.de_frequency[rbs], self.nde_frequency[rbs]])
        
    def gen_html(self, directory, align=50, alpha = 0.05):
        #fp = os.path.join(dir,outputname,title)
        link_d = {os.path.basename(directory):"fisher.html"}
        html = Html(name="Triplex", links_dict=link_d, fig_dir=os.path.join(os.path.dirname(directory),"fig"))
        #html.add_figure("projection_test.png", align="center")
        
        header_list = ["RNA Binding Site",
                       "DE genes<br>RBS binding", 
                       "DE genes<br>no binding",
                       "nonDE genes<br>RBS binding", 
                       "nonDE genes<br>no binding",
                       "Oddsratio",
                       "p-value"]
        
        type_list = 'sssssss'
        col_size_list = [10,10,10,10,10,10,10]
        data_table = []
        for rbs in self.txp.merged_dict.keys():
            if self.pvalue[rbs] < alpha:
                data_table.append([ rbs.region_str(), value2str(self.de_frequency[rbs][0]), value2str(self.de_frequency[rbs][1]), 
                                    value2str(self.nde_frequency[rbs][0]), value2str(self.nde_frequency[rbs][1]), 
                                    value2str(self.oddsratio[rbs]), "<font color=\"red\">"+value2str(self.pvalue[rbs])+"</font>" ])
            else:
                data_table.append([ rbs.region_str(), value2str(self.de_frequency[rbs][0]), value2str(self.de_frequency[rbs][1]), 
                                    value2str(self.nde_frequency[rbs][0]), value2str(self.nde_frequency[rbs][1]), 
                                    value2str(self.oddsratio[rbs]), value2str(self.pvalue[rbs])])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align)
        
        header_list=["Assumptions and hypothesis"]
        data_table = []
        
        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"fisher.html"))
    
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
