# Python Libraries
from __future__ import print_function
from collections import *
import os
import multiprocessing
# Local Libraries
from scipy import stats
import numpy
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# Distal Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from SequenceSet import Sequence, SequenceSet
from RNADNABindingSet import RNADNABinding, RNADNABindingSet
from Util import SequenceType, Html, OverlapType
from rgt.motifanalysis.Statistics import multiple_test_correction
#import multiprocessing

####################################################################################
####################################################################################
def print2(summary, string):
    """ Show the message on the console and also save in summary. """
    print(string)
    summary.append(string)

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

def run_triplexator(ss, ds, output, l=None, e=None, c=None, fr=None, fm=None, of=None, mf=None, rm=None):
    """Perform Triplexator"""
    path_triplexator = "~/Apps/triplexator/bin/triplexator"

    arguments = " "
    if ss: arguments += "-ss "+ss+" "
    if ds: arguments += "-ds "+ds+" "
    if l: arguments += "-l "+str(l)+" "
    if e: arguments += "-e "+str(e)+" "
    if c: arguments += "-c "+str(c)+" "
    if fr: arguments += "-fr "+fr+" "
    if fm: arguments += "-fm "+str(fm)+" "
    if of: arguments += "-of "+str(of)+" "
    if mf: arguments += "-mf "
    if rm: arguments += "-rm "+str(rm)+" "
    
    if output: arguments += "> "+output
    arguments += " 2>> "+os.path.join(os.path.dirname(output),"triplexator_errors.txt")
    os.system(path_triplexator+arguments)


def random_each(input):
    """Return the counts of DNA Binding sites with randomization
    For multiprocessing. 
    Input contains:
    number, rna, region, temp, remove_temp, organism
    """
    print("Processing: "+str(input[0]))

    random = input[2].random_regions(organism=input[5], multiply_factor=1, 
                                     overlap_result=True, overlap_input=True, 
                                     chrom_X=True, chrom_M=False)
    txp = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[5], prefix=str(input[0]), remove_temp=input[4])
    return [ len(dbss) for dbss in txp.merged_dict.values() ]

def find_triplex(rna_fasta, dna_region, temp, organism, prefix="", remove_temp=False):
    """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
    
    arguments_triplexator = "-l 15 -e 15 -c 2 -fr off -fm 0 -of 1 -mf"
    # Generate BED 
    dna_region.write_bed(os.path.join(temp,"dna_"+prefix+".bed"))
    # Bedtools
    os.system("bedtools getfasta -fi /data/genome/"+organism+"/"+organism+".fa -bed "+\
              os.path.join(temp,"dna_"+prefix+".bed")+" -fo "+os.path.join(temp,"dna_"+prefix+".fa"))
    # Triplexator
    run_triplexator(ss=rna_fasta, ds=os.path.join(temp,"dna_"+prefix+".fa"), 
                    output=os.path.join(temp, "dna_"+prefix+".txp"), 
                    l=15, e=15, c=2, fr="off", fm=0, of=1, mf=True)
    # Read txp
    txp = RNADNABindingSet("dna")
    txp.read_txp(os.path.join(temp, "dna_"+prefix+".txp"))
    #print(len(txp_de))
    txp.merge_rbs()

    if remove_temp:
        os.remove(os.path.join(temp,"dna_region"+prefix+".bed"))
        os.remove(os.path.join(temp,"dna_"+prefix+".fa"))
        os.remove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp

def read_ac(self, path, cut_off):
    """Read the RNA accessibility file and output its positions and values

    The file should be a simple table with two columns:
    The first column is the position and the second one is the value
    '#' will be skipped

    example:
        #pos      u4S            
        1        NA  
        2        NA  
        3        NA  
        4     1.740  
        5     2.785  
        6     3.367  
        7     2.727  
        8     2.315  
        9     3.005  
        10     2.679  
        11     3.803  
        12     4.267  
        13     1.096  
    """
    #pos = []
    values = []
    with open(path) as f:
        for line in f:
            line = line.split()
            if not line: continue
            elif line[0][0] == "#": continue
            elif len(line) < 2: continue
            else:
                #pos.append(int(line[0]))
                v = line[1]
                if v == "NA": v = 0
                else: v = float(v)
                v = 2**(-v)
                if v >= cut_off:
                    v = 0.7 # Make it brighter in colormap. 1 is too dark for colormap
                else:
                    v = 0
                values.append(v)
    return values
#####################################################################################

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

class PromoterTest:
    """Test the association between given triplex and differential expression genes"""
    def __init__(self, gene_list_file, bed, bg, organism, promoterLength, summary):
        """Initiation"""
        self.organism = organism
        self.de_regions = GenomicRegionSet("de genes")
        self.nde_regions = GenomicRegionSet("nde genes")

        if bed and bg:
            self.de_regions.read_bed(bed)
            self.nde_regions.read_bed(bg)
        else:
            # DE gene regions
            self.de_gene = GeneSet("de genes")
            self.de_gene.read(gene_list_file)
            gene_list_file = os.path.basename(gene_list_file)
            print2(summary, "   \t"+str(len(self.de_gene))+" genes are loaded from: "+gene_list_file)

            if len(self.de_gene) == 0:
                print("Error: No genes are loaded from: "+gene_list_file)
                print("Pleas check the format.")
                sys.exit(1)
            self.de_regions.get_promotors(gene_set=self.de_gene, organism=organism, 
                                          promoterLength=promoterLength)
            
            # nonDE gene regions
            self.nde_gene = GeneSet("nde genes")
            self.nde_gene.get_all_genes(organism=organism)
            print2(summary, "   \t"+str(len(self.nde_gene))+" genes are loaded from "+organism+" genome data")
            self.nde_gene.subtract(self.de_gene)
            print2(summary, "   \t"+str(len(self.nde_gene))+" genes are not in "+gene_list_file)
            self.nde_regions.get_promotors(gene_set=self.nde_gene, organism=organism, 
                                           promoterLength=promoterLength)
            
    def search_triplex(self, rna, temp, l, e, remove_temp=False):
        # DE
        self.de_regions.write_bed(os.path.join(temp,"de_regions.bed"))

        #os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
        #          "/home/joseph/Downloads/0122result/hotair_genes_promoters_up.bed"+" -fo "+os.path.join(temp,"de.fa"))
        os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
                  os.path.join(temp,"de_regions.bed")+" -fo "+os.path.join(temp,"de.fa"))
        run_triplexator(ss=rna, ds=os.path.join(temp,"de.fa"), 
                        output=os.path.join(temp, "de.txp"), 
                        l=l, e=e, c=2, fr="off", fm=0, of=1, mf=True)
        # non-DE
        self.nde_regions.write_bed(os.path.join(temp,"nde_regions.bed"))
        #os.system("bedtools getfasta -fi /data/genome/hg19/hg19.fa -bed "+\
        #          "/home/joseph/Downloads/0122result/hotair_promoters_bg_up.bed"+" -fo "+os.path.join(temp,"nde.fa"))
        os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
                  os.path.join(temp,"nde_regions.bed")+" -fo "+os.path.join(temp,"nde.fa"))
        run_triplexator(ss=rna, ds=os.path.join(temp,"nde.fa"), 
                        output=os.path.join(temp, "nde.txp"), 
                        l=l, e=e, c=2, fr="off", fm=0, of=1, mf=True)
        if remove_temp:
            os.remove(os.path.join(temp,"de_regions.bed"))
            os.remove(os.path.join(temp,"de.fa"))
            os.remove(os.path.join(temp,"nde_regions.bed"))
            os.remove(os.path.join(temp,"nde.fa"))
        
    def count_frequency(self, temp, remove_temp):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""
        
        # Read txp and merge RBS
        txp_de = RNADNABindingSet("DE")
        #txp_de.read_txp("/projects/lncRNA/data/fendrr/fendrr_genes.txp")
        #txp_de.read_txp("/projects/lncRNA/data/joseph_results/joseph.txp")
        txp_de.read_txp(os.path.join(temp, "de.txp"))
        txp_de.merge_rbs(rm_duplicate=True)

        txp_nde = RNADNABindingSet("non-DE")
        #txp_nde.read_txp("/projects/lncRNA/data/fendrr/fendrr_genes_bg.txp")
        #txp_nde.read_txp("/projects/lncRNA/data/joseph_results/joseph_bg.txp")
        txp_nde.read_txp(os.path.join(temp, "nde.txp"))
        #print(len(txp_nde))
        txp_nde.merge_rbs(rm_duplicate=True)

        self.de_frequency = OrderedDict()
        self.nde_frequency = OrderedDict()
        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)
        for rbs, regions in txp_de.merged_dict.iteritems():
            #if len(regions) < 10: continue
            # DE
            #inter = len(regions.intersect(self.de_regions, mode=OverlapType.ORIGINAL))
            #inter = len(self.de_regions.intersect(regions, mode=OverlapType.ORIGINAL))
            inter = len(regions)
            self.de_frequency[rbs] = [inter, len_de - inter]
            # non-DE
            #inter = len(regions.intersect(self.nde_regions, mode=OverlapType.ORIGINAL))
            #self.nde_frequency[rbs] = [0, len_nde]
            for rbs_n, regions_n in txp_nde.merged_dict.iteritems():
                
                if rbs.overlap(rbs_n):
                    #print("Overlap RBS: "+rbs.region_str_rna()+"   "+rbs_n.region_str_rna()+"   "+str(len(regions))+"   "+str(len(regions_n)))
                    #inter = len(self.nde_regions.intersect(regions_n, mode=OverlapType.ORIGINAL))
                    
                    inter = len(regions_n)
                    self.nde_frequency[rbs] = [inter, len_nde - inter]
                    #print(self.nde_frequency[rbs])
                
        self.txp_de = txp_de
        self.txp_nde = txp_nde
        
        if remove_temp:
            os.remove(os.path.join(temp,"de.txp"))
            os.remove(os.path.join(temp,"nde.txp"))

    def fisher_exact(self):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        pvalues = []
        self.sig_region = []
        for rbs in self.de_frequency.keys():
            #print(rbs.toString())
            #print(self.de_frequency[rbs])
            #print(self.nde_frequency[rbs])
            #self.oddsratio[rbs], p = stats.fisher_exact([self.de_frequency[rbs], self.nde_frequency[rbs]])
            table = numpy.array([self.de_frequency[rbs], self.nde_frequency[rbs]])
            #table = numpy.transpose(table)
            self.oddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
            pvalues.append(p)

        # correction
        if len(self.de_frequency.keys()) > 1:
            reject, pvals_corrected = multiple_test_correction(pvalues, alpha=0.05, method='indep')
        else:
            pvals_corrected = pvalues
        for i, rbs in enumerate(self.de_frequency.keys()):
            self.pvalue[rbs] = pvals_corrected[i]
            if pvals_corrected[i] < 0.05:
                self.sig_region.append(rbs)

    def plot_frequency_rna(self, rna, dir, cut_off, ac=None):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(rna)
        self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])

        rbsp = self.txp_de.get_rbs(orientation="P")
        rbsa = self.txp_de.get_rbs(orientation="A")
        
        # Extract data points
        x = range(self.rna_len)
        p_y = []
        a_y = []
        for i in range(self.rna_len):
            p_y.append(rbsp.count_rbs_position(i))
            a_y.append(rbsa.count_rbs_position(i))
        all_y = [sum(z) for z in zip(p_y, a_y)]
        max_y = float(max(all_y) * 1.05)
        if ac: min_y = float(max_y*(-0.1))
        else: min_y = 0

        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        for rbs in self.sig_region:
            rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="r", 
                                     edgecolor="none", alpha=0.1, lw=None)
            ax.add_patch(rect)
        ax.plot(x, all_y, color="b", alpha=.5, lw=1, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="g", alpha=.5, lw=1, label="Parallel")
        ax.plot(x, a_y, color="r", alpha=.5, lw=1, label="Anti-parallel")
        
        
        # RNA accessbility
        if ac:
            n_value = read_ac(ac, cut_off)
            ac = numpy.array([n_value])
            ax.imshow(ac, cmap='Oranges', interpolation='nearest', extent=[0, self.rna_len, min_y, 0],
                      aspect='auto', label="Accessibility")


        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.)
        ax.set_xlim(left=0, right=self.rna_len )
        ax.set_ylim( [min_y, max_y] ) 
        
        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        ax.set_ylabel("Frequency of RNA binding sites",fontsize=12, rotation=90)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "rna_plot.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)

    def plot_de(self, dir):
        """Demonstrate the difference between DE genes and non DE genes"""
        
        x = range(self.rna_len)
        de_y = [0] * self.rna_len 
        nde_y = [0] * self.rna_len
        zeros = [0] * self.rna_len
        
        max_de = max([ m[0] for m in self.de_frequency.values() ])
        max_nde = max([ m[0] for m in self.nde_frequency.values() ])
        
        for rbs in self.txp.merged_dict.keys():
            de_y[rbs.initial:rbs.final] = [ self.de_frequency[rbs][0]/float(max_de) ] * (rbs.final - rbs.initial) 
            nde_y[rbs.initial:rbs.final] = [ 0 - self.nde_frequency[rbs][0]/float(max_nde) ] * (rbs.final - rbs.initial) 

        f, ax = plt.subplots(1, 1, dpi=300)
        ax.plot(x, de_y, color="b", lw=1, label="DE genes")
        ax.plot(x, nde_y, color="g", lw=1, label="non-DE genes")

        ax.set_ylim([-1, 1])
        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        #ax.set_ylabel("DE genes with DBS",fontsize=12, rotation=90)
        ax.set_yticks([-1, -0.5, 0, 0.5, 1])
        ax.set_yticklabels(["100 %", "non-DE gene", "0 %", "DE gene", "100 %"],rotation=0, ha="right")
        ax.grid(b=True, which='major', color='gray', linestyle='-')
        ax.fill_between(x, de_y, zeros, where=de_y>zeros, interpolate=True, color='b')
        ax.fill_between(x, zeros, nde_y, where=zeros>nde_y, interpolate=True, color='g')

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "de_plot.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)

    def gen_html(self, directory, align=50, alpha = 0.05):
        def subtable(data):
            code = "<table>"
            for row in data:
                code += "<tr>"
                for col in data:
                    code += "<td>" + col + "</td>"
                code += "</tr>"
            code += "</table>"
            return code

        #fp = os.path.join(dir,outputname,title)
        link_d = {os.path.basename(directory):"promoter.html",
                  "All triplex binding sites":"all_triplex.html"}
        html = Html(name="Triplex", links_dict=link_d, fig_dir=os.path.join(directory,"fig"), fig_rpath="./fig")
        #html.add_figure("projection_test.png", align="center")
        
        html.add_figure("rna_plot.png", align="center")
        #html.add_figure("de_plot.png", align="center")
        # Table of merged TBS on promoters
        header_list = ["RBS",
                       "DE_genes<br>with DBS", 
                       "DE_genes<br>no DBS",
                       "nonDE_genes<br>with DBS", 
                       "nonDE_genes<br>no DBS",
                       "Oddsratio",
                       "P_value"]
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]
        data_table = []
        for rbs in self.de_frequency.keys():
            #dbss = []
            #for dbs in self.txp_de.merged_dict[rbs]:
            #    dbss.append(dbs.toString())

            if self.pvalue[rbs] < alpha:
                data_table.append([ rbs.region_str_rna(), #str(len(dbss)), 
                                    value2str(self.de_frequency[rbs][0]), value2str(self.de_frequency[rbs][1]), 
                                    value2str(self.nde_frequency[rbs][0]), value2str(self.nde_frequency[rbs][1]), 
                                    value2str(self.oddsratio[rbs]), "<font color=\"red\">"+value2str(self.pvalue[rbs])+"</font>" ])
            else:
                data_table.append([ rbs.region_str_rna(), #str(len(dbss)), 
                                    value2str(self.de_frequency[rbs][0]), value2str(self.de_frequency[rbs][1]), 
                                    value2str(self.nde_frequency[rbs][0]), value2str(self.nde_frequency[rbs][1]), 
                                    value2str(self.oddsratio[rbs]), value2str(self.pvalue[rbs])])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        header_list=["Notes"]
        data_table = [["RNA name: "+ self.rna_name ],
                      ["RBS stands for RNA Binding Site."],
                      ["DBS stands for DNA Binding Site."]]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        

        html.add_free_content(['<a href="summary.log" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"promoter.html"))
    
        # Table of all triplex
        html = Html(name="Triplex", links_dict=link_d, fig_dir=os.path.join(directory,"fig"), fig_rpath="./fig")
        header_list=["RBS", "DBS", "Score", "Motif", "Orientation", "Strand"]
        data_table = []
        for tbs in self.txp_de:
            data_table.append([ tbs.rna.region_str_rna(), tbs.dna.toString(), tbs.score, 
                                tbs.motif, tbs.orient, tbs.strand])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.write(os.path.join(directory,"all_triplex.html"))
####################################################################################
####################################################################################

class RandomTest:
    def __init__(self, rna_fasta, dna_region, organism):
    	self.organism = organism
        # RNA: Path to the FASTA file
        self.rna_fasta = rna_fasta
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(self.rna_fasta)
        self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])
        # DNA: GenomicRegionSet
        self.dna_region = GenomicRegionSet(name="target")
        self.dna_region.read_bed(dna_region)
        
    def target_dna(self, temp, remove_temp):
        """Calculate the true counts of triplexes on the given dna regions"""
        
        txp = find_triplex(rna_fasta=self.rna_fasta, dna_region=self.dna_region, 
                           temp=temp, organism=self.organism, remove_temp=remove_temp)
        self.rbss = txp.merged_dict.keys()
        
        self.counts_target = OrderedDict()
              
        for rbs in self.rbss:
            self.counts_target[rbs] = len(txp.merged_dict[rbs])

    def random_test(self, repeats, temp, remove_temp):
        """Perform randomization for the given times"""
        
        # Prepare the input lists for multiprocessing
        mp_input = []
        for i in range(repeats):
            mp_input.append([ str(i), self.rna_fasta, self.dna_region,
                              temp, remove_temp, self.organism ])
        # Multiprocessing
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        mp_output = pool.map(random_each, mp_input)
        print(mp_output)
        pool.close()
        pool.join()
        
        # Processing the result
        matrix = []
        self.p_random = []
        self.mean_random = []
        self.sig_region = []
        for i, rbs in enumerate(self.rbss):
            print(i)
            row = [ v[i] for v in mp_output ]
            print(row)
            matrix.append(row)
            p = float(len([ h for h in row if h > self.counts_target[rbs] ]))/repeats
            self.p_random.append(p)
            if p < 0.05: self.sig_region.append(rbs)
            self.mean_random.append(numpy.mean(row))
        self.counts_random = numpy.array(matrix)
        self.sd_random = numpy.std(self.counts_random, axis=0)

    def plot(self, dir, ac, cut_off):
        """Generate the visualized plot"""
        
        # Extract data points
        x = range(self.rna_len)
        y = []

        for i in range(self.rna_len):
            y.append(self.rbss.count_rbs_position(i))
        
        max_y = float(max(y) * 1.05)

        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        for rbs in self.sig_region:
            rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="r", 
                                     edgecolor="none", alpha=0.1, lw=None)
            ax.add_patch(rect)
        ax.plot(x, y, color="b", alpha=.5, lw=1)

        #ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, mode="expand", borderaxespad=0.)
        ax.set_xlim(left=0, right=self.rna_len )
        ax.set_ylim( [0, max_y] ) 

        if ac:
            n_value = read_ac(ac, cut_off)
            ac = numpy.array([n_value])
            ax.imshow(ac, cmap='Oranges', interpolation='nearest', extent=[0, self.rna_len, min_y, 0],
                      aspect='auto', label="Accessibility")

        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        ax.set_ylabel("Frequency of RNA binding sites",fontsize=12, rotation=90)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "randomtest.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)

    def gen_html(self, directory, align=50, alpha = 0.05):
        """Generate the HTML file"""
        link_d = {os.path.basename(directory):"randomtest.html"}
        html = Html(name="Triplex", links_dict=link_d, fig_dir=os.path.join(directory,"fig"), fig_rpath="./fig")
        
        html.add_figure("randomtest.png", align="center")
        header_list = ["RBS",
                       "Count of<br>target regions<br>with DBS",
                       "Average of<br>randomization",
                       "SD of<br>randomization",
                       "P value"]
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]
        data_table = []
        
        for i, rbs in enumerate(self.rbss):
            if self.p_random[i] < alpha:
                data_table.append([ rbs.region_str_rna(), value2str(self.counts_target[rbs]),
                                    value2str(self.mean_random[i]), value2str(self.sd_random[i]), 
                                    "<font color=\"red\">"+value2str(self.p_random[i])+"</font>" ])
            else:
                data_table.append([ rbs.region_str_rna(), value2str(self.counts_target[rbs]),
                                    value2str(self.mean_random[i]), value2str(self.sd_random[i]), 
                                    value2str(self.p_random[i]) ])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        header_list=["Notes"]
        data_table = [["RNA name: "+ self.rna_name ],
                      ["RBS stands for RNA Binding Site."],
                      ["DBS stands for DNA Binding Site."]]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        html.add_free_content(['<a href="summary.log" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"randomtest.html"))