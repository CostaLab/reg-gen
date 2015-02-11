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
from matplotlib import colors

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
    number, rna, region, temp, remove_temp, organism, rbss, step_interval
    """

    random = input[2].random_regions(organism=input[5], multiply_factor=1, 
                                     overlap_result=True, overlap_input=True, 
                                     chrom_X=True, chrom_M=False)
    txp = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[5], prefix=str(input[0]), remove_temp=input[4])
    txp.merge_by(rbss=input[6])

    if int(input[0])*41 % int(input[7]) == 0: print("=", end="")

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
    txp.read_txp(os.path.join(temp, "dna_"+prefix+".txp"), dna_fine_posi=True)
    #print(len(txp_de))

    if remove_temp:
        os.remove(os.path.join(temp,"dna_region"+prefix+".bed"))
        os.remove(os.path.join(temp,"dna_"+prefix+".fa"))
        os.remove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp

def read_ac(path, cut_off, rnalen):
    """Read the RNA accessibility file and output its positions and values

    The file should be a simple table with two columns:
    The first column is the position and the second one is the value
    '#' will be skipped

    """
    #pos = []
    access = []
    with open(path) as f:
        i = 0
        while i < rnalen:
            for line in f:
                line = line.split()
                if not line: continue
                elif line[0][0] == "#": continue
                elif len(line) < 2: continue
                else:
                    #pos.append(int(line[0]))
                    #print(line)
                    v = line[1]
                    if v == "NA": 
                        access.append(0)
                    else: 
                        try: v = 2**(-float(v))
                        except: continue
                        if v >= cut_off:
                            access.append(1)
                        else:
                            access.append(0)
                    i += 1
    return access

def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try:
        os.stat(path)
    except:
        os.mkdir(path)
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
    def __init__(self, gene_list_file, bed, bg, organism, promoterLength, rna_name, summary):
        """Initiation"""
        self.organism = organism
        self.rna_name = rna_name
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
        
        

        self.frequency = {}
        self.frequency["promoters"] = { "de": OrderedDict(), "nde": OrderedDict() }
        self.frequency["hits"] = { "de": OrderedDict(), "nde": OrderedDict() }
        
        ########################################################
        # Count the number of targeted promoters on each merged DNA Binding Domain
        txp_de = RNADNABindingSet("DE")
        txp_nde = RNADNABindingSet("non-DE")
        txp_de.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=False)
        txp_de.merge_rbs(rm_duplicate=True)
        txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        txp_nde.merge_by(rbss=txp_de.merged_dict.keys(), rm_duplicate=True)

        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)

        for rbs in txp_de.merged_dict.keys():
            #if len(regions) < 10: continue
            # DE
            self.frequency["promoters"]["de"][rbs] = [ len(txp_de.merged_dict[rbs]), len_de - len(txp_de.merged_dict[rbs]) ]
            #txp_de.merged_dict[rbs].write_bed("de.bed")
            # non-DE
            self.frequency["promoters"]["nde"][rbs] = [ len(txp_nde.merged_dict[rbs]), len_nde - len(txp_nde.merged_dict[rbs]) ]
            #txp_nde.merged_dict[rbs].write_bed("nde.bed")

        self.txp_de = txp_de
        self.txp_nde = txp_nde

        ########################################################
        # Count the number of hits on the promoters from each merged DBD 
        txp_def = RNADNABindingSet("DE")
        txp_ndef = RNADNABindingSet("non-DE")
        txp_def.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=True)
        txp_def.merge_rbs(rm_duplicate=True)
        txp_ndef.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=True)
        txp_ndef.merge_by(rbss=txp_def.merged_dict.keys(), rm_duplicate=True)

        for rbs in txp_def.merged_dict.keys():
            #if len(regions) < 10: continue
            # DE
            self.frequency["hits"]["de"][rbs] = len(txp_def.merged_dict[rbs])
            # non-DE
            self.frequency["hits"]["nde"][rbs] = len(txp_ndef.merged_dict[rbs])

        self.txp_def = txp_def
        self.txp_ndef = txp_ndef
        
        if remove_temp:
            os.remove(os.path.join(temp,"de.txp"))
            os.remove(os.path.join(temp,"nde.txp"))

    def fisher_exact(self):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        pvalues = []
        self.sig_region = []
        for rbs in self.frequency["promoters"]["de"]:
            table = numpy.array([self.frequency["promoters"]["de"][rbs], self.frequency["promoters"]["nde"][rbs]])
            self.oddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
            pvalues.append(p)

        # correction
        if len(self.frequency["promoters"]["de"].keys()) > 1:
            reject, pvals_corrected = multiple_test_correction(pvalues, alpha=0.05, method='indep')
        else:
            pvals_corrected = pvalues
        for i, rbs in enumerate(self.frequency["promoters"]["de"].keys()):
            self.pvalue[rbs] = pvals_corrected[i]
            if pvals_corrected[i] < 0.05:
                self.sig_region.append(rbs)

    def plot_promoter(self, rna, dir, cut_off, ac=None):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(rna)
        if not self.rna_name: self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])
        
        # Extract data points
        x = range(self.rna_len)
        p_y = [0] * self.rna_len
        a_y = [0] * self.rna_len

        for rbsm in self.txp_de.merged_dict.keys():
            p = 0
            a = 0
            for promoter in self.txp_de.merged_dict[rbsm]:
                if promoter.tri_orien == "P": p += 1
                if promoter.tri_orien == "A": a += 1
            for i in range(rbsm.initial, rbsm.final):
                p_y[i] += p
                a_y[i] += a

        all_y = [sum(z) for z in zip(p_y, a_y)]
        max_y = float(max(all_y) * 1.1)
        if ac: min_y = float(max_y*(-0.1))
        else: min_y = 0

        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        for rbs in self.sig_region:
            rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="lightgray", 
                                     edgecolor="none", alpha=0.5, lw=None)
            ax.add_patch(rect)
        ax.plot(x, all_y, color="r", alpha=.8, lw=1, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="g", alpha=.8, lw=1, label="Parallel")
        ax.plot(x, a_y, color="b", alpha=.8, lw=1, label="Anti-parallel")
        
        # RNA accessbility
        if ac:
            n_value = read_ac(ac, cut_off, rnalen=self.rna_len)
            #acn = numpy.array([n_value])
            #print(acn)
            ax.imshow([n_value], cmap=colors.ListedColormap(['white', 'orange']), interpolation='nearest', 
                      extent=[0, self.rna_len, min_y, 0], aspect='auto', label="Accessibility")

        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
                  prop={'size':9}, ncol=3)
        ax.set_xlim(left=0, right=self.rna_len )
        ax.set_ylim( [min_y, max_y] ) 
        
        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        ax.set_ylabel("Number of related promoters",fontsize=12, rotation=90)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "promoter.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)


    def plot_dbss(self, rna, dir, cut_off, ac=None):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        # Extract data points
        x = range(self.rna_len)
        p_y = [0] * self.rna_len
        a_y = [0] * self.rna_len

        for rbsm in self.txp_def.merged_dict.keys():
            p = 0
            a = 0
            for promoter in self.txp_def.merged_dict[rbsm]:
                if promoter.tri_orien == "P": p += 1
                if promoter.tri_orien == "A": a += 1
            for i in range(rbsm.initial, rbsm.final):
                p_y[i] += p
                a_y[i] += a

        all_y = [sum(z) for z in zip(p_y, a_y)]
        max_y = float(max(all_y) * 1.1)
        if ac: min_y = float(max_y*(-0.1))
        else: min_y = 0

        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        for rbs in self.sig_region:
            rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="lightgray", 
                                     edgecolor="none", alpha=0.5, lw=None)
            ax.add_patch(rect)
        ax.plot(x, all_y, color="r", alpha=.8, lw=1, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="g", alpha=.8, lw=1, label="Parallel")
        ax.plot(x, a_y, color="b", alpha=.8, lw=1, label="Anti-parallel")
        
        # RNA accessbility
        if ac:
            n_value = read_ac(ac, cut_off, rnalen=self.rna_len)
            #acn = numpy.array([n_value])
            #print(acn)
            ax.imshow([n_value], cmap=colors.ListedColormap(['white', 'orange']), interpolation='nearest', 
                      extent=[0, self.rna_len, min_y, 0], aspect='auto', label="Accessibility")

        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
                  prop={'size':9}, ncol=3)
        ax.set_xlim(left=0, right=self.rna_len )
        ax.set_ylim( [min_y, max_y] ) 
        
        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        ax.set_ylabel("Number of related DNA binding sites",fontsize=12, rotation=90)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "binding_sites.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)

    def plot_frequency_rna(self, rna, dir, cut_off, ac=None):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """

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
        max_y = float(max(all_y) * 1.1)
        if ac: min_y = float(max_y*(-0.1))
        else: min_y = 0

        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        for rbs in self.sig_region:
            rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="lightgray", 
                                     edgecolor="none", alpha=0.5, lw=None)
            ax.add_patch(rect)
        ax.plot(x, all_y, color="r", alpha=.8, lw=1, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="g", alpha=.8, lw=1, label="Parallel")
        ax.plot(x, a_y, color="b", alpha=.8, lw=1, label="Anti-parallel")
        
        
        # RNA accessbility
        if ac:
            n_value = read_ac(ac, cut_off, rnalen=self.rna_len)
            #acn = numpy.array([n_value])
            #print(acn)
            ax.imshow([n_value], cmap=colors.ListedColormap(['white', 'orange']), interpolation='nearest', 
                      extent=[0, self.rna_len, min_y, 0], aspect='auto', label="Accessibility")

        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
                  prop={'size':9}, ncol=3)
        ax.set_xlim(left=0, right=self.rna_len )
        ax.set_ylim( [min_y, max_y] ) 
        
        ax.set_xlabel("RNA sequence (bp)", fontsize=12)
        ax.set_ylabel("Number of DNA binding domains",fontsize=12, rotation=90)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "plot_rna.png"), facecolor='w', edgecolor='w',  
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

        check_dir(os.path.join(directory, "supplements"))
        html_header = "Triplex Domain Finder: Promoter Test"
        #fp = os.path.join(dir,outputname,title)
        self.link_d = OrderedDict()
        self.link_d["On RNA"] = "index.html"
        self.link_d["On DE Genes"] = "genes.html"
        self.link_d["On Background"] = "background.html"
        self.link_ds = OrderedDict()
        self.link_ds["On RNA"] = "../index.html"
        self.link_ds["On DE Genes"] = "../genes.html"
        self.link_ds["On Background"] = "../background.html"

        if self.organism == "hg19": self.ani = "human"
        elif self.organism == "mm9": self.ani = "mouse"

        #############################################################
        # DNA Binding Domains (on RNA)
        #############################################################
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        #html.add_figure("projection_test.png", align="center")
        html.add_figure("plot_rna.png", align="center")
        html.add_figure("promoter.png", align="center")
        html.add_figure("binding_sites.png", align="center")
        
        # Table of merged TBS on promoters
        header_list = ["DNA Binding Domain",
                       "Binding<br>DE_promoters", 
                       "Else<br>DE_promoters",
                       "Binding<br>Background", 
                       "Else<br>Background",
                       "Oddsratio",
                       "P_value",
                       "Unique_DBSs<br>DE_promoters",
                       "Unique_DBSs<br>Background" ]
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]
        data_table = []
        for rbs in self.frequency["promoters"]["de"].keys():
            #dbss = []
            #for dbs in self.txp_de.merged_dict[rbs]:
            #    dbss.append(dbs.toString())

            if self.pvalue[rbs] < alpha:
                data_table.append([ rbs.region_str_rna(), #str(len(dbss)), 
                                    '<a href="'+"supplements/pro_de_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["de"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["de"][rbs][1]), 
                                    '<a href="'+"supplements/pro_nde_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["nde"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["nde"][rbs][1]), 
                                    value2str(self.oddsratio[rbs]), 
                                    "<font color=\"red\">"+value2str(self.pvalue[rbs])+"</font>",
                                    value2str(self.frequency["hits"]["de"][rbs]),
                                    value2str(self.frequency["hits"]["nde"][rbs]) ])
            else:
                data_table.append([ rbs.region_str_rna(), #str(len(dbss)), 
                                    '<a href="'+"supplements/pro_de_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["de"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["de"][rbs][1]), 
                                    '<a href="'+"supplements/pro_nde_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["nde"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["nde"][rbs][1]), 
                                    value2str(self.oddsratio[rbs]), value2str(self.pvalue[rbs]),
                                    '<a href="'+"supplements/dbs_de_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["hits"]["de"][rbs])+'</a>', 
                                    '<a href="'+"supplements/dbs_nde_"+rbs.region_str_rna()+".html"+'" style="text-align:left">'+
                                    value2str(self.frequency["hits"]["nde"][rbs])+'</a>' ])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        
        html.add_heading("Notes")
        html.add_list([ "RNA name: "+ self.rna_name,
                        "DBD stands for DNA Binding Domain on RNA.",
                        "DBS stands for DNA Binding Site on DNA."])

        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"index.html"))

        #############################################################
        # Profile of targeted promoters for each merged DNA Binding Domain
        #############################################################
        
        header_list=["Targeted Promoter", "Gene", "Strand", "Width", "Motif", "Orientation"]
        # DE
        for rbsm in self.frequency["promoters"]["de"].keys():
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            data_table = []
            html.add_heading("Targeted promoters of DE genes relate to the merged domain: "+rbsm.region_str_rna())
            for dbs in self.txp_de.merged_dict[rbsm]:
                for promoter in self.de_regions:
                    if dbs.overlap(promoter):
                        proname = promoter.name
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+dbs.chrom+"%3A"+str(dbs.initial)+"-"+str(dbs.final)+
                                    '" style="text-align:left">'+dbs.toString()+'</a>', 
                                    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+self.ani+
                                    "&db="+self.organism+"&singleSearch=knownCanonical&position="+
                                    proname+'" style="text-align:left">'+proname+'</a>', 
                                    dbs.orientation, str(dbs.final-dbs.initial),
                                    dbs.motif, dbs.tri_orien ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "pro_de_"+rbsm.region_str_rna()+".html"))

        # non-DE
        for rbsm in self.frequency["promoters"]["nde"].keys():
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            data_table = []
            html.add_heading("Targeted promoters in Background relate to the merged domain: "+rbsm.region_str_rna())
            for dbs in self.txp_nde.merged_dict[rbsm]:
                for promoter in self.nde_regions:
                    if dbs.overlap(promoter):
                        proname = promoter.name
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+dbs.chrom+"%3A"+str(dbs.initial)+"-"+str(dbs.final)+
                                    '" style="text-align:left">'+dbs.toString()+'</a>', 
                                    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+self.ani+
                                    "&db="+self.organism+"&singleSearch=knownCanonical&position="+
                                    proname+'" style="text-align:left">'+proname+'</a>',
                                    dbs.orientation, str(dbs.final-dbs.initial),
                                    dbs.motif, dbs.tri_orien ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "pro_nde_"+rbsm.region_str_rna()+".html"))
        
        #############################################################
        # Profile of Hits for each merged DNA Binding Domain
        #############################################################
        
        header_list=["DNA_Binding_Site", "Strand", "Width", "Score", "Motif", "Orientation"]
        # DE
        for rbsm in self.frequency["hits"]["de"].keys():
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            data_table = []
            html.add_heading("DNA Binding Sites on DE promoters relate to the merged domain: "+rbsm.region_str_rna())
            for dbs in self.txp_def.merged_dict[rbsm]:
                data_table.append([ dbs.toString(), dbs.orientation, str(dbs.final-dbs.initial),
                                    dbs.score, dbs.motif, dbs.tri_orien ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "dbs_de_"+rbsm.region_str_rna()+".html"))

        # non-DE
        for rbsm in self.frequency["hits"]["nde"].keys():
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            data_table = []
            html.add_heading("DNA Binding Sites in Background relate to the merged domain: "+rbsm.region_str_rna())
            for dbs in self.txp_ndef.merged_dict[rbsm]:
                data_table.append([ dbs.toString(), dbs.orientation, str(dbs.final-dbs.initial),
                                    dbs.score, dbs.motif, dbs.tri_orien ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "dbs_nde_"+rbsm.region_str_rna()+".html"))


    def gen_html_genes(self, directory, align=50, alpha = 0.05, nonDE=False):

        html_header = "Triplex Domain Finder: Promoter Test"
        #fp = os.path.join(dir,outputname,title)
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]

        #############################################################
        # Gene centered
        #############################################################
        
        ##############################################################################################
        # DE
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        header_list=["Promoter", "Gene", "Total Bindings", "Parallel Bindings", "Antiparallel Bindings"]
        html.add_heading("The promoters of differential expression genes")
        data_table = []
        counts_p = self.de_regions.counts_per_region(self.txp_de.get_dbs(sort=True, orientation="P"))
        counts_a = self.de_regions.counts_per_region(self.txp_de.get_dbs(sort=True, orientation="A"))
        
        if not self.de_regions.sorted: self.de_regions.sort()
        # Iterate by each gene promoter
        for i, promoter in enumerate(self.de_regions):
            if counts_p[i] == 0 and counts_a[i] == 0:
                continue
            else:
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString()+'</a>',
                                    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+self.ani+
                                    "&db="+self.organism+"&singleSearch=knownCanonical&position="+
                                    promoter.name+'" style="text-align:left">'+promoter.name+'</a>',
                                    '<a href="'+os.path.join(directory,"supplements", "promoter_de"+str(i)+"_"+
                                    promoter.toString()+".html")+
                                    '" style="text-align:left">'+str(counts_p[i]+counts_a[i])+'</a>',
                                    str(counts_p[i]), str(counts_a[i])])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_heading("Notes")
        html.add_list(["All the genes without any bindings are ignored."])
        
        html.write(os.path.join(directory,"genes.html"))

        ############################
        header_list=["DNA_Binding_Domain", "DNA Binding Site", "Width", "Score", "Motif", "Orientation"]
        for i, promoter in enumerate(self.de_regions):
            if counts_p[i] == 0 and counts_a[i] == 0:
                continue
            else: 
                html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                            fig_rpath="../style", RGT_header=False)
                html.add_heading(promoter.name)
                data_table = []

                for rd in self.txp_de.sequences:
                    if rd.dna.overlap(promoter):
                        data_table.append([ rd.rna.region_str_rna(),
                                            '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                            "&position="+rd.dna.chrom+"%3A"+str(rd.dna.initial)+"-"+str(rd.dna.final)+'" style="text-align:left">'+
                                            rd.dna.toString()+'</a>',
                                            str(rd.dna.final-rd.dna.initial),
                                            rd.dna.score, 
                                            rd.rna.motif, rd.rna.orientation])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "promoter_de"+str(i)+"_"+promoter.toString()+".html"))

        ##############################################################################################
        # Non-DE
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        header_list=["Promoter", "Gene", "Total Bindings", "Parallel Bindings", "Antiparallel Bindings"]

        html.add_heading("The promoters of Background")
        data_table = []
        counts_p = self.nde_regions.counts_per_region(self.txp_nde.get_dbs(sort=True, orientation="P"))
        counts_a = self.nde_regions.counts_per_region(self.txp_nde.get_dbs(sort=True, orientation="A"))
        
        if not self.nde_regions.sorted: self.nde_regions.sort()
        # Iterate by each gene promoter
        for i, promoter in enumerate(self.nde_regions):
            if counts_p[i] == 0 and counts_a[i] == 0:
                continue
            else:
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString()+'</a>',
                                    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+self.ani+
                                    "&db="+self.organism+"&singleSearch=knownCanonical&position="+
                                    promoter.name+'" style="text-align:left">'+promoter.name+'</a>',
                                    '<a href="'+os.path.join(directory,"supplements", "promoter_bg"+str(i)+"_"+
                                    promoter.toString()+".html")+
                                    '" style="text-align:left">'+str(counts_p[i]+counts_a[i])+'</a>',
                                    str(counts_p[i]), str(counts_a[i])])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_heading("Notes")
        html.add_list(["All the genes without any bindings are ignored."])
        
        html.write(os.path.join(directory,"background.html"))
 
        ############################
        header_list=["DNA_Binding_Domain", "DNA Binding Site", "Width", "Score", "Motif", "Orientation"]
        for i, promoter in enumerate(self.nde_regions):
            if counts_p[i] == 0 and counts_a[i] == 0:
                continue
            else: 
                html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                            fig_rpath="../style", RGT_header=False)
                html.add_heading(promoter.name)
                data_table = []

                for rd in self.txp_nde.sequences:
                    if rd.dna.overlap(promoter):
                        data_table.append([ rd.rna.region_str_rna(),
                                            '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                            "&position="+rd.dna.chrom+"%3A"+str(rd.dna.initial)+"-"+str(rd.dna.final)+'" style="text-align:left">'+
                                            rd.dna.toString()+'</a>',
                                            str(rd.dna.final-rd.dna.initial),
                                            rd.dna.score, 
                                            rd.rna.motif, rd.rna.orientation])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "promoter_bg"+str(i)+"_"+promoter.toString()+".html"))


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
        txp.merge_rbs()
        self.rbss = txp.merged_dict.keys()
        self.counts_target = [ len(txp.merged_dict[rbs]) for rbs in self.rbss ]
        
    def random_test(self, repeats, temp, remove_temp):
        """Perform randomization for the given times"""
        
        # Prepare the input lists for multiprocessing
        mp_input = []
        for i in range(repeats):
            mp_input.append([ str(i), self.rna_fasta, self.dna_region,
                              temp, remove_temp, self.organism, self.rbss, str(repeats) ])
        # Multiprocessing
        print("Processing:    |0%                  |                100%|")
        print("               [", end="")
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        mp_output = pool.map(random_each, mp_input)
        #print(mp_output)
        pool.close()
        pool.join()
        print("]")
        
        # Processing the result
        matrix = []
        self.p_random = []
        self.mean_random = []
        self.sig_region = []
        for i, rbs in enumerate(self.rbss):
            #print(i)
            counts_rbs = [ v[i] for v in mp_output ]
            #print(len(counts_rbs))
            matrix.append(counts_rbs)
            p = float(len([ h for h in counts_rbs if h > self.counts_target[i] ]))/repeats
            self.p_random.append(p)
            if p < 0.05: self.sig_region.append(rbs)
            self.mean_random.append(numpy.mean(counts_rbs))
        self.counts_random = numpy.array(matrix)
        self.sd_random = numpy.std(self.counts_random, axis=0)

    def plot(self, dir, ac, cut_off):
        """Generate the visualized plot"""
        
        # Extract data points
        x = range(self.rna_len)
        y = []
        rbss = BindingSiteSet("RBSs")
        rbss.sequences = self.rbss
        for i in range(self.rna_len):
            y.append(rbss.count_rbs_position(i))
        
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
        
        html.add_figure("randomtest.png", align="left")
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

        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"randomtest.html"))