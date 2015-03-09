# Python Libraries
from __future__ import print_function
from collections import *
import os
import sys
import multiprocessing
import time, datetime
# Local Libraries
from scipy import stats
import numpy
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator
from matplotlib import colors

# Distal Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from rgt.SequenceSet import Sequence, SequenceSet
from RNADNABindingSet import RNADNABinding, RNADNABindingSet
from rgt.Util import SequenceType, Html, OverlapType
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
        elif 1 > value >= 0.05: r = "{:.2f}".format(value)
        elif 0.05 > value > 0.0001: r = "{:.4f}".format(value)
        else: r = "{:.1e}".format(value)
        return r

def uniq(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def random_each(input):
    """Return the counts of DNA Binding sites with randomization
    For multiprocessing. 
    Input contains:
    0       1               2                3     4              5          6                    
    str(i), self.rna_fasta, self.dna_region, temp, self.organism, self.rbss, str(marks.count(i)),
    number, rna,            region,          temp, organism,      rbss,      number of mark

    7  8  9  10  11  12  13  14  15
    l, e, c, fr, fm, of, mf, rm, filter_bed 
    """
    # Filter BED file
    if input[15]:
        random = input[2].random_regions(organism=input[4], multiply_factor=1, 
                                         overlap_result=True, overlap_input=True, 
                                         chrom_X=True, chrom_M=False, filter_path=input[15])
    else:
        random = input[2].random_regions(organism=input[4], multiply_factor=1, 
                                         overlap_result=True, overlap_input=True, 
                                         chrom_X=True, chrom_M=False)
        
    txp = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[4], prefix=str(input[0]), remove_temp=True, 
                       l=int(input[7]), e=int(input[8]),  c=input[9], fr=input[10], 
                       fm=input[11], of=input[12], mf=input[13], rm=input[14],
                       dna_fine_posi=True)
    
    txp.merge_by(rbss=input[5], rm_duplicate=True)
    #print("\t Randomization: \t"+input[0])
    sys.stdout.flush()
    print("".join(["="]*int(input[6])), end="")

    return [ len(dbss) for dbss in txp.merged_dict.values() ]

def find_triplex(rna_fasta, dna_region, temp, organism, l, e, dna_fine_posi, prefix="", remove_temp=False, 
                 c=None, fr=None, fm=None, of=None, mf=None, rm=None):
    """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
    
    # Generate BED 
    dna_region.write_bed(os.path.join(temp,"dna_"+prefix+".bed"))
    # Bedtools
    os.system("bedtools getfasta -fi /data/genome/"+organism+"/"+organism+".fa -bed "+\
              os.path.join(temp,"dna_"+prefix+".bed")+" -fo "+os.path.join(temp,"dna_"+prefix+".fa"))
    # Triplexator
    run_triplexator(ss=rna_fasta, ds=os.path.join(temp,"dna_"+prefix+".fa"), 
                    output=os.path.join(temp, "dna_"+prefix+".txp"), 
                    l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf)
    # Read txp
    txp = RNADNABindingSet("dna")
    txp.read_txp(os.path.join(temp, "dna_"+prefix+".txp"), dna_fine_posi=dna_fine_posi)
    txp.remove_duplicates()
    #print(len(txp_de))

    if remove_temp:
        os.remove(os.path.join(temp,"dna_"+prefix+".bed"))
        os.remove(os.path.join(temp,"dna_"+prefix+".fa"))
        os.remove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp

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
        filelist = [ f for f in os.listdir(path) ]
        for f in filelist:
            try: os.remove(os.path.join(path,f))
            except: pass
    except:
        os.mkdir(path)

def split_gene_name(gene_name, ani, org):
    p1 = '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+ani+\
         "&db="+org+"&singleSearch=knownCanonical&position="
    p2 = '" style="text-align:left">'
    p3 = '</a>'
    
    if ":" in gene_name:
        genes = gene_name.split(":")
        genes = set(genes)
        for i, g in enumerate(genes):
            if i == 0:
                result = p1+g+p2+g+p3
            else:
                result += ","+p1+g+p2+g+p3
    elif gene_name == ".":
        result = "none"

    else:
        result = p1+gene_name+p2+gene_name+p3
    return result


def lineplot(txp, rnalen, rnaname, dir, sig_region, cut_off, log, ylabel, linelabel, 
             filename, ac=None, showpa=False):
    # Plotting
    f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
    

    # Extract data points
    x = range(rna_len)
    if log:
        all_y = [1] * 
        p_y = [1] * rna_len
        a_y = [1] * rna_len
    else:
        all_y = [0] * rna_len
        p_y = [0] * rna_len
        a_y = [0] * rna_len

    txp.remove_duplicates_by_dbs()
    for rd in txp:
        if rd.rna.orientation == "P":
            for i in range(rd.rna.initial, rd.rna.final):
                p_y[i] += 1
                all_y[i] += 1
        if rd.rna.orientation == "A":
            for i in range(rd.rna.initial, rd.rna.final):
                a_y[i] += 1
                all_y[i] += 1
    # Log
    if log:
        all_y = numpy.log(all_y)
        p_y = numpy.log(p_y)
        a_y = numpy.log(a_y)
        max_y = max(all_y)+0.5
        min_y = 1
        ylabel += "(log10)"
    else:
        max_y = float(max(all_y) * 1.1)
        min_y = 0

    if ac:
        min_y = float(max_y*(-0.09))
    
    
    # Plotting
    for rbs in sig_region:
        rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor="powderblue", 
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant region")
        ax.add_patch(rect)
    
    lw = 1.5
    if showpa:
        ax.plot(x, all_y, color="mediumblue", alpha=1, lw=lw, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="purple", alpha=1, lw=lw, label="Parallel")
        ax.plot(x, a_y, color="dimgrey", alpha=.8, lw=lw, label="Anti-parallel")
    else:
        ax.plot(x, all_y, color="mediumblue", alpha=1, lw=lw, label=linelabel)
    # RNA accessbility
    if ac:
        n_value = read_ac(ac, cut_off, rnalen=rna_len)
        drawing = False
        for i in x:
            if n_value[i] > 0:
                if drawing:
                    continue
                else:
                    last_i = i
                    drawing = True
            elif drawing:
                pac = ax.add_patch(patches.Rectangle((last_i, min_y), i-last_i, -min_y,
                                   hatch='///', fill=False, snap=False, linewidth=0, label="RNA accessibility"))
                drawing = False
            else:
                continue

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    legend_h = []
    legend_l = []
    for uniqlabel in uniq(labels):
        legend_h.append(handles[labels.index(uniqlabel)])
        legend_l.append(uniqlabel)
    ax.legend(legend_h, legend_l, 
              bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
              prop={'size':9}, ncol=3)

    # XY axis
    ax.set_xlim(left=0, right=self.rna_len )
    ax.set_ylim( [min_y, max_y] ) 
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    ax.set_xlabel(self.rna_name+" sequence (bp)", fontsize=9)
    
    ax.set_ylabel(ylabel,fontsize=9, rotation=90)

    f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
    f.savefig(os.path.join(dir, filename), facecolor='w', edgecolor='w',  
              bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
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
            self.de_regions = self.de_regions.gene_association(organism=self.organism)
            self.nde_regions.read_bed(bg)
            self.nde_regions = self.nde_regions.gene_association(organism=self.organism)
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
            
    def search_triplex(self, rna, temp, l, e, c, fr, fm, of, mf, remove_temp=False):
        # DE
        self.de_regions.write_bed(os.path.join(temp,"targeted_promoters.bed"))

        #os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
        #          "/home/joseph/Downloads/0122result/hotair_genes_promoters_up.bed"+" -fo "+os.path.join(temp,"de.fa"))
        os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
                  os.path.join(temp,"targeted_promoters.bed")+" -fo "+os.path.join(temp,"de.fa"))
        run_triplexator(ss=rna, ds=os.path.join(temp,"de.fa"), 
                        output=os.path.join(temp, "de.txp"), 
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf)
        
        # non-DE
        self.nde_regions.write_bed(os.path.join(temp,"untargeted_promoters.bed"))
        #os.system("bedtools getfasta -fi /data/genome/hg19/hg19.fa -bed "+\
        #          "/home/joseph/Downloads/0122result/hotair_promoters_bg_up.bed"+" -fo "+os.path.join(temp,"nde.fa"))
        os.system("bedtools getfasta -fi /data/genome/"+self.organism+"/"+self.organism+".fa -bed "+\
                  os.path.join(temp,"untargeted_promoters.bed")+" -fo "+os.path.join(temp,"nde.fa"))
        run_triplexator(ss=rna, ds=os.path.join(temp,"nde.fa"), 
                        output=os.path.join(temp, "nde.txp"), 
                        l=l, e=e, c=2, fr="off", fm=0, of=1, mf=True)
        if remove_temp:
            os.remove(os.path.join(temp,"targeted_promoters.bed"))
            os.remove(os.path.join(temp,"de.fa"))
            os.remove(os.path.join(temp,"untargeted_promoters.bed"))
            os.remove(os.path.join(temp,"nde.fa"))
        
    def count_frequency(self, temp, remove_temp, cutoff, bed_output=False):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""
        
        self.frequency = {}
        self.frequency["promoters"] = { "de": OrderedDict(), "nde": OrderedDict() }
        self.frequency["hits"] = { "de": OrderedDict(), "nde": OrderedDict() }
        
        ########################################################
        # Count the number of targeted promoters on each merged DNA Binding Domain
        txp_de = RNADNABindingSet("DE")
        txp_nde = RNADNABindingSet("non-DE")
        
        txp_de.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=False)
        txp_de.remove_duplicates()
        txp_de.merge_rbs(rm_duplicate=True, asgene_organism=self.organism, cutoff=cutoff)
        
        self.rbss = txp_de.merged_dict.keys()

        txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        txp_nde.remove_duplicates()
        txp_nde.merge_by(rbss=self.rbss, rm_duplicate=True, asgene_organism=self.organism)

        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)

        for rbs in self.rbss:
            # DE
            self.frequency["promoters"]["de"][rbs] = [ len(txp_de.merged_dict[rbs]), 
                                                       len_de - len(txp_de.merged_dict[rbs]) ]
            # non-DE
            self.frequency["promoters"]["nde"][rbs] = [ len(txp_nde.merged_dict[rbs]), 
                                                        len_nde - len(txp_nde.merged_dict[rbs]) ]

        self.txp_de = txp_de
        self.txp_nde = txp_nde

        ########################################################
        # Count the number of hits on the promoters from each merged DBD 
        txp_def = RNADNABindingSet("DE")
        txp_def.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=True)

        txp_def.merge_by(rbss=self.rbss, rm_duplicate=True, asgene_organism=self.organism)

        txp_ndef = RNADNABindingSet("non-DE")
        txp_ndef.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=True)

        txp_ndef.merge_by(rbss=self.rbss, rm_duplicate=True, asgene_organism=self.organism)

        for rbs in self.rbss:
            #if len(txp_def.merged_dict[rbs]) < 10: continue
            # DE
            self.frequency["hits"]["de"][rbs] = len(txp_def.merged_dict[rbs])
            # non-DE
            self.frequency["hits"]["nde"][rbs] = len(txp_ndef.merged_dict[rbs])

        self.txp_def = txp_def
        self.txp_ndef = txp_ndef
        
        if remove_temp:
            os.remove(os.path.join(temp,"de.txp"))
            os.remove(os.path.join(temp,"nde.txp"))

        if bed_output:
            self.txp_de.write_bed(filename=os.path.join(temp,"targeted_binding_promoters.bed"), remove_duplicates=True)
            self.txp_nde.write_bed(filename=os.path.join(temp,"untargeted_binding_promoters.bed"), remove_duplicates=True)
            self.txp_def.write_bed(filename=os.path.join(temp,"dna_binding_sites_tar_pro.bed"), remove_duplicates=True)
            self.txp_ndef.write_bed(filename=os.path.join(temp,"dna_binding_sites_untar_pro.bed"), remove_duplicates=True)

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

    def plot_lines(self, txp, rna, dirp, cut_off, log, ylabel, linelabel, filename, ac=None, showpa=False):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(rna)
        if not self.rna_name: self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])
        
        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dir=dirp, sig_region=self.sig_region, 
                 cut_off=cut_off, log=log, ylabel=ylabel, linelabel=linelabel,  
                 filename=filename, ac=ac, showpa=showpa)

    def promoter_profile(self):
        """count the number of DBD and DBS regarding to each promoter"""
        self.promoter = { "de": {},
                          "nde": {} }

        #################################################
        # Each promoter to all related merged DBD
        print("\tEach promoter to all related merged DBD...")
        # Targeted promoters
        self.promoter["de"]["rd"] = self.txp_def.sort_rd_by_regions(regionset=self.de_regions)
        self.promoter["nde"]["rd"] = self.txp_ndef.sort_rd_by_regions(regionset=self.nde_regions)
        #################################################
        
        # Each promoter to all related DBS
        print("\tEach promoter to all related DBS...")
        # Targeted promoters
        
        self.promoter["de"]["merged_dbs"] = {}
        self.promoter["de"]["dbs"] = {}
        for k, rds in self.promoter["de"]["rd"].items():
            self.promoter["de"]["dbs"][k] = rds.get_dbs()
            self.promoter["de"]["merged_dbs"][k] = self.promoter["de"]["dbs"][k].merge(w_return=True) 
        # Untargeted promoters
        self.promoter["nde"]["merged_dbs"] = {}
        self.promoter["nde"]["dbs"] = {}
        for k, rds in self.promoter["nde"]["rd"].items():
            self.promoter["nde"]["dbs"][k] = rds.get_dbs()
            self.promoter["nde"]["merged_dbs"][k] = self.promoter["nde"]["dbs"][k].merge(w_return=True)
        
        # DBS coverage
        self.promoter["de"]["dbs_coverage"] = {}
        for promoter in self.de_regions:
            self.promoter["de"]["dbs_coverage"][promoter.toString()] = float(self.promoter["de"]["merged_dbs"][promoter.toString()].total_coverage()) / len(promoter)
        self.promoter["nde"]["dbs_coverage"] = {}
        for promoter in self.nde_regions:
            self.promoter["nde"]["dbs_coverage"][promoter.toString()] = float(self.promoter["nde"]["merged_dbs"][promoter.toString()].total_coverage()) / len(promoter)

    def gen_html(self, directory, align=50, alpha = 0.05):

        check_dir(os.path.join(directory, "supplements"))
        html_header = "Triplex Domain Finder: Promoter Test"
        #fp = os.path.join(dir,outputname,title)
        self.link_d = OrderedDict()
        self.link_d["RNA"] = "index.html"
        self.link_d["Target promoters"] = "promoters.html"
        self.link_d["Non-target promoters"] = "background.html"
        self.link_ds = OrderedDict()
        self.link_ds["RNA"] = "../index.html"
        self.link_ds["Target promoters"] = "../promoters.html"
        self.link_ds["Non-target promoters"] = "../background.html"

        if self.organism == "hg19": self.ani = "human"
        elif self.organism == "mm9": self.ani = "mouse"

        #############################################################
        # Index main page
        #############################################################
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        
        html.add_figure("plot_promoter.png", align="left")
        html.add_figure("plot_dbss.png", align="left")
        
        # Table of merged TBS on promoters
        header_list = [ ["DBD", "Target Promoter", None, "Non-target Promoter", None, "Statistics", None],
                        [" ", "with DBS", "without DNS", "with DBS", "without DBS", "Oddsratio", "P-value"] ]
                       #"Unique DBSs<br>Targeted promoters",
                       #"Unique DBSs<br>Untargeted promoters" ]
        
        type_list = 'sssssssss'
        col_size_list = [50,50,50,50,50,10,10,10,10]
        data_table = []
        for rbs in self.frequency["promoters"]["de"].keys():
            if self.frequency["promoters"]["de"][rbs][0] < 10: continue
            if self.pvalue[rbs] < alpha:
                data_table.append([ rbs.str_rna(pa=False), #str(len(dbss)), 
                                    '<a href="'+"supplements/dbds_promoters.html#"+rbs.str_rna()+
                                    '" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["de"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["de"][rbs][1]), 
                                    '<a href="'+"supplements/dbds_background.html#"+rbs.str_rna()+
                                    '" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["nde"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["nde"][rbs][1]), 
                                    value2str(self.oddsratio[rbs]), 
                                    "<font color=\"red\">"+value2str(self.pvalue[rbs])+"</font>" ])
                                    #'<a href="'+"supplements/dbds_dbss.html#"+rbs.str_rna()+
                                    #'" style="text-align:left">'+
                                    #value2str(self.frequency["hits"]["de"][rbs])+'</a>', 
                                    #'<a href="'+"supplements/dbds_bg_dbss.html#"+rbs.str_rna()+
                                    #'" style="text-align:left">'+
                                    #value2str(self.frequency["hits"]["nde"][rbs])+'</a>' ])
            else:
                data_table.append([ rbs.str_rna(pa=False), #str(len(dbss)), 
                                    '<a href="'+"supplements/dbds_promoters.html#"+rbs.str_rna()+
                                    '" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["de"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["de"][rbs][1]), 
                                    '<a href="'+"supplements/dbds_background.html#"+rbs.str_rna()+
                                    '" style="text-align:left">'+
                                    value2str(self.frequency["promoters"]["nde"][rbs][0])+'</a>', 
                                    value2str(self.frequency["promoters"]["nde"][rbs][1]), 
                                    value2str(self.oddsratio[rbs]), value2str(self.pvalue[rbs]) ])
                                    #'<a href="'+"supplements/dbds_dbss.html#"+rbs.str_rna()+
                                    #'" style="text-align:left">'+
                                    #value2str(self.frequency["hits"]["de"][rbs])+'</a>', 
                                    #'<a href="'+"supplements/dbds_bg_dbss.html#"+rbs.str_rna()+
                                    #'" style="text-align:left">'+
                                    #value2str(self.frequency["hits"]["nde"][rbs])+'</a>' ])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        
        html.add_heading("Notes")
        html.add_list([ "DBD stands for functional DNA Binding Domain on RNA.",
                        "RBS stands for RNA Binding Site on RNA.",
                        "DBS stands for DNA Binding Site on DNA."])

        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory,"index.html"))

        #############################################################
        # RNA subpage: Profile of targeted promoters for each merged DNA Binding Domain
        #############################################################
        
        header_list= ["Target Promoter", "Strand", "Gene", "No. DBSs", "DBS coverage"] 

        # dbds_promoters.html
        html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False)
     
        for rbsm in self.frequency["promoters"]["de"].keys():    
            html.add_heading("DNA Binding Domain: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for promoter in self.txp_de.merged_dict[rbsm]:
                # Add information
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString(space=True)+'</a>', 
                                    promoter.orientation,
                                    split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                                    str(len(self.promoter["de"]["rd"][promoter.toString()])),
                                    value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                                     ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.write(os.path.join(directory,"supplements", "dbds_promoters.html"))

        # dbds_background.html
        html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False)
      
        for rbsm in self.frequency["promoters"]["nde"].keys():
            
            html.add_heading("DNA Binding Domain: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for promoter in self.txp_nde.merged_dict[rbsm]:
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString(space=True)+'</a>', 
                                    promoter.orientation,
                                    split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                                    str(len(self.promoter["nde"]["rd"][promoter.toString()])),
                                    value2str(self.promoter["nde"]["dbs_coverage"][promoter.toString()])
                                     ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.write(os.path.join(directory,"supplements", "dbds_background.html"))
        
        #############################################################
        # Profile of Hits for each merged DNA Binding Domain
        #############################################################
        if None:
            header_list=["DNA_Binding_Site", "Gene", "Strand", "Width", "Score", "Motif", "Orientation"]
            # dbds_dbss.html
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            
            for rbsm in self.frequency["hits"]["de"].keys():
                
                html.add_heading("DBSs in the DBD: "+rbsm.str_rna(),
                                 idtag=rbsm.str_rna())
                data_table = []
                for dbs in self.txp_def.merged_dict[rbsm]:
                    data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+dbs.chrom+"%3A"+str(dbs.initial)+"-"+str(dbs.final)+
                                        '" style="text-align:left">'+dbs.toString(space=True)+'</a>',
                                        split_gene_name(gene_name=dbs.name, ani=self.ani, org=self.organism),
                                        dbs.orientation, str(dbs.final-dbs.initial),
                                        dbs.data[0], dbs.data[1], dbs.data[2]
                                        ])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "dbds_dbss.html"))

        
            # dbds_bg_dbss.html
            html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                        fig_rpath="../style", RGT_header=False)
            for rbsm in self.frequency["hits"]["nde"].keys():
                html.add_heading("DBSs in the DBD: "+rbsm.str_rna(),
                                 idtag=rbsm.str_rna())
                data_table = []
                for dbs in self.txp_ndef.merged_dict[rbsm]:
                    data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+dbs.chrom+"%3A"+str(dbs.initial)+"-"+str(dbs.final)+
                                        '" style="text-align:left">'+dbs.toString(space=True)+'</a>',
                                        split_gene_name(gene_name=dbs.name, ani=self.ani, org=self.organism),
                                        dbs.orientation, str(dbs.final-dbs.initial),
                                        dbs.data[0], dbs.data[1], dbs.data[2]
                                       ])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            html.write(os.path.join(directory,"supplements", "dbds_bg_dbss.html"))

    def gen_html_genes(self, directory, align=50, alpha = 0.05, nonDE=False):

        html_header = "Triplex Domain Finder: Promoter Test"
        #fp = os.path.join(dir,outputname,title)
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]

        #############################################################
        # Promoter centered
        #############################################################
        
        ##############################################################################################
        # promoters.html
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        header_promoter_centered=[ "Promoter", "Gene", "No. DBSs", "DBS coverage" ]
        header_list = header_promoter_centered
        html.add_heading("Target promoters")
        data_table = []
        #counts_p = self.de_regions.counts_per_region(self.txp_de.get_dbs(sort=True, orientation="P"))
        #counts_a = self.de_regions.counts_per_region(self.txp_de.get_dbs(sort=True, orientation="A"))
        
        if not self.de_regions.sorted: self.de_regions.sort()
        # Iterate by each gene promoter
        for i, promoter in enumerate(self.de_regions):
            if len(self.promoter["de"]["dbs"][promoter.toString()]) == 0:
            #counts_p[i] == 0 and counts_a[i] == 0:
                continue
            else:
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString(space=True)+'</a>',
                                    split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                                    '<a href="supplements/promoters_dbds.html#'+promoter.toString()+
                                    '" style="text-align:left">'+
                                    str(len(self.promoter["de"]["dbs"][promoter.toString()]))+'</a>',
                                    value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                                    ])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_heading("Notes")
        html.add_list(["All the prmoters without any bindings are ignored.",
                       "DBS stands for DNA Binding Site on DNA.",
                       "DBS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
        
        html.write(os.path.join(directory,"promoters.html"))

        ############################
        # Subpages for promoter centered page
        # promoters_dbds.html

        header_sub = ["RBS", "DBS", "Strand", "Score", "Motif", "Orientation"]
        header_list=header_sub
        html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False)
        
        for i, promoter in enumerate(self.de_regions):
            if len(self.promoter["de"]["dbs"][promoter.toString()]) == 0:
                continue
            else:         
                html.add_heading(split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism), idtag=promoter.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                        '" style="margin-left:50">'+
                                        promoter.toString(space=True)+'</a>'])
                data_table = []
                
                for rd in self.promoter["de"]["rd"][promoter.toString()]:
                    data_table.append([ rd.rna.str_rna(pa=False),
                                        rd.dna.toString(space=True),
                                        rd.dna.orientation,
                                        rd.score, 
                                        rd.motif, rd.orient])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
                
        html.write(os.path.join(directory,"supplements", "promoters_dbds.html"))
 
        # background.html
        html = Html(name=html_header, links_dict=self.link_d, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        header_list=header_promoter_centered
        html.add_heading("Non-target promoters")
        data_table = []
        #counts_p = self.nde_regions.counts_per_region(self.txp_nde.get_dbs(sort=True, orientation="P"))
        #counts_a = self.nde_regions.counts_per_region(self.txp_nde.get_dbs(sort=True, orientation="A"))
        
        if not self.nde_regions.sorted: self.nde_regions.sort()
        # Iterate by each gene promoter
        for i, promoter in enumerate(self.nde_regions):
            if len(self.promoter["nde"]["dbs"][promoter.toString()]) == 0:
                continue
            else:
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                    '" style="text-align:left">'+promoter.toString(space=True)+'</a>',
                                    split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                                    '<a href="supplements/background_dbds.html#'+promoter.toString()+
                                    '" style="text-align:left">'+
                                    str(len(self.promoter["nde"]["dbs"][promoter.toString()]))+'</a>',
                                    value2str(self.promoter["nde"]["dbs_coverage"][promoter.toString()])
                                    ])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_heading("Notes")
        html.add_list(["All the prmoters without any bindings are ignored.",
                       "DBS stands for DNA Binding Site on DNA.",
                       "DBS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
        
        html.write(os.path.join(directory,"background.html"))

        ############################
        # background_dbds.html
        html = Html(name=html_header, links_dict=self.link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False)
       
        header_list=header_sub
        for i, promoter in enumerate(self.nde_regions):
            if len(self.promoter["nde"]["dbs"][promoter.toString()]) == 0:
                continue
            else: 
                html.add_heading(split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism), idtag=promoter.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+
                                        '" style="margin-left:50">'+
                                        promoter.toString(space=True)+'</a>'])
                data_table = []
                for rd in self.promoter["nde"]["rd"][promoter.toString()]:
                    data_table.append([ rd.rna.str_rna(pa=False),
                                        rd.dna.toString(space=True),
                                        rd.dna.orientation,
                                        rd.score, 
                                        rd.motif, rd.orient])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
            
        html.write(os.path.join(directory,"supplements", "background_dbds.html"))


####################################################################################
####################################################################################

class RandomTest:
    def __init__(self, rna_fasta, rna_name, dna_region, organism):
    	self.organism = organism
        # RNA: Path to the FASTA file
        self.rna_fasta = rna_fasta
        

        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(self.rna_fasta)
        if rna_name: 
            self.rna_name = rna_name
        else:
            self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])
        # DNA: GenomicRegionSet
        self.dna_region = GenomicRegionSet(name="target")
        self.dna_region.read_bed(dna_region)
        self.dna_region = self.dna_region.gene_association(organism=self.organism)
        
    def target_dna(self, temp, remove_temp, l, e, c, fr, fm, of, mf, obed=False, obedname=None):
        """Calculate the true counts of triplexes on the given dna regions"""
        
        txp = find_triplex(rna_fasta=self.rna_fasta, dna_region=self.dna_region, 
                           temp=temp, organism=self.organism, remove_temp=remove_temp, 
                           l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf,
                           prefix="targeted_region", dna_fine_posi=False)
        txp.merge_rbs(rm_duplicate=False, asgene_organism=self.organism)
        self.txp = txp
        txp.remove_duplicates()
        self.rbss = txp.merged_dict.keys()

        txpf = find_triplex(rna_fasta=self.rna_fasta, dna_region=self.dna_region, 
                            temp=temp, organism=self.organism, remove_temp=remove_temp, 
                            l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf,
                            prefix="dbs", dna_fine_posi=True)
        txpf.remove_duplicates()
        txpf.merge_by(rbss=self.rbss, rm_duplicate=False, asgene_organism=self.organism)
        self.txpf = txpf
        
        self.counts_tr = OrderedDict()
        self.counts_dbs = OrderedDict()
        self.counts_dbsm = OrderedDict()
        self.mdbs = OrderedDict()

        for rbs in self.rbss:
            self.mdbs[rbs] = txpf.merged_dict[rbs].merge(w_return=True)
            tr = len(self.txp.merged_dict[rbs])
            self.counts_tr[rbs] = [tr, len(self.dna_region) - tr]
            self.counts_dbs[rbs] = len(self.txpf.merged_dict[rbs])
            self.counts_dbsm[rbs] = len(self.mdbs[rbs])
            
        self.region_dbd = self.txpf.sort_rbs_by_regions(self.dna_region)
        
        
        self.region_dbs = self.txpf.sort_rd_by_regions(regionset=self.dna_region)
        self.region_dbsm = {}
        dbss = self.txpf.get_dbs()
        for region in self.dna_region:
            self.region_dbsm[region.toString()] = self.region_dbs[region.toString()].get_dbs().merge(w_return=True)

        if obed:
            btr = txp.get_dbs()
            btr = dbss.gene_association(organism=self.organism)
            btr.write_bed(os.path.join(temp, "binding_region_"+obedname+".bed"))
            dbss = txpf.get_dbs()
            dbss = dbss.gene_association(organism=self.organism)
            dbss.write_bed(os.path.join(temp, "dbs_"+obedname+".bed"))

    def random_test(self, repeats, temp, remove_temp, l, e, c, fr, fm, of, mf, rm, filter_bed):
        """Perform randomization for the given times"""
        self.repeats = repeats
        marks = numpy.round(numpy.linspace(0, repeats-1, num=41)).tolist()
        
        # Prepare the input lists for multiprocessing
        
        mp_input = []
        for i in range(repeats):
            mp_input.append([ str(i), self.rna_fasta, self.dna_region,
                              temp, self.organism, self.rbss, str(marks.count(i)),
                              str(l), str(e), str(c), str(fr), str(fm), str(of), str(mf), str(rm),
                              filter_bed ])
        # Multiprocessing
        print("\t\t|0%                  |                100%|")
        print("\t\t[", end="")
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
        self.sig_boolean = []
        for i, rbs in enumerate(self.rbss):
            #print(i)
            counts_rbs = [ v[i] for v in mp_output ]
            #print(len(counts_rbs))
            matrix.append(counts_rbs)
            p = float(len([ h for h in counts_rbs if h > self.counts_dbs[rbs] ]))/repeats
            self.p_random.append(p)
            if p < 0.05: 
                self.sig_region.append(rbs)
                self.sig_boolean.append(True)
            else:
                self.sig_boolean.append(False)
            
            self.mean_random.append(numpy.mean(counts_rbs))
        self.counts_random = numpy.array(matrix)
        self.sd_random = numpy.std(self.counts_random, axis=1)

    def lineplot(self, txp, dirp, ac, cut_off, log, ylabel, linelabel, showpa, filename):
        """Generate lineplot for RNA"""
        
        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dir=dirp, sig_region=self.sig_region, 
                 cut_off=cut_off, log=log, ylabel=ylabel, linelabel=linelabel,  
                 filename=filename, ac=ac, showpa=showpa)

    def boxplot(self, dir):
        """Generate the visualized plot"""
        tick_size = 8
        label_size = 9
        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        max_y = int(self.counts_random.max()*1.1) + 1
        min_y = max(int(self.counts_random.min()*0.9) - 1, 0)
        bp = ax.boxplot(self.counts_random.transpose(), notch=False, sym='o', vert=True, 
                        whis=1.5, positions=None, widths=None, 
                        patch_artist=True, bootstrap=None)
        z = 10
        plt.setp(bp['boxes'], color='lightblue')
        plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8,zorder=z)
        plt.setp(bp['fliers'], markerfacecolor='gray',color='white',alpha=0.3,markersize=1.8,zorder=z)
        plt.setp(bp['caps'],color='white',zorder=-1)
        plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
        
        # Plot target regions
        ax.plot(range(1, len(self.rbss)+1), self.counts_dbs.values(), 'ro', markersize=4, zorder=z+2)
        ax.set_xlabel("Potential DNA Binding Domains", fontsize=label_size)
        ax.set_ylabel("Number of DNA Binding Sites",fontsize=label_size, rotation=90)
        

        for i, r in enumerate(self.sig_boolean):
            if r:
                rect = patches.Rectangle(xy=(i+0.6,min_y), width=0.8, height=max_y, facecolor="lightgray", 
                                         edgecolor="none", alpha=0.5, lw=None)
                ax.add_patch(rect)

        ax.set_ylim( [min_y, max_y] ) 
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        ax.set_xticklabels( [dbd.str_rna(pa=False) for dbd in self.rbss], rotation=35, 
                            ha="right", fontsize=tick_size)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(tick_size) 

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, "boxplot.png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)

    def gen_html(self, directory, align=50, alpha = 0.05):
        """Generate the HTML file"""

        html_header = "Triplex Domain Finder: Region Test"
        link_ds = {"RNA":"index.html",
                   "Target regions":"target_regions.html"}

        if self.organism == "hg19": self.ani = "human"
        elif self.organism == "mm9": self.ani = "mouse"

        ##################################################
        # index.html

        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        # Plots
        html.add_figure("lineplot_region.png", align="left")
        html.add_figure("lineplot_dbs.png", align="left")
        html.add_figure("boxplot.png", align="left")

        header_list = [ "DNA Binding Domain",
                        "Binding<br>Targeted regions", 
                        "No binding<br>Targeted regions",
                        "Randomization Average DBS",
                        "Randomization s.d.",
                        "P value",
                        "DBSs",
                        "Binding target regions" ]
        
        type_list = 'sssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10]
        data_table = []
        
        for i, rbs in enumerate(self.rbss):
            if self.p_random[i] < alpha:
                data_table.append([ rbs.str_rna(pa=False),
                                    '<a href="dbd_region.html#'+rbs.str_rna()+
                                    '" style="text-align:left">'+str(self.counts_tr[rbs][0])+'</a>',
                                    str(self.counts_tr[rbs][1]),
                                    value2str(self.mean_random[i]), 
                                    value2str(self.sd_random[i]), 
                                    "<font color=\"red\">"+value2str(self.p_random[i])+"</font>",
                                    '<a href="dbd_dbs.html#'+rbs.str_rna()+
                                    '" style="text-align:left">'+str(self.counts_dbs[rbs])+'</a>', 
                                    str(self.counts_dbsm[rbs]) ])
            else:
                data_table.append([ rbs.str_rna(pa=False),
                                    '<a href="dbd_region.html#'+rbs.str_rna()+
                                    '" style="text-align:left">'+str(self.counts_tr[rbs][0])+'</a>',
                                    str(self.counts_tr[rbs][1]),
                                    value2str(self.mean_random[i]), 
                                    value2str(self.sd_random[i]), 
                                    value2str(self.p_random[i]),
                                    '<a href="dbd_dbs.html#'+rbs.str_rna()+
                                    '" style="text-align:left">'+str(self.counts_dbs[rbs])+'</a>', 
                                    str(self.counts_dbsm[rbs]) ])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")


        html.add_heading("Notes")
        html.add_list([ "RNA name: "+ self.rna_name,
                        "Randomization is performed for "+str(self.repeats)+" times.",
                        "DBD stands for DNA Binding Domain on RNA.",
                        "DBS stands for DNA Binding Site on DNA."])

        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See summary</a>'])
        html.write(os.path.join(directory,"index.html"))


        #############################################################
        # RNA subpage: Profile of targeted regions for each merged DNA Binding Domain
        #############################################################
        
        header_list=[ "Targeted region", 
                      "Associated Gene", 
                      "Strand", 
                      "Binding of DBS" ]

        #########################################################
        # dbd_region.html
        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
     
        for rbsm in self.rbss:    
            html.add_heading("Targeted regions relate to DBD: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for region in self.txp.merged_dict[rbsm]:
                # Add information
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                    '" style="text-align:left">'+region.toString(space=True)+'</a>',
                                    split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism),
                                    region.orientation,
                                    str(len(self.region_dbs[region.toString()])) ])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.write(os.path.join(directory, "dbd_region.html"))


        #########################################################
        # dbd_dbs.html
        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        
        header_list=[ "DBS", 
                      "Associated Gene", 
                      "Strand",
                      "Score" ]

        for rbsm in self.rbss:    
            html.add_heading("Targeted regions relate to DBD: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for dbs in self.txpf.merged_dict[rbsm]:
                # Add information
                
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+dbs.chrom+"%3A"+str(dbs.initial)+"-"+str(dbs.final)+
                                    '" style="text-align:left">'+dbs.toString(space=True)+'</a>', 
                                    split_gene_name(gene_name=dbs.name, ani=self.ani, org=self.organism),
                                    dbs.orientation,
                                    dbs.data[0] ])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.write(os.path.join(directory, "dbd_dbs.html"))


        #############################################################
        # Targeted regions centered
        #############################################################
        
        ##############################################################################################
        # target_regions.html
        html = Html(name=html_header, links_dict=link_ds,fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        header_list=[ "Targeted region", 
                      "Associated Gene", 
                      "RBS",
                      "P-RBS", 
                      "A-RBS",
                      "DBS",
                      "Merged DBSs" ]

        html.add_heading("The targeted promoters")
        data_table = []
        
        #p = self.dna_region.counts_per_region(self.txp.get_dbs(sort=True, orientation="P"))
        #a = self.dna_region.counts_per_region(self.txp.get_dbs(sort=True, orientation="A"))

        if not self.dna_region.sorted: self.dna_region.sort()
        # Iterate by each gene promoter
        for i, region in enumerate(self.dna_region):
            
            if len(self.region_dbs[region.toString()]) == 0:
                continue
            else:
                p = len(self.region_dbd[region.toString()].get_bs(orientation="P"))
                a = len(self.region_dbd[region.toString()].get_bs(orientation="A"))
                data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                    '" style="text-align:left">'+region.toString(space=True)+'</a>',
                                    split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism),
                                    '<a href="region_dbd.html#'+region.toString()+
                                    '" style="text-align:left">'+str(p+a)+'</a>',
                                    str(p), str(a),
                                    '<a href="region_dbs.html#'+region.toString()+
                                    '" style="text-align:left">'+str(len(self.region_dbs[region.toString()]))+'</a>',
                                    '<a href="region_dbsm.html#'+region.toString()+
                                    '" style="text-align:left">'+str(len(self.region_dbsm[region.toString()]))+'</a>' ])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_heading("Notes")
        html.add_list(["All the genes without any bindings are ignored.",
                       "P-DBDs stands for Parallel DNA Binding Domains", 
                       "A-DBDs stands for Anti-parallel DNA Binding Domains"])
        
        html.write(os.path.join(directory,"target_regions.html"))

        ############################
        # Subpages for targeted region centered page
        # region_dbd.html
        
        header_list = ["RBS", "Motif", "Orientation"]

        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        html.add_heading("Notes")
        html.add_list(["Here we neglect the difference of DBSs, only the difference of DBDs are counted."])

        for i, region in enumerate(self.dna_region):
            if len(self.region_dbs[region.toString()]) == 0:
                continue
            else:         
                html.add_heading("Associated gene: "+split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism), idtag=region.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                        '" style="margin-left:50">'+
                                        region.toString(space=True)+'</a>'])
                data_table = []
                for dbd in self.region_dbd[region.toString()]:
                    data_table.append([ dbd.str_rna(pa=False),
                                        dbd.motif, dbd.orientation ])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
                
        html.write(os.path.join(directory, "region_dbd.html"))

        ############################
        # Subpages for targeted region centered page
        # region_dbs.html
        
        header_list = ["RBS", "DBS", "Strand", "Score", "Motif", "Orientation" ]

        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        
        for i, region in enumerate(self.dna_region):
            if len(self.region_dbs[region.toString()]) == 0:
                continue
            else:         
                html.add_heading("Associated gene: "+split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism), idtag=region.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                        '" style="margin-left:50">'+
                                        region.toString(space=True)+'</a>'])
                data_table = []
                for rd in self.region_dbs[region.toString()]:
                    data_table.append([ rd.rna.str_rna(pa=False),
                                        '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+rd.dna.chrom+"%3A"+str(rd.dna.initial)+"-"+str(rd.dna.final)+
                                        '" style="text-align:left">'+rd.dna.toString(space=True)+'</a>',
                                        rd.dna.orientation,
                                        rd.score,
                                        rd.motif, 
                                        rd.orient ])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
                
        html.write(os.path.join(directory, "region_dbs.html"))


        ############################
        # Subpages for targeted region centered page
        # region_dbsm.html
        
        header_list = ["Target region", "Strand", "Motif", "Orientation" ]

        html = Html(name=html_header, links_dict=link_ds, fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="./style", RGT_header=False)
        
        for i, region in enumerate(self.dna_region):
            if len(self.region_dbsm[region.toString()]) == 0:
                continue
            else:         
                html.add_heading("Associated gene: "+split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism), idtag=region.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                        '" style="margin-left:50">'+
                                        region.toString(space=True)+'</a>'])
                data_table = []
                for dbsm in self.region_dbsm[region.toString()]:
                    data_table.append([ '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                        "&position="+dbsm.chrom+"%3A"+str(dbsm.initial)+"-"+str(dbsm.final)+
                                        '" style="text-align:left">'+dbsm.toString(space=True)+'</a>',
                                        dbsm.orientation,
                                        dbsm.data[1], 
                                        dbsm.data[2] ])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
                
        html.write(os.path.join(directory, "region_dbsm.html"))

