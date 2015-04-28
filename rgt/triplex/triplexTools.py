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
import matplotlib
matplotlib.use('Agg')
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator, FuncFormatter
from matplotlib import colors
from matplotlib.backends.backend_pdf import PdfPages
import pysam
import pickle

# Distal Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from rgt.SequenceSet import Sequence, SequenceSet
from RNADNABindingSet import RNADNABinding, RNADNABindingSet
from rgt.Util import SequenceType, Html, OverlapType, ConfigurationFile, GenomeData
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.AnnotationSet import AnnotationSet

# Color code for all analysis
target_color = "mediumblue"
nontarget_color = "darkgrey"
sig_color = "powderblue"

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
        if abs(value) >= 1000: 
            try: r = "{}".format(int(value))
            except: r = "Inf"
        elif 1000 > abs(value) > 10: r = "{:.1f}".format(value)
        elif 10 > abs(value) >= 1: r = "{:.2f}".format(value)
        elif 1 > abs(value) >= 0.05: r = "{:.2f}".format(value)
        elif 0.05 > abs(value) > 0.0001: r = "{:.4f}".format(value)
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
                       fm=input[11], of=input[12], mf=input[13], rm=input[14], genome_path=input[16],
                       dna_fine_posi=False)
    
    txp.merge_by(rbss=input[5], rm_duplicate=True)

    txpf = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[4], prefix=str(input[0]), remove_temp=True, 
                       l=int(input[7]), e=int(input[8]),  c=input[9], fr=input[10], 
                       fm=input[11], of=input[12], mf=input[13], rm=input[14], genome_path=input[16],
                       dna_fine_posi=True)
    
    txpf.merge_by(rbss=input[5], rm_duplicate=True)
    sys.stdout.flush()
    print("".join(["="]*int(input[6])), end="")

    return [ [len(tr) for tr in txp.merged_dict.values() ], [len(dbss) for dbss in txpf.merged_dict.values()] ]

def get_sequence(dir, filename, regions, genome_path):
    """
    Fetch sequence into FASTA file according to the given BED file
    """
    genome = pysam.Fastafile(genome_path)
    with open(os.path.join(dir, filename), 'w') as output:
        for region in regions:
            print(">"+ region.toString(), file=output)
            print(genome.fetch(region.chrom, max(0, region.initial), region.final), file=output)
                    

def find_triplex(rna_fasta, dna_region, temp, organism, l, e, dna_fine_posi, genome_path, prefix="", remove_temp=True, 
                 c=None, fr=None, fm=None, of=None, mf=None, rm=None):
    """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
    
    # Generate BED 
    get_sequence(dir=temp, filename="dna_"+prefix+".fa", regions=dna_region, genome_path=genome_path)

    # Triplexator
    run_triplexator(ss=rna_fasta, ds=os.path.join(temp,"dna_"+prefix+".fa"), 
                    output=os.path.join(temp, "dna_"+prefix+".txp"), 
                    l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf)
    # Read txp
    txp = RNADNABindingSet("dna")
    txp.read_txp(os.path.join(temp, "dna_"+prefix+".txp"), dna_fine_posi=dna_fine_posi)
    txp.remove_duplicates()

    if remove_temp:
        os.remove(os.path.join(temp,"dna_"+prefix+".fa"))
        os.remove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp

def run_triplexator(ss, ds, output, l=None, e=None, c=None, fr=None, fm=None, of=None, mf=None, rm=None):
    """Perform Triplexator"""
    #triplexator_path = check_triplexator_path()

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
    #os.system(triplexator_path+arguments)
    os.system("triplexator"+arguments)

def read_ac(path, cut_off, rnalen):
    """Read the RNA accessibility file and output its positions and values

    The file should be a simple table with two columns:
    The first column is the position and the second one is the value
    '#' will be skipped

    """
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

    if not ani:
        if org == "tair10":
           p1 = '<a href="https://www.arabidopsis.org/servlets/TairObject?name='+ gene_name + '&type=locus"></a>'
           return p1
        else:
            return gene_name
    else:    
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


def lineplot(txp, rnalen, rnaname, dirp, sig_region, cut_off, log, ylabel, linelabel, 
             filename, ac=None, showpa=False):
    # Plotting
    f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
    
    # Extract data points
    x = range(rnalen)
    if log:
        all_y = [1] * rnalen
        p_y = [1] * rnalen
        a_y = [1] * rnalen
    else:
        all_y = [0] * rnalen
        p_y = [0] * rnalen
        a_y = [0] * rnalen

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
        rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor=sig_color, 
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        ax.add_patch(rect)
    
    lw = 1.5
    if showpa:
        ax.plot(x, all_y, color=target_color, alpha=1, lw=lw, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="purple", alpha=1, lw=lw, label="Parallel")
        ax.plot(x, a_y, color="dimgrey", alpha=.8, lw=lw, label="Anti-parallel")
    else:
        ax.plot(x, all_y, color="mediumblue", alpha=1, lw=lw, label=linelabel)

    # RNA accessbility
    if ac:
        n_value = read_ac(ac, cut_off, rnalen=rnalen)
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
    ax.set_xlim(left=0, right=rnalen )
    ax.set_ylim( [min_y, max_y] ) 
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    ax.set_xlabel(rnaname+" sequence (bp)", fontsize=9)
    
    ax.set_ylabel(ylabel,fontsize=9, rotation=90)
    
    f.tight_layout(pad=1.08, h_pad=None, w_pad=None)

    f.savefig(os.path.join(dirp, filename), facecolor='w', edgecolor='w',  
              bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
    # PDF
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14) 
    ax.legend(legend_h, legend_l, 
              bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
              prop={'size':12}, ncol=3)
    pp = PdfPages(os.path.splitext(os.path.join(dirp,filename))[0] +'.pdf')
    pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
    pp.close()

def load_dump(path, filename):
    file = open(os.path.join(path,filename),'r')
    object = pickle.load(file)
    file.close()
    print("\tLoading from file: "+filename)
    return object

def dump(object, path, filename):
    file = open(os.path.join(path,filename),'wb')
    pickle.dump(object,file)
    file.close()
    print("\tDump to file: "+filename)

def check_triplexator_path():
    try:
        cf = ConfigurationFile()
        with open(os.path.join(cf.data_dir,"triplexator_path.txt")) as f:
            l = f.readline()
            l = l.strip("\n")
            return l
    except:
        print("Please define the path to Triplexator by command: rgt-TDF triplexator -path <PATH>")
        sys.exit(1)
    

####################################################################################
####################################################################################

class PromoterTest:
    """Test the association between given triplex and differential expression genes"""
    def __init__(self, gene_list_file, bed, bg, organism, promoterLength, rna_name, 
                 summary, temp, showdbs=None, score=False, scoreh=False):
        """Initiation"""
        self.organism = organism
        genome = GenomeData(organism)
        self.rna_name = rna_name
        self.genome_path = genome.get_genome()
        self.showdbs = showdbs
        self.scores = None
        self.scoreh = scoreh

        # Input BED files
        if bed and bg:
            self.de_regions = GenomicRegionSet("de genes")
            self.nde_regions = GenomicRegionSet("nde genes")
            self.de_regions.read_bed(bed)
            self.de_regions.remove_duplicates()

            # score
            if score:
                self.scores = []
                for promoter in self.de_regions:
                    self.scores.append(promoter.data.split("\t")[0])

            try: self.de_regions = self.de_regions.gene_association(organism=self.organism)
            except: pass
            self.nde_regions.read_bed(bg)
            try: self.nde_regions = self.nde_regions.gene_association(organism=self.organism)
            except: pass
            self.nde_regions.remove_duplicates()

        # Input DE gene list
        else:
            try:
                ann = load_dump(path=temp, filename="annotation_"+organism)
            except:
                ann = AnnotationSet(organism, alias_source=organism)
                #print("DDD")
                #dump(object=ann, path=temp, filename="annotation_"+organism)

            # DE gene regions
            self.de_gene = GeneSet("de genes")
            self.de_gene.read(gene_list_file)
            if score:
                self.de_gene.read_expression(geneListFile=gene_list_file, header=scoreh)
            
            # When there is no genes in the list
            if len(self.de_gene) == 0:
                print("Error: No genes are loaded from: "+gene_list_file)
                print("Please check the format.")
                sys.exit(1)
            # Generate a dict for ID transfer


            #print("1")
            #print(len(self.de_gene.genes))

            de_ensembl, unmap_gs = ann.fix_gene_names(gene_set=self.de_gene)
            
            #print("2")
            #print(len(de_ensembl))

            # NonDE gene regions
            nde_ensembl = [ g for g in ann.symbol_dict.keys() if g not in de_ensembl ]
            
            self.de_gene.genes = de_ensembl
            
            self.nde_gene = GeneSet("nde genes")
            self.nde_gene.genes = nde_ensembl
            
            # Get promoters from de genes
            de_prom = ann.get_promoters(promoterLength=promoterLength, gene_set=self.de_gene)

            print("3")
            print(len(de_prom))

            de_prom.merge(namedistinct=True)

            print("4")
            print(len(de_prom))

            for promoter in de_prom:
                promoter.name = ann.get_official_symbol(gene_name_source=promoter.name)   
            self.de_regions = de_prom
            
            print2(summary, "   \t"+str(len(de_ensembl))+" unique target promoters are loaded")
            print2(summary, "   \t"+str(len(de_prom))+" unique target promoters are loaded")
            # Get promoters from nonDE gene
            nde_prom = ann.get_promoters(promoterLength=promoterLength, gene_set=self.nde_gene)
            nde_prom.merge(namedistinct=True)
            for promoter in nde_prom:
                promoter.name = ann.get_official_symbol(gene_name_source=promoter.name)
            self.nde_regions = nde_prom
            
            print2(summary, "   \t"+str(len(nde_ensembl))+" unique non-target promoters are loaded")
            print2(summary, "   \t"+str(len(nde_prom))+" unique non-target promoters are loaded")
            # Loading score
            if score:
                self.scores = []
                for promoter in self.de_regions:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
                    try: 
                        self.scores.append(abs(self.de_gene.values[promoter.name.upper()]))
                    except:
                        print("Warning: "+promoter.name+"\tcannot be mapped to get its score.")
                        self.scores.append(0)
                    

    def search_triplex(self, rna, temp, l, e, c, fr, fm, of, mf, remove_temp=False):
        print("    \tRunning Triplexator...")
        self.triplexator_p = [ l, e, c, fr, fm, of, mf ]
        # DE
        get_sequence(dir=temp, filename=os.path.join(temp,"de.fa"), regions=self.de_regions, 
                     genome_path=self.genome_path)
        run_triplexator(ss=rna, ds=os.path.join(temp,"de.fa"), 
                        output=os.path.join(temp, "de.txp"), 
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf)
        # non-DE
        get_sequence(dir=temp, filename=os.path.join(temp,"nde.fa"), regions=self.nde_regions, 
                     genome_path=self.genome_path)
        run_triplexator(ss=rna, ds=os.path.join(temp,"nde.fa"), 
                        output=os.path.join(temp, "nde.txp"), 
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf)
        if remove_temp:
            os.remove(os.path.join(temp,"de.fa"))
            os.remove(os.path.join(temp,"nde.fa"))
        
    def count_frequency(self, temp, remove_temp, cutoff, obedp=False):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""
        
        self.frequency = {}
        self.frequency["promoters"] = { "de": OrderedDict(), "nde": OrderedDict() }
        
        
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

        if self.showdbs:
            self.frequency["hits"] = { "de": OrderedDict(), "nde": OrderedDict() }
            numdbs_def = len(txp_def.get_dbs(rm_duplicate=True))
            numdbs_ndef = len(txp_ndef.get_dbs(rm_duplicate=True))
            for rbs in self.rbss:
                # DE
                self.frequency["hits"]["de"][rbs] = [ len(txp_def.merged_dict[rbs]), numdbs_def-len(txp_def.merged_dict[rbs]) ]
                # non-DE
                self.frequency["hits"]["nde"][rbs] = [ len(txp_ndef.merged_dict[rbs]), numdbs_ndef-len(txp_ndef.merged_dict[rbs]) ]

        self.txp_def = txp_def
        self.txp_ndef = txp_ndef
        
        if remove_temp:
            os.remove(os.path.join(temp,"de.txp"))
            os.remove(os.path.join(temp,"nde.txp"))

        if obedp:
            self.de_regions.write_bed(filename=os.path.join(temp,obedp+"_target_promoters.bed"))
            self.txp_de.write_bed(filename=os.path.join(temp,obedp+"_target_promoters_dbs.bed"), remove_duplicates=True)
            self.txp_def.write_bed(filename=os.path.join(temp,obedp+"_dbss.bed"), remove_duplicates=True)

    def fisher_exact(self, alpha):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        pvalues = []
        self.sig_region_promoter = []
        for rbs in self.frequency["promoters"]["de"]:
            table = numpy.array([self.frequency["promoters"]["de"][rbs], self.frequency["promoters"]["nde"][rbs]])
            self.oddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
            pvalues.append(p)

        # correction
        if len(self.frequency["promoters"]["de"].keys()) > 1:
            reject, pvals_corrected = multiple_test_correction(pvalues, alpha=alpha, method='indep')
        else:
            pvals_corrected = pvalues
        for i, rbs in enumerate(self.frequency["promoters"]["de"].keys()):
            self.pvalue[rbs] = pvals_corrected[i]
            if pvals_corrected[i] < alpha:
                self.sig_region_promoter.append(rbs)

        if self.showdbs:
            self.hoddsratio = {}
            self.hpvalue = {}
            pvalues = []
            self.sig_region_dbs = []
            for rbs in self.frequency["hits"]["de"]:
                table = numpy.array([self.frequency["hits"]["de"][rbs], self.frequency["hits"]["nde"][rbs]])
                self.hoddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
                pvalues.append(p)
            # correction
            if len(self.frequency["hits"]["de"].keys()) > 1:
                reject, pvals_corrected = multiple_test_correction(pvalues, alpha=alpha, method='indep')
            else:
                pvals_corrected = pvalues
            for i, rbs in enumerate(self.frequency["hits"]["de"].keys()):
                self.hpvalue[rbs] = pvals_corrected[i]
                if pvals_corrected[i] < alpha:
                    self.sig_region_dbs.append(rbs)
            
    def plot_lines(self, txp, rna, dirp, cut_off, log, ylabel, linelabel, filename, sig_region, ac=None, showpa=False):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(rna)
        if not self.rna_name: self.rna_name = rnas[0].name
        self.rna_len = len(rnas[0])

        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dirp=dirp, sig_region=sig_region, 
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

    def barplot(self, dirp, filename, sig_region, dbs=False):
        """Generate the barplot to show the difference between target promoters and non-target promoters"""
        def to_percent(y, position):
            # Ignore the passed in position. This has the effect of scaling the default
            # tick locations.
            s = str(100 * y)
            return s + '%'
                
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        ind = range(len(self.rbss))
        width = 0.35
        
        if not dbs:
            propor_de = [ float(b[0])/len(self.de_regions) for b in self.frequency["promoters"]["de"].values() ]
            propor_nde = [ float(b[0])/len(self.nde_regions) for b in self.frequency["promoters"]["nde"].values() ]
        else:
            propor_de = [ float(b[0])/(b[0]+b[1]) for b in self.frequency["hits"]["de"].values() ]
            propor_nde = [ float(b[0])/(b[0]+b[1]) for b in self.frequency["hits"]["nde"].values() ]
            
        max_y = max([ max(propor_de),max(propor_nde) ]) * 1.2
        
        # Plotting
        rect = patches.Rectangle(xy=(1,0), width=0.8, height=max_y, facecolor=sig_color, 
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        for i, rbs in enumerate(self.rbss):
            if rbs in sig_region:
                rect = patches.Rectangle(xy=(i+0.05 ,0), width=0.9, height=max_y, facecolor=sig_color, 
                                         edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
                ax.add_patch(rect)
        
        rects_de = ax.bar([i+0.15 for i in ind], propor_de, width, color=target_color, 
                          edgecolor = "none", label="Target promoters")
        rects_nde = ax.bar([i+0.15+width for i in ind], propor_nde, width, color=nontarget_color, 
                           edgecolor = "none", label="Non-target promoters")
        
        # Legend
        tr_legend, = plt.plot([1,1], color=target_color, linewidth=6, alpha=1)
        ntr_legend, = plt.plot([1,1], color=nontarget_color, linewidth=6, alpha=1)
        ax.legend( [tr_legend, ntr_legend, rect], ["Target promoters", "Non-target promoters", "Significant DBD"], 
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
                  prop={'size':9}, ncol=3)
        tr_legend.set_visible(False)
        ntr_legend.set_visible(False)

        tick_size=8
        # Y axis
        ax.set_ylim( [ 0, max_y ] ) 
        formatter = FuncFormatter(to_percent)
        # Set the formatter
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9) 
        ax.set_ylabel("Proportion of binding promoters (%)",fontsize=9, rotation=90)
        
        # X axis
        ax.set_xlim( [ 0, len(self.rbss) ] )
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax.set_xticks([i + 0.5 for i in range(len(self.rbss))])
        ax.set_xticklabels( [dbd.str_rna(pa=False) for dbd in self.rbss], rotation=35, 
                            ha="right", fontsize=tick_size)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        
        for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9) 
        ax.set_xlabel(self.rna_name+" DNA Binding Domains", fontsize=9)
    
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dirp, filename), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14) 
        pp = PdfPages(os.path.splitext(os.path.join(dirp,filename))[0] +'.pdf')
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()
        
    def gen_html(self, directory, parameters, ccf, align=50, alpha = 0.05):
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = "Promoter Test: "+dir_name
        self.link_d = OrderedDict()
        self.link_d["RNA"] = "index.html"
        self.link_d["Target promoters"] = "promoters.html"
        self.link_d["Parameters"] = "parameters.html"
        

        if self.organism == "hg19": self.ani = "human"
        elif self.organism == "mm9": self.ani = "mouse"
        else: 
            self.ani = None

        #############################################################
        # Index main page
        #############################################################
        html = Html(name=html_header, links_dict=self.link_d, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        
        html.add_figure("plot_promoter.png", align="left", width="45%", more_images=["bar_promoter.png"])
        
        if self.showdbs:
            html.add_figure("plot_dbss.png", align="left", width="45%", more_images=["bar_dbss.png"])
            
        # Table of merged TBS on promoters
        if self.showdbs:
            header_list = [ [ "", "", "Promoters", None, None, None, None, None, "DBSs", None, None, None, None, None ],
                            [ "#", "DBD", 
                              "Target Promoter", None, "Non-target Promoter", None, "Statistics", None, 
                              "Target Promoter", None, "Non-target Promoter", None, "Statistics", None ],
                            [ " ", " ", 
                              "with DBS", "without DBS", "with DBS", "without DBS", "OR", "<i>p</i>-value", 
                              "No. DBSs", "Other DBSs", "No. DBSs", "Other DBSs","OR", "<i>p</i>-value"] ]
            header_titles = [ [ "", "", "Statistics on promoter level", None, None, None, None, None, 
                                "Statistics on DBS level", None, None, None, None, None ],
                              [ "Rank of the talbe",
                                "DNA Binding Domain which is the functional region on RNA.",
                                "Promoters of the differential expression genes.", None,
                                "Promoters of the non-differential expression genes.", None,
                                "Statistics based on promoters", None,
                                "Promoters of the differential expression genes.", None,
                                "Promoters of the non-differential expression genes.", None,
                                "Statistics based on DBSs", None ], 
                              [ "",
                                "",
                                "Number of target promoters which contain DBSs (DNA Binding Sites).",
                                "Number of target promoters which don't contain DBSs (DNA Binding Sites).",
                                "Number of non-target promoters which contain DBSs (DNA Binding Sites).",
                                "Number of non-target promoters which don't contain DBSs (DNA Binding Sites).",
                                "Odds Ratio", "P-value",
                                "Number of DBSs found in the target promoters.",
                                "Number of DBSs not found in the target promoters.",
                                "Number of DBSs found in the non-target promoters.",
                                "Number of DBSs not found in the non-target promoters.",
                                "Odds Ratio", "P-value" ] 
                                ]
            border_list = [ " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:2pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"" ]
        else:
            header_list = [ ["#", "DBD", "Target Promoter", None, "Non-target Promoter", None, "Statistics", None ],
                            [" ", " ", "with DBS", "without DBS", "with DBS", "without DBS", "OR", "<i>p</i>" ] ]
            header_titles = [ [ "Rank of the talbe",
                                "DNA Binding Domain which is the functional region on RNA.",
                                "Promoters of the differential expression genes.", None,
                                "Promoters of the non-differential expression genes.", None,
                                "Statistics based on promoters", None ], 
                              [ "",
                                "",
                                "Number of target promoters which contain DBSs (DNA Binding Sites).",
                                "Number of target promoters which don't contain DBSs (DNA Binding Sites).",
                                "Number of non-target promoters which contain DBSs (DNA Binding Sites).",
                                "Number of non-target promoters which don't contain DBSs (DNA Binding Sites).",
                                "Odds Ratio", "P-value" ] 
                                ]
            border_list = [ "style=\"border-right:1pt solid gray\"",
                            "style=\"border-right:1pt solid gray\"", "",
                            "style=\"border-right:1pt solid gray\"", "",
                            "style=\"border-right:1pt solid gray\"", "",
                            "style=\"border-right:1pt solid gray\"" ]
        
        type_list = 'ssssssssssssssssssss'
        col_size_list = [50,50,50,50,50,50,50,50,50,50,50,50,50,50,50]
        data_table = []
        rank=0
        for rbs in self.frequency["promoters"]["de"].keys():
            if self.frequency["promoters"]["de"][rbs][0] < ccf: continue
            rank += 1
            if self.pvalue[rbs] < alpha:
                p_promoter = "<font color=\"red\">"+value2str(self.pvalue[rbs])+"</font>"
            else:
                p_promoter = value2str(self.pvalue[rbs])

            if self.showdbs:
                if self.hpvalue[rbs] < alpha:
                    p_hit = "<font color=\"red\">"+value2str(self.hpvalue[rbs])+"</font>"
                else:
                    p_hit = value2str(self.hpvalue[rbs])
            
            new_row = [ str(rank),
                        rbs.str_rna(pa=False),
                        '<a href="'+"dbds_promoters.html#"+rbs.str_rna()+
                        '" style="text-align:left">'+
                        value2str(self.frequency["promoters"]["de"][rbs][0])+'</a>', 
                        value2str(self.frequency["promoters"]["de"][rbs][1]), 
                        value2str(self.frequency["promoters"]["nde"][rbs][0]), 
                        value2str(self.frequency["promoters"]["nde"][rbs][1]), 
                        value2str(self.oddsratio[rbs]), 
                        p_promoter ]
            if self.showdbs:
                new_row += [ value2str(self.frequency["hits"]["de"][rbs][0]),
                             value2str(self.frequency["hits"]["de"][rbs][1]),
                             value2str(self.frequency["hits"]["nde"][rbs][0]),
                             value2str(self.frequency["hits"]["nde"][rbs][1]),
                             value2str(self.hoddsratio[rbs]),
                             p_hit ]
    
            data_table.append(new_row)
            
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titles, border_list=border_list, sortable=True)
        
        html.add_heading("Notes")
        html.add_list([ "DBD stands for functional DNA Binding Domain on RNA.",
                        "RBS stands for RNA Binding Site on RNA.",
                        "DBS stands for DNA Binding Site on DNA."])
        
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory,"index.html"))

        
        #############################################################
        # RNA subpage: Profile of targeted promoters for each merged DNA Binding Domain
        #############################################################
        header_list= ["#", "Target Promoter", "Gene", "DBSs counts", "DBS coverage"] 
        header_title= ["Rank", "The given target promoter", "Gene symbol", 
                       "Number of DNA Binding Sites binding to the promoter", 
                       "The proportion of the promoter covered by DBS binding"]
        # dbds_promoters.html
        html = Html(name=html_header, links_dict=self.link_d, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
     
        for rbsm in self.frequency["promoters"]["de"].keys():    
            html.add_heading("DNA Binding Domain: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for i, promoter in enumerate(self.txp_de.merged_dict[rbsm]):
                # Add information
                if self.ani:
                    pr = '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+"&position="+promoter.chrom+"%3A"+str(promoter.initial)+"-"+str(promoter.final)+'" style="text-align:left">'+promoter.toString(space=True)+'</a>'
                else:
                    pr = promoter.toString(space=True)

                data_table.append([ str(i+1),
                                    pr, 
                                    split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                                    str(len(self.promoter["de"]["rd"][promoter.toString()])),
                                    value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                                     ])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_title, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "dbds_promoters.html"))

        ################################################################
        ############# Parameters
        html = Html(name=html_header, links_dict=self.link_d, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_heading("Parameters")
        header_list = ["Description", "Arguments","Value"]

        if parameters.de:
            de = os.path.basename(parameters.de)
            bed = "False"
            bg = "False"
        else:
            de = "False"
            bed = os.path.basename(parameters.bed)
            bg = os.path.basename(parameters.bg)

        data_table = [ ["RNA sequence name", "-rn", parameters.rn ],
                       ["Input RNA sequence file", "-r", os.path.basename(parameters.r)],
                       ["Input file for defferentially expression gene list", "-de", de ],
                       ["Input BED file as promoters", "-bed", bed ],
                       ["Input BED file as backgrounds", "-bg", bg ],
                       ["Output directory", "-o", os.path.basename(parameters.o) ],
                       ["Organism", "-organism", parameters.organism ],
                       ["Promoter length", "-pl", str(parameters.pl) ],
                       ["Alpha level for rejection p value", "-a", str(parameters.a) ],
                       ["Cut off value for filtering out the low counts of DBSs", "-ccf", str(parameters.ccf) ],
                       ["Remove temporary files", "-rt", str(parameters.rt) ],
                       ["Input file for RNA accecibility", "-ac", str(parameters.ac) ],
                       ["Cut off value for RNA accecibility", "-accf", str(parameters.accf) ],
                       ["Output the BED files for DNA binding sites.", "-obed", str(parameters.obed) ],
                       ["Show parallel and antiparallel bindings in the plot separately.", "-showpa", str(parameters.showpa) ],
                       ["Minimum length", "-l", str(self.triplexator_p[0])],
                       ["Maximum error rate", "-e", str(self.triplexator_p[1])],
                       ["Tolerated number of consecutive errors", "-c", str(self.triplexator_p[2])],
                       ["Filtering repeats", "-fr", str(self.triplexator_p[3])],
                       ["Filtering mode", "-fm", str(self.triplexator_p[4])],
                       ["Output format", "-of", str(self.triplexator_p[5])],
                       ["Merge features", "-mf", str(self.triplexator_p[6])] ]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True)      
        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory,"parameters.html"))

            
    def gen_html_genes(self, directory, align=50, alpha = 0.05, nonDE=False):
        dir_name = os.path.basename(directory)
        html_header = "Promoter Test: "+dir_name
        
        type_list = 'sssssssssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10,10,10,10,10]

        #############################################################
        # Promoter centered
        #############################################################
        
        ##############################################################################################
        # promoters.html
        html = Html(name=html_header, links_dict=self.link_d, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if self.scores:
            if self.scoreh and self.de_gene:
                score_header = [self.de_gene.cond[0]]
            else:
                if "(" in self.scores[0]:
                    scs = self.scores[0].replace("(","")
                    scs = scs.replace(")","")
                    scs = scs.split(",")
                    score_header = ["Score"] * len(scs)
                else:        
                    score_header = ["Score"]
            header_list = [ "#", "Promoter", "Gene", "DBSs Count", "DBS coverage" ]
            header_list += score_header
            header_list += [ "Sum of Ranks" ]
        
            header_titles = [ "", "Target promoters", "Gene symbol", 
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites" ]
            header_titles += [ "Scores loaded by their absolute values from gene list or BED input. If there is annotation error for the gene names, it shows zero instead."] * len(score_header)
            header_titles += [ "Sum up the ranks from left-hand side columns" ]

        else:
            header_list=[ "#", "Promoter", "Gene", "DBSs Count", "DBS coverage", "Sum of Ranks" ]
        
            header_titles = [ "", "Target promoters", "Gene symbol", 
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites",
                              "Sum up the ranks from left-hand side columns"]
            
        html.add_heading("Target promoters")
        data_table = [] 
        
        if not self.de_regions.sorted: self.de_regions.sort()
        # Iterate by each gene promoter
        
        # Calculate the ranking
        rank_count = len(self.de_regions)-stats.rankdata([ len(self.promoter["de"]["dbs"][p.toString()]) for p in self.de_regions ], method='max')+1
        rank_coverage = len(self.de_regions)-stats.rankdata([ self.promoter["de"]["dbs_coverage"][p.toString()] for p in self.de_regions ], method='max')+1
        
        if self.scores:
            multiple_scores = False
            if isinstance(self.scores[0], str):
                if "(" in self.scores[0]:
                    def ranking(scores):
                        rank_score = len(self.de_regions)-stats.rankdata(scores, method='max')+1
                        return rank_score

                    multiple_scores = True
                    scs = []
                    for s in self.scores:
                        s = s.replace("(","")
                        s = s.replace(")","")
                        s = s.split(",")
                        scs.append([float(ss) for ss in s])
                    ar = numpy.array(scs)
                    #ar.transpose()
                    #print(ar)
                    rank_score = numpy.apply_along_axis(ranking, axis=0, arr=ar)
                    rank_score = rank_score.transpose()
                    rank_sum = numpy.sum(rank_score, axis=0).tolist()
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_sum)]
                    rank_score = rank_score.tolist()

                else:
                    self.scores = [ float(s) for s in self.scores ]     
                    rank_score = len(self.de_regions)-stats.rankdata(self.scores, method='max')+1
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]
    
            else: 
                rank_score = len(self.de_regions)-stats.rankdata(self.scores, method='max')+1
                rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

        for i, promoter in enumerate(self.de_regions):
            
            if len(self.promoter["de"]["dbs"][promoter.toString()]) == 0:
                dbssount = str(0)
            else:
                dbssount = '<a href="promoters_dbds.html#'+promoter.toString()+'" style="text-align:left">'+str(len(self.promoter["de"]["dbs"][promoter.toString()]))+'</a>'
            
            if self.ani:
                region_link = "".join(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=', self.organism,
                                       "&position=", promoter.chrom, "%3A", str(promoter.initial), "-",
                                       str(promoter.final), '" style="text-align:left">',
                                       promoter.toString(space=True), '</a>'])
            else:
                if self.organism == "tair10":
                    region_link = "".join(['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=', 
                                           promoter.chrom, "%3A", str(promoter.initial), "..", str(promoter.final),'">',
                                           promoter.toString(space=True), '</a>'])
                else: region_link = promoter.toString(space=True)

            newline = [ str(i+1),
                        region_link,
                        split_gene_name(gene_name=promoter.name, ani=self.ani, org=self.organism),
                        dbssount,
                        value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                        ]

            if self.scores: 
                if multiple_scores:
                    for j in range(len(rank_score)):
                        newline += value2str(rank_score[j][i])
                else:
                    newline += [value2str(self.scores[i])]

            newline += ["<i>"+str(int(rank_sum[i]))+"</i>"]
            print(newline)
            data_table.append(newline)

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titles, border_list=None, sortable=True)
        html.add_heading("Notes")
        html.add_list(["DBS stands for DNA Binding Site on DNA.",
                       "DBS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory,"promoters.html"))

        ############################
        # Subpages for promoter centered page
        # promoters_dbds.html
        header_sub = ["#", "RBS", "DBS", "Strand", "Score", "Motif", "Orientation"]
        header_titles = ["", "RNA Binding Site", "DNA Binding Site", "Strand of DBS on DNA",
                         "Score of binding event", "Motif of binding by triple helix rule",
                         "Orientation of interaction between DNA and RNA. 'P'- Parallel; 'A'-Antiparallel"]
        header_list=header_sub
        html = Html(name=html_header, links_dict=self.link_d, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        
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
                
                for j, rd in enumerate(self.promoter["de"]["rd"][promoter.toString()]):
                    data_table.append([ str(j+1),
                                        rd.rna.str_rna(pa=False),
                                        rd.dna.toString(space=True),
                                        rd.dna.orientation,
                                        rd.score, 
                                        rd.motif, rd.orient])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     header_titles=header_titles, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "promoters_dbds.html"))



####################################################################################
####################################################################################

class RandomTest:
    def __init__(self, rna_fasta, rna_name, dna_region, organism, showdbs=False):
    	self.organism = organism
        genome = GenomeData(organism)
        self.genome_path = genome.get_genome()
        # RNA: Path to the FASTA file
        self.rna_fasta = rna_fasta
        self.showdbs = showdbs

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
        
    def target_dna(self, temp, remove_temp, cutoff, l, e, c, fr, fm, of, mf, obed=False):
        """Calculate the true counts of triplexes on the given dna regions"""
        self.triplexator_p = [ l, e, c, fr, fm, of, mf ]

        txp = find_triplex(rna_fasta=self.rna_fasta, dna_region=self.dna_region, 
                           temp=temp, organism=self.organism, remove_temp=remove_temp, 
                           l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, genome_path=self.genome_path,
                           prefix="targeted_region", dna_fine_posi=False)
        txp.merge_rbs(rm_duplicate=True, asgene_organism=self.organism, cutoff=cutoff)
        self.txp = txp
        txp.remove_duplicates()
        self.rbss = txp.merged_dict.keys()
        if len(self.rbss) == 0:
            print("ERROR: No potential binding event. Please change the parameters.")
            sys.exit(1)

        txpf = find_triplex(rna_fasta=self.rna_fasta, dna_region=self.dna_region, 
                            temp=temp, organism=self.organism, remove_temp=remove_temp, 
                            l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, genome_path=self.genome_path,
                            prefix="dbs", dna_fine_posi=True)
        txpf.remove_duplicates()
        txpf.merge_by(rbss=self.rbss, rm_duplicate=True, asgene_organism=self.organism)
        self.txpf = txpf
        
        self.counts_tr = OrderedDict()
        self.counts_dbs = OrderedDict()

        for rbs in self.rbss:
            tr = len(self.txp.merged_dict[rbs])
            self.counts_tr[rbs] = [tr, len(self.dna_region) - tr]
            self.counts_dbs[rbs] = len(self.txpf.merged_dict[rbs])
            
        self.region_dbd = self.txpf.sort_rbs_by_regions(self.dna_region)
        
        self.region_dbs = self.txpf.sort_rd_by_regions(regionset=self.dna_region)
        self.region_dbsm = {}
        self.region_coverage = {}
        dbss = self.txpf.get_dbs()
        for region in self.dna_region:
            self.region_dbsm[region.toString()] = self.region_dbs[region.toString()].get_dbs().merge(w_return=True)
            self.region_coverage[region.toString()] = float(self.region_dbsm[region.toString()].total_coverage()) / len(region)
        
        if obed:
            btr = txp.get_dbs()
            btr = dbss.gene_association(organism=self.organism)
            btr.write_bed(os.path.join(temp, obed+"_target_region_dbs.bed"))
            dbss = txpf.get_dbs()
            dbss = dbss.gene_association(organism=self.organism)
            dbss.write_bed(os.path.join(temp, obed+"_dbss.bed"))

    def random_test(self, repeats, temp, remove_temp, l, e, c, fr, fm, of, mf, rm, filter_bed, alpha):
        """Perform randomization for the given times"""
        self.repeats = repeats
        marks = numpy.round(numpy.linspace(0, repeats-1, num=41)).tolist()
        
        # Prepare the input lists for multiprocessing
        mp_input = []
        for i in range(repeats):
            mp_input.append([ str(i), self.rna_fasta, self.dna_region,
                              temp, self.organism, self.rbss, str(marks.count(i)),
                              str(l), str(e), str(c), str(fr), str(fm), str(of), str(mf), str(rm),
                              filter_bed, self.genome_path ])
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
        self.region_matrix = []
        self.dbss_matrix = []
        self.data = {"region": {"ave": [],
                                "sd": [],
                                "p": [],
                                "sig_region": [],
                                "sig_boolean": [] },
                     "dbs": {"ave": [],
                             "sd": [],
                             "p": [],
                             "sig_region": [],
                             "sig_boolean": [] }   }

        region_counts = [ v[0] for v in mp_output ]
        dbss_counts = [ v[1] for v in mp_output ]

        for i, rbs in enumerate(self.rbss):
            
            counts_regions = [ v[i] for v in region_counts ]
            
            self.data["region"]["ave"].append(numpy.mean(counts_regions))
            self.data["region"]["sd"].append(numpy.std(counts_regions))
            num_sig = len([ h for h in counts_regions if h > self.counts_tr[rbs][0] ])
            p_region = float(num_sig)/repeats
            self.data["region"]["p"].append(p_region)
            self.region_matrix.append(counts_regions)

            if p_region < alpha: 
                self.data["region"]["sig_region"].append(rbs)
                self.data["region"]["sig_boolean"].append(True)
            else:
                self.data["region"]["sig_boolean"].append(False)
                
            # Analysis based on DBSs
            if self.showdbs:
                counts_dbss = [ v[i] for v in dbss_counts ]
            
                self.data["dbs"]["ave"].append(numpy.mean(counts_dbss))
                self.data["dbs"]["sd"].append(numpy.std(counts_dbss))
                num_sig = len([ h for h in counts_dbss if h > self.counts_dbs[rbs] ])
                p_dbs = float(num_sig)/repeats
                self.data["dbs"]["p"].append(p_dbs)
                self.dbss_matrix.append(counts_dbss)
         
                if p_dbs < alpha: 
                    self.data["dbs"]["sig_region"].append(rbs)
                    self.data["dbs"]["sig_boolean"].append(True)
                else:
                    self.data["dbs"]["sig_boolean"].append(False)
            
        self.region_matrix = numpy.array(self.region_matrix)
        if self.showdbs: self.dbss_matrix = numpy.array(self.dbss_matrix)
        
    def lineplot(self, txp, dirp, ac, cut_off, log, ylabel, linelabel, showpa, sig_region, filename):
        """Generate lineplot for RNA"""
        
        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dirp=dirp, sig_region=sig_region, 
                 cut_off=cut_off, log=log, ylabel=ylabel, linelabel=linelabel,  
                 filename=filename, ac=ac, showpa=showpa)

    def boxplot(self, dir, matrix, sig_region, truecounts, sig_boolean, ylabel, filename):
        """Generate the visualized plot"""
        tick_size = 8
        label_size = 9

        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
        max_y = int(max([matrix.max()] + truecounts) *1.1) + 1
        min_y = max(int(matrix.min()*0.9) - 1, 0)
        
        # Significant DBD
        rect = patches.Rectangle(xy=(1,0), width=0.8, height=max_y, facecolor=sig_color, 
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        for i, r in enumerate(sig_boolean):
            if r:
                rect = patches.Rectangle(xy=(i+0.6,min_y), width=0.8, height=max_y, facecolor=sig_color, 
                                         edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
                ax.add_patch(rect)
        
        # Plotting
        
        bp = ax.boxplot(matrix.transpose(), notch=False, sym='o', vert=True, 
                        whis=1.5, positions=None, widths=None, 
                        patch_artist=True, bootstrap=None)
        z = 10
        plt.setp(bp['boxes'], color=nontarget_color, alpha=1, edgecolor="none")
        plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=1,zorder=z, alpha=1)
        plt.setp(bp['fliers'], markerfacecolor='gray',color='white',alpha=0.3, markersize=1.8,zorder=z)
        plt.setp(bp['caps'],color='white',zorder=-1)
        plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
        
        # Plot target regions
        plt.plot(range(1, len(self.rbss)+1), truecounts, markerfacecolor=target_color,
                 marker='o', markersize=5, linestyle='None', markeredgecolor="white", zorder=z+5)
        
        ax.set_xlabel(self.rna_name+" DNA Binding Domains", fontsize=label_size)
        ax.set_ylabel(ylabel,fontsize=label_size, rotation=90)

        ax.set_ylim( [min_y, max_y] ) 
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        
        ax.set_xticklabels( [dbd.str_rna(pa=False) for dbd in self.rbss], rotation=35, 
                            ha="right", fontsize=tick_size)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(tick_size) 

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')

        # Legend
        dot_legend, = plt.plot([1,1], color=target_color, marker='o', markersize=5, markeredgecolor="white", linestyle='None')
        bp_legend, = plt.plot([1,1], color=nontarget_color, linewidth=6, alpha=1)
        

        ax.legend( [dot_legend, bp_legend, rect], ["Target Regions", "Non-target regions", "Significant DBD"], 
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
                  prop={'size':9}, ncol=3, numpoints=1)
        bp_legend.set_visible(False)
        dot_legend.set_visible(False)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, filename+".png"), facecolor='w', edgecolor='w',  
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12) 
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14) 
        
        pp = PdfPages(os.path.join(dir,filename+'.pdf'))
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def gen_html(self, directory, parameters, align=50, alpha=0.05, score=False):
        """Generate the HTML file"""
        dir_name = os.path.basename(directory)
        html_header = "Genomic Region Test: "+dir_name
        link_ds = OrderedDict()
        link_ds["RNA"] = "index.html"
        link_ds["Target regions"] = "target_regions.html"
        link_ds["Parameters"] = "parameters.html"

        if self.organism == "hg19": self.ani = "human"
        elif self.organism == "mm9": self.ani = "mouse"

        ##################################################
        # index.html

        html = Html(name=html_header, links_dict=link_ds, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        # Plots
        html.add_figure("lineplot_region.png", align="left", width="45%", more_images=["boxplot_regions.png"])
        if self.showdbs:
            html.add_figure("lineplot_dbs.png", align="left", width="45%", more_images=["boxplot_dbs.png"])

        if self.showdbs:
            header_list = [ ["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics", 
                             "Target Regions", "Non-target Regions", None, "Statistics"],
                            ["", "", "with DBS", "without DBS", "with DBS (average)", "s.d.", "<i>p</i>-value", 
                             "NO. DBSs", "NO. DBSs (average)", "s.d.", "<i>p</i>-value"] ]
            header_titles = [ ["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                               "Regions from randomization", None, "Statistics based on target regions",
                               "Given target regions on DNA","Regions from randomization", None, 
                               "Statistics based on DNA Binding Sites"],
                               ["", "", 
                                "Number of target regions with DBS binding",
                                "Number of target regions without DBS binding",
                                "Average number of regions from randomization with DBS binding",
                                "Standard deviation", "P value",
                                "Number of related DNA Binding Sites binding to target regions",
                                "Average number of DNA Binding Sites binding to random regions",
                                "Standard deviation", "P-value"]]
            border_list = [ " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:2pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"" ]
        else:
            header_list = [ ["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics" ],
                            ["", "", "with DBS", "without DBS", "with DBS (average)", "s.d.", "<i>p</i>-value" ] ]
            header_titles = [ ["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                               "Regions from randomization", None, "Statistics based on target regions"],
                               ["", "", 
                                "Number of target regions with DBS binding",
                                "Number of target regions without DBS binding",
                                "Average number of regions from randomization with DBS binding",
                                "Standard deviation", "P value"]]
            border_list = [ " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"", "",
                            " style=\"border-right:1pt solid gray\"",
                            " style=\"border-right:1pt solid gray\"" ]
            
        type_list = 'ssssssssssssssss'
        col_size_list = [50,50,50,50,50,50,50,50,50,50,50,50,50,50,50]
        data_table = []
        
        for i, rbs in enumerate(self.rbss):
            if self.data["region"]["p"][i] < alpha:
                p_region = "<font color=\"red\">"+value2str(self.data["region"]["p"][i])+"</font>"
            else:
                p_region = value2str(self.data["region"]["p"][i])
            
            new_line = [str(i+1),
                        rbs.str_rna(pa=False),
                        '<a href="dbd_region.html#'+rbs.str_rna()+
                        '" style="text-align:left">'+str(self.counts_tr[rbs][0])+'</a>',
                        str(self.counts_tr[rbs][1]),
                        value2str(self.data["region"]["ave"][i]), 
                        value2str(self.data["region"]["sd"][i]), 
                        p_region ]
            if self.showdbs:
                if self.data["dbs"]["p"][i] < alpha:
                    p_dbs = "<font color=\"red\">"+value2str(self.data["dbs"]["p"][i])+"</font>"
                else:
                    p_dbs = value2str(self.data["dbs"]["p"][i])

                new_line += [str(self.counts_dbs[rbs]),
                             value2str(self.data["dbs"]["ave"][i]), 
                             value2str(self.data["dbs"]["sd"][i]),
                             p_dbs]
            data_table.append(new_line)

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, border_list=border_list, sortable=True)

        html.add_heading("Notes")
        html.add_list([ "RNA name: "+ self.rna_name,
                        "Randomization is performed for "+str(self.repeats)+" times.",
                        "DBD stands for DNA Binding Domain on RNA.",
                        "DBS stands for DNA Binding Site on DNA."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory,"index.html"))

        

        #############################################################
        # RNA subpage: Profile of targeted regions for each merged DNA Binding Domain
        #############################################################
        
        header_list=[ "#", "Target Region", 
                      "Associated Gene",  
                      "No. of DBSs",
                      "DBS coverage" ]
        header_titles=["Rank", "Given target regions from BED files",
                       "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                       "Number of DNA Binding Sites locate within the region",
                       "The proportion of the region covered by DBS binding"]
                      
        #########################################################
        # dbd_region.html
        html = Html(name=html_header, links_dict=link_ds, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
     
        for rbsm in self.rbss:    
            html.add_heading("DNA Binding Domain: "+rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for i, region in enumerate(self.txp.merged_dict[rbsm]):
                # Add information
                data_table.append([ str(i+1),
                                    '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                                    "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                                    '" style="text-align:left">'+region.toString(space=True)+'</a>',
                                    split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism),
                                    str(len(self.region_dbs[region.toString()])),
                                    value2str(self.region_coverage[region.toString()])
                                     ])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 auto_width=True, header_titles=header_titles, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "dbd_region.html"))

        #############################################################
        # Targeted regions centered
        #############################################################
        
        ##############################################################################################
        # target_regions.html
        html = Html(name=html_header, links_dict=link_ds,#fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        
        if score:
            header_list = [ "#", "Target region", "Associated Gene", "DBSs Count", 
                            "DBS coverage", "Score", "Sum of ranks" ]
            header_titles = [ "Rank", 
                              "Target regions loaded from the given BED file",
                              "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                              "Number of DNA Binding Sites within the region",
                              "The proportion of the region covered by DBS binding",
                              "Scores from BED file",
                              "Sum of all the left-hand-side ranks"]
        else:
            header_list = [ "#", "Target region", "Associated Gene", "DBSs Count", 
                            "DBS coverage", "Sum of ranks" ]
            header_titles = [ "Rank", 
                              "Target regions loaded from the given BED file",
                              "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                              "Number of DNA Binding Sites within the region",
                              "The proportion of the region covered by DBS binding",
                              "Sum of all the left-hand-side ranks"]
        html.add_heading("Target Regions")
        data_table = []
        
        if not self.dna_region.sorted: self.dna_region.sort()
        # Iterate by each gene promoter

        nz_promoters = []
        for promoter in self.dna_region:
            if len(self.region_dbs[promoter.toString()]) > 0:
                nz_promoters.append(promoter) 

        # Calculate the ranking
        rank_count = len(nz_promoters)-stats.rankdata([ len(self.region_dbs[p.toString()]) for p in nz_promoters ], method='max')
        rank_coverage = len(nz_promoters)-stats.rankdata([ self.region_coverage[p.toString()] for p in nz_promoters ], method='max')
        
        if score:
            rank_score = len(nz_promoters)-stats.rankdata([ float(p.data.split("\t")[0]) for p in nz_promoters ], method='max')
            rank_sum = [x + y +z for x, y, z in zip(rank_count, rank_coverage, rank_score)]
            sum_rank = stats.rankdata(rank_sum, method='min')
        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]
            sum_rank = stats.rankdata(rank_sum, method='min')

        for i, region in enumerate(nz_promoters):
            newline = [ str(i+1),
                        '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db='+self.organism+
                        "&position="+region.chrom+"%3A"+str(region.initial)+"-"+str(region.final)+
                        '" style="text-align:left">'+region.toString(space=True)+'</a>',
                        split_gene_name(gene_name=region.name, ani=self.ani, org=self.organism),
                        '<a href="region_dbs.html#'+region.toString()+
                        '" style="text-align:left">'+str(len(self.region_dbs[region.toString()]))+'</a>',
                        value2str(self.region_coverage[region.toString()]) ]
            if score:
                newline += [ str(region.data.split("\t")[0]),
                             str(int(rank_sum[i]+3)) ]
            else:
                newline += [ str(int(rank_sum[i]+2)) ]
            data_table.append(newline)
                                
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, sortable=True)
        html.add_heading("Notes")
        html.add_list(["All target regions without any bindings are ignored." ])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory,"target_regions.html"))

        ############################
        # Subpages for targeted region centered page
        # region_dbs.html
        header_list = ["RBS", "DBS", "Strand", "Score", "Motif", "Orientation" ]

        html = Html(name=html_header, links_dict=link_ds, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        
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
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     auto_width=True)
        html.write(os.path.join(directory, "region_dbs.html"))

        ###############################################################################33
        ################ Parameters.html

        html = Html(name=html_header, links_dict=link_ds, #fig_dir=os.path.join(directory,"style"), 
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        html.add_heading("Parameters")
        header_list = ["Description", "Arguments","Value"]

        data_table = [ ["RNA sequence name", "-rn", parameters.rn ],
                       ["Input RNA sequence file", "-r", os.path.basename(parameters.r)],
                       ["Input BED file", "-bed", os.path.basename(parameters.bed) ],
                       ["Output directory", "-o", os.path.basename(parameters.o) ],
                       ["Organism", "-organism", parameters.organism ],
                       ["Number of repitetion of andomization", "-n", str(parameters.n) ],
                       ["Alpha level for rejection p value", "-a", str(parameters.a) ],
                       ["Cut off value for filtering out the low counts of DBSs", "-ccf", str(parameters.ccf) ],
                       ["Remove temporary files", "-rt", str(parameters.rt) ],
                       ["Input BED file for masking in randomization", "-f", str(parameters.f) ],
                       ["Input file for RNA accecibility", "-ac", str(parameters.ac) ],
                       ["Cut off value for RNA accecibility", "-accf", str(parameters.accf) ],
                       ["Output the BED files for DNA binding sites.", "-obed", str(parameters.obed) ],
                       ["Show parallel and antiparallel bindings in the plot separately.", "-showpa", str(parameters.showpa) ],
                       ["Minimum length", "-l", str(self.triplexator_p[0])],
                       ["Maximum error rate", "-e", str(self.triplexator_p[1])],
                       ["Tolerated number of consecutive errors", "-c", str(self.triplexator_p[2])],
                       ["Filtering repeats", "-fr", str(self.triplexator_p[3])],
                       ["Filtering mode", "-fm", str(self.triplexator_p[4])],
                       ["Output format", "-of", str(self.triplexator_p[5])],
                       ["Merge features", "-mf", str(self.triplexator_p[6])] ]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True)
        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory,"parameters.html"))