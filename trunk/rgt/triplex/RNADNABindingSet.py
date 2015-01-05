# Python Libraries
from __future__ import print_function
from collections import *
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from BindingSiteSet import BindingSite, BindingSiteSet

class RNADNABinding:
    """Describe a binding event of RNA and DNA to form a triplex 
    
    Attributes:
        self.rna        A BindingSite object which includes RNA name, \
                        region (start position and end position), \
                        score, error and its Sequence
        self.dna        A BindingSite object which includes DNA name, \
                        region (start position and end position), \
                        score, error and its Sequence
        self.score      The score from comparing RNA region and DNA region 
        self.err_rate   The ration of error base pair between RNA and DNA
        self.err        *************
        self.motif      The motif of triplex formation:
                        R - the purine motif that permit guanines (G) and adenines (A)
                        Y - the pyrimidine motif that permit cytosines (C) and thymines (T)
                        M - the mixed motif, purine-pyrimdine, that permit guanines (G) and thymines (T)
                        P - parallel binding, i.e. motifs facilitating Hoogsten bonds; 
                            this covers the pyrimidine motif and the purine-pyrimidine motif in parallel configuration
                        A - anti-parallel binding, i.e. motifs facilitating reverse Hoogsten bonds; 
                            this covers the purine motif and the purine-pyrimidine motif in anti-parallel configuration
        self.strand     The strand of DNA
        self.orient     The orientation of RNA regarding to DNA
        self.guan_rate  *************
        self.rna_seq    A string of RNA extended by gaps to match the position of DNA
        self.dna_seq    A string of DNA extended by gaps to match the position of RNA
        self.match      A string with '|' for perfect match, '*' for mismatch
    
    """

    def __init__(self, rna, dna, score, err_rate, err, guan_rate, rna_seq, dna_seq, match):
        """Initiation"""
        self.rna = rna
        self.dna = dna
        self.motif = rna.motif
        self.strand = dna.strand
        self.orient = rna.orientation
        self.score = score
        self.err_rate = err_rate
        self.err = err
        self.guan_rate = guan_rate
        self.rna_seq = rna_seq
        self.dna_seq = dna_seq
        self.match = match
    
        
class RNADNABindingSet:
    """A framework to map DNA binding sites to corresponding RNA binding sites 
    
    Attributes:
        self.rna    A BindingSequence-ID   
        [1] TFO start   
        [2] TFO end 
        [3] Duplex-ID   
        [4] TTS start   
        [5] TTS end 
        [6] Score   
        [7] Error-rate  
        [8] Errors  
        [9] Motif   
        [10] Strand  
        [11] Orientation 
        [12] Guanine-rate
    
    """

    def __init__(self, name):
        """Initialize"""
        self.name = name       # RNA name
        self.sequences = []    # A list to contain all RNA/DNA interactions and extra information 
        self.sorted_dna = False
        self.sorted_rna = False
    
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
                rna = BindingSite(chrm=line[0], initial=rna_start, final=rna_end, score=line[6], 
                                  errors_bp=line[8], motif=line[9], orientation=line[11], 
                                  guanine_rate=line[12])
                
                # DNA binding site
                dna_start = int(line[3].split(":")[1].split("-")[0]) + int(line[4])
                dna_end = int(line[3].split(":")[1].split("-")[0]) + int(line[5])
                dna = GenomicRegion(chrom=line[3].split(":")[0], 
                                    initial=dna_start, final=dna_end, 
                                    name=line[3], orientation=line[10])
                
                # Map RNA binding site to DNA binding site
                self.add(RNADNABinding(self, rna, dna, score, err_rate, err, guan_rate, rna_seq, dna_seq, match))

    def __len__(self):
        """Return the number of interactions between DNA and RNA """
        return len(self.dict.keys())

    def __iter__(self):
        return self.dict.iteritems()
    
    def add(self, rna, dna, score, match):
        """Add one pair of BindingSite on RNA and GenomicRegion on DNA"""
        self.dict[rna] = [dna
        
    def get_rna(self, sort=False):
        """Return a BindingSiteSet which contains all RNA binding sites"""
        rna_set = BindingSiteSet(name=self.name)
        rna_set.sequences =  self.dict.keys()
        if sort: rna_set.sort()
        return rna_set

    def get_dna(self, sort=False):
        """Return GenomicRegionSet which contains all DNA binding sites"""
        dna_set = GenomicRegionSet(name="DNA_binding_sites")
        dna_set.sequences = self.dict.values()
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

if __name__ == '__main__':
    a = BindingSiteSet(name="a")
    a.add(BindingSite("a",1,5))
    a.add(BindingSite("a",10,15))
    
    b = BindingSiteSet(name="b")
    b.add(BindingSite("b",4,8))
    
    file_txp = "/projects/lncRNA/data/sample.txp"
    rd = RNADNABindingSet(name="test")
    rd.get_rna()
    rd.get_dna()
    rd.sort_dna()
    rd.sort_rna()
    
    for i,j in rd.dict.iteritems():
    	print(str(i))
    	print("   "+str(j))
