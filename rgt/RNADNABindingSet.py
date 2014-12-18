# Python Libraries
from __future__ import print_function
from collections import *
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from RNABinding import *
from RNABindingSet import *

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
    

if __name__ == '__main__':
    a = RNABindingSet(name="a")
    a.add(RNABinding("a",1,5))
    a.add(RNABinding("a",10,15))
    
    b = RNABindingSet(name="b")
    b.add(RNABinding("b",4,8))
    
    file_txp = "/projects/lncRNA/data/sample.txp"
    rd = RNADNABindingSet(organism="hg19", filename=file_txp)
    rd.get_rna()
    rd.get_dna()
    rd.sort_dna()
    rd.sort_rna()
    
    for i,j in rd.dict.iteritems():
    	print(str(i))
    	print("   "+str(j))
