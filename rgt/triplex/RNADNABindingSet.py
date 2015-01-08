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

    def __init__(self, rna, dna, score, err_rate, err, guan_rate, rna_seq=None, dna_seq=None, match=None):
        """Initiation"""
        self.rna = rna  # BindingSite
        self.dna = dna  # GenomicRegion
        self.motif = rna.motif
        self.strand = dna.orientation
        self.orient = rna.orientation
        self.score = score
        self.err_rate = err_rate
        self.err = err
        self.guan_rate = guan_rate
        self.rna_seq = rna_seq
        self.dna_seq = dna_seq
        self.match = match
    
    def __str__(self):
        return "\t".join( [ self.rna.chrom, str(self.rna.initial), str(self.rna.final), self.dna.toString(), "0", 
                            str(len(self.dna)), self.score, self.err_rate, self.err, self.motif,
                            self.strand, self.orient, self.guan_rate ] )
        
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
 

    def __len__(self):
        """Return the number of interactions between DNA and RNA """
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)
    
    def add(self, rnadnabinding):
        """Add one pair of BindingSite on RNA and GenomicRegion on DNA"""
        self.sequences.append(rnadnabinding)
        
    def get_rbs(self, sort=False):
        """Return a BindingSiteSet which contains all RNA binding sites"""
        rna_set = BindingSiteSet(name=self.name)
        for rd in self.sequences:
            rna_set.add(rd.rna)
        if sort: rna_set.sort()
        return rna_set

    def get_dbs(self, sort=False):
        """Return GenomicRegionSet which contains all DNA binding sites"""
        dna_set = GenomicRegionSet(name="DNA_binding_sites")
        for rd in self.sequences:
            dna_set.add(rd.dna)
        if sort: dna_set.sort()
        return dna_set

    def sort_rbs(self):
        """Sort the dictionary by RNA"""
        self.sequences = sorted(self.sequences, key=lambda x: x.rna, cmp=GenomicRegion.__cmp__)
        
    def sort_dbs(self):
        """Sort the dictionary by DNA"""
        self.sequences = sorted(self.sequences, key=lambda x: x.dna, cmp=GenomicRegion.__cmp__)
    
    def merge_rbs(self):
        """Merge the RNA binding regions which have overlap to each other and 
           combine their corresponding DNA binding regions.
        
        extend -> Define the extending length in basepair of each RNA binding regions
        perfect_match -> Merge only the exactly same RNA binding regions
        """
        # Merge RBS
        rna_merged = self.get_rbs()
        rna_merged.merge()
        print(len(rna_merged))
        # A dict: RBS as key, and GenomicRegionSet as its value
        new_dict = OrderedDict()
        
        for rbsm in rna_merged:
            for rd in self:
                if rbsm.overlap(rd.rna):
                # Add to new dictionary
                    try:
                        new_dict[rbsm].add(rd.dna)
                    except:
                        new_dict[rbsm] = GenomicRegionSet(rbsm.toString())
                        new_dict[rbsm].add(rd.dna)

        self.merged_dict = new_dict

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
                rna = BindingSite(chrom=line[0], initial=rna_start, final=rna_end, score=line[6], 
                                  errors_bp=line[8], motif=line[9], orientation=line[11], 
                                  guanine_rate=line[12])
                # DNA binding site
                dna_start = int(line[3].split(":")[1].split("-")[0]) + int(line[4])
                dna_end = int(line[3].split(":")[1].split("-")[0]) + int(line[5])
                dna = GenomicRegion(chrom=line[3].split(":")[0], initial=dna_start, final=dna_end, 
                                    name=line[3], orientation=line[10])
                
                # Map RNA binding site to DNA binding site
                self.add(RNADNABinding(rna=rna, dna=dna, score=line[6], err_rate=line[7], err=line[8], guan_rate=line[12]))

    def write_txp(self, filename):
        """Write RNADNABindingSet into the file"""
        try: os.stat(os.path.dirname(filename))
        except: os.mkdir(os.path.dirname(filename))   
        with open(filename, "w") as f:
            print("# RNA-ID\tRBS-start\tRBS-end\tDNA-ID\tDBS-start\tDBS-end\tScore\tError-rate\tErrors\
                   \tMotif\tStrand\tOrientation\tGuanine-rate", file=f)
            for rd in self:
                print(str(rd), file=f)

if __name__ == '__main__':
    
    file_txp = "/projects/lncRNA/data/sample.txp"
    rd = RNADNABindingSet(name="test")
    rd.read_txp(file_txp)
    rd.get_rnas()
    rd.get_dnas()
    rd.sort_dna()
    rd.sort_rna()
    
    for i in rd:
    	print(str(i))

    rd.write_txp("/projects/lncRNA/data/write.txp")
    	