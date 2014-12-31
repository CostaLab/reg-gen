from Util import SequenceType
"""
Describe the Sequence with ATCG alphabets as well as its types

Author: JosephKuo
"""

class Sequence():
    def __init__(self, seq_type=SequenceType.DNA, name="", seq="", strand="+"):
        """Initiation"""
        self.seq_type = seq_type # DNA or RNA
        self.name = name # Extra information about the sequence
        self.seq = seq.upper() # Convert all alphabets into upper case
        self.strand = strand
        
    def __len__(self):
        """Return the length of the sequence"""
        return len(self.seq)
        
    def dna_to_rna(self):
        """Convert the sequence from DNA sequence to RNA sequence"""
        self.seq = self.seq.replace("T","U")
        self.seq_type = SequenceType.DNA
        
    def rna_to_dna(self):
        """Convert the sequence from RNA sequence to DNA sequence"""
        self.seq = self.seq.replace("U", "T")
        self.seq_type = SequenceType.RNA
        
    def GC_content(self):
        """Return the ratio of GC content in this sequence"""
        gc = self.seq.count("G") + self.seq.count("C")
        return gc/float(len(self))
    
    def find_TTS():
        """Return a triplex target sites"""

    def find_triplex_region(self, min_len, max_len, parallel):
        """Return a dictionary of possible which takes RNABinding (or GenomicRegion when the Sequence is DNA) as the keys and its correspoinding """
    

if __name__ == '__main__':
    a = Sequence(name="xx", seq="AGGCCTT")
    print(len(a))
    print(a.GC_content())
    
    print(a.seq)
        
        
        