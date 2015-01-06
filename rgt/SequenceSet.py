# Python Libraries
from __future__ import print_function
import os
import copy
# Local Libraries

# Distal Libraries
from Util import SequenceType
from rgt.GenomicRegionSet import GenomicRegionSet

####################################################################################
####################################################################################

class Sequence():
    """Describe the Sequence with ATCG as alphabets as well as its types

    Author: JosephKuo
    """

    def __init__(self, seq, strand, name=None):
        """Initiation"""
        self.name = name # Extra information about the sequence
        self.seq = seq.upper() # Convert all alphabets into upper case
        self.strand = strand
        
    def __len__(self):
        """Return the length of the sequence"""
        return len(self.seq)

    def __str__(self):
        return ":".join([str(x) for x in [self.name, self.seq] if x])
        
    def dna_to_rna(self):
        """Convert the sequence from DNA sequence to RNA sequence"""
        self.seq = self.seq.replace("T","U")
        
    def rna_to_dna(self):
        """Convert the sequence from RNA sequence to DNA sequence"""
        self.seq = self.seq.replace("U", "T")
        
    def GC_content(self):
        """Return the ratio of GC content in this sequence"""
        gc = self.seq.count("G") + self.seq.count("C")
        return gc/float(len(self))
        
    def complement(self):
        """Return another Sequence which is the complement to original sequence"""
        if self.strand == "RNA":
            print("No complement strand for RNA "+self.name)
            return
        elif self.strand == "+": strand = "-"
        else: strand = "-"
        
        seq = copy.deepcopy(self.seq)
        seq.replace("A","G")
        seq.replace("G","A")
        seq.replace("C","T")
        seq.replace("T","C")
        
        s = Sequence(seq, name=self.name+"_complement",  strand=strand)
        return s

####################################################################################
####################################################################################
 
class SequenceSet:           
    """Define the collection of sequences and their functions
    Author: joseph Kuo
    """

    def __init__(self, name, seq_type):
        """Initiation"""
        self.name = name
        self.sequences = []
        self.seq_type = seq_type
        
    def __len__(self):
        return len(self.sequences)
        
    def __iter__(self):
        return iter(self.sequences)

    def __getitem__(self, key):
        return self.sequences[key]

    def add(self, seq):
        """Add one more sequence into the object"""
        self.sequences.append(seq)

    def read_fasta(self, fasta_file):
        """Read all the sequences in the given fasta file"""
        pre_seq = False
        with open(fasta_file) as f:
            for line in f:
                line = line.strip()
                if not line: pass
                elif line[0] == ">":
                    if pre_seq:
                        self.sequences.append(Sequence(seq=seq, strand=strand, name=info))
                        pre_seq = False
                    info = line.split()[0][1:]
                    seq = ""
                    try: strand = line[line.index("strand")+7]
                    except: strand = "+"
                else:
                    seq = seq + line
                    pre_seq = True
            self.sequences.append(Sequence(seq=seq, strand=strand, name=info))

    def read_bed(self, bedfile, genome_file_dir):
        """Read the sequences defined by BED file on the given genomce"""

        # Read BED into GenomicRegionSet
        bed = GenomicRegionSet(os.path.basename(bedfile))
        bed.read_bed(bedfile)
        
        # Parse each chromosome and fetch the defined region in this chromosome
        chroms = list(set(bed.get_chrom()))

        chro_files = [x.split(".")[0] for x in os.listdir(genome_file_dir)]

        for ch in chroms:
            if ch not in chro_files: print(" *** There is no genome FASTA file for: "+ch)

            # Read genome in FASTA according to the given chromosome
            ch_seq = SequenceSet(name=ch, seq_type=SequenceType.DNA)
            try: 
                ch_seq.read_fasta(os.path.join(genome_file_dir, ch+".fa"))
            except:
                continue
            
            # Regions in given chromosome
            beds = bed.any_chrom(chrom=ch)

            for s in beds:
                seq = ch_seq[0].seq[s.initial:s.final]
                try: strand = s.strand
                except: strand = "+"
                self.sequences.append(Sequence(seq=seq, name=s.__repr__(), 
                                               strand=strand))

####################################################################################
####################################################################################

if __name__ == '__main__':
    a = SequenceSet(name="test")
    a.read_bed("/projects/stemaging/data/beds/f6_hyper.bed", "/media/joseph/931ef578-eebe-4ee8-ac0b-0f4690f126ed/data/genome")
    print(len(a))
    print(str(a[0]))