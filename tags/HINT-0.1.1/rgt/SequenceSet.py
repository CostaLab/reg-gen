# Python Libraries
from __future__ import print_function
import os
# Local Libraries

# Distal Libraries
from Sequence import Sequence
from Util import SequenceType
from rgt.GenomicRegionSet import GenomicRegionSet


"""
Define the collection of sequences and their functions
Author: joseph Kuo
"""

class SequenceSet:
    def __init__(self, name, seq_type):
        """Initiation"""
        self.name = name
        self.sequences = []
        self.seq_type = seq_type
        
    def __len__(self):
        return len(self.sequences)
        
    def __iter__(self):
        return iter(self.sequences)
        
    def read_fasta(self, fasta_file):
        """Read all the sequences in the given fasta file"""
        with open(fasta_file) as f:
            for line in f:
                line.strip()
                print(line)
                if not line: pass
                elif line[0] == ">":
                    try: # Save the previous sequence
                        self.sequences.append(s)
                    except: # Start to parse a new sequence
                        print("new s")
                        s = Sequence(seq_type=self.seq_type, info=line[1:])
                else:
                    print("s con")
                    s.seq = s.seq + line
                
    def read_bed(self, bedfile, genome_file):
        """Read the sequences defined by BED file on the given genomce"""
        bed = GenomicRegionSet(os.path.basename(bedfile))
        bed.read_bed(bedfile)
        
        # Read genome by chromosome

if __name__ == '__main__':
    a = SequenceSet(name="test", seq_type=SequenceType.DNA)
    a.read_fasta("/Users/jovesus/data/sample.fasta")
    print(len(a))