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

    def __getitem__(self, key):
        return self.sequences[key]

    def read_fasta(self, fasta_file):
        """Read all the sequences in the given fasta file"""
        with open(fasta_file) as f:
            for line in f:
                line.strip()
                if not line: pass
                elif line[0] == ">":
                    try: # Save the previous sequence
                        self.sequences.append(s)
                    except: # Start to parse a new sequence
                        s = Sequence(seq_type=self.seq_type, name=line[1:])
                else:
                    s.seq = s.seq + line
            try: # Save the previous sequence
                self.sequences.append(s)
            except: # Start to parse a new sequence
                pass


    def read_bed(self, bedfile, genome_file_dir):
        """Read the sequences defined by BED file on the given genomce"""
        # Read BED into GenomicRegionSet
        bed = GenomicRegionSet(os.path.basename(bedfile))
        bed.read_bed(bedfile)
        
        # Parse each chromosome and fetch the defined region in this chromosome
        for s in bed:
            
        
        
        
        # Read genome by chromosome
        

if __name__ == '__main__':
    a = SequenceSet(name="test", seq_type=SequenceType.DNA)
    a.read_fasta("/Users/jovesus/data/sample.fasta")
    print(len(a))