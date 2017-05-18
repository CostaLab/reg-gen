# Python Libraries
from __future__ import print_function
import os
import copy
import pysam
# Local Libraries

# Distal Libraries
from Util import SequenceType

####################################################################################
####################################################################################
"""
Sequence
===================
Sequence describes the sequence with ATCG as alphabets as well as its types.

"""
class Sequence():

    def __init__(self, seq, strand, name=None):
        """*Keyword arguments:*

            - seq -- The sequence composef of A, T, C, G, or U.
            - strand -- Orientation of the sequence, '+' or '-'.
        """
        self.name = name # Extra information about the sequence
        self.seq = seq.upper() # Convert all alphabets into upper case
        self.strand = strand
        
    def __len__(self):
        """Return the length of the sequence."""
        return len(self.seq)

    def __str__(self):
        return ":".join([str(x) for x in [self.name, self.seq] if x])
        
    def dna_to_rna(self):
        """Convert the sequence from DNA sequence to RNA sequence."""
        self.seq = self.seq.replace("T","U")
        
    def rna_to_dna(self):
        """Convert the sequence from RNA sequence to DNA sequence."""
        self.seq = self.seq.replace("U", "T")
        
    def GC_content(self):
        """Return the ratio of GC content in this sequence."""
        gc = self.seq.count("G") + self.seq.count("C")
        return gc/float(len(self))

    def methylated_sequence(self, CpG):
#ALPHABET "DNA with covalent modifications" DNA-LIKE

# Core symbols
# A "Adenine" 8510A8 ~ T "Thymine" A89610
# C "Cytosine" A50026 ~ G "Guanine" 313695
# m "5-Methylcytosine" D73027 ~ 1 "Guanine:5-Methylcytosine" 4575B4
        aux=list(self.seq)
        for pos in CpG:
          if self.seq[pos]=="C":
            aux[pos]="m"
          elif self.seq[pos]=="G":
            aux[pos]="1"
          else:
            raise Exception("Position "+str(pos)+"-"+self.seq[pos]+" from "+self.seq+" can not be methylated")         
        self.seq="".join(aux)
        
    def complement(self):
        """Return another Sequence which is the complement to original sequence."""
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

"""
SequenceSet
===================
SequenceSet represents the collection of sequences and their functions.

"""

class SequenceSet:

    def __init__(self, name, seq_type):
        """*Keyword arguments:*

            - name -- The name of this collection of sequences.
            - seq_type -- DNA or RNA.
        """
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
        """Add one more sequence into the object.

        *Keyword arguments:*

            - seq -- A Sequence object.
        """
        self.sequences.append(seq)

    def read_fasta(self, fasta_file):
        """Read all the sequences in the given fasta file.

        *Keyword arguments:*
            -fasta_file -- The path to the FASTA file.
        """
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
                    try: seq = seq + line
                    except: seq = line
                    pre_seq = True
            self.sequences.append(Sequence(seq=seq, strand=strand, name=info))


    def read_regions(self, regionset, genome_fasta, ex=0):
        genome = pysam.Fastafile(genome_fasta)

        for region in regionset:
            seq = genome.fetch(region.chrom, region.initial - ex, region.final + ex)
            seq = seq.upper()

            if region.orientation == "-":
                t = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
                ns = ""
                for s in seq: ns += t[s]
                seq = ns[::-1]

            self.add(Sequence(seq=seq, strand=region.orientation, name=region.name))

    # def read_genomic_set(self, genomic_set, genome_file_dir):
    #     """Read the sequences defined by a given genomic set.s
    #     *Keyword arguments:*
    #
    #         - genomic_set - genomic set with regions to obtain the fasta file
    #         - genome_file_dir -- A directory which contains the FASTA files for each chromosome.
    #     """
    #
    #     bed = genomic_set
    #     # Parse each chromosome and fetch the defined region in this chromosome
    #     chroms = list(set(bed.get_chrom()))
    #
    #     chro_files = [x.split(".")[0] for x in os.listdir(genome_file_dir)]
    #
    #     for ch in chroms:
    #         if ch not in chro_files: print(" *** There is no genome FASTA file for: "+ch)
    #
    #         # Read genome in FASTA according to the given chromosome
    #         ch_seq = SequenceSet(name=ch, seq_type=SequenceType.DNA)
    #         try:
    #             ch_seq.read_fasta(os.path.join(genome_file_dir, ch+".fa"))
    #         except:
    #             continue
    #
    #         # Regions in given chromosome
    #         beds = bed.any_chrom(chrom=ch)
    #
    #         for s in beds:
    #             seq = ch_seq[0].seq[s.initial:s.final]
    #             try: strand = s.strand
    #             except: strand = "+"
    #             self.sequences.append(Sequence(seq=seq, name=s.__repr__(), strand=strand))



    # def read_bed(self, bedfile, genome_file_dir):
    #     """Read the sequences defined by BED file on the given genomce.
    #
    #     *Keyword arguments:*
    #
    #         - bedfile -- The path to the BED file which defines the regions.
    #         - genome_file_dir -- A directory which contains the FASTA files for each chromosome.
    #     """
    #
    #     # Read BED into GenomicRegionSet
    #     from rgt.GenomicRegionSet import GenomicRegionSet
    #     bed = GenomicRegionSet(os.path.basename(bedfile))
    #     bed.read_bed(bedfile)
    #     self.read_genomic_set(bed, genome_file_dir)
        


    def total_len(self):
        """Retrun the total length of the given SequenceSet."""
        tol = 0
        for s in self:
            tol += len(s)
        return tol
 
    def cal_motif_statistics(self):
        self.motif_statistics_1 = {"A":0, "T":0, "C":0, "G":0}
        motif_statistics_2 = {"AA": 0, "AT": 0, "AC": 0, "AG": 0,
                              "TA": 0, "TT": 0, "TC": 0, "TG": 0,
                              "CA": 0, "CT": 0, "CC": 0, "CG": 0,
                              "GA": 0, "GT": 0, "GC": 0, "GG": 0 }
        for s in self:
            for nt1 in self.motif_statistics_1.keys():
                self.motif_statistics_1[nt1] += s.seq.count(nt1)
            for nt2 in motif_statistics_2.keys():
                motif_statistics_2[nt2] += s.seq.count(nt2)

        # motif2 = {"AA": motif_statistics_2["AA"],
        #           "AT/TA": motif_statistics_2["AT"] + motif_statistics_2["TA"],
        #           "AC/CA": motif_statistics_2["AC"] + motif_statistics_2["CA"],
        #           "AG/GA": motif_statistics_2["AG"] + motif_statistics_2["GA"],
        #           "TT": motif_statistics_2["TT"],
        #           "TC/CT": motif_statistics_2["TC"] + motif_statistics_2["CT"],
        #           "TG/GT": motif_statistics_2["TG"] + motif_statistics_2["GT"],
        #           "CC": motif_statistics_2["CC"],
        #           "CG/GC": motif_statistics_2["CG"] + motif_statistics_2["GC"],
        #           "GG": motif_statistics_2["GG"] }

        # self.motif_statistics_2 = motif2
        self.motif_statistics_2 = motif_statistics_2
