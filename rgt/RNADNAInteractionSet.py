# Python Libraries
from __future__ import print_function

# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
#from rgt.ExperimentalMatrix import *
#from rgt.AnnotationSet import *
from rgt.Util import GenomeData, OverlapType, Html
#from rgt.CoverageSet import *
#from rgt.motifanalysis.Statistics import multiple_test_correction

# Local test
#dir = os.getcwd()

"""
Represent list of RNA/DNA interactions (e.g. triplex formation) from txp format.

Authors: Joseph Kuo

Methods:


"""

class RNADNAInteractionSet():

    def __init__(self, organism, filename):
        """Initialize a set of RNA/DNA interactions """
        
        #self.name = name          # RNA name
        self.interactions = []    # All RNA/DNA interactions and extra information 
        self.sorted_tfo = False
        self.sorted_tts = False
        self.organism = organism
        self.read_txp(filename)
        self.fileName = filename
    
    def read_txp(self, filename):
        """Read txp file to load all interactions. """
        
        with open(filename) as f:
            for line in f:
                if line[0] == "#": continue # skip the comment line
                
                line = line.strip("\n")
                line = line.split()
                
                if len(line) < 12: continue # skip the unimportant lines in txp
                
                # Single strand ID (RNA)
                if not self.name: self.name = line[0]
                
                # Single strand region (RNA)
                TFO_start, TFO_end = int(line[1]), int(line[2])
                if TFO_start > TFO_end: TFO_start, TFO_end =  TFO_end, TFO_start
                
                # Double strand (DNA)
                TTS_start = int(line[3].split(":")[1].split("-")[0]) + int(line[4])
                TTS_end = int(line[3].split(":")[1].split("-")[0]) + int(line[5])
                ds = GenomicRegion(chrom=line[3].split(":")[0], 
                                   start=TTS_start, 
                                   end=TTS_end, 
                                   name=line[3])
                
                # Orientation
                orientation = line[11]
                
                # Extra information: Score, Error-rate, Errors, Motif, Strand, Guanine-rate
                data = line[6,7,8,9,10,12]
                
                # Save information of each interaction into a list
                self.interactions.append([TFO_start, TFO_end, ds, orientation, data])
                
            
    def __len__(self):
        """Return the number of triplex-forming oligonucleotides (TFO) on RNA """
        
        return len(self.interactions)
    
    def sort_tfo(self):
        """Sort the interactions by TFO (triplex-forming oligonucleotides)"""
        
        self.interactions.sort(key=lambda x: x[1])
        self.interactions.sort(key=lambda x: x[0])
        self.sorted_tfo = True
        
    def sort_tts(self):
        """Sort the interactions by TTS (triplex target sites)"""
        
        self.interactions.sort(key=lambda x: x[2])
        self.sorted_tts = True