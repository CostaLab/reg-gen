# Python Libraries
from __future__ import print_function
from collections import *
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

    def __init__(self, organism, name=None, filename=None):
        """Initialize a set of RNA/DNA interactions """
        
        self.name = name          # RNA name
        self.interactions = []    # All RNA/DNA interactions and extra information 
        self.sorted_tfo = False
        self.sorted_tts = False
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
                
                # Single strand ID (RNA)
                if not self.name: self.name = line[0]
                
                # Single strand region (RNA)
                TFO_start, TFO_end = int(line[1]), int(line[2])
                if TFO_start > TFO_end: TFO_start, TFO_end =  TFO_end, TFO_start
                
                # Double strand (DNA)
                TTS_start = int(line[3].split(":")[1].split("-")[0]) + int(line[4])
                TTS_end = int(line[3].split(":")[1].split("-")[0]) + int(line[5])
                ds = GenomicRegion(chrom=line[3].split(":")[0], 
                                   initial=TTS_start, 
                                   final=TTS_end, 
                                   name=line[3])
                
                # Orientation
                orientation = line[11]
                
                # Extra information: 
                # Score, Error-rate, Errors, Motif, Strand, Orientation, Guanine-rate
                data = line[6:]

                # Save information of each interaction into a list
                self.interactions.append([(TFO_start, TFO_end), ds, orientation, data])
                
            
    def __len__(self):
        """Return the number of triplex-forming oligonucleotides (TFO) on RNA """
        return len(self.interactions)

    def __iter__(self):
        return iter(self.interactions)

    def __getitem__(self, key):
        return self.interactions[key]

    def sort_tfo(self):
        """Sort the interactions by TFO (triplex-forming oligonucleotides)"""
        
        self.interactions.sort(key=lambda x: x[0][1])
        self.interactions.sort(key=lambda x: x[0][0])
        self.sorted_tfo = True
        
    def sort_tts(self):
        """Sort the interactions by TTS (triplex target sites)"""
        
        self.interactions.sort(key=lambda x: x[1]) #cmp=GenomicRegion.__cmp__
        self.sorted_tts = True

    def filter_tfo(self, start, end, output=True, adverse=False):
        """Return the filtered RNADNAInteractionSet by given boundary on RNA

        output:    True  --> return a new RNADNAInteractionSet
                   False --> change the self
        adverse:   False --> return the interaction within the given limit
                   True  --> return the interaction outside the given limit 
        """
        #if not self.sorted_tfo: self.sort_tfo()
        
        result = []
        for s in self:
            if adverse:
                if end <= s[0][0] or s[0][1] <= start:
                    result.append(s)
            else:
                if start <= s[0][0] and s[0][1] <= end:
                    result.append(s)

        if output: 
            new = RNADNAInteractionSet(name=self.name, organism=self.organism)
            new.interactions = result
            return new

        else: self.interactions = result

    def filter_tts(self, regionset, bedfile=False, output=True, adverse=False):
        """Return the filtered RNADNAInteractionSet which includes the TTS in the given regionset 
        
        bedfile    False --> regionset should be a GenomicRegionSet
                   True  --> regionset should be the path to the BED file
        output:    True  --> return a new RNADNAInteractionSet
                   False --> change the self
        adverse:   False --> return the interaction within the given limit
                   True  --> return the interaction outside the given limit 
        """
        #if not self.sorted_tfo: self.sort_tfo()
        

        if bedfile: 
            bed = GenomicRegionSet("mask")
            bed.read_bed(regionset)
            regionset = bed

        if adverse: regionset = regionset.complement(self.organism)

        result = []
        for s in self:
                if regionset.include(s[1]): 
                    result.append(s)

        if output: 
            new = RNADNAInteractionSet(name="filtered_"+self.name, organism=self.organism)
            new.interactions = result
            return new

        else: self.interactions = result

    def count_tfo(self, start, end):
        """Count the number of TFOs within the given boundary """
        select = self.filter_tfo(start, end)
        return len(select)

    def count_tts(self, regionset, bedfile=False):
        """Count the number of TTSs within the given GenomicRegionSet """
        select = self.filter_tts(regionset, bedfile=bedfile)
        return len(select)
    
    def get_tfos(self):
        """Return a list of all TFOs as tuples (strat position and end position) """
        return [t[0] for t in self]

    def get_ttss(self):
        """Return a GenomicRegionSet of all TTSs """
        g = GenomicRegionSet("TTS")
        for t in self:
            g.add(t[1])
        return g

    def unique_tfo(self):
        """Return a list of unique tfos """
        tfos = self.get_tfos()
        self.unique_tfos = list(set(tfos))
        return self.unique_tfos

    def merge_tfo(self):
        """Return a dictionary of interactions for merged TFOs 
        (TFO start, TFO end) as the key, 
        a GenomicRegionSet as the value
        
        TFO start and TFO end are position after merging
        GenomicRegionSet contains all the TTSs of the merged TFOs

        """
        # Empty or has only one tfo, do nothing
        if len(self) < 2: return
        
        if not self.sorted_tfo: self.sort_tfo()
        new = OrderedDict()
        ms = None
        me = None
        g = None
        for t in self:
            if me < t[0][0]:
                if ms: new[(ms,me)] = g
                ms = t[0][0]
                me = t[0][1]
                g = GenomicRegionSet("TTS")
                g.add(t[1])
                continue
            elif me < t[0][1]:
                me = t[0][1]
                g.add(t[1])
                continue
        new[(ms,me)] = g # Add the last region
        
        self.merged = new









