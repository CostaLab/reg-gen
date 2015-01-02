# Python Libraries
from __future__ import print_function
from collections import *
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from BindingSiteSet import BindingSite, BindingSiteSet
from SequenceSet import Sequence, SequenceSet

####################################################################################
####################################################################################

class RegionSequenceDict:
    """A structure to connect a region (GenomicRegion or RNABinding) to its sequence (Sequence)

    Author: Joseph Kuo 
    """
    def __init__(self, name, seq_type):
        """Initiation"""
        self.name = name
        self.seq_type = seq_type
        self.dict = OrderedDict()

    def add(self, a_region, its_sequence):
        """Add a region with its corresponding sequence"""
        self.dict[a_region] = its_sequence

    def __len__(self):
        """Return the number of region-sequece pair """
        return len(self.dict.keys())

    def __iter__(self):
        """Iteration for the dictionary"""
        return self.dict.iteritems()

    def get_binding_sites(self):
        """Return all the regions by GenomicRegionSet or BindingSiteSet"""
        regions = BindingSiteSet("all binding sites")
        for s in self.dict.keys():
            regions.add(s)
        return regions

    def get_sequences(self):
        """Return all sequences in a SequenceSet"""
        sequences = SequenceSet("all sequences", self.seq_type)
        for s in self.values():
            sequences.add(s)
        return sequences

