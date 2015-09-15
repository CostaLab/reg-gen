"""
Author: Manuel Allhoff (allhoff@aices.rwth-aachen.de)

"""

from __future__ import print_function
import sys

from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from os.path import splitext, basename

class RegionGiver:
    regionset = GenomicRegionSet('')
    chrom_sizes_dict = {}
    
    def __init__(self, chrom_sizes, regions=None):
        self.counter = 0
        if regions is not None:
            print("Call DPs on specified regions.", file=sys.stderr)
            with open(regions) as f:
                for line in f:
                    line = line.strip()
                    line = line.split('\t')
                    c, s, e = line[0], int(line[1]), int(line[2])
            #if c in contained_chrom:                
                    self.regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
                    self.chrom_sizes_dict[c] = e
        else:
            print("Call DPs on whole genome.", file=sys.stderr)
            with open(chrom_sizes) as f:
                for line in f:
                    line = line.strip()
                    line = line.split('\t')
                    chrom, end = line[0], int(line[1])
            #if chrom in contained_chrom:
                    self.regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
                    self.chrom_sizes_dict[chrom] = end
        
        if not self.regionset.sequences:
            print('something wrong here', file=sys.stderr)
            sys.exit(2)
    
    def __len__(self):
        return len(self.regionset)
    
    def __iter__(self):
        for el in self.regionset:
            tmp = GenomicRegionSet('')
            tmp.add(el)
            yield tmp
        #return iter(self.regionset)
    
    def get_regionset(self):
        return self.regionset
    
    def get_chrom_dict(self):
        return self.chrom_sizes_dict
    
    
    def get_training_regionset(self):
        r = GenomicRegionSet('')
        r.add(self.regionset[self.counter])
        self.counter += 1
        if self.counter == len(self.chrom_sizes_dict):
            return None
        else:
            return r
            
            
#if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
#     if test:
#         contained_chrom = ['chr1', 'chr2']
#     else:
#         #contained_chrom = get_all_chrom(bamfiles)
#         contained_chrom = ['chr1', 'chr2']