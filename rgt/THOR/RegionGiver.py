"""
THOR detects differential peaks in multiple ChIP-seq profiles associated
with two distinct biological conditions.

Copyright (C) 2014-2016 Manuel Allhoff (allhoff@aices.rwth-aachen.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@author: Manuel Allhoff
"""


import sys

from ..GenomicRegion import GenomicRegion
from ..GenomicRegionSet import GenomicRegionSet

class RegionGiver:
    regionset = GenomicRegionSet('')
    chrom_sizes_dict = {}
    
    def __init__(self, chrom_sizes, regions=None):
        self.counter = 0
        if regions is not None:
            print("Call DPs on specified regions.", file=sys.stderr)
            with open(regions) as f:
                for line in f:
                    if line:
                        line = line.strip()
                        line = line.split()
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
        
        if self.counter == len(self.chrom_sizes_dict):
            return None
        else:
            self.counter += 1
            return r
        
            
            
#if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
#     if test:
#         contained_chrom = ['chr1', 'chr2']
#     else:
#         #contained_chrom = get_all_chrom(bamfiles)
#         contained_chrom = ['chr1', 'chr2']