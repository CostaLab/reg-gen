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

from __future__ import print_function
import sys

# from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
# from os.path import splitext, basename

class RegionGiver:
    """ provide Region to analyse or masked
    regionset = region_giver.region
    regionset.sequences.sort()
    chrom_sizes = region_giver.get_chrom_sizes
    chrom_sizes_dict = chrom_sizes.get_chrom_sizes_dict()

    region_giver = RegionGiver(chrom_sizes, options.regions)

    For this regions we should have valid_regions which is chrom region with non_zeros data
    For zero chroms regions we don't consider it ...
    we have statistical information in format ['chrm1',size, read_nums,unmapped_read_nums]

    """
    # why here to use class paramter, it belong to whole class not instance.
    # However we don't need to create new object form it
    # from simple ones begin !!!
    def __init__(self, chrom_file, regions_file=None, file_statics=None, ignored_regions_file=None):
        self.counter = 0

        self.chrom_sizes_file = chrom_file
        self.regionset = GenomicRegionSet('')

        self.ignored_region_file = ignored_regions_file
        self.ignored_regionset = GenomicRegionSet('Ignored_regions')

        self.chrom_sizes_dict = {}
        ## stats_data in dictionary format
        # valid_chroms are union of chroms from stats_data, then we get data from them and initalize data for training
        valid_chroms = set(tmp[0] for each_statics in file_statics for tmp in each_statics)
        if regions_file:
            with open(regions_file) as f:
                for line in f:
                    if line:
                        line = line.strip()
                        line = line.split()
                        chrom, start, end = line[0], int(line[1]), int(line[2])
                        if chrom in valid_chroms:
                            self.regionset.add(GenomicRegion(chrom=chrom, initial=start, final=end))
                            self.chrom_sizes_dict[chrom] = end
        else:
            print("Call DPs on whole genome.", file=sys.stderr)
            with open(chrom_file) as f:
                for line in f:
                    line = line.strip()
                    line = line.split('\t')
                    chrom, end = line[0], int(line[1])
                    if chrom in valid_chroms:
                        self.regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
                        self.chrom_sizes_dict[chrom] = end


        """
        if  ignored_regions_file: 
            with open( ignored_regions_file) as f:
                for line in f:
                    if line:
                        line = line.strip()
                        line = line.split()
                        c, s, e = line[0], int(line[1]), int(line[2])
                        self.ignored_regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
        """
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