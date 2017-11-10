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
    # use class attributes to create common regions
    regionset = None
    chrom_sizes_file = None
    mask_file = None
    chrom_sizes_dict = {}

    def __init__(self, chrom_file, regions_file=None, mask_file=None):
        self.counter = 0

        self.chrom_sizes_file = chrom_file
        self.regionset = GenomicRegionSet('')

        self.mask_file = mask_file
        self.ignored_regionset = GenomicRegionSet('Ignored_regions')

        ## instance attribute
        self.valid_regionset = self.regionset
        self.valid_chrom_sizes = self.chrom_sizes_dict

        if regions_file and not self.regionset:
            self.regionset = GenomicRegionSet('')
            with open(regions_file) as f:
                for line in f:
                    if line:
                        line = line.strip()
                        line = line.split()
                        chrom, start, end = line[0], int(line[1]), int(line[2])
                        self.regionset.add(GenomicRegion(chrom=chrom, initial=start, final=end))
                        self.chrom_sizes_dict[chrom] = end
        elif not self.regionset:
            self.regionset = GenomicRegionSet('')
            print("Call DPs on whole genome.", file=sys.stderr)
            with open(chrom_file) as f:
                for line in f:
                    line = line.strip()
                    line = line.split('\t')
                    chrom, end = line[0], int(line[1])
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
        return len(self.valid_regionset)
    
    def __iter__(self):
        for el in self.valid_regionset:
            tmp = GenomicRegionSet('')
            tmp.add(el)
            yield tmp
        #return iter(self.regionset)
    
    def get_regionset(self):
        return self.valid_regionset
    
    def get_chrom_dict(self):
        return self.valid_chrom_sizes
    
    
    def get_training_regionset(self):
        r = GenomicRegionSet('')
        r.add(self.valid_regionset[self.counter])
        
        if self.counter == len(self.valid_chrom_sizes):
            return None
        else:
            self.counter += 1
            return r

    def update_regions(self, signal_statics):
        """this function update regionset from signal_statics and only consider the valid regionset"""
        valid_chroms = set()
        self.valid_regionset = GenomicRegionSet('valid_regionset')
        self.valid_chrom_sizes = {}
        for i in range(signal_statics['dim'][0]):
            for j in range(signal_statics['dim'][1]):
                valid_chroms |= set(tmp[0] for tmp in  signal_statics['data'][i][j]['stats_data'])

        for region in self.regionset:
            if region.chrom in valid_chroms:
                self.valid_regionset.add(GenomicRegion(chrom=region.chrom, initial=region.initial, final=region.final))
                self.valid_chrom_sizes[region.chrom] = region.final
