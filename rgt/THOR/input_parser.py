#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

def get_data_block(filepath, feature):
    with open(filepath) as file:
        data = []
        read = False
        for line in file:
            line = line.strip()
            if line.startswith("#") and line == "#" + str(feature):
                read = True
            
            if line.startswith("#") and line != "#" + str(feature):
                read = False
            
            if not line.startswith("#") and read:
                data.append(line)
                
    if len(data) == 1 and not (feature == "rep1" or feature == "rep2" or feature == "inputs1" or feature == "inputs2"):
        return data[0]
    else:
        return data

def input_parser(filepath):
    dim_1, dim_2 = 0, 0
    
    bamfiles_1 = get_data_block(filepath, "rep1")
    bamfiles_2 = get_data_block(filepath, "rep2")
    #regions = get_data_block(filepath, "regions")
    genome = get_data_block(filepath, "genome")
    chrom_sizes = get_data_block(filepath, "chrom_sizes")
    inputs1 = get_data_block(filepath, "inputs1")
    inputs2 = get_data_block(filepath, "inputs2")
        
    dims = [len(bamfiles_1), len(bamfiles_2)]
    
    if not inputs1 and not inputs2:
        inputs = None
    else:
        inputs = inputs1 + inputs2
    if not genome:
        genome = None
    #if not regions:
    #    regions = None
    
    return bamfiles_1 + bamfiles_2, genome, chrom_sizes, inputs, dims

if __name__ == '__main__':
    print(input_parser('/home/manuel/blueprint.h3k4me1.cluster.config')[1])


