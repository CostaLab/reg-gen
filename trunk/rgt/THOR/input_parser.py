#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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
    if len(data) == 1:
        return data[0]
    else:
        return data

def input_parser(filepath):
    dim_1, dim_2 = 0, 0
    
    bamfiles_1 = get_data_block(filepath, "rep1")
    bamfiles_2 = get_data_block(filepath, "rep2")
    regions = get_data_block(filepath, "regions")
    genome = get_data_block(filepath, "genome")
    chrom_sizes = get_data_block(filepath, "chrom_sizes")
    inputs1 = get_data_block(filepath, "inputs1")
    inputs2 = get_data_block(filepath, "inputs2")
        
    dims = [len(bamfiles_1), len(bamfiles_2)]
    
    if not inputs:
        inputs = None
    if not regions:
        regions = None
    
    return bamfiles_1 + bamfiles_2, regions, genome, chrom_sizes, inputs1 + inputs2, dims

if __name__ == '__main__':
    print(input_parser('/home/manuel/workspace/eclipse/office_share/blueprint/playground/input_test'))


