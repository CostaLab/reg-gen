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
                
    if len(data) == 1 and not (feature == "rep1" or feature == "rep2" or feature == "inputs1" or feature == "inputs2"):
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
    
    if not inputs1 and not inputs2:
        inputs = None
    else:
        inputs = inputs1 + inputs2
    if not regions:
        regions = None
    
    return bamfiles_1 + bamfiles_2, regions, genome, chrom_sizes, inputs, dims

if __name__ == '__main__':
    print(input_parser('/home/manuel/workspace/eclipse/office_share/blueprint/playground/blueprint.h3k4me1.cluster.config'))


