#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Manuel Allhoff

"""
from __future__ import print_function
import sys

def get_data_block(file):
    data = []
    for line in file:
        line = line.strip()
        if not line.startswith("#"):
            data.append(line)
        else:
            break
    
    if len(data) == 1:
        return data[0]
    else:
        return data

def input_parser(filepath):
    dim_1, dim_2 = 0, 0
    
    with open(filepath) as file:
        _ = file.readline() #read first header
        bamfiles_1 = get_data_block(file)
        bamfiles_2 = get_data_block(file)
        regions = get_data_block(file)
        genome = get_data_block(file)
        chrom_sizes = get_data_block(file)
        inputs = get_data_block(file)
        
    dims = [len(bamfiles_1), len(bamfiles_2)]
    
    if not inputs:
        inputs = None
    
    return bamfiles_1 + bamfiles_2, regions, genome, chrom_sizes, inputs, dims

if __name__ == '__main__':
    print(input_parser('/home/manuel/workspace/eclipse/office_share/blueprint/playground/input_test'))