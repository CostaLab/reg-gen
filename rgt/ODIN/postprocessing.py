#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ODIN detects differential peaks in multiple ChIP-seq profiles associated
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
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
import sys

def merge_data(regions):
    for el in regions:
        tmp = el.data
        d = tmp.split('_$_')
        c1, c2 = 0, 0
        logpvalue = []
        for tmp in d:
            tmp = tmp.split(',')
            c1 += int(tmp[0].replace("(", ""))
            c2 += int(tmp[1])
            logpvalue.append(float(tmp[2].replace(")", "")))

        el.data = str((c1, c2, max(logpvalue)))


def output(name, regions):
    f = open(name + '-diffpeaks.bed', 'w')
    color = {'+': '255,0,0', '-': '0,255,0'}
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        c1 = int(tmp[0].replace("(", ""))
        c2 = int(tmp[1])
        logpvalue = float(tmp[2].replace(")", ""))
        
        print(el.chrom, el.initial, el.final, 'Peak'+str(i), 1000, el.orientation, el.initial, el.final, \
              color[el.orientation], 0, str(c1)+','+str(c2)+','+str(logpvalue), sep='\t', file=f)
    f.close()
    

def merge_delete(ext_size, merge, peak_list, pvalue_list, name):
#     peaks_gain = read_diffpeaks(path)
    
    regions_plus = GenomicRegionSet('regions') #pot. mergeable
    regions_minus = GenomicRegionSet('regions') #pot. mergeable
    regions_unmergable = GenomicRegionSet('regions')
    last_orientation = ""
    
    for i, t in enumerate(peak_list):
        chrom, start, end, c1, c2, strand = t[0], t[1], t[2], t[3], t[4], t[5] 
        r = GenomicRegion(chrom = chrom, initial = start, final = end, name = '', \
                          orientation = strand, data = str((c1, c2, pvalue_list[i])))
        if end - start > ext_size:
            if strand == '+':
                if last_orientation == '+':
                    region_plus.add(r)
                else:
                    regions_unmergable.add(r)
            elif strand == '-':
                if last_orientation == '-':
                    region_mins.add(r)
                else:
                    regions_unmergable.add(r)
                    
                    
    if merge:
        regions_plus.extend(ext_size/2, ext_size/2)
        regions_plus.merge()
        regions_plus.extend(-ext_size/2, -ext_size/2)
        merge_data(regions_plus)
        
        regions_minus.extend(ext_size/2, ext_size/2)
        regions_minus.merge()
        regions_minus.extend(-ext_size/2, -ext_size/2)
        merge_data(regions_minus)
    
    results = GenomicRegionSet('regions')
    for el in regions_plus:
        results.add(el)
    for el in regions_minus:
        results.add(el)
    for el in regions_unmergable:
        results.add(el)
    results.sort()
    output(name, results)

if __name__ == '__main__':
    ext_size1 = int(sys.argv[1]) #100
    ext_size2 = int(sys.argv[2]) #100
    path = sys.argv[3] # '/home/manuel/merge_test.data'
    merge = sys.argv[4] #True #for histones
    
    regions = merge_delete(path, ext_size1, '+', merge)
    #regions_minus = merge_delete(path, ext_size2, '-', merge)
    
    i = 0
    i = output(regions, i)
    #i = output(regions_minus, i, '-')
    
    