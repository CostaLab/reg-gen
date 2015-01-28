#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Created on Mar 6, 2013

@author: manuel
'''
from __future__ import print_function
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
import sys
import re

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
    color = {'+': '255,0,0', '-': '0,255,0'}
    output = []
    
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        #counts = ",".join(map(lambda x: re.sub("\D", "", x), tmp[:len(tmp)-1]))
        main_sep = ':' #sep <counts> main_sep <counts> main_sep <pvalue>
        int_sep = ';' #sep counts in <counts>
        counts = ",".join(tmp).replace('], [', ';').replace('], ', int_sep).replace('([', '').replace(')', '').replace(', ', main_sep)
        pvalue = float(tmp[len(tmp)-1].replace(")", "").strip())
        
        output.append("\t".join([el.chrom, el.initial, el.final, 'Peak'+str(i), 1000, el.orientation, el.initial, el.final, \
              color[el.orientation], 0, counts]), pvalue)
    
    return output

def merge_delete(ext_size, merge, peak_list, pvalue_list):
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
    
    return results

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
    
    
