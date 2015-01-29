#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: Manuel Allhoff
'''
from __future__ import print_function
import numpy as np
import pysam
import sys
from HTSeq import FastaReader

class help_content():
    def __init__(self):
        self.content = [[] for _ in range(101)]
        self.g_gc = []
        self.g = -1
        
    def add(self, d, c):
        self.content[int(c*100)].append(round(float(d),2))
        
    def _compute(self):
        for l in self.content:
            r = sum(l) / float(len(l)) if len(l) > 0 else 0
            self.g_gc.append(r)

        self.g = sum(self.g_gc) / float(len(self.g_gc))

    def _map(self, x):
        return self.g_gc[x]

def get_gc_context(stepsize, binsize, genome_path, cov_list, chrom_sizes_dict):
    chromosomes = []
    #get first chromosome, typically chr1
    #for s in FastaReader(genome_path):
    #    chromosomes.append(s.name)
    tmp = chrom_sizes_dict.keys()
    #tmp = map(lambda x: x.replace('chr',''), chromosomes)
    tmp.sort()
    genome_fasta = pysam.Fastafile(genome_path)
    genome = genome_fasta.fetch(reference = tmp[0])
    
    content = help_content()
    gc_content_cov = []

    for cov in cov_list:
        cur_gc_value = []
        
        for i in range(len(cov)):
            s = i * stepsize
            e = i * stepsize + binsize
            seq = genome[s : e+1]
            seq = seq.upper()
            count = seq.count("C") + seq.count("G")
            if len(seq) > 0:
                gc_content = float(count) / len(seq)
                #print(float(count), len(seq), float(count) / len(seq), cov[i], file=sys.stderr)
                content.add(cov[i], gc_content)
                cur_gc_value.append(int(gc_content*100))
            else:
#                print("Bins exceed genome (%s - %s), add pseudo counts" %(s, e), file=sys.stderr)
#                     content.add(cov[i], 0.5)
                cur_gc_value.append(0)
        
        gc_content_cov.append(cur_gc_value)
    
    content._compute()
    r = []
    for l in gc_content_cov:
        tmp = map(content._map, l)
        r.append(tmp)
    
    return r, content.g, content.g_gc
    