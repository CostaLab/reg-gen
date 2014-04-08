#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
@author: manuel
'''
from __future__ import print_function
import numpy as np
import pysam
import sys

class help_content():
    def __init__(self):
        self.content = [[] for _ in range(100)]
        self.g_gc = []
        self.g = -1
        
    def add(self, d, c):
#        print(c*100, round(float(d),2), file=sys.stderr)
        self.content[int(c*100)].append(round(float(d),2))
        
    def _compute(self):
        for l in self.content:
            r = sum(l) / float(len(l)) if len(l) > 0 else 0
            self.g_gc.append(r)

        self.g = sum(self.g_gc) / float(len(self.g_gc))

    def _map(self, x):
        return self.g_gc[x]

def get_gc_context(stepsize, binsize, genome_path, cov_list):
    genome_fasta = pysam.Fastafile(genome_path)
    genome = genome_fasta.fetch(reference="chr1") #TODO: only chr1
    content = help_content()
    gc_content_cov = []

    for cov in cov_list:
        cur_gc_value = []
        
        for i in range(len(cov)):
#             if cov[i] < 1e-300:
#                 cur_gc_value.append(0)
#                 continue
            s = i * stepsize
            e = i * stepsize + binsize
            seq = genome[s : e+1]
            seq = seq.upper()
            count = seq.count("C") + seq.count("G")
            if len(seq) > 0:
                gc_content = float(count) / len(seq)
                content.add(cov[i], gc_content)
                cur_gc_value.append(int(gc_content*100))
            else:
#                print("Bins exceed genome (%s - %s), add pseudo counts" %(s, e), file=sys.stderr)
#                     content.add(cov[i], 0.5)
                cur_gc_value.append(0)
        
        gc_content_cov.append(cur_gc_value)
    
#    gc_content_cov = np.matrix(gc_content_cov)
    
    content._compute()
    r = []
    for l in gc_content_cov:
        tmp = map(content._map, l)
        r.append(tmp)
    
    return r, content.g, content.g_gc
    
if __name__ == '__main__':
    
    stepsize=50
    binsize=100
    genome_path='/home/manuel/g.fa'
    c = np.random.uniform(0, 10, size=10000/stepsize)
    c = np.around(c, 2)
    a,b=get_gc_context(stepsize, binsize, genome_path, [c])
#     print(c)
#     print(b)
#     print(len(a[0]))
#     print(a[0])
    
    for i in a[0]:
        print(i)
    
    
    
    