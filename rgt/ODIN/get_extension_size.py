#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Estimate extension size of reads for an ChIP-seq experiment.

Author: Manuel Allhoff

Calculate the cross correlation between coverage of forward and backward reads
on first chromosome.
Methods is based on:

Kharchenko et al., Design and analysis of ChIP-seq experiments for DNA-binding
proteins (2008)

Methods:

get_extension_size(bamfile)
Return shift/extension size of reads descriebed by BAM file.

"""

from __future__ import print_function
import pysam, sys

cov_f = {}
cov_r = {}

def get_value(d, x):
    """Return d[x] of dictionary d or 0"""
    return d[x] if d.has_key(x) else 0 

def init_cov(filename):
    file = pysam.Samfile(filename, "rb")
    
    for read in file.fetch(file.references[0]):
        if not read.is_unmapped:
            if not read.seq:
                h = 0
            else:
                h = len(read.seq)
            pos = read.pos + read.rlen - h if read.is_reverse else read.pos
            if read.is_reverse:
                if not cov_r.has_key(pos):
                    cov_r[pos] = 1
            else:
                if not cov_f.has_key(pos):
                    cov_f[pos] = 1

def ccf(k):
    """Return value of cross-correlation function"""
    s = 0
    forward_keys = set(cov_f.keys())
    reverse_keys = set(map(lambda x: x - k, cov_r.keys()))
    keys = forward_keys & reverse_keys
    
    for p in keys:
        s += get_value(cov_f, p) & get_value(cov_r, p + k)
    
    return s, k

def get_extension_size(filename, start=0, end=600, stepsize=5):
    """Return extension/shift size of reads and all computed values of the convolution. 
    Search value with a resolution of <stepsize> from start to end."""
    init_cov(filename)
    
    r = map(ccf, range(start, end, stepsize) )
    
    return max(r[1:])[1], r[1:]


if __name__ == '__main__':
    #a, b = get_extension_size('/home/manuel/workspace/cluster_p/blueprint/raw/input/C000S5H1.Input.bwa_filtered.20130415.bam')
    a, b = get_extension_size('/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/HEL_hg_Rux_H3K9ac_rep1.bam')
    print(a, b)
    
    
    
    