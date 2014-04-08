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
import multiprocessing


cov_f = {}
cov_r = {}

def get_value(d, x):
    """Return d[x] of dictionary d or 0"""
    return d[x] if d.has_key(x) else 0 

# def init_cov(filename):
#     """Initialize coverage data, consider only first chromosome"""
#     file = pysam.Samfile(filename, "rb")
#     i = 1
#     for column in file.pileup(file.references[0]):
#         if i % 10000000 == 0:
#             print("%sm pileups for initializing considered" % (i/1000000), file=sys.stderr)
# #             break
#         i += 1
#         pos = column.pos
#         read_strands = [ read.alignment.is_reverse for read in column.pileups ]
#         cov_f[pos] = read_strands.count(False)
#         cov_r[pos] = read_strands.count(True)

def init_cov(filename):
    file = pysam.Samfile(filename, "rb")
    
    for read in file.fetch(file.references[0]):
        if not read.is_unmapped:
            pos = read.pos + read.rlen - len(read.seq) if read.is_reverse else read.pos
            if read.is_reverse:
                if not cov_r.has_key(pos):
                    cov_r[pos] = 1
            else:
                if not cov_f.has_key(pos):
                    cov_f[pos] = 1

def ccf(k):
    """Return value of cross-correlation function"""
    s = 0
#     print(k, file=sys.stderr)
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
    
    #pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()/2)
    r = map(ccf, range(start, end, stepsize) )
    
    return max(r)[1], r

if __name__ == '__main__':
    filename = sys.argv[1]
#    filename='/home/manuel/data/project_chipseq_norm/data/PU1_CDP_1000000.bam'
    
    init_cov(filename)
    #pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()/2)
    r = map(ccf, range(0,int(sys.argv[2]),5) )

    for i,k in r:
        print(i,k)
    