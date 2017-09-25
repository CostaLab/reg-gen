#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
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

@author: Manuel Allhoff

"""

from __future__ import print_function
import pysam
import numpy as np
from rgt.THOR.RegionGiver import RegionGiver

cov_f = {}
cov_r = {}


def get_value(d, x):
    """Return d[x] of dictionary d or 0"""
    return d[x] if d.has_key(x) else 0


def get_read_size(filename):
    f = pysam.Samfile(filename, "rb")

    s = []
    i = 0
    for read in f.fetch():
        i += 1
        if i == 1000:
            break
        if not read.is_unmapped:
            s.append(read.rlen)

    # print('read length is %s' %str(sum(s)/len(s)))

    return sum(s) / len(s)


def init_cov(filename):
    file = pysam.Samfile(filename, "rb")
    # error happens here, if we choose the first one
    ## file.references is ('chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',....] all chromosomes in files,
    # so we need to choose them according to different reference part, which is not right.
    # if extension size is for each bam files, then we need to get all the data, so at least 1000 data to simulate them..
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

    if not cov_f and not cov_r:
        return False
    else:
        return True


def new_init(bam_filename, chrom_file_name):
    """:return boolean for reading cov_f and cov_r for each bam files; True, if cov_f reading right, else False
        and the read length for each bamfiles
       Now there are some differences between old and new init.. I would use the old firstly and see how it's.
    """
    f = pysam.Samfile(bam_filename, "rb")

    region_giver = RegionGiver(chrom_file_name, None)
    chroms = np.intersect1d(region_giver.get_chrom_dict().keys(), f.references)

    s = []
    #chrom_sizes = [5000] * len(f.references)
    # here we only need to consider chromosomes in chromosomes files, so we need to change it..
    for chrom_idx in range(len(chroms)):
        i = 0
        print(chroms[chrom_idx])
        for read in f.fetch(chroms[chrom_idx]):
            i += 1
            #if i == chrom_sizes[chrom_idx]:
            #    break
            if not read.is_unmapped:
                s.append(read.rlen)
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
        print(i)
    if not cov_f and not cov_r:
        return False, -1
    else:
        read_len = sum(s) / len(s)
        return True, read_len


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
    read_length = get_read_size(filename)
    start -= read_length

    init_cov(filename)

    # if init_cov is empty, then we return empty.
    if not init_cov(filename):
        return None, None


    r = map(ccf, range(start, end, stepsize))

    r = map(lambda x: (x[0], x[1]), r)

    # print('extension size is %s' %max(r[read_length/stepsize*2:])[1])

    return max(r[read_length / stepsize * 2:])[1], r

def new_get_extension_size(bam_filename, chrom_file_name, start=0, end=600, stepsize=5):

    read_suc, read_length = new_init(bam_filename, chrom_file_name)
    print(read_length)
    start = max(0, start - read_length)


    r = map(ccf, range(start, end, stepsize))

    return max(r)[1], r


if __name__ == '__main__':
    # a, b = get_extension_size('/home/manuel/workspace/cluster_p/blueprint/raw/input/C000S5H1.Input.bwa_filtered.20130415.bam')
    a, b = new_get_extension_size('/home/kefang/programs/THOR_example_data/bug/test/FL8_H3K27ac.100k.bam','/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes')
    print(a, b)

    #for el in b:
    #    print(el[1], el[0], sep='\t')


