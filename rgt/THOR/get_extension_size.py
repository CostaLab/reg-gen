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

import math
import pysam

cov_f = {}
cov_r = {}


def get_value(d, x):
    """Return d[x] of dictionary d or 0"""
    return d[x] if x in d else 0


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

    for read in file.fetch(file.references[0]):
        if not read.is_unmapped:
            if not read.seq:
                h = 0
            else:
                h = len(read.seq)
            pos = read.pos + read.rlen - h if read.is_reverse else read.pos
            if read.is_reverse:
                if pos not in cov_r:
                    cov_r[pos] = 1
            else:
                if pos not in cov_f:
                    cov_f[pos] = 1


def ccf(k):
    """Return value of cross-correlation function"""
    s = 0
    forward_keys = set(cov_f.keys())
    reverse_keys = set([x - k for x in list(cov_r.keys())])
    keys = forward_keys & reverse_keys

    for p in keys:
        s += get_value(cov_f, p) & get_value(cov_r, p + k)

    return s, k


def get_extension_size(filename, start=0, end=600, stepsize=5):
    """Return extension/shift size of reads and all computed values of the convolution. 
    Search value with a resolution of <stepsize> from start to end."""
    read_length = math.ceil(get_read_size(filename))
    start -= read_length

    init_cov(filename)

    r = list(map(ccf, list(range(start, end, stepsize))))

    r = [(x[0], x[1]) for x in r]

    # print('extension size is %s' % max(r[read_length//stepsize*2:])[1])

    return max(r[read_length // stepsize * 2:])[1], r


if __name__ == '__main__':
    # a, b = get_extension_size('/home/manuel/workspace/cluster_p/blueprint/raw/input/C000S5H1.Input.bwa_filtered.20130415.bam')
    a, b = get_extension_size('/home/manuel/workspace/cluster_p/allhoff/project_THOR/data/payton/FL1_H3K27ac_input.bam')
    print(a, b)

    for el in b:
        print(el[1], el[0], sep='\t')


