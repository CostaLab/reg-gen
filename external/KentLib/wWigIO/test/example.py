#!/usr/bin/python
#Last-modified: 04 Oct 2013 02:27:00 PM

#         Module/Scripts Description
# 
# Copyright (c) 2008 Yunfei Wang <Yunfei.Wang1@utdallas.edu>
# 
# This code is free software; you can redistribute it and/or modify it
# under the terms of the BSD License (see the file COPYING included with
# the distribution).
# 
# @status:  experimental
# @version: 1.0.0
# @author:  Yunfei Wang
# @contact: yfwang0405@gmail.com

# ------------------------------------
# python modules
# ------------------------------------

import os,sys
import string
import wWigIO

# ------------------------------------
# constants
# ------------------------------------

# ------------------------------------
# Misc functions
# ------------------------------------

# ------------------------------------
# Classes
# ------------------------------------

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    # get information from bigwig file
    wWigIO.open('test.bw')
    chroms = wWigIO.getChromSize('test.bw')
    wigs = wWigIO.getIntervals('test.bw', 'chr1', 10, 200)
    wWigIO.close('test.bw')
    print wigs

    # bigwig -> wig
    wWigIO.bigWigToWig('test.bw','test.wig')

    # write the chrom sizes into test.sizes
    with open('test.sizes','w') as fh:
        for chrom in chroms:
            print >>fh, chrom+"\t"+str(chroms[chrom])
    
    # wig -> bigwig
    wWigIO.wigToBigWig('test.wig','test.sizes','test2.bw')
