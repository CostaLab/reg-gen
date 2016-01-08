#!/usr/bin/python
#Last-modified: 03 Oct 2013 09:52:14 PM

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

import sys
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

class BigWigFile(object):
    '''  
    Fast reader of BigWig format file.
    Usage:
        Open file:
            fh=BigWigFile("test.bw")
        Get chromosome sizes:
            chroms=fh.chromSizes() # {key:value = chrom:size}
        Fetch regions:
            wigs=fh.fetch(chrom="chr1",start=100,stop=200)
            for start,stop,depth in wigs:
                #do some thing with wig
                print start,stop,depth
        Close file:
            fh.close()
    
    Parameters:
        chrom=None: return empty list.
        start=None: start at first position.
        stop=None:  stop at the end of chromosome.
    '''
    def __init__(self,fname):
        ''' Open BigWig file. '''
        self.fname=fname
        wWigIO.open(self.fname)
    def chromSizes(self):
        ''' Get chromosome sizes.'''
        chroms=wWigIO.getChromSize(self.fname)
        return chroms
    def fetch(self,**kwargs):
        ''' Fetch intervals in a given region. '''
        if kwargs.has_key('chrom'):
            wigs=wWigIO.getIntervals(self.fname,kwargs['chrom'],kwargs.get('start',0),kwargs.get('stop',0))
            if isinstance(wigs,basestring): # bad value
                raise ValueError("Couldn't get intervals.")
            return wigs
        raise ValueError("Chromosome not provided.")
        return
    def close(self):
        ''' Close BigWig file. '''
        wWigIO.close(self.fname)
    def __del__(self): 
        ''' Close BigWig file.  Avoid memory leaks.'''
        wWigIO.close(self.fname)

# ------------------------------------
# Main
# ------------------------------------

if __name__=="__main__":
    fh=BigWigFile("test.bw")
    chroms=fh.chromSizes() # {key:value = chrom:size}
    print chroms
    wigs=fh.fetch(chrom="chr1",start=100,stop=200)
    for start,stop,depth in wigs:
        #do some thing with wig
        print start,stop,depth

