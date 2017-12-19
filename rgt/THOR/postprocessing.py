#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
THOR detects differential peaks in multiple ChIP-seq profiles associated
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
import os
from scipy.stats.mstats import zscore
import numpy as np
import matplotlib.pyplot as plt

from rgt.Util import npath
from rgt.motifanalysis.Statistics import multiple_test_correction

import configuration


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


def merge_delete(ext_size, merge, peak_list, pvalue_list):
    regions_plus = GenomicRegionSet('regions') #pot. mergeable
    regions_minus = GenomicRegionSet('regions') #pot. mergeable
    regions_unmergable = GenomicRegionSet('regions')
    last_orientation = ""
    
    for i, t in enumerate(peak_list):
        chrom, start, end, c1, c2, strand, ratio = t[0], t[1], t[2], t[3], t[4], t[5], t[6]
        # change c1, c2 into integer
        # c1 = map(int, c1)
        # c2 = map(int, c2)
        r = GenomicRegion(chrom = chrom, initial = start, final = end, name = '', \
                          orientation = strand, data = str((c1, c2, pvalue_list[i], ratio)))
        if end - start > ext_size:
            if strand == '+':
                if last_orientation == '+':
                    regions_plus.add(r)
                else:
                    regions_unmergable.add(r)
            elif strand == '-':
                if last_orientation == '-':
                    regions_minus.add(r)
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


def filter_by_pvalue_strand_lag(ratios, pcutoff, pvalues, output, no_correction, name, singlestrand):
    """Filter DPs by strang lag and pvalue"""
    if not singlestrand:
        zscore_ratios = zscore(ratios)
        ratios_pass = np.where(np.bitwise_and(zscore_ratios > -2, zscore_ratios < 2) == True, True, False)
    if not no_correction:
        pv_pass = [True] * len(pvalues)
        pvalues = map(lambda x: 10**-x, pvalues)
        
        output_BED(name + '-uncor', output, pvalues, pv_pass)
        output_narrowPeak(name + '-uncor', output, pvalues, pv_pass)
        
        pv_pass, pvalues = multiple_test_correction(pvalues, alpha=pcutoff)
    else:
        pv_pass = np.where(np.asarray(pvalues) >= -np.log10(pcutoff), True, False)
    
    if not singlestrand:
        filter_pass = np.bitwise_and(ratios_pass, pv_pass)
        assert len(pv_pass) == len(ratios_pass)
    else:
        filter_pass = pv_pass
    
    return output, pvalues, filter_pass


def output_BED(name, output, pvalues, filter):

    f = open(name, 'w')
     
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, strand, counts = output[i]
        p_tmp = -np.log10(pvalues[i]) if pvalues[i] > 0 else sys.maxint
        counts = ';'.join(counts.split(';')[:2] + [str(p_tmp)])
        
        if filter[i]:
            print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, colors[strand], 0, counts, sep='\t', file=f)
    
    f.close()

def output_peaks(name, output,pvalues, filter, output_narrow=True,separate_peaks=False):
    """output peaks files.
    output_narrow=True -- Also output narrow peaks
    separate_peaks=True -- it will output peaks into gain and lose files
    """
    if not separate_peaks:
        output_BED(name + '_diffpeaks', output, pvalues,filter)
        if output_narrow:
            output_narrowPeak(name + '-diffpeaks.narrow', output, pvalues, filter)
    else:
        # first to create gain or lose data list
        gains, loses = separate_peaks(output, pvalues, filter)
        # pass gain / lose list to output_BED
        output_BED(name + '_diffpeaks_gain', gains[0], gains[1], gains[2])
        output_BED(name + '_diffpeaks_lose', loses[0], loses[1], loses[2])
        if output_narrow:
            output_narrowPeak(name + '_diffpeaks_gain.narrow', gains[0], gains[1], gains[2])
            output_narrowPeak(name + '_diffpeaks_lose.narrow', loses[0], loses[1], loses[2])



def separate_peaks(output, pvalues, filter):
    """accept both output and pvalues, maybe also filter; not clean codes"""
    gain_peaks, lose_peaks = [], []
    gain_pvalues, lose_pvalues = [], []
    gain_filter, lose_filter = [], []
    for i in range(len(pvalues)):  ## add name parameter here
        strand = output[i].orientation
        output[i].orientation = '.'
        if filter[i]:
            if strand == '+':
                gain_peaks.append(output[i])
                gain_pvalues.append(pvalues[i])
                gain_filter.append(filter[i])
            elif strand == '-':
                lose_peaks.append(output[i])
                lose_pvalues.append(pvalues[i])
                lose_filter.append(filter[i])
    return [gain_peaks, gain_pvalues, gain_filter],  [lose_peaks, lose_pvalues, lose_filter]

def output_narrowPeak(name, output, pvalues, filter):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name, 'w')
    for i in range(len(pvalues)):
        c, s, e, strand, _ = output[i]
        p_tmp = -np.log10(pvalues[i]) if pvalues[i] > 0 else sys.maxint
        if filter[i]:
            print(c, s, e, 'Peak' + str(i), 0, strand, 0, p_tmp, 0, -1, sep='\t', file=f)
    f.close()


def merge_bw_output(signal_statics, options, no_bw_files, chrom_sizes):
    dim = signal_statics['dim']

    for i in range(dim[0]):
        for j in range(dim[1]):
            ## here are to output bed files for each signal file and output files
            temp_bed = npath(options.name + '-s%s-rep%s_temp.bed'% (i, j))

            files = [options.name + '-' + str(num) + '-s%s-rep%s.bw'%(i, j) for num in no_bw_files]
            if len(files) > 1:
                files = filter(lambda x: os.path.isfile(x), files)
                t = ['bigWigMerge'] + files + [temp_bed]
                c = " ".join(t)
                os.system(c)

                os.system("LC_COLLATE=C sort -k1,1 -k2,2n " + temp_bed + ' > ' + temp_bed +'.sort')

                t = ['bedGraphToBigWig', temp_bed + '.sort', chrom_sizes, options.name + '-s%s-rep%s.bw' % (i, j)]
                c = " ".join(t)
                os.system(c)

                for f in files:
                    os.remove(f)
                os.remove(temp_bed)
                os.remove(temp_bed + ".sort")
            else:
                ftarget = [options.name + '-s%s-rep%s.bw' %(i, j)]
                for k in range(len(ftarget)):
                    c = ['mv', files[k], ftarget[k]]
                    c = " ".join(c)
                    os.system(c)


def _output_ext_data(ext_data_list, bamfiles):
    """Output textfile and png file of read size estimation"""
    names = [os.path.splitext(os.path.basename(bamfile))[0] for bamfile in bamfiles]

    for k, ext_data in enumerate(ext_data_list):
        f = open(configuration.FOLDER_REPORT_DATA + 'fragment_size_estimate_' + names[k] + '.data', 'w')
        for d in ext_data:
            print(d[0], d[1], sep='\t', file=f)
        f.close()

    for i, ext_data in enumerate(ext_data_list):
        d1 = map(lambda x: x[0], ext_data)
        d2 = map(lambda x: x[1], ext_data)
        ax = plt.subplot(111)
        plt.xlabel('shift')
        plt.ylabel('convolution')
        plt.title('Fragment Size Estimation')
        plt.plot(d2, d1, label=names[i])

    ax.legend()
    plt.savefig(configuration.FOLDER_REPORT_PICS + 'fragment_size_estimate.png')
    plt.close()


if __name__ == '__main__':
    ext_size1 = int(sys.argv[1]) #100
    ext_size2 = int(sys.argv[2]) #100
    path = sys.argv[3] # '/home/manuel/merge_test.data'
    merge = sys.argv[4] #True #for histones
    
    regions = merge_delete(path, ext_size1, '+', merge)
    #regions_minus = merge_delete(path, ext_size2, '-', merge)
    
    i = 0

