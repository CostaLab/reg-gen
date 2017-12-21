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
import sys
import os
from scipy.stats.mstats import zscore
import numpy as np
import matplotlib.pyplot as plt

from rgt.Util import npath
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet

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


def merge_delete(ext_size, merge, peak_list):
    regions_plus = GenomicRegionSet('regions_plus') #pot. mergeable
    regions_minus = GenomicRegionSet('regions_minus') #pot. mergeable

    for i, r in enumerate(peak_list):
        if r.final - r.initial > ext_size:
            if r.orientation == '+':
                regions_plus.add(r)
            elif r.orientation == '-':
                regions_minus.add(r)

    if merge:
        regions_plus.extend(ext_size/2, ext_size/2)
        regions_plus.merge()
        regions_plus.extend(-ext_size/2, -ext_size/2)
        merge_data(regions_plus)

        regions_minus.extend(ext_size/2, ext_size/2)
        regions_minus.merge()
        regions_minus.extend(-ext_size/2, -ext_size/2)
        merge_data(regions_minus)
    """no idea what this means"""
    results = GenomicRegionSet('regions')
    for el in regions_plus:
        results.add(el)
    for el in regions_minus:
        results.add(el)
    results.sort()

    return results


def separate_peaks(output, filter):
    """accept both output and pvalues, maybe also filter; not clean codes"""
    gain_peaks, lose_peaks = GenomicRegionSet("gain"), GenomicRegionSet("lose")
    gain_filter, lose_filter = [], []
    for i in range(len(filter)):  ## add name parameter here
        strand = output[i].orientation
        if filter[i]:
            if strand == '+':# cause output is not a region, so we can't output
                gain_peaks.add(output[i])
                gain_filter.append(filter[i])
            elif strand == '-':
                lose_peaks.add(output[i])
                lose_filter.append(filter[i])
    return [gain_peaks, gain_filter],  [lose_peaks, lose_filter]


def add_gene_name(peaks_set, organism_name):
    """add genes names to output peaks"""
    peaks = peaks_set.gene_association(organism=organism_name, promoterLength=1000,
                                         threshDist=500000, show_dis=True)
    return peaks


def get_peak_column(peaks_set, col_name, in_data=True):
    """extract one column data from peaks"""
    cols = []
    for peak in peaks_set:
        if in_data:
            cols.append(peak.data[col_name])
        else:
            cols.append(peak.col_name)
    return cols


def filter_by_pvalue_strand_lag(peaks_set, pcutoff, no_correction, name, singlestrand, separate=False):
    """Filter DPs by strang lag and pvalue,
    singlestrand: Allow single strand BAM file as input."""
    # first make ratios and pvalues to be extracted from peaks_set
    ratios = get_peak_column(peaks_set, col_name='ratio')
    pvalues = get_peak_column(peaks_set, col_name='pvalue')

    if not no_correction:
        pv_pass = [True] * len(pvalues)
        # we only create pvalues and then
        pvalues = map(lambda x: 10**-x, pvalues)
        
        output_peaks(name + '-uncor', peaks_set, pv_pass, separate=separate)
        
        pv_pass, pvalues = multiple_test_correction(pvalues, alpha=pcutoff)
    else:
        pv_pass = np.where(np.asarray(pvalues) >= -np.log10(pcutoff), True, False)
    
    if not singlestrand:
        zscore_ratios = zscore(ratios)
        ratios_pass = np.where(np.bitwise_and(zscore_ratios > -2, zscore_ratios < 2) == True, True, False)
        filter_pass = np.bitwise_and(ratios_pass, pv_pass)
    else:
        filter_pass = pv_pass
    
    return peaks_set, filter_pass


def transfrom_2string(dict_data, mask_col, item_sep = ';',data_sep=':'):
    """change data in dictionary into string, return string"""
    result_str= ''
    for keys, values in dict_data.items():
        if keys in mask_col:# without no inputs
            continue
        if isinstance(values, list):
            result_str += data_sep.join(map(str,values)) + item_sep
        else:
            result_str += str(values) + item_sep
    return result_str.rstrip(item_sep)


def output_BED(name, output, filter):

    f = open(name, 'w')
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000

    for i in range(len(filter)):
        if filter[i]:
            c, s, e, name, state, counts = output[i].chrom, output[i].initial, output[i].final,  output[i].name, output[i].orientation, output[i].data
            counts['pvalue'] = -np.log10(counts['pvalue']) if counts['pvalue'] > 0 else sys.maxint
            counts_str = transfrom_2string(counts, mask_col=['ratio'])
            print(c, s, e, 'Peak' + str(i) + '_' + name, bedscore, state, s, e, colors[state], 0, counts_str, sep='\t', file=f)

    f.close()


def output_peaks(name, output, filter, output_narrow=True, separate=False):
    """output peaks files.
    output_narrow=True -- Also output narrow peaks
    separate_peaks=True -- it will output peaks into gain and lose files
    What we should do ?? Draw a graph to see relationship of ratios and pvalues, and output
    If we all use the same data, we could get it easily, right ??
    """
    if not separate:
        output_BED(name + '_diffpeaks', output, filter)
        if output_narrow:
            output_narrowPeak(name + '-diffpeaks.narrow', output, filter)
    else:
        # first to create gain or lose data list
        gains, loses = separate_peaks(output,  filter)
        # pass gain / lose list to output_BED
        output_BED(name + '_diffpeaks_gain', gains[0], gains[1])
        output_BED(name + '_diffpeaks_lose', loses[0], loses[1])
        if output_narrow:
            output_narrowPeak(name + '_diffpeaks_gain.narrow', gains[0], gains[1])
            output_narrowPeak(name + '_diffpeaks_lose.narrow', loses[0], loses[1])


def output_narrowPeak(name, output,filter):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name, 'w')
    for i in range(len(filter)):
        if filter[i]:
            c, s, e, name, state, counts = output[i].chrom, output[i].initial, output[i].final, output[i].name, output[
                i].orientation, output[i].data
            counts['pvalue'] = -np.log10(counts['pvalue']) if counts['pvalue'] > 0 else sys.maxint
            print(c, s, e, 'Peak' + str(i) + '_' + name, 0, state, 0, counts['pvalue'], 0, -1, sep='\t', file=f)
    f.close()


def merge_bw_output(name, dim, no_bw_files, chrom_sizes):
    for i in range(dim[0]):
        for j in range(dim[1]):
            ## here are to output bed files for each signal file and output files
            temp_bed = npath(name + '-s%s-rep%s_temp.bed'% (i, j))

            files = [name + '-' + str(num) + '-s%s-rep%s.bw'%(i, j) for num in no_bw_files]
            if len(files) > 1:
                files = filter(lambda x: os.path.isfile(x), files)
                t = ['bigWigMerge'] + files + [temp_bed]
                c = " ".join(t)
                os.system(c)

                os.system("LC_COLLATE=C sort -k1,1 -k2,2n " + temp_bed + ' > ' + temp_bed +'.sort')

                t = ['bedGraphToBigWig', temp_bed + '.sort', chrom_sizes, name + '-s%s-rep%s.bw' % (i, j)]
                c = " ".join(t)
                os.system(c)

                for f in files:
                    os.remove(f)
                os.remove(temp_bed)
                os.remove(temp_bed + ".sort")
            else:
                ftarget = [name + '-s%s-rep%s.bw' %(i, j)]
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
    """
    ext_size1 = int(sys.argv[1]) #100
    ext_size2 = int(sys.argv[2]) #100
    path = sys.argv[3] # '/home/manuel/merge_test.data'
    merge = sys.argv[4] #True #for histones
    
    regions = merge_delete(path, ext_size1, '+', merge)
    #regions_minus = merge_delete(path, ext_size2, '-', merge)
    
    i = 0
    """
    counts = {'cov1':[6,9], 'cov2':[2,5], 'pvalue':3.6, 'ratio':5}
    result= transfrom_2string(counts,mask_col=['ratio'])
    print(result)

