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


from ..GenomicRegion import GenomicRegion
from ..GenomicRegionSet import GenomicRegionSet
import sys
# import re
from scipy.stats.mstats import zscore
from ..motifanalysis.Statistics import multiple_test_correction
import numpy as np
from numpy import log10


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


def output(name, regions):
    color = {'+': '255,0,0', '-': '0,255,0'}
    output = []
    
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        #counts = ",".join(map(lambda x: re.sub("\D", "", x), tmp[:len(tmp)-1]))
        main_sep = ':' #sep <counts> main_sep <counts> main_sep <pvalue>
        int_sep = ';' #sep counts in <counts>
        counts = ",".join(tmp).replace('], [', ';').replace('], ', int_sep).replace('([', '').replace(')', '').replace(', ', main_sep)
        pvalue = float(tmp[len(tmp)-1].replace(")", "").strip())
        
        output.append("\t".join([el.chrom, el.initial, el.final, 'Peak' + str(i), 1000,
                                 el.orientation, el.initial, el.final, color[el.orientation], 0, counts, pvalue]))
    
    return output


def merge_delete(ext_size, merge, peak_list, pvalue_list):
#     peaks_gain = read_diffpeaks(path)
    
    regions_plus = GenomicRegionSet('regions') #pot. mergeable
    regions_minus = GenomicRegionSet('regions') #pot. mergeable
    regions_unmergable = GenomicRegionSet('regions')
    last_orientation = ""
    
    for i, t in enumerate(peak_list):
        chrom, start, end, c1, c2, strand, ratio = t[0], t[1], t[2], t[3], t[4], t[5], t[6]
        r = GenomicRegion(chrom=chrom, initial=int(start), final=int(end), name='',
                          orientation=strand, data=str((c1, c2, pvalue_list[i], ratio)))
        # if end - start > ext_size:
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


def filter_by_pvalue_strand_lag(ratios, pcutoff, logpvalues, output, no_correction, name, singlestrand):
    """Filter DPs by strang lag and pvalue"""

    # convert the logpvalues to original pvalues
    pvalues = [10 ** -x for x in logpvalues]

    # create filter using the pvalue cutoff
    pv_pass = np.where(np.asarray(pvalues) <= pcutoff, True, False)

    if not no_correction:
        # save the data with the uncorrected pvalues and without filtering by cutoff
        # (this is only provided to the user in case they desire to perform their own
        # correction and/or filtering)
        # NB: the two _output_XXX functions take original p-values
        pv_pass = [True] * len(pvalues)
        _output_BED(name + '-uncor', output, pvalues, pv_pass)
        _output_narrowPeak(name + '-uncor', output, pvalues, pv_pass)

        # apply multiple test correction (takes original p-values) and filter using the cutoff
        pv_pass, pvalues = multiple_test_correction(pvalues, alpha=pcutoff)
    
    if not singlestrand:
        zscore_ratios = zscore(ratios)
        ratios_pass = np.where(np.bitwise_and(zscore_ratios > -2, zscore_ratios < 2) == True, True, False)

        filter_pass = np.bitwise_and(ratios_pass, pv_pass)

        assert len(pv_pass) == len(ratios_pass)
    else:
        filter_pass = pv_pass
    
    assert len(output) == len(pvalues)
    assert len(filter_pass) == len(pvalues)

    # we return the filtered (and maybe corrected) pvalues in original format (NOT log-pvalues)

    return output, pvalues, filter_pass


def _output_BED(name, output, pvalues, filter):
    f = open(name + '-diffpeaks.bed', 'w')
     
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, strand, counts = output[i]
        p_tmp = -log10(pvalues[i]) if pvalues[i] > 0 else sys.maxsize
        counts = ';'.join(counts.split(';')[:2] + [str(p_tmp)])
        
        if filter[i]:
            print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, colors[strand], 0, counts, sep='\t', file=f)
    
    f.close()


def _output_narrowPeak(name, output, pvalues, filter):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name + '-diffpeaks.narrowPeak', 'w')
    for i in range(len(pvalues)):
        c, s, e, strand, _ = output[i]
        p_tmp = -log10(pvalues[i]) if pvalues[i] > 0 else sys.maxsize
        if filter[i]:
            print(c, s, e, 'Peak' + str(i), 0, strand, 0, p_tmp, 0, -1, sep='\t', file=f)
    f.close()


def filter_deadzones(bed_deadzones, peak_regions):
    """Filter by peaklist by deadzones"""
    deadzones = GenomicRegionSet('deadzones')
    deadzones.read(bed_deadzones)
    peak_regions = peak_regions.subtract(deadzones, whole_region=True)
    
    return peak_regions


if __name__ == '__main__':
    ext_size1 = int(sys.argv[1]) #100
    ext_size2 = int(sys.argv[2]) #100
    path = sys.argv[3] # '/home/manuel/merge_test.data'
    merge = sys.argv[4] #True #for histones
    
    regions = merge_delete(path, ext_size1, '+', merge)
    #regions_minus = merge_delete(path, ext_size2, '-', merge)
    
    i = 0
    i = output(regions, i)
    #i = output(regions_minus, i, '-')

