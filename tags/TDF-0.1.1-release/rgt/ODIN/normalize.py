#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
%prog <BED,BAM> <BED,BAM> <NAME> <STEP_WIDTH> <Zero-counts 1/0> <hg19/mm9>

Computes p-q list, k and a from Diaz et al., 2012.

Features must have same length in both input files.
'''
from __future__ import print_function
from optparse import OptionParser
from HTSeq import GenomicPosition, GenomicArray, GenomicInterval
from math import fabs
import pysam, sys, operator, os.path
from itertools import chain
import cProfile
from time import time

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def _accumulate(iterable):
    """Return running totals"""
    it = iter(iterable)
    total = next(it)
    yield [total[0], total[1]]
    for e1, e2 in it:
        total[0] = total[0] + e1
        total[1] = total[1] + e2
        yield [total[0], total[1]]
              
def get_feature_len(path):
    """Return length of reads"""
    filename, fileextension = os.path.splitext(path)
    
    if fileextension == '.bed':
        with open(path, 'r') as f:
            for line in f:
                tmp = line.split('\t')
                return int(tmp[2]) - int(tmp[1])
    elif fileextension == '.bam':
        samfile = pysam.Samfile(path, "rb")
        for read in samfile.fetch():
            if not read.is_unmapped:
                return read.alen
            
def _get_read_info(path):
    """Yield tuple (chr, pos) for different fileformats."""
    filename, fileextension = os.path.splitext(path)
    
    if fileextension == '.bed':
        with open(path, 'r') as f:
            for line in f:
                tmp = line.split('\t')
                yield tmp[0], int(tmp[1])
    elif fileextension == '.bam':
        samfile = pysam.Samfile(path, "rb")
        for read in samfile.fetch():
            if not read.is_unmapped:
                yield samfile.getrname(read.tid), read.pos

def get_count_list(CHROM_LEN, path, stop=False):
    """Compute list of read's starting positions based on HTSeq's genomic array."""
    genomic_array = GenomicArray(CHROM_LEN, stranded=False, typecode='i', storage='step')
    i = 0
    chromosomes = set()
    for chrom, pos in _get_read_info(path):
        i += 1
        
        if stop and i == 5000000:
            break
        
        if not chrom.startswith("chr"):
            chrom = "chr" + chrom

        if chrom not in chromosomes:
            chromosomes.add(chrom)

        #ignore reads that fall out of chromosome borders, should not be necassary!
        if chrom not in CHROM_LEN.keys() or pos >= CHROM_LEN[chrom] or pos < 0:
            #print("illegal read on chromosome %s at positions %s not"%(chrom, pos), file=sys.stderr)
            continue
        
        genomic_array[ GenomicPosition(chrom, pos) ] += 1
        
    #print('count_list', i, file=sys.stderr)
    
    return genomic_array, chromosomes

def _get_overrun(chrom, i, end, step_width, count_list, feature_len):
    """Return overrun of reads that fall in two bins"""
    help_c1 = filter( lambda x: x[0].start + feature_len > end and x[1] is not 0, list(count_list[ GenomicInterval(chrom, i, end) ].steps()) )
    overrun = 0 if not help_c1 else sum(map(lambda x: x[1], help_c1))
    
    return overrun

def write_bedgraph(chrom_len, chromosomes, countstable, filename, step_width):
    """Write countstable to file"""
    with open(filename, 'w') as f:
        for chrom in chromosomes:
            if not chrom_len.has_key(chrom):
                print("Warning: %s not found, ignore" %chrom, file=sys.stderr)
            else:
                assert chrom_len[chrom]/step_width+1 == len(countstable[chrom]) #maybe errorprone if chrom len is 100000
                for i in range(0, chrom_len[chrom]/step_width+1):
                    start = i * step_width
                    end = min( (i+1) * step_width, chrom_len[chrom] )
                    print(chrom, start, end, countstable[chrom][i], sep='\t', file=f)

def write_pq_list(pq_list, max_index, max_value, factor1, factor2, filename):
    """Write p,q-list to file"""
    if pq_list:
        with open(filename, 'w') as f:
            print('#max index', 'max value', 'factor1', 'factor2', sep='\t', file=f)
            print('#' + str(max_index), str(max_value), str(factor1), str(factor2), sep='\t', file=f)
            for p, q in pq_list:
                print(p, q, file=f)
 
def get_bins(chrom_len, chromosomes, count_list, step_width, feature_len):
    """Creates list of bins of length <step_width> with values describing 
    the number of reads that fall into a bin.
    <count_list> has to be created with 'get_count_list' 
    It returns a dict like: {'chr1' [0,10,2,0,...], 'chr2' ...} 
    where the first list entry gives the first bin and so on."""
    result = {}
    for chrom in chromosomes:
        overrun = 0
        if not chrom_len.has_key(chrom):
#             print("Warning: %s not found, do not consider" %chrom, file=sys.stderr)
            pass
        else:
            #print("... considering %s..."%chrom, file=sys.stderr)
            for i in range(0, chrom_len[chrom], step_width):
                end = min(i + step_width, chrom_len[chrom])
                counts = reduce(lambda x, y: x + y, count_list[ GenomicInterval(chrom, i, end) ])
                count_list[ GenomicInterval(chrom, i, end) ] = 0
                counts += overrun
    
                if chrom in result.keys():
                    result[chrom].append(counts)
                else:
                    result[chrom] = [counts] 
    
                overrun = _get_overrun(chrom, i, end, step_width, count_list, feature_len)

    return result

# def get_normalization_factor_directly(counts_dict_x, counts_dict_y, zero_counts):
#     counts_dict_x = counts_dict_x[:197196432] #only chr1
#     counts_dict_y = counts_dict_y[:197196432]
#     count_list = zip(counts_dict_x, counts_dict_y)
#     count_list = map(lambda x: [float(x[0]), float(x[1])], count_list)
#     
#     pq_list, max_index, max_value, factor1, factor2 = _get_lists(count_list, zero_counts)
#     
#     return factor1

def _get_lists(count_list, zero_counts, two_sample=False):
    if two_sample:
        count_list.sort(key = lambda x: x[0] + x[1])
        if not zero_counts:
            count_list = filter(lambda x: x[0] != 0.0 or x[1] != 0.0, count_list)
        pre_pq_list = list(_accumulate(count_list))
        pq_list = pre_pq_list
        max_index = int(len(pq_list)*0.5)
        max_value = fabs(pq_list[max_index][1] - pq_list[max_index][0])
    else:
        count_list.sort(key = lambda a: a[0])
        if not zero_counts:
            count_list = filter(lambda x: x[0] != 0.0, count_list)
        pre_pq_list = list(_accumulate(count_list))
        #get k, a from Diaz et al., 2012, compute p, q from Diaz et al., 2012
        pq_list = map(lambda x: [ x[0]/float(pre_pq_list[-1][0]), x[1]/float(pre_pq_list[-1][1]) ], pre_pq_list)
        max_index, max_value = max(enumerate(map(lambda x: fabs(x[0] - x[1]), pq_list)), key=operator.itemgetter(1))

    return pq_list, max_index, max_value, \
        float(pq_list[max_index][0]) / pq_list[max_index][1], \
        float(pq_list[max_index][1]) / pq_list[max_index][0]

def get_binstats(chrom_len, count_list_1, count_list_2, feature_len, chromosomes, step_width=1000, zero_counts=True, two_sample=False):
    """Compute p,q-list from Diaz et al., 2012 with k, a and coressponding factors
    for normalization.
    Also give list of tuples (x,y) describing bins where x,y are the values for each bin."""
    #create dict of counts
    #print("Dividing genome into bins...", file=sys.stderr)
    if set(chrom_len.keys()) & chromosomes == set(): #if chrom in input does not overlap chrom in IP channel
        return None, None, None, None, None
        
    counts_dict_2 = get_bins(chrom_len, chromosomes, count_list_2, step_width, feature_len)
    counts_dict_1 = get_bins(chrom_len, chromosomes, count_list_1, step_width, feature_len)
    #merge values with zip, obtain [ [(),()...], [(),(),..] ]
    tmp_list = [zip( counts_dict_1[chrom], counts_dict_2[chrom]) for chrom in counts_dict_2.keys()]
    #convert tmp_list to [[],[],[]...]
    count_list = map(lambda x: list(x), list(chain.from_iterable(tmp_list)) )
    
    return _get_lists(count_list, zero_counts, two_sample)

def work(first_path, second_path, step_width, zero_counts, two_sample, chrom_sizes_dict, stop):
    """work"""
    CHROM_LEN = chrom_sizes_dict #CHROM_LEN_HUMAN if genome == 'hg19' else CHROM_LEN_MOUSE
    #counts_1, counts_1 is a genomicarray (HTSeq) describing chr und pos of reads
    #print("Reading first input file...", file=sys.stderr)
    counts_1, chromosomes1 = get_count_list(CHROM_LEN, first_path, stop=stop)
    
    #print("Reading second input file...", file=sys.stderr)
    counts_2, chromosomes2 = get_count_list(CHROM_LEN, second_path, stop=stop)
#     print("...done", file=sys.stderr)
    chromosomes = chromosomes1 | chromosomes2
    
    #count_list1, count_list2 are dict like {'chr1': [10,2,0...], ...} describing the bins
    pq_list, max_index, max_value, factor1, factor2 \
    = get_binstats(CHROM_LEN, counts_1, counts_2, \
                               get_feature_len(first_path), chromosomes, step_width, zero_counts, two_sample)
    
    return pq_list, max_index, max_value, factor1, factor2, chromosomes

def get_normalization_factor(first_path, second_path, step_width, zero_counts, filename, debug, chrom_sizes_dict, two_sample=False, stop=False):
    """Return normalization factor (see Diaz et al) for the input
    if two_sample is True: compare sample with index of 0.15"""
    #take two largest chromosomes for analysis:
    #tmp = chrom_sizes_dict.items()
    #tmp.sort(key=lambda x: x[1], reverse=True)
    #tmp = dict(tmp[:min(len(tmp), 5)])
    tmp = {}
    for k in chrom_sizes_dict.keys():
        new_k = k if k.startswith('chr') else 'chr' + k
        tmp[new_k] = chrom_sizes_dict[k]
    #tmp = chrom_sizes_dict
    #print("For input normalization consider chromosomes: %s" %(", ".join(tmp.keys())), file=sys.stderr)
    
    pq_list, max_index, max_value, factor1, factor2, chromosomes =\
    work(first_path, second_path, step_width, zero_counts, two_sample, tmp, stop)
    
    if debug:
        write_pq_list(pq_list, max_index, max_value, factor1, factor2, filename + '-pqlist')
#    write_bedgraph(chrom_len, chromosomes, countstable, filename, step_width)
    
    #print("norm sums 1 : ", sum([i for (i,j) in pq_list]), file=sys.stderr)
    #print("norm sums 2 : ", sum([j for (i,j) in pq_list]), file=sys.stderr)
    
    if pq_list == None:
        return None, None
    
    if two_sample:
        l = 0.5
        s1 = sum([i for (i,j) in pq_list[:int(len(pq_list)*l)]])
        s2 = sum([j for (i,j) in pq_list[:int(len(pq_list)*l)]])
        
        #print("norm sums half: ", s1, s2, file=sys.stderr)
        
        if s1 > s2:
            return 2, factor1
        else:
            return 1, factor1
    else:
        return -1, factor1

def test_run():
    p1 = '/home/manuel/data/testdata/PU1_CDP.chr1-2.bam'
    p2 = '/home/manuel/data/testdata/PU1_Input4mix.chr1-2.bam'
    
#     p1 = '/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/MPP_WT_H3K27ac_1.bam'
#     p2 = '/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/CDP_WT_H3K27ac_1.bam'
    
    chrom_sizes_dict = { 'chr1' : 197195432, 'chr2':181748087,'chr3': 159599783,'chr4': 155630120,
             'chr5': 152537259,'chr6': 149517037,'chr7': 152524553,'chr8': 131738871,
             'chr9': 124076172,'chr10': 129993255,'chr11': 121843856,'chr12': 121257530,
             'chr13': 120284312,'chr14': 125194864,'chr15': 103494974,'chr16': 98319150,
             'chr17': 95272651,'chr18': 90772031,'chr19': 61342430,'chrX': 166650296,
             'chrY': 15902555} #, 'chrM': 16299}
    
    chrom_sizes_dict_mm10_enseml = { 'chr1' : 197195432, 'chr2':181748087,'chr3': 159599783,'chr4': 155630120,
             'chr5': 152537259,'chr6': 149517037,'chr7': 152524553,'chr8': 131738871,
             'chr9': 124076172,'chr10': 129993255,'chr11': 121843856,'chr12': 121257530,
             'chr13': 120284312,'chr14': 125194864,'chr15': 103494974,'chr16': 98319150,
             'chr17': 95272651,'chr18': 90772031,'chr19': 61342430,'chrX': 166650296,
             'chrY': 15902555} #, 'chrM': 16299}
    
    chrom_sizes_dict_mm10_enseml = {'10': 130694993}
    
    print('run')
    
    
    p1 = '/home/manuel/Downloads/odin/TC1-H3K4-A-D1.GRCm38.p3.q30_chr10.bam'
    p2 = '/home/manuel/Downloads/odin/TC1-I-A-D5.GRCm38.p3.q30_chr10.bam'
    
    a = get_normalization_factor(p1, p2, step_width=1000, zero_counts=0, \
                                                filename='test', debug=False, chrom_sizes_dict=chrom_sizes_dict_mm10_enseml, two_sample=False)
    print(a)
    
if __name__ == '__main__':
    #cProfile.run("test_run()")
    a = time()
    test_run()
    print(time() -a)
    