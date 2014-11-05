#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <EXP MATRIX> <PATH>

Determine experimental matrix M with columns name, type and file.
Give path to organism-specific rgtgen data folder.

Choose between modes:
mode 1 (default):
Output all genes that are associated to all regions given by M.

mode 2:
Output for each region of M the associated genes.
Create *.data file for each row in M.

mode 3:
Assign to each gene a list of peaks.
Create *.data file for each row in M.

@Author: Ivan Costa, Manuel Allhoff

"""
from __future__ import print_function
import sys
from optparse import OptionParser
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
import scipy.stats
from rgt.Util import GenomeData

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def mode_1(exp_matrix):
    for region in exp_matrix.get_regionsets():
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, _ = region_set.filter_by_gene_association(region.fileName, None, gene_file, genome_file, threshDist=50000)
        print('#number of mapped genes:', mappedGenes)
        print(region.name+"\t"+("\t".join(region_set.genes)))

def mode_2(exp_matrix):
    
    #remember value of bedgraph, ugly way
    value = {}
    for regions in exp_matrix.get_regionsets():
        for region in regions:
            value[(region.chrom, region.initial, region.final)] = region.data
    
    for region in exp_matrix.get_regionsets():
        f = open("region_" + str(region.name) + ".data", 'w')
        
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association(region.fileName, None, gene_file, genome_file, threshDist=2000)
        
        for k in gene_peaks_mapping.keys():
            chr, raw_positions = k.split(':')
            start, end = map(lambda x: int(x), raw_positions.split('-'))
            
            #if peak is not assigned, an empty string occurs
            if "" in gene_peaks_mapping[k]:
                gene_peaks_mapping[k].remove("")
            
            list = 'NA' if not gene_peaks_mapping[k] else ','.join(gene_peaks_mapping[k])
            
            print(chr, start, end, value[(chr, start, end)], list, sep='\t', file = f)
        
        f.close()

def mode_3(exp_matrix):
    #remember value of bedgraph, ugly way
    score = {}
    for regions in exp_matrix.get_regionsets():
        for region in regions:
            score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = region.data
    
    
    for region in exp_matrix.get_regionsets():
        f = open("region_" + str(region.name) + ".data", 'w')
        
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association(region.fileName, None, gene_file, genome_file, threshDist=2000)
        
        avg_score = {} #score per peak
        genes = {}
        
        for peak, gene_list in gene_peaks_mapping.items():
            for gen in gene_list: #reverse mapping peak -> gene to gene -> peak
                if not gen:
                    continue
                genes[gen] = genes.get(gen, set())
                genes[gen].add(peak)
                avg_score[gen] = avg_score.get(gen, [])
                avg_score[gen].append(score[peak]) #join all scores of peaks assigned to a gen
        
        for gen in genes.keys():
            avg = sum(map(lambda x: float(x), avg_score[gen]))/ float(len(avg_score[gen]))
            print(gen, avg, ", ".join(str(t) for t in genes[gen]), sep='\t', file = f)
               
        f.close()
     

if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    parser.add_option("--mode", "-m", dest="mode", default=1, help="choose mode", type="int")
    (options, args) = parser.parse_args()
    
    i = 2
    if len(args) != i:
        parser.error("Exactly %s parameters are needed" %i)
      
    path_exp_matrix = args[0]
    path_annotation = args[1]
    
    #options.mode = 3
    #path_exp_matrix = '/workspace/cluster_p/hematology/exp/exp03_rerun_chipseq/assign_peak/exp_matrix_peak_assign_chipseq'
    #path_annotation = '/home/manuel/data/rgt-data/mm9/'
    
    genome_file = os.path.join(path_annotation, "chrom.sizes")
    gene_file = os.path.join(path_annotation, "association_file.bed")
    
    exp_matrix = ExperimentalMatrix()
    exp_matrix.read(path_exp_matrix, is_bedgraph=True)
    
    if options.mode is 1:
        mode_1(exp_matrix)
    elif options.mode is 2:
        mode_2(exp_matrix)
    elif options.mode is 3:
        mode_3(exp_matrix)

