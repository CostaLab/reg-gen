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
    for region in exp_matrix.get_regionsets():
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association(region.fileName, None, gene_file, genome_file, threshDist=2000)
        for k in gene_peaks_mapping.keys():
            chr, raw_positions = k.split(':')
            start, end = map(lambda x: int(x), raw_positions.split('-'))
            
            #if peak is not assigned, an empty string occurs
            if "" in gene_peaks_mapping[k]:
                gene_peaks_mapping[k].remove("")
            
            list = 'NA' if not gene_peaks_mapping[k] else ','.join(gene_peaks_mapping[k])
            
            print(chr, start, end, list, sep='\t')
        

if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    parser.add_option("--mode", "-m", dest="mode", default=1, help="choose mode", type="int")
    (options, args) = parser.parse_args()
    path_exp_matrix = args[0]
    path_annotation = args[1]
    
    i = 2
    if len(args) != i:
        parser.error("Exactly %s parameters are needed" %i)
    
    #options.mode = 2
    #path_exp_matrix = '/home/manuel/test_exp_matrix'
    #path_annotation = '/home/manuel/data/rgtdata/mm9'
    
    genome_file = os.path.join(path_annotation, "chrom.sizes")
    gene_file = os.path.join(path_annotation, "association_file.bed")
    
    exp_matrix = ExperimentalMatrix()
    exp_matrix.read(path_exp_matrix)
    
    if options.mode is 1:
        mode_1(exp_matrix)
    elif options.mode is 2:
        mode_2(exp_matrix)

