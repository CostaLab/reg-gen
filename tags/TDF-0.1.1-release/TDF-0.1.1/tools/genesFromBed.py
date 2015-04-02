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
import numpy as np

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def mode_1(exp_matrix,thresh):
    for region in exp_matrix.get_regionsets():
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, _ = region_set.filter_by_gene_association_old(region.fileName, None, gene_file, genome_file, threshDist=thresh)
        print('#number of mapped genes:', mappedGenes)
        print(region.name+"\t"+("\t".join(region_set.genes)))

def mode_2(exp_matrix,thresh):
    #remember value of bedgraph, ugly way
    value = {}
    for regions in exp_matrix.get_regionsets():
        for region in regions:
            value[(region.chrom, region.initial, region.final)] = region.data
    
    for region in exp_matrix.get_regionsets():
        f = open("region_" + str(region.name) + ".data", 'w')
        
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association_old(region.fileName, None, gene_file, genome_file, threshDist=thresh)

        print(mappedGenes)
        
        for k in gene_peaks_mapping.keys():
            chr, raw_positions = k.split(':')
            start, end = map(lambda x: int(x), raw_positions.split('-'))
            
            #if peak is not assigned, an empty string occurs
            if "" in gene_peaks_mapping[k]:
                gene_peaks_mapping[k].remove("")
            
            list = 'NA' if not gene_peaks_mapping[k] else ','.join(gene_peaks_mapping[k])
            
            print(chr, start, end, value[(chr, start, end)], list, sep='\t', file = f)
        
        f.close()

def mode_3(exp_matrix, thresh, type_file):
    #remember value of bedgraph, ugly way
    score = {}
    for regions in exp_matrix.get_regionsets():
        for region in regions:
            if type_file=="ODIN":
              aux=(region.data).split("\t")
              aux=aux[-1].split(";")
              score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = float(region.data[-1])
            if type_file=="THOR":
              aux=(region.data).split(";")
              score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = float(aux[-1])
            else:
               score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = region.data
    
    
    for i, region in enumerate(exp_matrix.get_regionsets()):
        f = open("region_" + str(region.name) + ".data", 'w')
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association_old(region.fileName, None, gene_file, genome_file, threshDist=thresh)

        avg_score = {} #score per peak
        genes = {}
        
        print('Consider row %s of exp. matrix, number of mapped genes is %s' %(i, mappedGenes), file=sys.stderr)
        for peak, gene_list in gene_peaks_mapping.items():            
            for gen in gene_list: #reverse mapping peak -> gene to gene -> peak
                if not gen:
                    continue
                genes[gen] = genes.get(gen, set())
                genes[gen].add(peak)
                
                avg_score[gen] = avg_score.get(gen, [])
                avg_score[gen].append(score[peak]) #join all scores of peaks assigned to a gen
        
        for gen in genes.keys():
            if options.metric == 'mean':
                avg = np.mean(avg_score[gen])
            elif options.metric == 'max':
                avg = np.max(avg_score[gen])
            print(gen, avg, ", ".join(str(t) for t in genes[gen]), sep='\t', file = f)
        
        f.close()
     

def mode_4(exp_matrix,thresh,type_file,geneexp_file):
    #remember value of bedgraph, ugly way
        
    gene_set = GeneSet("")    
    gene_set.read_expression(geneexp_file)

    score = {}
    for regions in exp_matrix.get_regionsets():
        for region in regions:
            if type_file=="ODIN":
              aux=(region.data).split("\t")
              aux=aux[-1].split(";")
              score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = aux[-1]
            else:
              score[(region.chrom + ':' + str(region.initial) + '-' + str(region.final))] = region.data
    
    
    for region in exp_matrix.get_regionsets():
        f = open("region_" + str(region.name) + ".data", 'w')
        
        region_set = GenomicRegionSet("")
        _, _, mappedGenes, _, gene_peaks_mapping = region_set.filter_by_gene_association_old(region.fileName, gene_set.genes, gene_file, genome_file, threshDist=thresh)

        print(mappedGenes)

        #region.filter_by_gene_association(organism=organism,threshDist=thresh)
        # _, _, mappedGenes, _, gene_peaks_mapping

         
        
        avg_score = {} #score per peak
        genes = {}
        
        print(region)
        for peak, gene_list in gene_peaks_mapping.items():            
            for gen in gene_list: #reverse mapping peak -> gene to gene -> peak
                if not gen:
                    continue
                genes[gen] = genes.get(gen, set())
                genes[gen].add(peak)
                    
                
                avg_score[gen] = avg_score.get(gen, [])
                avg_score[gen].append(score[peak]) #join all scores of peaks assigned to a gen

        print(avg_score)
        
        for gen in gene_set.genes:
            try:
              avg = sum(map(lambda x: float(x), avg_score[gen]))/ float(len(avg_score[gen]))
              peaks = ", ".join(str(t) for t in genes[gen])
              siz=avg*len(avg_score[gen])
            except:
              avg = 0.0 
              siz=0
              peaks = "_"           
            try:
              print(gen, "\t".join([str(t) for t in gene_set.values[gen.upper()]]),  avg, siz,peaks , sep='\t', file = f)
            except:
              pass
               
        f.close()


if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    parser.add_option("--mode", "-m", dest="mode", default=1, help="choose mode", type="int")
    parser.add_option("--distance", "-d", dest="distance", default=50000, help="distance from peak to gene", type="int")
    parser.add_option("--type", "-t", dest="type", default="bed", help="type of bed file (<bed>, <ODIN>, <THOR>)", type="str")
    parser.add_option("--metric", dest="metric", default="max", help="metric to merge peaks' scores (mean, max)", type="str")
    (options, args) = parser.parse_args()
     
    i = 3
    if len(args) > i:
        parser.error("Exactly %s parameters are needed" %i)
        
    path_exp_matrix = args[0]
    path_annotation = args[1]
    
#     options.mode = 3
#     options.distance = 50000
#     options.type='THOR'
#     options.metric = 'max'
#     path_exp_matrix = '/home/manuel/workspace/cluster_p/hematology/exp/exp12_peak_gene_assignment/assign_peaks_mm.expmatrix'
#     path_annotation = '/home/manuel/rgt-data/mm9/'
    
    genome_file = os.path.join(path_annotation, "chrom.sizes")
    gene_file = os.path.join(path_annotation, "association_file.bed")
    
    exp_matrix = ExperimentalMatrix()
    exp_matrix.read(path_exp_matrix, is_bedgraph=False)
    
    print("Use metric %s to merge peaks' score." %options.metric, file=sys.stderr)
    
    if options.mode is 1:
        mode_1(exp_matrix,options.distance)
    elif options.mode is 2:
        mode_2(exp_matrix,options.distance)
    elif options.mode is 3:
        mode_3(exp_matrix,options.distance,options.type)
    elif options.mode is 4:
        geneexp_file = args[2]
        mode_4(exp_matrix,options.distance,options.type,geneexp_file)

