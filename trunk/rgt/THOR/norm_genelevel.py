#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
from optparse import OptionParser
import sys
from rgt.ExperimentalMatrix import ExperimentalMatrix
from rgt.CoverageSet import CoverageSet
import numpy as np
import os
from rgt.GenomicRegionSet import GenomicRegionSet

def get_experimental_matrix(bams, bed):
    """Load artificially experimental matrix. Only genes in BED file are needed."""
    m = ExperimentalMatrix()
    
    m.fields = ['name', 'type', 'file']
    m.fieldsDict = {}
    
    names = []
    for bam in bams:
        n, _ = os.path.splitext(os.path.basename(bam))
        m.files[n] = bam
        names.append(n) 
    m.names = np.array(['housekeep'] + names)
    m.types = np.array(['regions'] + ['reads']*len(names))
    g = GenomicRegionSet('RegionSet')
    g.read_bed(bed)
    m.objectsDict['housekeep'] = g
    
    return m

def get_factor_matrix(d, colnames):
    """Give matrix describing factors between genes. Idealy, factors in a column should be approx. the same."""
    res = []
    if d.shape[0] > 1 and d.shape[1] > 1:
        for i in range(d.shape[0]):
            for j in range(d.shape[1]-1):
                res.append(d[i,j] / d[i,j+1])
    
                res = np.matrix(res).reshape((d.shape[0], d.shape[1]-1))
                colnames = map(lambda x: x[0] + "-" + x[1], zip(colnames[:len(colnames)-1], colnames[1:]))
    else:
        res, colnames = d, colnames
    
    return res, colnames

def output_R_file(name, res, colnames):
    """"Write R code to file to check whether genes give same signal among the samples"""
    f = open(name + 'norm.R', 'w')
    #output for R
    #everthing in one vector
    
    #if res.shape[1] > 0:
    l = reduce(lambda x, y: x+y, [map(lambda x: str(x), list(np.array(res[:,i]).reshape(-1,))) for i in range(res.shape[1])])
    #else:
    #    l = list(np.array(res.reshape(-1,)))
    
    print('d = c(', ', '.join(l), ')', sep='', file=f)
    print('d = matrix(d, %s)' %res.shape[0], file=f)
    print('names = c("', '", "'.join(colnames), '")', sep='', file=f)
    print('par(mar=c(15,5,5,5))', file=f)
    print('barplot(d, beside = TRUE, names.arg = names, las=2, main="Housekeeping Genes Ratio", ylab="Signal")', file=f)
    

def norm_gene_level(bams, bed, name, verbose, promotor_length = 1000):
    """Normalize bam files on a gene level. Give out list of normalization factors."""
    m = get_experimental_matrix(bams, bed)
    
    d = zip(m.types, m.names)
    d = map(lambda x: x[1], filter(lambda x: x[0] == 'reads', d)) #list of names which are reads
    
    regions = m.objectsDict['housekeep'] #GenomicRegionSet containing housekeeping genes
    for el in regions:
        if el.data == '-':
            el.extend(promotor_length, 0)
        elif el.data == '+':
            el.extend(0, promotor_length)
    
    covs = []
    
    for cond in d:
        bam_path = m.files[cond]
        c = CoverageSet(cond, regions) 
        c.coverage_from_bam(bam_file=bam_path)
        c.genomicRegions.sort()
        covs.append(c)
    
    #create matrix sample x gene for signal
    signals = [[sum(covs[k].coverage[i]) for i in range(len(covs[k].genomicRegions))] for k in range(len(covs))]
    
    assert len(covs) > 0
    
    gene_names = [covs[0].genomicRegions[i].name for i in range(len(covs[0].genomicRegions))]
    
    colnames = gene_names
    d = np.matrix(signals, dtype=float)
    
    if verbose:
        #output R code to check wether gene give same signal
        res, colnames = get_factor_matrix(d, colnames)
        output_R_file(name, res, colnames)
    
    #normalize: increase values to highest value  
    colmax = np.amax(d, axis=0)

    for i in range(d.shape[0]):
        for j in range(d.shape[1]):
            d[i,j] = colmax[:,j][0,0]/d[i,j]
    
    return list(np.array(np.mean(d, axis=1)).reshape(-1))
    

if __name__ == '__main__':
    #bams = ['/home/manuel/test1.bam', '/home/manuel/test2.bam']
    bams = ['/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_BCRABL_H3K9ac_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_BCRABL_H3K9ac_rep2.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_BCRABL_IM_H3K9ac_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_BCRABL_IM_H3K9ac_rep2.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_JAK2VF_H3K9ac_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_JAK2VF_H3K9ac_rep2.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_JAK2VF_Rux_H3K9ac_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_JAK2VF_Rux_H3K9ac_rep2.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_LV_H3K9ac_forBCRABL_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_LV_H3K9ac_forBCRABL_rep2.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_LV_H3K9ac_forJAK2VF_rep1.bam', '/home/manuel/workspace/cluster_p/hematology/local/new_run/bam/32D_mm_LV_H3K9ac_forJAK2VF_rep2.bam']
    bed = '/home/manuel/workspace/cluster_p/hematology/exp/exp16_check_housekeeping_genes/pot_housekeeping_genes_mm9.bed'
    
    bams = ['/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/MPP_WT_H3K27ac_1.bam', '/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/MPP_WT_H3K27ac_2.bam', '/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/CDP_WT_H3K27ac_1.bam', '/home/manuel/workspace/cluster_p/dendriticcells/local/zenke_histones/bam/CDP_WT_H3K27ac_2.bam']
    bed = '/home/manuel/hk.bed'
    
    
    print(norm_gene_level(bams, bed, 'testname', True))
    
    
    
    
    
    