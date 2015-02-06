#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <experimental matrix file> <gene expression file> <annotation path> <output dir>

A intersectGenomicSets perform intersection statistics for a set of bed files.

Authors: Ivan G. Costa, Manuel Allhoff

It recieves as input a experimental matrix with a list of bed files and outputs simple overlap statistics.

"""

from __future__ import print_function
from optparse import OptionParser
import sys
from rgt.ExperimentalMatrix import *
from rgt.GenomicRegionSet import *
from rgt.CoverageSet import *
import rgt.GeneSet
import numpy


def averageExpression(region, expression, regionsToGenes):
    """Compute average gene expression"""
    values = numpy.zeros((len(region), len(expression.cond)))
    labels = []
    for i, b in enumerate(region):
        region_string = b.toString()
        labels.append(region_string)
        try:
            genes = regionsToGenes[region_string]

            for g in genes:
                values[i,] += expression.values[g.upper()]

            values[i,] = values[i,]/float(len(genes))
        except:
            values[i,] = 'NaN'
    
    return values, labels

def output(namesCol, namesLines, table, file_name):
    """Output table"""
    f = open(file_name,"w")
    f.write("\t"+("\t".join(namesCol))+"\n")
    for i,line in enumerate(table):
        f.write(namesLines[i]+"\t"+("\t".join([str(j) for j in line]))+"\n")

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

if __name__ == '__main__':
    parser = HelpfulOptionParser(usage=__doc__)
    (options, args) = parser.parse_args()
     
    if len(args) != 4:
        parser.error("Exactly three parameters are needed: experimental matrix, gene expression, annotation path and prefix for output")
    
    #map arguments
    experimental_matrix_file = args[0]
    gene_exp = args[1]
    annotation_path = args[2]
    outputdir = args[3]
    
    
#     experimental_matrix_file = "/home/manuel/workspace/cluster_p/THOR/exp/exp23_macs2_payton/1"
#     gene_exp = "/home/manuel/workspace/cluster_p/allhoff/project_THOR/data/payton/gene_expression/CCmean.data"
#     annotation_path = "/home/manuel/workspace/cluster_h/rgtdata/hg19/"
#     outputdir = "/home/manuel/test/"
    
    exps = ExperimentalMatrix()
    exps.read(experimental_matrix_file)
    regionsets = exps.get_regionsets()
    
    genome_file = annotation_path + "/chrom.sizes"
    gene_file = annotation_path + "/association_file.bed"
    
    genes = GeneSet("Expression")
    genes.read_expression(gene_exp)
    
    for region in regionsets:
        bedNew = GenomicRegionSet("")
        [degenes, de_peak_genes, mappedGenes, totalPeaks, regionsToGenes] \
        = bedNew.filter_by_gene_association_old(region.fileName, genes.genes, gene_file, genome_file)
        
        [ct, labels] = averageExpression(region, genes, regionsToGenes)
        aux = region.fileName.split("/")
        fileName = aux[-1]
        fileName = fileName.split(".")
        output(genes.cond, labels, ct, outputdir + "/" + fileName[0] + ".txt")
        
        

