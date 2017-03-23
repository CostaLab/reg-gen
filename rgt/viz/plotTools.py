# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import time, datetime, argparse
from collections import *
import copy
import itertools
import pickle
import multiprocessing
import multiprocessing.pool
import urllib2
import re
import numpy
from scipy.stats import mstats, wilcoxon, mannwhitneyu, rankdata
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib import cm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib_venn import venn2, venn3

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.AnnotationSet import *
from rgt.Util import GenomeData, OverlapType, Html
from rgt.CoverageSet import *
from shared_function import *
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.helper import pretty

# Local test
dir = os.getcwd()

###########################################################################################
#                    Venn Diagram
###########################################################################################

class Venn:
    def __init__(self, sets, organism):
        self.sets = sets
        self.num = len(sets)
        self.organism = organism
        self.process_genesets()

    def process_genesets(self):
        self.gene_sets = []
        for gs in self.sets:
            if gs.endswith(".bed"):
                gregions = GenomicRegionSet(gs)
                gregions.read_bed(gs)
                associated_gr = gregions.gene_association(organism=self.organism, 
                    promoterLength=1000, threshDist=0, show_dis=False)
                genes = []
                for r in associated_gr: 
                    if ":" in r.name:
                        rgg = r.name.upper().split(":")
                        genes += rgg
                    else: genes.append(r.name.upper())
                self.gene_sets.append(set(genes))
                
            elif gs.endswith(".txt"):
                genes = []
                with open(gs) as f:
                    for line in f:
                        g = line.strip().split()
                        genes.append(g[0].upper())
                self.gene_sets.append(set(genes))

    def venn_diagram(self, directory, title, labels):
        
        if len(self.sets) == 2: 
            from matplotlib_venn import venn2

            f = plt.figure(figsize=(10,10))
            inter = len(self.gene_sets[0].intersection(self.gene_sets[1]))
            fig_venn = venn2(subsets = (len(self.gene_sets[0]) - inter, inter, len(self.gene_sets[1]) - inter),
                             set_labels = (self.sets[0].partition("/")[2].partition(".")[0], 
                                           self.sets[1].partition("/")[2].partition(".")[0]))
            plt.title("Sample Venn diagram")
            return f

        elif len(self.sets) == 3:
            from matplotlib_venn import venn3
            def write_genes(filename, geneset):
                with open(filename, "w") as g:
                    for gene in geneset: print(gene, file=g)

            f = plt.figure(figsize=(10,10))

            s100 = self.gene_sets[0] - self.gene_sets[1] - self.gene_sets[2]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[0]+".txt"),
                        geneset=s100)
            s010 = self.gene_sets[1] - self.gene_sets[0] - self.gene_sets[2]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[1]+".txt"),
                        geneset=s010)
            s001 = self.gene_sets[2] - self.gene_sets[0] - self.gene_sets[1]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[2]+".txt"),
                        geneset=s001)
            s111 = self.gene_sets[0].intersection(self.gene_sets[1].intersection(self.gene_sets[2]))
            write_genes(filename=os.path.join(directory, title,"list_"+labels[0]+"_"+labels[1]+"_"+labels[2]+".txt"),
                        geneset=s111)
            s110 = self.gene_sets[0].intersection(self.gene_sets[1]) - self.gene_sets[2]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[0]+"_"+labels[1]+".txt"),
                        geneset=s110)
            s011 = self.gene_sets[1].intersection(self.gene_sets[2]) - self.gene_sets[0]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[1]+"_"+labels[2]+".txt"),
                        geneset=s011)
            s101 = self.gene_sets[0].intersection(self.gene_sets[2]) - self.gene_sets[1]
            write_genes(filename=os.path.join(directory, title,"list_"+labels[0]+"_"+labels[2]+".txt"),
                        geneset=s101)
            
            fig_venn = venn3(subsets = (len(s100),len(s010),len(s110),len(s001),len(s101),len(s011),len(s111)),
                             set_labels = labels)
            return f