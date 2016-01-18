#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""GeneSet describes genes and their expression.

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

@author: Ivan G. Costa, Manuel Allhoff, Joseph Kuo
"""

###############################################################################
# Libraries
###############################################################################
# Python
# Internal
from Util import GenomeData
# External

###############################################################################
# Class
###############################################################################

class GeneSet:
    """Describe genes and their expression."""
    
    def __init__(self,name):
        """Initiate a GeneSet"""
        self.name = name
        self.genes = [] #list of genes to consider
        self.values = {} #keys: gene, value: expression data as a list
        self.cond = []

    def __len__(self):
        """Return the number of genes"""
        return len(self.genes)
    
    def __iter__(self):
        """Iterate this GeneSet"""
        return iter(self.genes)

    def read(self, geneListFile):
        """Read genes from the file.

        Keyword arguments:
        geneListFile -- Path to the file which contains a list of genes
        """
        with open(geneListFile) as f:
            for line in f:            
                line = line.strip()
                if line:
                    l = line.split()
                    if l[0] != "":
                        self.genes.append(l[0].upper())
            
    def read_expression(self, geneListFile, header=False, valuestr=False):
        """Read gene expression data

        Keyword arguments:
        geneListFile -- Path to the file which contains genes and expression value
        header -- Read first line as header
        valuestr -- Keep the value as a string, otherwise convert to float number
        """
        with open(geneListFile) as f:
            if header:
                l = f.readline()
                l = l.strip("\n")
                l = l.split("\t")
                self.cond = l[1:len(l)]
            else:
                l = f.readline()
                l = l.strip("\n")
                l = l.split("\t")
                self.cond = [str(e) for e in range(len(l)-1)]
            for line in f.readlines():
                line = line.strip("\n")
                l = line.split("\t")
                if l[0] != "":
                    try:
                        self.genes.append(l[0].upper())
                        #self.values[l[0].upper()] = [float(v) for v in l[1:len(l)]]
                        if not valuestr:
                            self.values[l[0].upper()] = float(l[1])
                        else:
                            self.values[l[0].upper()] = l[1]
                    except:
                        print("*** error in loading gene: "+line)

    def get_all_genes(self, organism):
        """Get all gene names for a given organism
        
        Keyword arguments:
        organism -- Define the organism
        """
        genome = GenomeData(organism=organism)
        self.genes = []
        f = open(genome.get_association_file())
        for l in f.readlines():
            l = l.strip("\n")
            l = l.split("\t")
            self.genes.append(l[3].upper())
        f.close()

        self.genes = list(set(self.genes))

    def subtract(self, gene_set):
        """Subtract another GeneSet
        
        Keyword arguments:
        gene_set -- Another GeneSet for subtracting with
        """
        self.genes = [gene for gene in self.genes if gene not in gene_set.genes]
