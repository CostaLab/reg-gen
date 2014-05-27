

import os
import sys
from copy import deepcopy
from GenomicRegion import GenomicRegion
from Util import GenomeData, MotifData, AuxiliaryFunctions

"""
Standardization of annotation.

Authors: Eduardo G. Gusmao.

"""

"""
To do:

- connect motif with file and biopython object
- keep track of databases (???)

"""

class Motif:  

    def __init__(self, tf_id, name, database, tf_class, genes, genes_suffix):
        """Initialize GenomicRegion"""
        self.id = tf_id
        self.name = name
        self.database = database
        self.tf_class=tf_class
        self.genes = genes
        self.genes_suffix = genes_suffix


class MotifSet:

    def __init__(self):
        """Initialize GenomicRegion"""
        self.motifs_map = {}
        self.genes_map={}
        self.genes_suffix_map={}
        

    def add(self, new_motif):
        #new_motif=Motif(tf_id, name, database, tf_class, genes, genes_suffix)
        try:
            self.motifs_map[new_motif.name]
            print tf_id,"repeated motif"
        except:
            self.motifs_map[new_motif.name]=new_motif
            for g in new_motif.genes:
              g=g.upper()
              try:
                motifs_aux=self.genes_map[g]
                motifs_aux.append(new_motif)
              except:  
                self.genes_map[g]=[new_motif]
            for g in new_motif.genes_suffix:
              try:
                motifs_aux=self.genes_suffix_map[g]
                motifs_aux.append(new_motif)
              except:  
                self.genes_suffix_map[g]=[new_motif]

    def match_suffix(self,gene_name):
        ''' XXX - this could be optmized'''
        res=[]
        gene_name=gene_name.upper()
        for s in self.genes_suffix_map.keys():
          if gene_name.startswith(s):
            res.append(s)
        return res   

    def filter_by_motifs(self,motifs):
        motif_set=MotifSet()
        for m in motifs:
          try:
            motif_set.add(self.motifs_map[m])
          except:
            print "motif not found", m
        return motif_set


    def filter_by_genes(self,genes,search="exact"):
        ''' This method returns motifs associated to genes. T
             The search has three modes: 
               exact - exact match only
               inexact - genes with no exact match are searched for inexact matcth
               all - all genes are applied to an inexact match
        '''
        motif_set=MotifSet()
        motifs_genes={}
        genes_motifs={}
        not_found_genes= [] # keep genes for inexact search
        genes=genes.genes
        for g in genes:
           g=g.upper()
           try:
             motifs=self.genes_map[g]
             for m in motifs:
                motif_set.add(m)
                try:
                  genes_motifs[g].append(m.name)   
                except:
                  genes_motifs[g]=[m.name]
                try:
                  motifs_genes[m.name].append(g)   
                except:
                  motifs_genes[m.name]=[g]

           except Exception, e:
             #print "not found" ,g
             not_found_genes.append(g) # keep genes for inexact search 
        if search == "inexact":
            genes = not_found_genes;
        elif search == "exact":
            genes=[];                  
        for g in genes:
             suffs=self.match_suffix(g)
             for s in suffs:
               motifs=self.genes_suffix_map[s]
               for m in motifs:
                  motif_set.add(m)
                  try:
                    genes_motifs[g].append(m.name)   
                  except:
                    genes_motifs[g]=[m.name]
                  try:
                    motifs_genes[m.name].append(g)   
                  except:
                    motifs_genes[m.name]=[g]
        return motif_set,genes_motifs,motifs_genes
                     

    def read_file(self, file_name_list):
        """
        Reads TF annotation in mtf (internal -- check manual) format.
        
        Keyword arguments:
        file_name_list -- A list with .mtf files.
        
        Return: void.
        """

        # Iterating over the file name list
        for file_name in file_name_list:

            # Opening MTF file
            try: mtf_file = open(file_name,"r")
            except Exception: pass # TODO

            # Reading file
            for line in mtf_file:
                # Processing line
                line_list = line.strip().split("\t")
                tf_id=line_list[0]
                name=line_list[1]
                database=line_list[2]
                tf_class=int(line_list[3])
                genes=line_list[4].split(";")
                genes_suffix=line_list[5].split(";")

                self.add(Motif(tf_id, name, database, tf_class, genes, genes_suffix))


            # Termination
            mtf_file.close()



