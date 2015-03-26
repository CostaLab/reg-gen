import os
import sys
import sets
from copy import deepcopy
from GenomicRegion import GenomicRegion
from Util import GenomeData, MotifData, AuxiliaryFunctions
import glob

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
        self.motifs_map={}
        self.genes_map={}
        self.genes_suffix_map={}
        self.networks={}
        self.motifs_enrichment={}
        self.conditions=[]
        

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
        motif_sets=MotifSet()
        for m in motifs:
          try:
            motif_sets.add(self.motifs_map[m])
          except:
            print "motif not found", m
        return motif_sets


    def filter_by_genes(self,genes,search="exact"):
        ''' This method returns motifs associated to genes. T
             The search has three modes: 
               exact - exact match only
               inexact - genes with no exact match are searched for inexact matcth
               all - all genes are applied to an inexact match
        '''
        motif_sets=MotifSet()
        motifs_genes={}
        genes_motifs={}
        not_found_genes= [] # keep genes for inexact search
        genes=genes.genes
        print genes
        for g in genes:
           g=g.upper()
           try:
             motifs=self.genes_map[g]
             for m in motifs:
                motif_sets.add(m)
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
                  motif_sets.add(m)
                  try:
                    genes_motifs[g].append(m.name)   
                  except:
                    genes_motifs[g]=[m.name]
                  try:
                    motifs_genes[m.name].append(g)   
                  except:
                    motifs_genes[m.name]=[g]
        return motif_sets,genes_motifs,motifs_genes
                     

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
            #try: 
            mtf_file = open(file_name,"r")
            #except Exception: pass # TODO

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


    def  read_motif_targets_enrichment(self,enrichment_files,pvalue_threshold):
        ''' reads current output of motif enrichment analysis to get gene targets
            XXX - This function could be made internal by integration with motif enrichment search
        '''
        # reading networks
        for f in glob.glob(enrichment_files): 
          # use last dir name as name for condition
          condition=os.path.dirname(f)
          condition=condition.split("/")[-1]
          self.conditions.append(condition)
          network={}
          for line in open(f):
            line=line.strip("\n")
            values=line.split("\t")
            motif=values[0]
            try:
                self.motifs_map[motif]  
                p_value=float(values[2])
                genes=values[9].split(",")
                if (pvalue_threshold>=p_value):
                    network[motif]=genes
                try:
                  p_values=self.motifs_enrichment[motif]
                  p_values[condition]=p_value
                except:
                  p_values={}                   
                  p_values[condition]=p_value
                  self.motifs_enrichment[motif]=p_values                     
            except Exception, e:
                print("motif not found: "+motif+str(e))

          self.networks[condition]=network


    def write_enrichment_table(self,threshold,out_file,motifs_map):
      f = open(out_file,"w")
      f.write("\t"+("\t".join(self.conditions))+"\n")
      motifs=motifs_map.keys()       
      for v in self.motifs_enrichment.keys():
        values=self.motifs_enrichment[v]
        filter_p=False
        p_values=[]
        for c in self.conditions:
           try:
             pvalue=values[c]
             p_values.append(str(pvalue))    
             if pvalue<=threshold:
               filter_p=True
           except:
             p_values.append("1")
        if ((filter_p) & (v in motifs)):
          genes="|".join(motifs_map[v])
          f.write(v+"|"+genes+"\t"+("\t".join(p_values))+"\n")


#    def read_motif_targets(self,file_name,condition,pvalue_threshold):
#        ''' reads current output of motif enrichment analysis to get gene targets
#            XXX - This function could be made internal by integration with motif enrichment search
#        '''
#        network={}
#        for line in open(file_name):
#            line=line.strip("\n")
#            values=line.split("\t")
#            motif=values[0]
#            try:
#                self.motifs_map[motif]  
#                pvalue=float(values[2])
#                genes=values[9].split(",")
#                if (pvalue_threshold>=pvalue):
#                    network[motif]=genes
#            except Exception, e:
#                print("motif not found: "+motif+str(e))
#        self.networks[condition]=network

#    def get_motifs_networks():
#        '''
#        Getting all motifs presents in at least one network
#        '''
#        motifs={}
#        for net in self.networks:
#            for m in net.keys():
#              motifs[m]=[]
#        return motifs.keys()

    def write_cytoscape_network(self,genes,gene_mapping_search,out_path,targets,threshold):
        '''
        Write files to be used as input for cytoscape.
        it recieves a list of genes to map to, a mapping search strategy and path for outputting files
        '''
        #getting mapping of genes to motifs
        [motif_sets,genes_motifs,motifs_genes]=self.filter_by_genes(genes,gene_mapping_search)
        #print genes_motifs
        print "mapping of motifs to factors"

        f=open(out_path+"/mapping_tf_genes.txt","w")
        motifs_all={}
        for gene in genes_motifs.keys():
          motifs = genes_motifs[gene]
          for m in motifs:
            try:
              g=motifs_all[m]
              if gene not in g:
                motifs_all[m].append(gene)
            except:
              motifs_all[m]=[gene]
            f.write(gene+"\t"+m+"\n")
        f.close()

        import sets 

        self.write_enrichment_table(threshold,out_path+"/pvalue_table_"+str(threshold*100)+".txt",motifs_all)


        net_pairs={}
        net_tfs={}
        all_pairs=sets.Set()
        all_tfs=sets.Set()
        all_genes=sets.Set()
        if targets == None:
          filter_targets=False
        else:
          filter_targets=True
          targets=[g.upper() for g in targets.genes]
          #print targets

        # using genes to motif mapping to get network in all conditions
        for net_name in self.networks.keys():
            print net_name
            net=self.networks[net_name]
            pairs=sets.Set()
            tfs=sets.Set()
            net_pairs[net_name]=pairs
            net_tfs[net_name]=tfs
            for tf in genes_motifs.keys():
                motifs=genes_motifs[tf]
                for m in motifs:
                    try:
                      for target in net[m]:
                        if (not filter_targets or (target in targets)):
                          pairs.add((tf,target))                          
                          tfs.add(tf)
                          all_genes.add(tf)
                          all_genes.add(target)
                    except Exception, e:
                        print "motif not in network",m, tf, e
            all_pairs=all_pairs.union(pairs)
            all_tfs=all_tfs.union(tfs)
         
        #printing out network
        for net_name in net_pairs.keys():
            f=open(out_path+"/"+net_name+"_targets.txt","w")
            pairs_aux=net_pairs[net_name]
            for pair in all_pairs:
              #check if pair is active in the network  
              if pair in pairs_aux:
                f.write(pair[0]+"\t"+pair[1]+"\tactive\n")
              else:
                f.write(pair[0]+"\t"+pair[1]+"\tinactive\n")
            f.close()          
            f=open(out_path+"/"+net_name+"_genes.txt","w")
            for gene in all_genes:
              #check if gene is tf active in network
              if gene in net_tfs[net_name]:
                f.write(gene+"\ttf_active\n")
              elif gene in all_tfs:
                f.write(gene+"\ttf_inactive\n")
              else:
                f.write(gene+"\ttarget\n")
            f.close()       
              
