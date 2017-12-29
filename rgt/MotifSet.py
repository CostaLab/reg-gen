#!/usr/bin/env python

"""
MotifSet
===================
Represents a transcription factor motif and the standardization of motif annotation.

"""

# Python
from __future__ import print_function

import glob
import os


class MotifAnnotation:
    """Represents a transcription factor with all available annotation.

    *Keyword arguments:*

      - tf_id -- Transcription factor ID.
      - name -- Transcription factor name (symbol).
      - database -- Database/repository in which this motif was obtained from.
      - tf_class -- Class of transcription factor motif (can be any string).
      - gene_names -- List of gene names for this transcription factor (usually only one, sometimes two)
      - uniprot_ids -- List of UniProt accession IDs for this transcription factor (like above)
    """

    # TODO: merge with Motif and/or Thresholds for fpr information?

    def __init__(self, tf_id, name, database, version, gene_names, tf_class, uniprot_ids):
        self.id = tf_id
        self.name = name
        self.database = database
        self.version = version
        self.gene_names = gene_names
        self.tf_class = tf_class
        self.uniprot_ids = uniprot_ids


class MotifSet:
    """Represents a set of motifs."""

    def __init__(self):
        self.motifs_map = {}
        self.genes_map = {}
        self.genes_suffix_map = {}
        self.networks = {}
        self.motifs_enrichment = {}
        self.conditions = []

    def add(self, motif):
        """Adds a new motif to this set.

        *Keyword arguments:*

          - new_motif -- New motif to be added.
        """

        if motif.name not in self.motifs_map:
            self.motifs_map[motif.name] = motif

            # adding genes
            for g in motif.genes:
                if g not in self.genes_map:
                    self.genes_map[g] = []

                self.genes_map[g].append(motif)

            # adding gene suffixes
            for g in motif.genes_suffix:
                if g not in self.genes_suffix_map:
                    self.genes_suffix_map[g] = []

                self.genes_suffix_map[g].append(motif)

    def match_suffix(self, gene_name):
        """Match with gene suffix

        *Keyword arguments:*

          - gene_name -- Gene name to perform the match.

        *Return:*

          - res -- ID of mapped genes.
        """

        res = []
        gene_name = gene_name.upper()
        for s in self.genes_suffix_map.keys():
            if gene_name.startswith(s):
                res.append(s)
        return res

    def filter_by_motifs(self, motifs):
        """Filter this motif set by defined motifs.

        *Keyword arguments:*

          - motifs -- Motifs in which to filter this set.

        *Return:*

          - motif_set -- Set of filtered motifs.
        """

        motif_set = MotifSet()

        for m in motifs:
            if m in self.motifs_map:
                motif_set.add(self.motifs_map[m])
            else:
                print("motif not found: " + str(m))

        return motif_set

    def filter_by_genes(self, genes, search="exact"):
        """This method returns motifs associated to genes. The search has three modes: 
            1. 'exact' - exact match only
            2. 'inexact' - genes with no exact match are searched for inexact matcth
            3. 'all' - all genes are applied to an inexact match

        *Keyword arguments:*

          - genes -- Gene set to perform the filtering.
          - search -- Search mode (default = 'exact').

        *Return:*

          - motif_set -- set of filtered motifs.
          - genes_motifs -- Dictionary of genes to motifs.
          - motifs_genes -- Dictionary of motifs to genes.
        """
        motif_set = MotifSet()
        motifs_genes = {}
        genes_motifs = {}
        not_found_genes = []  # keep genes for inexact search
        genes = genes.genes
        for g in genes:
            g = g.upper()

            if g in self.genes_map:
                motifs = self.genes_map[g]
                for m in motifs:
                    motif_set.add(m)

                    if g in genes_motifs:
                        genes_motifs[g].append(m.name)
                    else:
                        genes_motifs[g] = [m.name]

                    if m.name in motifs_genes:
                        motifs_genes[m.name].append(g)
                    else:
                        motifs_genes[m.name] = [g]
            else:
                not_found_genes.append(g)  # keep genes for inexact search

        if search == "inexact":
            genes = not_found_genes
        elif search == "exact":
            genes = []
        for g in genes:
            suffs = self.match_suffix(g)
            for s in suffs:
                motifs = self.genes_suffix_map[s]
                for m in motifs:
                    motif_set.add(m)

                    if g in genes_motifs:
                        genes_motifs[g].append(m.name)
                    else:
                        genes_motifs[g] = [m.name]

                    if m.name in motifs_genes:
                        motifs_genes[m.name].append(g)
                    else:
                        motifs_genes[m.name] = [g]

        return motif_set, genes_motifs, motifs_genes

    def read_file(self, file_name_list):
        """Reads TF annotation in mtf (internal format; check manual) format.
        
        *Keyword arguments:*

          - file_name_list -- A list with .mtf files.
        """

        # Iterating over the file name list
        for file_name in file_name_list:

            # Opening MTF file
            mtf_file = open(file_name, "r")

            # Reading file
            for line in mtf_file:
                # Processing line
                line_list = line.strip().split("\t")
                tf_id = line_list[0].strip()
                name = line_list[1].strip()
                database = line_list[2].strip()
                version = int(line_list[3].strip())
                gene_names = line_list[4].strip().split(";")
                tf_class = line_list[5].strip()
                uniprot_ids = line_list[6].strip().split(";")

                self.add(MotifAnnotation(tf_id, name, database, version, gene_names, tf_class, uniprot_ids))

            # Termination
            mtf_file.close()

    def read_motif_targets_enrichment(self, enrichment_files, pvalue_threshold):
        """Reads current output of motif enrichment analysis to get gene targets.

        *Keyword arguments:*

          - enrichment_files -- Enrichment files to read.
          - pvalue_threshold -- P-value threshold for motif acceptance.
        """

        # reading networks
        for f in glob.glob(enrichment_files):
            # use last dir name as name for condition
            condition = os.path.dirname(f)
            condition = condition.split("/")[-1]
            self.conditions.append(condition)
            network = {}
            for line in open(f):
                line = line.strip("\n")
                values = line.split("\t")
                motif = values[0]
                if motif in self.motifs_map:
                    p_value = float(values[2])
                    genes = values[9].split(",")
                    if pvalue_threshold >= p_value:
                        network[motif] = genes
                    p_values = self.motifs_enrichment[motif]
                    if condition not in p_values:
                        p_values = {condition: p_value}
                        self.motifs_enrichment[motif] = p_values
                    p_values[condition] = p_value
                else:
                    print("motif not found: " + motif)

            self.networks[condition] = network

    def write_enrichment_table(self, threshold, out_file, motifs_map):
        """Writes enrichment table for network generation.

        *Keyword arguments:*

          - threshold -- P-value threshold for motif acceptance.
          - out_file -- Output file name.
          - motifs_map -- Mapping of motifs.
        """

        f = open(out_file, "w")
        f.write("\t" + ("\t".join(self.conditions)) + "\n")
        motifs = motifs_map.keys()
        for v in self.motifs_enrichment.keys():
            values = self.motifs_enrichment[v]
            filter_p = False
            p_values = []
            for c in self.conditions:
                try:
                    pvalue = values[c]
                    p_values.append(str(pvalue))
                    if pvalue <= threshold:
                        filter_p = True
                except:
                    p_values.append("1")
            if filter_p and (v in motifs):
                genes = "|".join(motifs_map[v])
                f.write(v + "|" + genes + "\t" + ("\t".join(p_values)) + "\n")

    def write_cytoscape_network(self, genes, gene_mapping_search, out_path, targets, threshold):
        """
        Write files to be used as input for cytoscape. It receives a list of genes to map to, a mapping search
        strategy and path for outputting files.

        *Keyword arguments:*

          - genes -- Gene set.
          - gene_mapping_search -- Gene mapping.
          - out_path -- Output path.
          - targets -- Gene targets.
          - threshold -- Threshold for motif acceptance.
        """

        # getting mapping of genes to motifs
        [_, genes_motifs, _] = self.filter_by_genes(genes, gene_mapping_search)

        f = open(out_path + "/mapping_tf_genes.txt", "w")
        motifs_all = {}
        for gene in genes_motifs.keys():
            motifs = genes_motifs[gene]
            for m in motifs:
                if m in motifs_all:
                    if gene not in motifs_all[m]:
                        motifs_all[m].append(gene)
                else:
                    motifs_all[m] = [gene]

                f.write(gene + "\t" + m + "\n")
        f.close()

        self.write_enrichment_table(threshold, out_path + "/pvalue_table_" + str(threshold * 100) + ".txt", motifs_all)

        net_pairs = {}
        net_tfs = {}
        all_pairs = set()
        all_tfs = set()
        all_genes = set()
        if not targets:
            filter_targets = False
        else:
            filter_targets = True
            targets = [g.upper() for g in targets.genes]

        # using genes to motif mapping to get network in all conditions
        for net_name in self.networks.keys():
            net = self.networks[net_name]
            pairs = set()
            tfs = set()
            net_pairs[net_name] = pairs
            net_tfs[net_name] = tfs
            for tf in genes_motifs.keys():
                motifs = genes_motifs[tf]
                for m in motifs:
                    if m in net:
                        for target in net[m]:
                            if not filter_targets or (target in targets):
                                pairs.add((tf, target))
                                tfs.add(tf)
                                all_genes.add(tf)
                                all_genes.add(target)
                    else:
                        print("motif not in network: " + str(m) + " " + str(tf) + " ")
            all_pairs = all_pairs.union(pairs)
            all_tfs = all_tfs.union(tfs)

        # printing out network
        for net_name in net_pairs.keys():
            f = open(out_path + "/" + net_name + "_targets.txt", "w")
            pairs_aux = net_pairs[net_name]
            for pair in all_pairs:
                # check if pair is active in the network
                if pair in pairs_aux:
                    f.write(pair[0] + "\t" + pair[1] + "\tactive\n")
                else:
                    f.write(pair[0] + "\t" + pair[1] + "\tinactive\n")
            f.close()
            f = open(out_path + "/" + net_name + "_genes.txt", "w")
            for gene in all_genes:
                # check if gene is tf active in network
                if gene in net_tfs[net_name]:
                    f.write(gene + "\ttf_active\n")
                elif gene in all_tfs:
                    f.write(gene + "\ttf_inactive\n")
                else:
                    f.write(gene + "\ttarget\n")
            f.close()
