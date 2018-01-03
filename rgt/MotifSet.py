"""
MotifSet
===================
Represents a transcription factor motif and the standardization of motif annotation.

"""

# Python 3 compatibility
from __future__ import print_function

# Python
import glob
import os

# Internal
from rgt.GeneSet import GeneSet


class MotifAnnotation:
    """
    Represents a transcription factor with all available annotation (from MTF files).

    *Keyword arguments:*

      - tf_id -- Transcription factor ID.
      - name -- Transcription factor name (symbol).
      - database -- Database/repository in which this motif was obtained from.
      - family -- Class/Family of transcription factor motif (can be any string).
      - gene_names -- List of gene names for this transcription factor (usually only one, sometimes two)
      - uniprot_ids -- List of UniProt accession IDs for this transcription factor (like above)
    """

    def __init__(self, tf_id, name, database, version, gene_names, family, uniprot_ids):
        self.id = tf_id
        self.name = name
        self.database = database
        self.version = version
        self.gene_names = gene_names
        self.family = family
        self.uniprot_ids = uniprot_ids

    def __str__(self):
        return "MotifAnnotation[id={},name={},database={},version={},gene_names={},family={},uniprot_ids={}]" \
            .format(self.id, self.name, self.database, self.version, self.gene_names, self.family, self.uniprot_ids)


class MotifSet:
    """
    Represents a set of motifs. It contains MotifAnnotation instances.
    """

    def __init__(self):
        self.motifs_map = {}
        self.networks = {}
        self.motifs_enrichment = {}
        self.conditions = []

    def __len__(self):
        return len(self.motifs_map)

    def __iter__(self):
        return iter(self.motifs_map)

    def __getitem__(self, key):
        return self.motifs_map[key]

    def add(self, motif):
        """
        Adds a new motif to this set.

        *Keyword arguments:*

          - motif -- Must be an instance of MotifAnnotation
        """

        self.motifs_map[motif.name] = motif

    def filter(self, values, key_type="name"):
        """
        Returns a new MotifSet containing all matching motifs. By default, it expects a list of motif names, but this
        can be configured via the key_type parameter.

        *Keyword arguments:*

          - values -- List of strings representing the motif to filter this set on. Actual values depend on key_type.
          - key_type -- "name" for exact matching on the motif name; "family" for (partial) matching in the motif family
            description; "uniprot_id" for exact matching on UniProt IDs; "gene_name" for matching on the exact
            gene name of the motif (see also filter_by_genes for a more flexible alternative).

        *Return:*

          - motif_set -- Set of filtered motifs.
        """

        valid_keys = ["name", "gene_name", "family", "uniprot_id"]

        if key_type not in valid_keys:
            raise ValueError("key_type must be one of", valid_keys)

        # no need to iterate here if we are trying to do gene-name matching,
        # let's delegate
        if key_type == "gene_name":
            gene_set = GeneSet("")
            gene_set.genes = values
            motif_set, _, _ = self.filter_by_genes(gene_set)

            return motif_set

        # in all other cases, we create a new motif set and populate it with matching motifs

        motif_set = MotifSet()

        for v in values:
            if key_type == "name" and v in self.motifs_map:
                motif_set.add(self.motifs_map[v])
            else:
                for m in self.motifs_map.values():
                    if v in getattr(m, key_type)():
                        motif_set.add(m)

        return motif_set

    def filter_by_genes(self, genes, search="exact"):
        """
        This method returns a new MotifSet of motifs associated to those genes. The search has three modes:

            1. 'exact' - exact match only
            2. 'inexact' - genes with no exact match are searched for inexact match
            3. 'all' - all genes are applied to an inexact match

        *Keyword arguments:*

          - genes -- GeneSet to perform the filtering.
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

        # get the actual list of genes inside the GeneSet
        genes = genes.genes

        # exact matching
        for motif in self.motifs_map.values():
            for g in motif.gene_names:
                # TODO: this can be replaced with a "match" function, to avoid re-doing
                # inexact matching later on (depending on search parameter)
                if g in genes:
                    motif_set.add(motif)

                    # add this motif to gene2motif dict
                    if g in genes_motifs:
                        genes_motifs[g].append(motif.name)
                    else:
                        genes_motifs[g] = [motif.name]

                    # add this gene to motif2gene dict
                    if motif.name in motifs_genes:
                        motifs_genes[motif.name].append(g)
                    else:
                        motifs_genes[motif.name] = [g]
                else:
                    # keep genes for inexact search
                    not_found_genes.append(g)

        if search == "inexact":
            genes = not_found_genes
        elif search == "exact":
            genes = []

        # TODO: inexact matching - how to? We could use some string distance metric, or
        # a "suffix approach" like there was before

        return motif_set, genes_motifs, motifs_genes

    def read_mtf(self, mtf_filenames):
        """
        Reads TF annotation in mtf (internal format; check manual) format.

        *Keyword arguments:*

          - file_name_list -- A list with .mtf files.
        """

        # Iterating over the file name list
        for filename in mtf_filenames:

            # Opening MTF file
            mtf_file = open(filename, "r")

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

    def read_enrichment(self, enrichment_files, pvalue_threshold=1):
        """
        Reads current output of motif enrichment analysis to get gene targets.

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

                    if motif in self.motifs_enrichment:
                        self.motifs_enrichment[motif][condition] = p_value
                    else:
                        self.motifs_enrichment[motif] = {condition: p_value}
                else:
                    print("motif not found: " + motif)

            self.networks[condition] = network

    def write_enrichment(self, threshold, out_file, motifs_map):
        """
        Writes enrichment table for network generation.

        *Keyword arguments:*

          - threshold -- P-value threshold for motif acceptance.
          - out_file -- Output file name.
          - motifs_map -- Mapping of motifs.
        """

        f = open(out_file, "w")
        f.write("\t" + ("\t".join(self.conditions)) + "\n")

        for v in self.motifs_enrichment:
            values = self.motifs_enrichment[v]
            filter_p = False
            p_values = []

            for c in self.conditions:
                if c in values:
                    pvalue = values[c]
                    p_values.append(str(pvalue))

                    if pvalue <= threshold:
                        filter_p = True
                else:
                    p_values.append("1")

            if filter_p and (v in motifs_map):
                genes = "|".join(motifs_map[v])
                f.write(v + "|" + genes + "\t" + ("\t".join(p_values)) + "\n")

    def write_network(self, genes, gene_mapping_search, out_path, targets, threshold):
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
        [filtered_motifs, genes_motifs, _] = self.filter_by_genes(genes, gene_mapping_search)

        f = open(out_path + "/mapping_tf_genes.txt", "w")
        motifs_all = {}
        for gene in genes_motifs:
            motifs = genes_motifs[gene]
            for m in motifs:
                if m in motifs_all:
                    if gene not in motifs_all[m]:
                        motifs_all[m].append(gene)
                else:
                    motifs_all[m] = [gene]

                f.write(gene + "\t" + m + "\n")
        f.close()

        filtered_motifs.write_enrichment(threshold, out_path + "/pvalue_table_" + str(threshold * 100) + ".txt", motifs_all)

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
        for net_name in self.networks:
            net = self.networks[net_name]
            pairs = set()
            tfs = set()
            net_pairs[net_name] = pairs
            net_tfs[net_name] = tfs
            for tf in genes_motifs:
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
                        print("motif not in network: " + m + " " + str(tf) + " ")
            all_pairs = all_pairs.union(pairs)
            all_tfs = all_tfs.union(tfs)

        # printing out network
        for net_name in net_pairs:
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
