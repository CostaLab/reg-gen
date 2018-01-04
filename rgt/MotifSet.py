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
from rgt.Util import npath, MotifData, strmatch


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
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)


class MotifSet:
    """
    Represents a set of motifs. It contains MotifAnnotation instances.
    """

    def __init__(self, preload_motifs=False):
        self.motifs_map = {}
        self.networks = {}
        self.motifs_enrichment = {}
        self.conditions = []

        if preload_motifs:
            motif_data = MotifData()

            self.read_mtf(motif_data.mtf_list)

    def __len__(self):
        return len(self.motifs_map)

    def __iter__(self):
        return iter(self.motifs_map.values())

    def __getitem__(self, key):
        return self.motifs_map[key]

    def __str__(self):
        return "MotifSet:"+str(self.__dict__)

    def __repr__(self):
        return self.__str__()

    def add(self, motif):
        """
        Adds a new motif to this set.

        *Keyword arguments:*

          - motif -- Must be an instance of MotifAnnotation
        """

        self.motifs_map[motif.name] = motif

    def filter(self, keys, key_type="name", search="exact"):
        """
        Returns a new MotifSet containing all matching motifs. By default, it expects a list of motif names, but this
        can be configured via the key_type parameter.

        *Keyword arguments:*

          - values -- List of strings representing the motif to filter this set on. Actual values depend on key_type.
          - key_type -- "name" for matching on the motif name;
                        "family" for motif family/description;
                        "uniprot_ids" for matching on UniProt IDs (might be more than one);
                        "gene_names" for matching on the gene names (symbols).
          - search -- Search mode (default = 'exact'). If "exact", only perfect matches will be accepted. If
            "inexact", key inclusion will be considered a match. For example, is keys is ["ARNT"] and
            key_type="gene_names" and search="inexact", all motifs corresponding to the gene names "ARNT", "ARNT2",
            "ARNTL", etc will be included. If search="exact", only the motifs corresponding to the gene name "ARNT"
            are included.

        *Return:*

          - motif_set -- Set of filtered motifs.
          - motif2keys -- Dictionary mapping all found motif names to their corresponding attribute values. For example,
            if filtering by genes, it will provide a quick mapping between motif names and all the matching gene names.
          - key2motifs -- Inverse of motif2keys. It mapes the found values
        """

        if not isinstance(keys, list):
            raise ValueError("keys must be a list")

        valid_keys = ["name", "gene_names", "family", "uniprot_ids"]

        if key_type not in valid_keys:
            raise ValueError("key_type must be one of {}".format(valid_keys))

        motif_set = MotifSet()
        motif2keys = {}
        key2motifs = {}

        for key in keys:
            for m in self.motifs_map.values():
                attr_vals = getattr(m, key_type)
                # this is to avoid duplicating code for string-attributes and list-attributes
                if not isinstance(attr_vals, list):
                    attr_vals = [attr_vals]

                for attr_val in attr_vals:
                    if strmatch(key, attr_val, search=search):
                        motif_set.add(m)

                        # add this motif+attribute to motif2key dict
                        if m.name in motif2keys:
                            motif2keys[m.name].append(attr_val)
                        else:
                            motif2keys[m.name] = [attr_val]

                        # add this attribute+motif to key2motif dict
                        if attr_val in key2motifs:
                            key2motifs[attr_val].append(m.name)
                        else:
                            key2motifs[attr_val] = [m.name]

        return motif_set, motif2keys, key2motifs

    def read_mtf(self, mtf_filenames):
        """
        Reads TF annotation in mtf (internal format; check manual) format.

        *Keyword arguments:*

          - file_name_list -- A list with .mtf files.
        """

        # Iterating over the file name list
        for filename in mtf_filenames:

            # Opening MTF file
            mtf_file = open(npath(filename), "r")

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
            f = npath(f)
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
          - gene_mapping_search -- Gene mapping mode. See "filter/3" for more information.
          - out_path -- Output path.
          - targets -- Gene targets.
          - threshold -- Threshold for motif acceptance.
        """

        # getting mapping of genes to motifs
        [filtered_motifs, _, genes_motifs] = self.filter(genes.genes, key_type="gene_names", search=gene_mapping_search)

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


if __name__ == "__main__":
    ms = MotifSet(preload_motifs=True)
    print("gene_name ARN:", ms.filter(["ARN"], key_type="gene_names", search="inexact"))
    print("family GC:", ms.filter(["GC"], key_type="family", search="inexact"))
    print("name CEBP:", ms.filter(["CEBP"], key_type="name", search="inexact"))
    print("uniprot Q999:", ms.filter(["Q999"], key_type="uniprot_ids", search="inexact"))
