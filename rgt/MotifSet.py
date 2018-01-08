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
      - data_source -- A string representing the 'Source' this transcription factor was generated from,
        eg ChiP-Seq, SELEX..
    """

    def __init__(self, tf_id, name, database, version, gene_names, family, uniprot_ids, data_source):
        self.tf_id = tf_id
        self.name = name
        self.database = database
        self.version = version
        self.gene_names = gene_names
        self.family = family
        self.uniprot_ids = uniprot_ids
        self.data_source = data_source

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
        return "MotifSet:" + str(self.__dict__)

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
                        "gene_names" for matching on the gene names (symbols);
                        "data_source" for Chip-Seq, SELEX, etc.
          - search -- Search mode (default = 'exact'). If "exact", only perfect matches will be accepted. If
            "inexact", key inclusion will be considered a match. For example, if keys=["ARNT"] and
            key_type="gene_names" and search="inexact", all motifs corresponding to the gene names "ARNT", "ARNT2",
            "ARNTL", etc will be included. If search="exact", only the motifs corresponding to the gene name "ARNT"
            are included.

        *Return:*

          - motif_set -- Set of filtered motifs.
        """

        if not isinstance(keys, list):
            raise ValueError("keys must be a list")

        valid_keys = ["name", "gene_names", "family", "uniprot_ids", "data_source"]

        if key_type not in valid_keys:
            raise ValueError("key_type must be one of {}".format(valid_keys))

        motif_set = MotifSet(preload_motifs=False)

        for key in keys:
            for m in self.motifs_map.values():
                attr_vals = getattr(m, key_type)
                # this is to avoid duplicating code for string-attributes and list-attributes
                if not isinstance(attr_vals, list):
                    attr_vals = [attr_vals]

                for attr_val in attr_vals:
                    if strmatch(key, attr_val, search=search):
                        motif_set.add(m)

        return motif_set

    def get_mappings(self, key_type="gene_names"):
        """
        Returns two dictionaries, the first mapping each motif to its specified "keys" (see filter/3 for more details),
        the second mapping each "key" to corresponding motifs.

        *Keyword arguments:*

          - key_type -- "name" for matching on the motif name;
                        "family" for motif family/description;
                        "uniprot_ids" for matching on UniProt IDs (might be more than one);
                        "gene_names" for matching on the gene names (symbols);
                        "data_source" for Chip-Seq, SELEX, etc.

        *Return:*

          - motif2keys -- Dictionary mapping all motif names to their corresponding attribute values. For example,
            if filtering by genes, it will provide a quick mapping between motif names and all the matching gene names.
          - key2motifs -- Inverse of motif2keys. It maps the key values to their corresponding motifs.
        """

        valid_keys = ["gene_names", "family", "uniprot_ids", "data_source"]

        if key_type not in valid_keys:
            raise ValueError("key_type must be one of {}".format(valid_keys))

        motif2keys = {}
        key2motifs = {}

        for m in self.motifs_map.values():
            attr_vals = getattr(m, key_type)
            # this is to avoid duplicating code for string-attributes and list-attributes
            if not isinstance(attr_vals, list):
                attr_vals = [attr_vals]

            for attr_val in attr_vals:
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

        return motif2keys, key2motifs

    def read_mtf(self, mtf_filenames):
        """
        Reads TF annotation in mtf (internal format; check manual) format.

        *Keyword arguments:*

          - file_name_list -- A string, or a list of strings, representing .mtf file paths.
        """

        if isinstance(mtf_filenames, list):
            file_list = [filename for pattern in mtf_filenames for filename in glob.glob(npath(pattern))]
        else:
            file_list = glob.glob(npath(mtf_filenames))

        # Iterating over the file name list
        for filename in file_list:

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
                data_source = line_list[7].strip()

                self.add(MotifAnnotation(tf_id, name, database, version, gene_names, tf_class, uniprot_ids, data_source))

            # Termination
            mtf_file.close()

    def read_enrichment(self, enrichment_files, threshold=1):
        """
        Reads current output of motif enrichment analysis to get gene targets.

        *Keyword arguments:*

          - enrichment_files -- One string, or a list of strings, representing enrichment file paths.
          - threshold -- P-value threshold for motif acceptance.
        """

        if isinstance(enrichment_files, list):
            file_list = [filename for pattern in enrichment_files for filename in glob.glob(npath(pattern))]
        else:
            file_list = glob.glob(npath(enrichment_files))

        # reading networks
        for filename in file_list:
            # use last dir name as name for condition
            condition = os.path.dirname(filename)
            condition = condition.split("/")[-1]
            self.conditions.append(condition)

            network = {}

            f = open(filename, "r")

            # skip header
            next(f)

            for line in f:
                line = line.strip("\n")
                values = line.split("\t")
                motif = values[0]

                if motif in self.motifs_map:
                    p_value = float(values[2])
                    genes = values[9].split(",")

                    if threshold >= p_value:
                        network[motif] = genes

                    if motif in self.motifs_enrichment:
                        self.motifs_enrichment[motif][condition] = p_value
                    else:
                        self.motifs_enrichment[motif] = {condition: p_value}
                else:
                    print("motif not found: " + motif)

            self.networks[condition] = network

            f.close()

    def write_enrichment(self, out_file, threshold=1):
        """
        Writes enrichment table for network generation.

        *Keyword arguments:*

          - out_file -- Output file name.
          - threshold -- P-value threshold for motif acceptance.
        """

        f = open(npath(out_file), "w")
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

            if filter_p and (v in self.motifs_map) and self.motifs_map[v].gene_names:
                genes = "|".join(self.motifs_map[v].gene_names)
                f.write(v + "|" + genes + "\t" + ("\t".join(p_values)) + "\n")

    def write_network(self, targets, out_path, threshold=1):
        """
        If enrichment information has been loaded before (via read_enrichment), this function creates
        a cytoscape-compatible network into the output folder.

        *Keyword arguments:*

          - targets -- Gene targets.
          - out_path -- Output path.
          - threshold -- Threshold for motif acceptance.
        """

        self.write_enrichment(out_path + "/pvalue_table_" + str(threshold * 100) + ".txt", threshold)

        out_path = npath(out_path)

        _, genes_motifs = self.get_mappings(key_type="gene_names")

        net_pairs = {}
        net_tfs = {}
        all_pairs = set()
        all_tfs = set()
        all_genes = set()

        if targets:
            filter_targets = True
        else:
            filter_targets = False

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
        for net_name, pairs_aux in net_pairs.items():
            f = open(out_path + "/" + net_name + "_targets.txt", "w")

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
