"""
MotifSet
===================
Represents a transcription factor motif and the standardization of motif annotation.

"""
# Python 3 compatibility


# Python
from glob import glob
import glob
import os

# Internal
from rgt.Util import npath, MotifData, strmatch
from rgt.motifanalysis.Motif import Motif

from MOODS import tools, parsers


class MotifAnnotation:
    """
    Represents a transcription factor with all available annotation (from MTF files).

    *Keyword arguments:*

      - tf_id -- Transcription factor ID.
      - name -- Transcription factor name (symbol).
      - database -- Database/repository in which this motif was obtained from.
      - family -- Class/Family of transcription factor motif (can be any string).
      - version -- A string representing the version of the motif.
      - gene_names -- List of gene names for this transcription factor (usually only one, sometimes two)
      - uniprot_ids -- List of UniProt accession IDs for this transcription factor (like above)
      - data_source -- A string representing the 'Source' this transcription factor was generated from,
        eg ChiP-Seq, SELEX..
      - tax_group -- A string representing the taxonomic group of the organism this transcription factor was found in (vertebrates, plants, ...)
      - species -- A string representing the species of the organism this transcription factor was found in (Homo sapiens, Mus musculus,...)
      - threshold -- A dictionary of motif matching thresholds using corresponding fpr values as keys

    """

    def __init__(self, tf_id, name, database, version, gene_names, family, uniprot_ids, data_source, tax_group, species,
                 thresholds):
        self.tf_id = tf_id
        self.name = name
        self.database = database
        self.version = version
        self.gene_names = gene_names
        self.family = family
        self.uniprot_ids = uniprot_ids
        self.data_source = data_source
        self.tax_group = tax_group
        self.species = species
        self.thresholds = thresholds

    def __str__(self):
        return str(self.__dict__)

    def __repr__(self):
        return str(self.__dict__)


class MotifSet:
    """
    Represents a set of motifs. It contains MotifAnnotation instances.

     *Keyword arguments:*

          - preload_motifs -- Must be a list of repositories from which the motifs should be taken
          if preload_motifs == None an empty MotifSet is created
          if preload_motifs == "default" all available(determined in config file) repositories are used to load motifs

          - motif_dbs -- if True, preload motifs is not a list of repositories but a list of paths to pwm files
    """

    def __init__(self, preload_motifs=None, motif_dbs=False):
        self.motifs_map = {}
        self.networks = {}
        self.motifs_enrichment = {}
        self.conditions = []

        if preload_motifs:
            if motif_dbs:
                # create empty MotifData and set attributes manually (preload_motifs is a list of paths to pwm files)
                self.motif_data = MotifData()
                self.motif_data.set_custom(preload_motifs)
                # add motifs to self.motifs_map (repositories_list = preload motifs)
                self.load_directory(self.motif_data.pwm_list)
            else:
                self.motif_data = MotifData(repositories=preload_motifs)
                # add motifs to self.motifs_map
                self.read_mtf(self.motif_data.mtf_list)
        else:
            self.motif_data = MotifData()

    def __len__(self):
        return len(self.motifs_map)

    def __iter__(self):
        return iter(list(self.motifs_map.values()))

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

    def filter(self, values, search="exact"):
        """
        Returns a new MotifSet containing all matching motifs. By default, it expects a list of motif names, but this
        can be configured via the key_type parameter.

        *Keyword arguments:*
          - values -- dictionary whose values are a list of strings representing the motif to filter this set on.
                        Actual meaning of the items' values depends on the respective key.
                        valid keys:
                        "name" for the motif name;
                        "family" for motif family/description;
                        "uniprot_ids" for UniProt IDs (might be more than one);
                        "gene_names" for the gene names (symbols);
                        "data_source" for Chip-Seq, SELEX, etc.
                        "tax_group" for taxonomic group (vertebrates, plants, ...)
                        "species" for species (Homo sapiens, Mus musculus,...)
                        "database" where the TF data is taken from
          - search -- Search mode (default = 'exact'). If "exact", only perfect matches will be accepted. If
            "inexact", key inclusion will be considered a match. For example, if keys=["ARNT"] and
            key_type="gene_names" and search="inexact", all motifs corresponding to the gene names "ARNT", "ARNT2",
            "ARNTL", etc will be included. If search="exact", only the motifs corresponding to the gene name "ARNT"
            are included.

        *Return:*

          - motif_set -- Set of filtered motifs.
        """

        if not isinstance(values, dict):
            raise ValueError("values must be a dictionary")

        valid_keys = ["name", "gene_names", "family", "uniprot_ids", "data_source", "tax_group", "species", "database"]

        for key_type in list(values.keys()):
            if not key_type in valid_keys:
                raise ValueError("wrong key-type: key_type must be one of {}".format(valid_keys))

        # TODO: maybe just ignore invalid keys instead of raising an error

        current = list(self.motifs_map.values())

        motif_set = MotifSet(preload_motifs=None)

        if not values:
            for m in current:
                motif_set.add(m)

        for key_type in list(values.keys()):

            motif_set = MotifSet(preload_motifs=None)

            for key in values[key_type]:
                for m in current:
                    attr_vals = getattr(m, key_type)
                    # this is to avoid duplicating code for string-attributes and list-attributes
                    if not isinstance(attr_vals, list):
                        attr_vals = [attr_vals]

                    for attr_val in attr_vals:
                        if strmatch(key, attr_val, search=search):
                            motif_set.add(m)
            current = list(motif_set.motifs_map.values())

        motif_set.motif_data = self.motif_data  # contains data from more motifs than those in motif_set.motifs_map

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
                        "tax_group" for matching taxonomic group (vertebrates, plants, ...)
                        "species" for matching species (Homo sapiens, Mus musculus,...)

        *Return:*

          - motif2keys -- Dictionary mapping all motif names to their corresponding attribute values. For example,
            if filtering by genes, it will provide a quick mapping between motif names and all the matching gene names.
          - key2motifs -- Inverse of motif2keys. It maps the key values to their corresponding motifs.
        """

        valid_keys = ["gene_names", "family", "uniprot_ids", "data_source", "tax_group", "species", "database"]

        if key_type not in valid_keys:
            raise ValueError("get_mappings key_type must be one of {}".format(valid_keys))

        motif2keys = {}
        key2motifs = {}

        for m in iter(self):
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

          - mtf_filenames -- A string, or a list of strings, representing .mtf file paths.
        """

        if not isinstance(mtf_filenames, list):
            mtf_filenames = [mtf_filenames]

        file_list = [filename for pattern in mtf_filenames for filename in glob.glob(npath(pattern))]

        # Iterating over the file name list
        for filename in file_list:

            database = os.path.splitext(os.path.basename(filename))[0]

            # Opening MTF file
            mtf_file = open(filename, "r")

            # Reading file
            for line in mtf_file:
                # Processing line
                line_list = line.strip().split("\t")
                tf_id = line_list[0].strip()
                name = line_list[1].strip()
                version = line_list[2].strip()
                gene_names = line_list[3].strip().split("+")
                tf_class = line_list[4].strip()
                uniprot_ids = line_list[5].strip().split(";")
                data_source = line_list[6].strip()
                tax_group = line_list[7].strip()
                species = line_list[8].strip()
                threshold_list = line_list[9].strip().split(",")
                fpr_list = [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
                thresholds = {}
                for i in range(0, 6):
                    thresholds[fpr_list[i]] = float(threshold_list[i])

                self.add(MotifAnnotation(tf_id, name, database, version, gene_names, tf_class, uniprot_ids, data_source,
                                         tax_group, species, thresholds))

            # Termination
            mtf_file.close()

    def load_directory(self, db_list):

        for directory in db_list:
            for file_name in glob.glob(directory + "/*.pwm"):
                tf_id = os.path.splitext(os.path.basename(file_name))[0]
                name = tf_id
                database = os.path.basename(directory)
                version = "0"
                gene_names = None
                tf_class = None
                uniprot_ids = None
                data_source = None
                tax_group = None
                species = None
                thresholds = {}

                self.add(MotifAnnotation(tf_id, name, database, version, gene_names, tf_class, uniprot_ids, data_source,
                                         tax_group, species, thresholds))

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
        for net_name, pairs_aux in list(net_pairs.items()):
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

    def get_motif_list(self, pseudocounts=1.0, fpr=0.0001):

        motif_list = []

        # iterate over all available PWM files
        for motif_dir_path in self.motif_data.pwm_list:

            # iterate over all motif elements in this set
            for motif_name, ma in list(self.motifs_map.items()):
                motif_file_name = os.path.join(motif_dir_path, motif_name + ".pwm")

                # if the motif annotation has a corresponding PWM file, add to return list
                if os.path.isfile(motif_file_name):
                    # check whether ma provides the motif matching threshold for the given fpr
                    # recalculate (and store) it otherwise
                    if fpr in ma.thresholds and ma.thresholds[fpr]:
                        threshold = ma.thresholds[fpr]
                    else:
                        pfm = parsers.pfm(str(motif_file_name))
                        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
                        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
                        threshold = tools.threshold_from_p(pssm, bg, fpr)
                        ma.thresholds[fpr] = threshold

                    motif_list.append(Motif(motif_file_name, pseudocounts, threshold))

        return motif_list
