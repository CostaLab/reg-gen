#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
AnnotationSet
===================
AnnotationSet represent genomic annotation from genes.

"""

# Python
import os

# Internal
from rgt.GenomicRegionSet import *
from rgt.Util import GenomeData, MotifData, AuxiliaryFunctions, cmp


class AnnotationSet:
    """This class represents genomic annotation from genes.

    *Keyword arguments:*

        - gene_source -- Gene source annotation. It will be used to create the gene_list element. It can be:
            - A matrix (list of lists): An AnnotationSet will be created based on such matrix.
            - A string representing a gtf file: An AnnotationSet will be created based on such gtf file.
            - A string representing an organism: An AnnotationSet will be created based on the gtf file for that organism in data.config file.

        - tf_source -- TF source annotation. After initialization, this object is mapped with gene_list. It can be:
            - A matrix (list of lists): Represents a final tf_list element.
            - A list of mtf files: The tf_list will be created based on all mtf files.
            - A list of repositories: The tf_list will be created based on the mtf files associated with such repositories in data.config.

        - alias_source -- Alias dictionary source annotation. It can be:
            - A dictionary: An alias dictionary will be created based on such dictionary.
            - A string representing a alias (txt) file: An alias dictionary will be created based on such txt file.
            - A string representing an organism: An alias dictionary will be created based on the txt file for that organism in data.config file.
    """

    class GeneField:
        """Gtf fields constants.

        *Constants:*

            - GENOMIC_REGION.
            - ANNOTATION_SOURCE.
            - FEATURE_TYPE.
            - GENOMIC_PHASE.
            - GENE_ID.
            - TRANSCRIPT_ID.
            - GENE_TYPE.
            - GENE_STATUS.
            - GENE_NAMES.
            - TRANSCRIPT_TYPE.
            - TRANSCRIPT_STATUS.
            - TRANSCRIPT_NAME.
            - LEVEL.
            - EXACT_GENE_MATCHES.
            - INEXACT_GENE_MATCHES.
        """

        # Gff Fields 
        def __init__(self):
            pass

        GENOMIC_REGION = 0
        ANNOTATION_SOURCE = 1
        FEATURE_TYPE = 2
        GENOMIC_PHASE = 3

        # Gtf Fields
        GENE_ID = 4
        TRANSCRIPT_ID = 5
        GENE_TYPE = 6
        GENE_STATUS = 7
        GENE_NAMES = 8
        TRANSCRIPT_TYPE = 9
        TRANSCRIPT_STATUS = 10
        TRANSCRIPT_NAME = 11
        LEVEL = 12

        # Internal Fields
        EXACT_GENE_MATCHES = 13
        INEXACT_GENE_MATCHES = 14

    class TfField:
        """Mtf fields constants.

        *Constants:*

            - MATRIX_ID.
            - SOURCE.
            - VERSION.
            - GENE_NAMES.
            - GROUP.
            - EXACT_GENE_MATCHES.
            - INEXACT_GENE_MATCHES.
        """

        # Mtf Fields
        def __init__(self):
            pass

        MATRIX_ID = 0
        SOURCE = 1
        VERSION = 2
        GENE_NAMES = 3
        GROUP = 4

        # Internal Fields
        EXACT_GENE_MATCHES = 5
        INEXACT_GENE_MATCHES = 6

    class DataType:
        """Data type constants.

        *Constants:*

            - GENE_LIST.
            - TF_LIST.
        """

        def __init__(self):
            pass

        GENE_LIST = 0
        TF_LIST = 1

    class ReturnType:
        """Return type constants.

        *Constants:*

            - ANNOTATION_SET.
            - LIST.
        """

        def __init__(self):
            pass

        ANNOTATION_SET = 0
        LIST = 1

    def __init__(self, gene_source, tf_source=None, alias_source=None,
                 filter_havana=False, protein_coding=False, known_only=False):

        # Class Objects
        self.gene_list = []  # Represents gene annotation.
        self.tf_list = []  # Represents TF PWM annotation.
        self.alias_dict = dict()  # Gene Symbol or other IDs -> ENSEMBL ID
        self.symbol_dict = dict()  # ENSEMBL ID -> Official gene symbol

        # Initializing Required Field - Gene List
        if isinstance(gene_source, list):  # It can be a matrix - Used by internal methods.
            self.gene_list = gene_source
        if isinstance(gene_source, str):  # It can be a string.
            if os.path.isfile(gene_source):  # The string may represent a path to a gtf file.
                # FTT for TDF True
                # filter_havana = False
                protein_coding = False
                known_only = False
                self.load_gene_list(gene_source,
                                    filter_havana=filter_havana,
                                    protein_coding=protein_coding,
                                    known_only=known_only)
            else:  # The string may represent an organism which points to a gtf file within data.config.
                genome_data = GenomeData(gene_source)
                self.load_gene_list(genome_data.get_annotation(),
                                    filter_havana=filter_havana,
                                    protein_coding=protein_coding,
                                    known_only=known_only)

        # Initializing Optional Field - TF List
        if tf_source:
            if isinstance(tf_source, list):
                if isinstance(tf_source[0], list):  # It can be a matrix
                    self.tf_list = tf_source
                else:
                    mtf_file_list = []
                    motif_data = MotifData()
                    for e in tf_source:
                        if os.path.isfile(e):  # It can be a path to a mtf file.
                            mtf_file_list.append(e)
                        else:  # It can represent an organism which points to an mtf file within data.config.
                            mtf_file = motif_data.get_mtf_path(e)
                            mtf_file_list.append(mtf_file)
                    self.load_tf_list(mtf_file_list)
            else:
                pass  # TODO Throw error.

        # Initializing Optional Field - Alias Dictionary
        if alias_source:
            if isinstance(alias_source, dict):  # It can be a dictionary - Used by internal methods.
                self.alias_dict = alias_source
            if isinstance(alias_source, str):  # It can be a string.
                if os.path.isfile(alias_source):  # The string may represent a path to a txt alias file.
                    self.load_alias_dict(alias_source)
                else:  # The string may represent an organism which points to a txt alias file within data.config.
                    genome_data = GenomeData(alias_source)
                    self.load_alias_dict(genome_data.get_gene_alias())
            else:
                pass  # TODO Throw error

    def load_gene_list(self, file_name, filter_havana=True, protein_coding=False, known_only=False):
        """Reads gene annotation in gtf (gencode) format. It populates self.gene_list with such entries.
        
        *Keyword arguments:*

            - file_name -- The gencode .gtf file name.
        """
        # Opening GTF file
        try:
            gtf_file = open(file_name, "r")
        except Exception:
            print("Error: Cannot find the annotation file: " + file_name)
            print("Please check the path in ~/rgtdata/data.config")
            sys.exit(1)

        # Reading GTF file
        for line in gtf_file:

            # Processing line
            line = line.strip()
            if line[0] == "#": continue
            line_list = line.split("\t")
            try:
                if filter_havana and line_list[1] == "HAVANA": continue
            except:
                pass

            addt_list = line_list[8].split(";")
            addt_list = [_f for _f in addt_list if _f]

            # Processing additional list of options
            addt_dict = dict()
            for addt_element in addt_list:
                addt_element_list = addt_element.split(" ")
                addt_element_list = [_f for _f in addt_element_list if _f]
                # Removing " symbol from string options
                addt_element_list[1] = addt_element_list[1].replace("\"", "")
                addt_dict[addt_element_list[0]] = addt_element_list[1]

            # filter non-protein-coding sequences, if required
            if protein_coding:
                if "gene_type" not in addt_dict or addt_dict["gene_type"] != "protein_coding":
                    continue
                if "transcript_type" in addt_dict and addt_dict["transcript_type"] != "protein_coding":
                    continue

            # filter unknown sequences, if required
            if known_only:
                if "gene_status" not in addt_dict or addt_dict["gene_status"] != "KNOWN":
                    continue
                if "transcript_status" in addt_dict and addt_dict["transcript_status"] != "KNOWN":
                    continue

            # Removing dot from IDs
            addt_dict["gene_id"] = addt_dict["gene_id"].split(".")[0]
            try:
                addt_dict["transcript_id"] = addt_dict["transcript_id"].split(".")[0]
            except:
                pass

            # Creating final version of additional arguments
            final_addt_list = []
            for addt_key in ["gene_id", "transcript_id", "gene_type", "gene_status", "gene_name",
                             "transcript_type", "transcript_status", "transcript_name", "level"]:
                try:
                    final_addt_list.append(addt_dict[addt_key])
                except Exception:
                    final_addt_list.append(None)

            # Handling score
            current_score = 0
            if AuxiliaryFunctions.string_is_int(line_list[5]):
                current_score = AuxiliaryFunctions.correct_standard_bed_score(line_list[5])

            # Creating GenomicRegion
            genomic_region = GenomicRegion(chrom=line_list[0],
                                           initial=int(line_list[3]) - 1,
                                           final=int(line_list[4]),
                                           orientation=line_list[6],
                                           data=current_score)

            # Creating final vector
            extra_index_elements = [[], []]  # One list for each: EXACT_GENE_MATCHES, INEXACT_GENE_MATCHES
            final_vector = [genomic_region, line_list[1], line_list[2],
                            line_list[7]] + final_addt_list + extra_index_elements
            self.gene_list.append(final_vector)

        # Termination
        gtf_file.close()

    def load_alias_dict(self, file_name):
        """Reads an alias.txt file and creates a dictionary to translate gene symbols/alternative IDs to ensembl gene ID
        
        *Keyword arguments:*

            - file_name -- Alias file name.
        """

        # Opening alias file
        alias_file = open(file_name, "r")

        # Iterating over alias file entries
        for line in alias_file:
            ll = line.strip().split("\t")
            ensembl_id = ll[0]
            official_name = ll[1]
            alias_vec = ll[2].split("&")
            self.symbol_dict[ensembl_id] = official_name
            for e in alias_vec:
                try:
                    self.alias_dict[e].append(ensembl_id)
                except Exception:
                    self.alias_dict[e] = [ensembl_id]

        # Termination
        alias_file.close()

    def load_tf_list(self, file_name_list):
        """Reads TF annotation in mtf (internal -- check manual) format. It populates self.tf_list with such entries. Everytime a TF annotation is loaded, a mapping with gene list is performed.
        
        *Keyword arguments:*

            - file_name_list -- A list with .mtf files.
        """

        # Iterating over the file name list
        for file_name in file_name_list:

            # Opening MTF file
            try:
                mtf_file = open(file_name, "r")
            except Exception:
                pass  # TODO

            # Reading file
            for line in mtf_file:

                # Processing line
                line_list = line.strip().split("\t")
                while len(line_list) < 5: line_list.append(".")

                # Creating final vector
                extra_index_elements = [[], []]  # One list for each: EXACT_GENE_MATCHES, INEXACT_GENE_MATCHES
                final_vector = line_list + extra_index_elements
                self.tf_list.append(final_vector)

            # Termination
            mtf_file.close()

        # Mapping gene and tf lists
        self.map_lists()

    def map_lists(self):
        """Maps self.gene_list with self.tf_list in various ways."""
        self.exact_mapping()
        # self.inexact_mapping()

    def exact_mapping(self):
        """Maps (O(n log n)) exact entries of self.gene_list's gene names with self.tf_list's gene names."""

        # Sorting lists
        self.gene_list = sorted(self.gene_list, key=lambda k: k[self.GeneField.GENE_NAMES])
        tf_assist = []
        curr_tf_index = 0
        for e in self.tf_list:
            gene_name_list = ";".join(e[self.TfField.GENE_NAMES].split("+")).split(";")
            for k in gene_name_list:
                tf_assist.append([k, curr_tf_index])
            curr_tf_index += 1
        tf_assist = sorted(tf_assist, key=lambda k: k[0])

        # Initializing index
        curr_gene_index = 0
        curr_tf_index = 0

        # Linear while comparison loop
        while True:
            gene_genename = self.gene_list[curr_gene_index][self.GeneField.GENE_NAMES]
            tf_genename = tf_assist[0]
            gene_vec = self.gene_list[curr_gene_index]
            tf_vec = self.tf_list[tf_assist[1]]
            try:
                cmp_res = cmp(gene_genename, tf_genename)
            except IndexError:
                break
            if cmp_res == 0:
                gene_vec[self.GeneField.EXACT_GENE_MATCHES] = tf_vec
                tf_vec[self.TfField.EXACT_GENE_MATCHES] = gene_vec
                curr_gene_index += 1
            elif cmp_res == -1:
                curr_gene_index += 1
            elif cmp_res == 1:
                curr_tf_index += 1

    def inexact_mapping(self):
        """Comming soon!"""
        pass  # TODO

    def fix_gene_names(self, gene_set, output_dict=False, mute_warn=True):
        """
        Checks if all gene names in gene_set are ensembl IDs. If a gene is not in ensembl format, it will be
        converted using alias_dict. If the gene name cannot be found then it is reported in a separate gene_set.
        
        *Keyword arguments:*

            - gene_set -- A GeneSet object.
            - output_dict -- Also output the mapping dictionary (default = False).
            - mute_warn -- Do not print warnings regarding genes that mapped to multiple entries (default = True).
        
        *Return:*

            - mapped_gene_list -- A list of ensembl IDs
            - unmapped_gene_list -- A list of unmapped gene symbols/IDs
        """

        # Creating resulting lists
        mapped_gene_list = []
        unmapped_gene_list = []

        if output_dict: mapping_dict = {}
        # Iterating on gene names
        for gene_name in gene_set.genes:

            # Verifying if it is ensembl name
            flag_ensembl = False
            try:
                if len(gene_name) >= 15 and gene_name[:3] == "ENS": flag_ensembl = True
            except Exception:
                flag_ensembl = False
            if flag_ensembl:
                mapped_gene_list.append(gene_name)
                sym = self.get_official_symbol(gene_name)
                if output_dict and sym and gene_name:
                    mapping_dict[gene_name] = sym
            else:
                try:
                    alias_list = self.alias_dict[gene_name.upper()]
                    if len(alias_list) > 1:
                        if not mute_warn:
                            print("Warning: The gene " + gene_name +
                                  " contains more than one matching IDs, both will be used.")
                    for e in alias_list:
                        mapped_gene_list.append(e)
                        if output_dict and e and gene_name:
                            mapping_dict[e] = gene_name

                except Exception:
                    unmapped_gene_list.append(gene_name)
        if output_dict:
            return mapped_gene_list, unmapped_gene_list, mapping_dict
        else:
            return mapped_gene_list, unmapped_gene_list

    def get(self, query=None, list_type=DataType.GENE_LIST, return_type=ReturnType.ANNOTATION_SET):
        """Gets subsets of either self objects and returns different types.

        *Keyword arguments:*

            - query -- A parameter that allows for subsets of self to be fetched. It can be:
                - None: All fields/values are going to be returned.
                - A dictionary: Subsets the desired list according to this structure. Each
                      key must be a field (please refer to AnnotationSet.GeneField or AnnotationSet.TfField)
                      that must point to a single value or a list of values.
            - list_type -- Indicates which list should be subsetted/returned. Please refer to AnnotationSet.DataType.
            - return_type -- Indicates what should be returned. Please refer to AnnotationSet.ReturnType.

        *Return:*

            - result_list -- A <return_type> containing the requested <list_type> subsetted according to <query>.
        """

        # Fetching local copies of the lists
        gene_list = deepcopy(self.gene_list)
        tf_list = deepcopy(self.tf_list)
        current_list = None

        # Deciding which list to be subsetted/returned
        if list_type == self.DataType.GENE_LIST:
            current_list = gene_list
        elif list_type == self.DataType.TF_LIST:
            current_list = tf_list

        # Subsetting
        if isinstance(query, dict):

            # Iterating over query elements
            for query_key in list(query.keys()):

                # O(n) operation if a single value is being queried.
                if not isinstance(query[query_key], list):
                    current_list = [k for k in current_list if k[query_key] == query[query_key]]

                # O(n log n) operation if multiple values are being queried.
                else:

                    # Sorting input lists
                    current_list = sorted(current_list, key=lambda k: k[query_key])
                    values_list = sorted(query[query_key])

                    # Initializing index and output structures
                    result_list = []
                    curr_value_index = 0
                    curr_list_index = 0

                    # Linear while comparison loop
                    while True:
                        try:
                            cmp_res = cmp(current_list[curr_list_index][query_key], values_list[curr_value_index])
                        except IndexError:
                            break
                        if cmp_res == 0:
                            result_list.append(current_list[curr_list_index])
                            curr_list_index += 1
                        elif cmp_res == -1:
                            curr_list_index += 1
                        elif cmp_res == 1:
                            curr_value_index += 1

                    # Resolving reference
                    current_list = result_list

        # Deciding which list to be subsetted/returned
        if list_type == self.DataType.GENE_LIST:
            gene_list = current_list
        elif list_type == self.DataType.TF_LIST:
            tf_list = current_list

        # Return
        if return_type == self.ReturnType.ANNOTATION_SET:
            # FIXME is deepcopy necessary here? seems overkill
            return AnnotationSet(gene_list, tf_source=tf_list, alias_source=deepcopy(self.alias_dict))
        elif return_type == self.ReturnType.LIST:
            return current_list

    def get_official_symbol(self, gene_name_source):
        """Returns the official symbol(s) from gene_name_source.

        *Keyword arguments:*

            - gene_source -- It can be a string (single gene name) or a GeneSet (multiple genes).

        *Return:*

            - if gene_source is string then returns the converted string gene name or None if gene name could not be converted.
            - if gene_source is list then returns two lists containing, respectively, converted and not-converted gene names.
        """

        if isinstance(gene_name_source, str):
            try:
                return self.symbol_dict[gene_name_source]
            except Exception:
                return None
        else:
            if isinstance(gene_name_source, list):
                curr_list = gene_name_source
            else:
                curr_list = gene_name_source.genes
            mapped_list = []
            unmapped_list = []
            for e in curr_list:
                try:
                    mapped_list.append(self.symbol_dict[e])
                except:
                    try:
                        mapped_list.append(self.symbol_dict[e.split(".")[0]])
                    except Exception:
                        unmapped_list.append(e)
            return mapped_list, unmapped_list

    def get_promoters(self, promoter_length=1000, tss=0, gene_set=None, unmaplist=False, variants=False, gene_id=False,
                      regiondata=False):
        """
        Gets promoters of genes given a specific promoter length. It returns a GenomicRegionSet with such promoters.
        The ID of each gene will be put in the NAME field of each GenomicRegion.
        Each promoter includes also the coordinate of the 5' base pair, therefore each promoter actual
        length is promoter_length+1.

        *Keyword arguments:*

            - promoter_length -- The length of the promoter region.
            - gene_set -- A set of genes to narrow the search.
            - unmaplist -- If True than also return the unmappable genes list (default = False).

        *Return:*

            - result_grs -- A GenomicRegionSet containing the promoters.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None

        if gene_set:
            mapped_gene_list, unmapped_gene_list, mapping_dict = self.fix_gene_names(gene_set, output_dict=True)

        # Fetching genes

        if not variants:
            target = "gene"
        else:
            target = "transcript"
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: target,
                                self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: target}

        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("promoters")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            if gr.orientation == "+":
                gr.final = gr.initial + 1 + tss
                gr.initial = gr.initial - promoter_length
            else:
                gr.initial = gr.final - 1 - tss
                gr.final = gr.initial + promoter_length + 1

            if gene_set:
                try:
                    gr.name = mapping_dict[e[self.GeneField.GENE_ID]]
                except:
                    gr.name = e[self.GeneField.GENE_ID]
            elif gene_id:
                gr.name = e[self.GeneField.GENE_ID]
            else:
                gr.name = e[self.GeneField.GENE_NAMES]

            if gene_set and regiondata:
                gr.data = gene_set.values[gr.name]
            result_grs.add(gr)

        if unmaplist:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_tss(self, gene_set=None):
        """Gets TSS(Transcription start site) of genes. It returns a GenomicRegionSet with such TSS. The ID of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - gene_set -- A set of genes to narrow the search.
        
        *Return:*

            - result_grs -- A GenomicRegionSet containing TSS.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None
        if gene_set: mapped_gene_list, unmapped_gene_list = self.fix_gene_names(gene_set)

        # Fetching genes
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene", self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("TSS")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            if gr.orientation == "+":
                gr.initial = gr.initial
                gr.final = gr.initial + 1
            else:
                gr.initial = gr.final - 1
                gr.final = gr.final
            gr.name = e[self.GeneField.GENE_ID]
            result_grs.add(gr)
        result_grs.merge()
        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_tts(self, gene_set=None):
        """Gets TTS(Transcription termination site) of genes. It returns a GenomicRegionSet with such TTS. The ID of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - gene_set -- A set of genes to narrow the search.
        
        *Return:*

            - result_grs -- A GenomicRegionSet containing TTS.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None
        if gene_set: mapped_gene_list, unmapped_gene_list = self.fix_gene_names(gene_set)

        # Fetching genes
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene", self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("TTS")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            if gr.orientation == "+":
                gr.initial = gr.initial
                gr.final = gr.initial + 1
            else:
                gr.initial = gr.final - 1
                gr.final = gr.final
            gr.name = e[self.GeneField.GENE_ID]
            result_grs.add(gr)
        result_grs.merge()
        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_exons(self, start_site=False, end_site=False, gene_set=None, merge=True):
        """Gets exons of genes. It returns a GenomicRegionSet with such exons. The id of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - start_site -- Whether to relocate the start sites.
            - end_site -- Whether to relocate the end sites.
            - gene_set -- A set of genes to narrow the search.

        *Return:*

            - result_grs -- A GenomicRegionSet containing the exons.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None
        if gene_set: mapped_gene_list, unmapped_gene_list = self.fix_gene_names(gene_set)

        # Fetching exons
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "exon", self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "exon"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("exon")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            # gr.name = e[self.GeneField.GENE_ID]
            gr.name = e[self.GeneField.TRANSCRIPT_ID]
            result_grs.add(gr)
        if start_site:
            result_grs.relocate_regions("leftend", left_length=1, right_length=1)
        elif end_site:
            result_grs.relocate_regions("rightend", left_length=1, right_length=1)
        if merge: result_grs.merge()
        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_genes(self, gene_set=None):
        """Gets regions of genes. It returns a GenomicRegionSet with such genes. The id of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - gene_set -- A set of genes to narrow the search.

        *Return:*

            - result_grs -- A GenomicRegionSet containing the genes.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None
        if gene_set: mapped_gene_list, unmapped_gene_list = self.fix_gene_names(gene_set)

        # Fetching genes
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene", self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "gene"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("genes")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            gr.name = e[self.GeneField.GENE_ID]
            result_grs.add(gr)
        result_grs.merge()
        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_introns(self, start_site=False, end_site=False, gene_set=None):
        """Gets introns of genes. It returns a GenomicRegionSet with such introns. The id of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - start_site -- Whether to relocate the start sites.
            - end_site -- Whether to relocate the end sites.
            - gene_set -- A set of genes to narrow the search.

        *Return:*

            - result_grs -- A GenomicRegionSet containing the introns.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        if gene_set:
            genes, unmapped_gene_list = self.get_genes(gene_set=gene_set)
            exons, unmapped_gene_list_e = self.get_exons(gene_set=gene_set)
        else:
            genes = self.get_genes()
            exons = self.get_exons()

        result_grs = genes.subtract(exons)
        if start_site:
            result_grs.relocate_regions("leftend", left_length=1, right_length=1)
        elif end_site:
            result_grs.relocate_regions("rightend", left_length=1, right_length=1)
        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_transcripts(self, gene_set=None):
        """Gets transcripts of genes. It returns a GenomicRegionSet with such transcripts. The id of each gene will be put in the NAME field of each GenomicRegion.

        *Keyword arguments:*

            - gene_set -- A set of genes to narrow the search.

        *Return:*

            - result_grs -- A GenomicRegionSet containing the exons.
            - unmapped_gene_list -- A list of genes that could not be mapped to an ENSEMBL ID.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None
        if gene_set: mapped_gene_list, unmapped_gene_list = self.fix_gene_names(gene_set)

        # Fetching exons
        if gene_set:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "exon", self.GeneField.GENE_ID: mapped_gene_list}
        else:
            query_dictionary = {self.GeneField.FEATURE_TYPE: "exon"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("exon")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            gr.name = e[self.GeneField.TRANSCRIPT_ID]
            result_grs.add(gr)

        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs

    def get_biotypes(self, gene_set=None):
        """Get the region sets of different Biotypes.

        *Keyword arguments:*

        *Return:*

            - result_grs -- A list of GenomicRegionSets containing the regions for each Biotype.
        """

        # Fetching gene names
        mapped_gene_list = None
        unmapped_gene_list = None

        # Fetching exons
        query_dictionary = {self.GeneField.FEATURE_TYPE: "exon"}
        query_annset = self.get(query_dictionary)

        # Creating GenomicRegionSet
        result_grs = GenomicRegionSet("exon")
        for e in query_annset.gene_list:
            gr = e[self.GeneField.GENOMIC_REGION]
            gr.name = e[self.GeneField.TRANSCRIPT_ID]
            result_grs.add(gr)

        if gene_set:
            return result_grs, unmapped_gene_list
        else:
            return result_grs
