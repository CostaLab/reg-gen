import os
import sys
from copy import deepcopy
from GenomicRegion import GenomicRegion
from Util import GenomeData, MotifData, AuxiliaryFunctions

"""
Standardization of annotation.

Authors: Eduardo G. Gusmao.

"""

class AnnotationSet:
    """
    Annotation of genes and TFs' PWMs.
    """

    def __init__(self, gene_source, tf_source=None):
        """
        Initializes AnnotationSet.

        Keyword arguments:
        gene_source -- Gene source annotation. It will be used to create the gene_list
                       element. It can be:
            * A matrix (list of lists): An AnnotationSet will be created based on such
                 matrix.
            * A string representing a gtf file: An AnnotationSet will be created based
                 on such gtf file.
            * A string representing an organism: An AnnotationSet will be created based
                 on the gtf file for that organism in data.config file.

        tf_source -- TF source annotation. After initialization, this object is mapped with 
                     gene_list. It can be:
            * A matrix (list of lists): Represents a final tf_list element.
            * A list of mtf files: The tf_list will be created based on all mtf files.
            * A list of repositories: The tf_list will be created based on the mtf files
                associated with such repositories in data.config.
        """

        # Class Objects
        self.gene_list = [] # Represents gene annotation.
        self.tf_list = [] # Represents TF PWM annotation.

        # Initializing Required Field - Gene List
        if(isinstance(gene_source,list)): # It can be a matrix - Used by internal methods.
            self.gene_list = gene_source
        if(isinstance(gene_source,str)): # It can be a string.
            if(os.path.isfile(gene_source)): # The string may represent a path to a gtf file.
                self.load_gene_list(gene_source)
            else: # The string may represent an organism which points to a gtf file within data.config.
                genome_data = GenomeData(gene_source)
                self.load_gene_list(genome_data.get_gencode_annotation())

        # Initializing Optional Field - TF List
        if(tf_source):
            if(isinstance(tf_source,list)):
                if(isinstance(tf_source[0],list)): # It can be a matrix
                    self.tf_list = tf_source
                else:
                    mtf_file_list = []
                    motif_data = MotifData()
                    for e in tf_source:
                        if(os.path.isfile(e)): # It can be a path to a mtf file.
                            mtf_file_list.append(e)
                        else: # It can represent an organism which points to an mtf file within data.config.
                            mtf_file = motif_data.get_mtf_path(e)
                            mtf_file_list.append(mtf_file)
                    self.tf_list = self.load_tf_list(mtf_file_list)
            else: pass # TODO Throw error.

    def load_gene_list(self, file_name):
        """
        Reads gene annotation in gtf (gencode) format. It populates self.gene_list with such entries.

        Keyword arguments:
        file_name -- The gencode .gtf file name.
        
        Return: void.
        """
        
        # Opening GTF file
        try: gtf_file = open(file_name,"r")
        except Exception: pass # TODO

        # Reading GTF file
        for line in gtf_file:

            # Processing line
            line = line.strip()
            if(line[0] == "#"): continue
            line_list = line.split("\t")
            addt_list = line_list[8].split(";")
            addt_list = filter(None,addt_list)

            # Processing additional list of options
            addt_dict = dict()
            for addt_element in addt_list:
                addt_element_list = addt_element.split(" ")
                addt_element_list = filter(None,addt_element_list)
                addt_element_list[1] = addt_element_list[1].replace("\"","") # Removing " symbol from string options
                addt_dict[addt_element_list[0]] = addt_element_list[1]

            # Creating final version of additional arguments
            final_addt_list = []
            for addt_key in ["gene_id", "transcript_id", "gene_type", "gene_status", "gene_name", 
                             "transcript_type", "transcript_status", "transcript_name", "level"]:
                try: final_addt_list.append(addt_dict[addt_key])
                except Exception: final_addt_list.append(None)

            # Handling score
            current_score = 0
            if(AuxiliaryFunctions.string_is_int(line_list[5])):
                current_score = AuxiliaryFunctions.correct_standard_bed_score(line_list[5])

            # Creating GenomicRegion
            genomic_region = GenomicRegion(chrom = line_list[0], 
                                           initial = int(line_list[3])-1, 
                                           final = int(line_list[4]), 
                                           orientation = line_list[6], 
                                           data = current_score)

            # Creating final vector
            extra_index_elements = [[],[]] # One list for each: EXACT_GENE_MATCHES, INEXACT_GENE_MATCHES
            final_vector = [genomic_region,line_list[1],line_list[2],line_list[7]] + final_addt_list + extra_index_elements
            self.gene_list.append(final_vector)
            
        # Termination
        gtf_file.close()

    def load_tf_list(self, file_name_list):
        """
        Reads TF annotation in mtf (internal -- check manual) format. It populates self.tf_list with such entries.
        Everytime a TF annotation is loaded, a mapping with gene list is performed.
        
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
                while(len(line_list) < 5): line_list.append(".")
                
                # Creating final vector
                extra_index_elements = [[],[]] # One list for each: EXACT_GENE_MATCHES, INEXACT_GENE_MATCHES
                final_vector = line_list + extra_index_elements
                self.tf_list.append(final_vector)

            # Termination
            mtf_file.close()

        # Mapping gene and tf lists
        self.map_lists()

    def map_lists(self):
        """
        Maps self.gene_list with self.tf_list in various ways.
        """
        self.exact_mapping(caps=True)
        #self.inexact_mapping()

    def exact_mapping(self, caps=True):
        """
        Maps (O(n log n)) exact entries of self.gene_list's gene names with self.tf_list's gene names.
        The mapping populates self.mapping_list with
        """
        
        # Sorting lists
        self.gene_list = sorted(self.gene_list, key=lambda k: k[self.GeneField.GENE_NAMES])
        tf_assist = []
        curr_tf_index = 0
        for e in self.tf_list:
            gene_name_list = ";".join(e[self.TfField.GENE_NAMES].split("+")).split(";")
            for k in gene_name_list:
                tf_assist.append([k,curr_tf_index])
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
            try: cmp_res = cmp(gene_genename,tf_genename)
            except IndexError: break
            if(cmp_res == 0):
                gene_vec[self.GeneField.EXACT_GENE_MATCHES] = tf_vec
                tf_vec[self.TfField.EXACT_GENE_MATCHES] = gene_vec
                curr_gene_index += 1
            elif(cmp_res == -1):
                curr_gene_index += 1
            elif(cmp_res == 1):
                curr_tf_index += 1

    def inexact_mapping(self):
        """
        TODO
        """
        pass

    def get(self, query=None, list_type=DataType.GENE_LIST, return_type=ReturnType.ANNOTATION_SET):
        """
        Gets subsets of either self objects and returns different types.

        Keyword arguments:
        query -- A parameter that allows for subsets of self to be fetched. It can be:
            * None: All fields/values are going to be returned.
            * A dictionary: Subsets the desired list according to this structure. Each
                  key must be a field (please refer to AnnotationSet.GeneField or AnnotationSet.TfField)
                  that must point to a single value or a list of values.
        list_type -- Indicates which list should be subsetted/returned. Please refer to AnnotationSet.DataType.
        return_type -- Indicates what should be returned. Please refer to AnnotationSet.ReturnType.

        Return:
        result_list -- A <return_type> containing the requested <list_type> subsetted according to <query>.
        """

        # Fetching local copies of the lists
        gene_list = deepcopy(self.gene_list)
        tf_list = deepcopy(self.tf_list)

        # Deciding which list to be subsetted/returned
        if(list_type == DataType.GENE_LIST): current_list = gene_list
        elif(list_type == DataType.TF_LIST): current_list = tf_list

        # Subsetting
        if(isinstance(query,dict)):
        
            # Iterating over query elements
            for query_key in query.keys():

                # O(n) operation if a single value is being queried.
                if(not isinstance(query[query_key],list)):
                    current_list = filter(lambda k: k[query_key] == query[query_key], current_list)

                # O(n log n) operation if multiple values are being queried.
                else:

                    # Sorting input lists
                    current_list = sorted(current_list, key=lambda k: k[query_key])
                    values_list = sorted(query[query_key])

                    # Initializing index and output structures
                    result_list = []
                    curr_value_index = 0;
                    curr_list_index = 0;

                    # Linear while comparison loop
                    while True:
                        try: cmp_res = cmp(current_list[curr_list_index][query_key],values_list[curr_value_index])
                        except IndexError: break
                        if(cmp_res == 0):
                            result_list.append(current_list[curr_list_index])
                            curr_list_index += 1
                        elif(cmp_res == -1):
                            curr_list_index += 1
                        elif(cmp_res == 1):
                            curr_value_index += 1

                    # Resolving reference
                    current_list = result_list

        # Return
        if(return_type == ReturnType.ANNOTATION_SET): 
            return AnnotationSet(gene_list, tf_list)
        elif(return_type == ReturnType.LIST): 
            return current_list

    class GeneField:
        """
        Gtf fields constants.
        """

        # Gff Fields 
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
        """
        Mtf fields constants.
        """

        # Mtf Fields
        MATRIX_ID = 0
        SOURCE = 1
        VERSION = 2
        GENE_NAMES = 3
        GROUP = 4

        # Internal Fields
        EXACT_GENE_MATCHES = 5
        INEXACT_GENE_MATCHES = 6

    class DataType:
        """
        Data type constants.
        """
        GENE_LIST = 0
        TF_LIST = 1

    class ReturnType:
        """
        Return type constants.
        """
        ANNOTATION_SET = 0
        LIST = 1


