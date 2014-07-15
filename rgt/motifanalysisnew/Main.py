
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
from glob import glob
from multiprocessing import Pool
from pickle import load, dump

# Internal
from .. Util import PassThroughOptionParser, ErrorHandler, MotifData, GenomeData
from .. ExperimentalMatrix import ExperimentalMatrix
from .. GeneSet import GeneSet
from .. GenomicRegion import GenomicRegion
from .. GenomicRegionSet import GenomicRegionSet
from Motif import Motif
from Match import match
from Statistics import multiple_test_correction

# External
from Bio import motifs
from Bio.Seq import Seq
from pysam import Fastafile
from fisher import pvalue

"""
Contains functions to common motif analyses.

Basic Input:
- XXX TODO

Dependencies:
- XXX TODO
- bedToBigBed, XXXX scripts in $PATH (if the option is used)

Authors: Eduardo G. Gusmao.
"""

def main():
    """
    Main function that redirects tool usage.

    Keyword arguments: None
        
    Return: None
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    current_version = "0.0.1"
    usage_message = ("\n--------------------------------------------------\n"
                     "The motif analysis program performs various motif-based analyses. "
                     "In order to use these tools, please type: \n\n"
                     "%prog [analysis type] [options]\n\n"
                     "Where [analysis type] refers to the type of the motif analysis performed "
                     "and [options] are the analysis-specific arguments.\n\n"
                     "Bellow you can find all current available analysis types. "
                     "To check the analyses specific options, please use:\n\n"
                     "%prog [analysis type] -h\n\n"
                     "For more information, please refer to our wiki:\n\n"
                     "https://code.google.com/p/reg-gen/wiki/RegGen\n\n"
                     "--------------------------------------------------\n\n"
                     "Options:\n"
                     "--version     show program's version number and exit.\n"
                     "-h, --help    show this help message and exit.\n"
                     "--matching    Performs motif matching analysis.\n"
                     "--enrichment  Performs motif enrichment analysis.\n")
    version_message = "Motif Analysis - Regulatory Analysis Toolbox (RGT). Version: "+str(current_version)

    # Processing Help/Version Options
    if(len(sys.argv) <= 1 or sys.argv[1] == "-h" or sys.argv[1] == "--help"):
        print(usage_message)
        sys.exit(0)
    elif(sys.argv[1] == "--version"):
        print(version_message)
        sys.exit(0)

    # Initializing Error Handler
    main_error_handler = ErrorHandler()

    ###################################################################################################
    # Redirecting to Specific Functions
    ###################################################################################################

    # Redirecting Loop
    if(sys.argv[1] == "--matching"):
        main_matching()
    elif(sys.argv[1] == "--enrichment"):
        main_enrichment()
    else:
        main_error_handler.throw_error("MOTIF_ANALYSIS_OPTION_ERROR")

def main_matching():
    """
    Performs motif matching.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    usage_message = "%prog --matching [options] <experiment_matrix>"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage = usage_message)

    # Parameters Options
    parser.add_option("--organism", dest = "organism", type = "string", metavar="STRING", default = "hg19",
                      help = ("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file."))
    parser.add_option("--fpr", dest = "fpr", type = "float", metavar="FLOAT", default = 0.0001,
                      help = ("False positive rate cutoff for motif matching."))
    parser.add_option("--precision", dest = "precision", type = "int", metavar="INT", default = 10000,
                      help = ("Score distribution precision for motif matching."))
    parser.add_option("--pseudocounts", dest = "pseudocounts", type = "float", metavar="FLOAT", default = 0.1,
                      help = ("Pseudocounts to be added to raw counts of each PWM."))
    parser.add_option("--rand-proportion", dest = "rand_proportion", type = "float", metavar="FLOAT", default = 10.0,
                      help = ("If random coordinates need to be created, then it will be created a number of "
                              "coordinates that equals this parameter x the number of input regions. If zero (0) "
                              "is passed, then no random coordinates are created."))
    parser.add_option("--processes", dest = "processes", type = "int", metavar="INT", default = 1,
                      help = ("Number of processes for multi-CPU based machines."))

    # Output Options
    parser.add_option("--output-location", dest = "output_location", type = "string", metavar="PATH", default = os.getcwd(),
                      help = ("Path where the output files will be written."))
    parser.add_option("--bigbed", dest = "bigbed", action = "store_true", default = False,
                      help = ("If this option is used, all bed files will be written as bigbed."))

    # Processing Options
    options, arguments = parser.parse_args()

    # Additional Parameters
    matching_folder_name = "Match"
    random_region_name = "random_regions"
    dump_file_name = "dump.p"

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    matching_output_location = os.path.join(options.output_location,matching_folder_name)
    try:
        if(not os.path.isdir(matching_output_location)): os.makedirs(matching_output_location)
    except Exception: pass # TODO ERROR

    # Default genomic data
    genome_data = GenomeData(options.organism)

    # Default motif data
    motif_data = MotifData()

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading arguments
    try:
        input_matrix = arguments[0]
        if(len(arguments) > 1): pass # TODO WARNING
    except Exception: pass # TODO ERROR

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception: pass # TODO ERROR

    ###################################################################################################
    # Reading Regions
    ###################################################################################################

    # Initialization
    max_region_len = 0
    max_region = None
    input_regions = []

    try:
        exp_matrix_objects_dict = exp_matrix.objectsDict
    except Exception: pass # TODO ERROR

    # Iterating on experimental matrix objects
    for k in exp_matrix_objects_dict.keys():

        curr_genomic_region = exp_matrix_objects_dict[k]

        # If the object is a GenomicRegionSet
        if(isinstance(curr_genomic_region,GenomicRegionSet)):

            # Sorting input region
            curr_genomic_region.sort()

            # Append label and GenomicRegionSet
            input_regions.append(curr_genomic_region)

            # Verifying max_region_len for random region generation
            curr_len = len(curr_genomic_region)
            if(curr_len > max_region_len):
                max_region_len = curr_len
                max_region = exp_matrix_objects_dict[k]

    ###################################################################################################
    # Creating random region
    ###################################################################################################

    # Create random coordinates
    rand_region = None
    if(options.rand_proportion > 0):

        # Create random coordinates and name it random_regions
        rand_region = max_region.random_regions(options.organism, multiply_factor=options.rand_proportion, chrom_X=True)
        rand_region.sort()
        rand_region.name = random_region_name

        # Put random regions in the end of the input regions
        input_regions.append(rand_region)

        # Writing random regions
        output_file_name = os.path.join(matching_output_location, random_region_name)
        rand_bed_file_name = output_file_name+".bed"
        rand_region.write_bed(rand_bed_file_name)

        # Verifying condition to write bb
        if(options.bigbed):

            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            rand_bb_file_name = output_file_name+".bb"
            try:
                os.system(" ".join(["bedToBigBed", rand_bed_file_name, chrom_sizes_file, rand_bb_file_name, "-verbose=0"]))
                os.remove(rand_bed_file_name)
            except Exception: pass # WARNING

    else: pass # TODO WARNING

    ###################################################################################################
    # Creating PWMs
    ###################################################################################################

    # TODO XXX Alterar setup.py para colocar os repositorios corretos de volta.

    # Initialization
    motif_list = []

    # Fetching list with all motif file names
    motif_file_names = []
    for motif_repository in motif_data.get_pwm_list():
        for motif_file_name in glob(os.path.join(motif_repository,"*.pwm")):
            motif_file_names.append(motif_file_name)

    # Grouping motif file names by the number of processes requested
    if(options.processes <= 0): pass # TODO ERROR
    elif(options.processes == 1): motif_file_names = [[e] for e in motif_file_names]
    else: motif_file_names = map(None, *(iter(motif_file_names),) * options.processes)

    # Iterating on grouped file name list
    for file_group in  motif_file_names:

        # Creating input data for multiprocessing
        curr_data_input = []
        curr_proc_nb = 0
        for motif_file_name in file_group:
            if(motif_file_name):
                curr_data_input.append([motif_file_name, options.pseudocounts, options.precision, options.fpr])
                curr_proc_nb += 1

        # Applying the Motif creation function with multiprocessing
        pool = Pool(curr_proc_nb)
        curr_motif_list = pool.map(Motif,curr_data_input)
        pool.close()
        pool.join()
        
        # Append curr_motif_list (motif_group -- group of Motif objects) to motif_list
        motif_list.append(curr_motif_list)

    ###################################################################################################
    # Motif Matching
    ###################################################################################################

    # Creating output structure
    # GENOMIC REGION -> FACTOR -> GENOMIC REGION SET (containing all the mpbs of FACTOR for GENOMIC REGION)
    mpbs_output_dict = dict()

    # Creating genome file
    genome_file = Fastafile(genome_data.get_genome())

    # Iterating on list of genomic regions
    for genomic_region_set in input_regions:

        # Initializing factor dictionary
        mpbs_output_dict[genomic_region_set.name] = dict()
    
        # Iterating on genomic regions
        for genomic_region in genomic_region_set.sequences:

            # Reading sequence associated to genomic_region
            sequence = str(genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

            # Iterating on motif group list
            for motif_group in motif_list:

                # Creating dataset for multiprocessing
                curr_data_input = [[m,sequence,genomic_region] for m in motif_group]
                curr_proc_nb = len(curr_data_input)

                # Applying motif matching function with multiprocessing
                pool = Pool(curr_proc_nb)
                curr_mpbs_list = pool.map(match,curr_data_input)
                pool.close()
                pool.join()
                for i in range(len(curr_mpbs_list)):
                    gr_list = curr_mpbs_list[i]
                    curr_motif_name = motif_group[i].name
                    for gr in gr_list:
                        try: mpbs_output_dict[genomic_region_set.name][curr_motif_name].add(gr)
                        except Exception:
                            mpbs_output_dict[genomic_region_set.name][curr_motif_name] = GenomicRegionSet(curr_motif_name)
                            mpbs_output_dict[genomic_region_set.name][curr_motif_name].add(gr)
        
    ###################################################################################################
    # Writing output
    ###################################################################################################

    # Dumping list of GenomicRegionSet for fast access by --enrichment operation
    output_file_name = os.path.join(matching_output_location, dump_file_name)
    output_file = open(output_file_name, "wb")
    dump(mpbs_output_dict, output_file)
    output_file.close()

    # Iterating on MPBS output list
    for k in mpbs_output_dict.keys():

        # Initializations
        mpbs_dict = mpbs_output_dict[k]
        output_file_name = os.path.join(matching_output_location, k+"_mpbs")

        # Writing bed file
        bed_file_name = output_file_name+".bed"
        bed_file = open(bed_file_name,"w")
        for kk in mpbs_dict.keys():
            for e in mpbs_dict[kk]:
                bed_file.write("\t".join([e.chrom,str(e.initial),str(e.final),e.name,str(e.data),e.orientation])+"\n")
        bed_file.close()

        # Verifying condition to write bb
        if(options.bigbed):

            # Fetching file with chromosome sizes
            chrom_sizes_file = genome_data.get_chromosome_sizes()

            # Converting to big bed
            bb_file_name = output_file_name+".bb"
            try:
                os.system(" ".join(["bedToBigBed", bed_file_name, chrom_sizes_file, bb_file_name, "-verbose=0"]))
                os.remove(bed_file_name)
            except Exception: pass # WARNING


def main_enrichment():
    """
    Performs motif enrichment.

    Authors: Eduardo G. Gusmao.
    """

    ###################################################################################################
    # Processing Input Arguments
    ###################################################################################################

    # Parameters
    usage_message = "%prog --enrichment [options] <experiment_matrix> <input_path>"

    # Initializing Option Parser
    parser = PassThroughOptionParser(usage = usage_message)

    # Optional Input Options
    #parser.add_option("--xxxxxx", dest = "xxxxxx", type = "xxxxxx", metavar="xxxxxx", default = None,
    #                  help = ("xxxxxx "))

    # Parameters Options
    parser.add_option("--organism", dest = "organism", type = "string", metavar="STRING", default = "hg19",
                      help = ("Organism considered on the analysis. Check our full documentation for all available "
                              "options. All default files such as genomes will be based on the chosen organism "
                              "and the data.config file."))
    parser.add_option("--promoter-length", dest = "promoter_length", type = "int", metavar="INT", default = 1000,
                      help = ("Length of the promoter region (in bp) considered on the creation of the regions-gene association."))
    parser.add_option("--maximum-association-length", dest = "maximum_association_length", type = "int", metavar="INT", default = 50000,
                      help = ("Maximum distance between a coordinate and a gene (in bp) in order for the former to "
                              "be considered associated with the latter."))
    parser.add_option("--multiple-test-alpha", dest = "multiple_test_alpha", type = "float", metavar="FLOAT", default = 0.05,
                      help = ("Alpha value for multiple test."))

    # Output Options
    parser.add_option("--output-location", dest = "output_location", type = "string", metavar="PATH", default = None,
                      help = ("Path where the output files will be written. Default is the input PATH."))
    parser.add_option("--print-thresh", dest = "print_thresh", type = "float", metavar="FLOAT", default = 0.05,
                      help = ("Only MPBSs whose factor's enrichment corrected p-value are less than equal "
                              "this option are print. Use 1.0 to print all MPBSs."))
    parser.add_option("--bigbed", dest = "bigbed", action = "store_true", default = False,
                      help = ("If this option is used, all bed files will be written as bigbed."))

    # Processing Options
    options, arguments = parser.parse_args()

    # Additional Parameters
    matching_folder_name = "Match"
    random_region_name = "random_regions"
    dump_file_name = "dump.p"
    gene_column_name = "genegroup"
    output_association_name = "coord_association"
    output_mpbs_filtered = "mpbs.bed"
    output_mpbs_filtered_ev = "mpbs_ev.bed"
    output_mpbs_filtered_ev = "mpbs_nev.bed"
    output_stat_genetest = "genetest_statistics"
    output_stat_randtest = "randtest_statistics"

    ###################################################################################################
    # Initializations
    ###################################################################################################

    # Output folder
    if(options.output_location): output_location = options.output_location
    else:
        try:
            matrix_name_without_ext = ".".join(arguments[0].split(".")[:-1])
            output_location = os.path.join(arguments[1],os.path.basename(matrix_name_without_ext))
        except Exception: pass # TODO ERROR 

    # Default genomic data
    genome_data = GenomeData(options.organism)

    # Default motif data
    motif_data = MotifData()

    ###################################################################################################
    # Reading Input Matrix
    ###################################################################################################

    # Reading arguments
    try:
        input_matrix = arguments[0]
        input_location = arguments[1]
        if(len(arguments) > 2): pass # TODO WARNING
    except Exception: pass # TODO ERROR

    # Create experimental matrix
    try:
        exp_matrix = ExperimentalMatrix()
        exp_matrix.read(input_matrix)
    except Exception: pass # TODO ERROR

    ###################################################################################################
    # Reading Regions & Gene Lists
    ###################################################################################################

    # Input Class
    class Input:
        def __init__(self,gene_set,region_list):
            self.gene_set = gene_set
            self.region_list = region_list

    # Initializations
    input_list = []

    # Reading dictionary grouped by fields
    flag_gene = True
    try:
        exp_matrix_fields_dict = exp_matrix.fieldsDict[gene_column_name]
    except Exception: flag_gene = False

    # Reading dictionary of objects
    try:
        exp_matrix_objects_dict = exp_matrix.objectsDict
    except Exception: pass # TODO ERROR

    if(flag_gene): # Genelist and randomic analysis will be performed

        # Iterating on experimental matrix fields
        for g in exp_matrix_fields_dict.keys():

            # Create input which will contain all regions associated with such gene group
            curr_input = Input(None, [])

            # This flag will be used to see if there are two gene files associated with the same gene label on genegroup column
            flag_foundgeneset = False

            # Iterating on experimental matrix objects
            for k in exp_matrix_fields_dict[g]:

                curr_object = exp_matrix_objects_dict[k]

                # If the current object is a GenomicRegionSet
                if(isinstance(curr_object,GenomicRegionSet)):

                    # Sorting input region
                    curr_object.sort()

                    # Updating Input object
                    curr_input.region_list.append(curr_object)

                # If the current object is a GeneSet
                if(isinstance(curr_object,GeneSet)):

                    # Updating Input object
                    curr_object.name = g # The name in gene_group column will be used. The 'name' column for genes are not used.
                    if(not flag_foundgeneset):
                        curr_input.gene_set = curr_object
                        flag_foundgeneset = True
                    else: pass # TODO WARNING - Only first geneset will be used
                    
            if(not flag_foundgeneset): pass # TODO ERROR - There must be a valid geneset associated with each region

            # Update input list
            input_list.append(curr_input)

    else: # Only randomic analysis will be performed

        # Create single input which will contain all regions
        single_input = Input(None, [])

        # Iterating on experimental matrix objects
        for k in exp_matrix_objects_dict.keys():

            curr_object = exp_matrix_objects_dict[k]

            # If the current object is a GenomicRegionSet
            if(isinstance(curr_genomic_region,GenomicRegionSet)):

                # Sorting input region
                curr_object.sort()

                # Updating Input object
                single_input.region_list.append(curr_object)

        # Updating input list with single input (only randomic analysis will be performed)
        input_list = [single_input]

    ###################################################################################################
    # Reading Motif Matching
    ###################################################################################################

    # Verifying if file exists
    curr_dump_file_name = os.path.join(input_location,matching_folder_name,dump_file_name)
    if(not os.path.exists(curr_dump_file_name)): pass # TODO ERROR
    
    # Opening dump
    dump_file = open(curr_dump_file_name, "rb")
    mpbs_dict = load(dump_file)
    dump_file.close()

    ###################################################################################################
    # Reading Random Coordinates
    ###################################################################################################

    # Verifying if file exists
    random_region_glob = glob(os.path.join(input_location,matching_folder_name,random_region_name+".*"))
    try: random_region_file_name = random_region_glob[0]
    except Exception: pass # TODO ERROR

    # Creating GenomicRegionSet
    random_regions = GenomicRegionSet(random_region_name)

    # Reading bed file
    if(random_region_file_name.split(".")[-1] == "bed"):
        random_regions.read_bed(random_region_file_name)

    # Reading bigbed file
    elif(random_region_file_name.split(".")[-1] == "bb"):
        try:
            bed_file_name = ".".join(random_region_file_name.split(".")[:-1])
            os.system(" ".join(["bigBedToBed", random_region_file_name, bed_file_name]))
        except Exception: pass # TODO ERROR
        random_regions.read_bed(bed_file_name)
        os.remove(bed_file_name)

    else: pass # TODO ERROR

    # Evaluating random statistics
    #rand_c_dict, rand_d_dict = random_regions.fishertable(mpbs_dict[random_regions.name])

    ###################################################################################################
    # Enrichment Statistics
    ###################################################################################################

    """

    # TODO - APAGAR
    output_association_name = "coord_association"
    output_mpbs_filtered = "mpbs.bed"
    output_mpbs_filtered_ev = "mpbs_ev.bed"
    output_mpbs_filtered_ev = "mpbs_nev.bed"
    output_stat_genetest = "genetest_statistics"
    output_stat_randtest = "randtest_statistics"

    # Iterating on each input object
    for curr_input in input_list:

        # Iterating on each input genomic region set
        for grs in curr_input.region_list:

            # Initialization
            original_name = grs.name

            # Creating output folder
            if(curr_input.gene_set): curr_output_folder_name = os.path.join(output_location,grs.name+"__"+curr_input.gene_set.name)
            else: curr_output_folder_name = os.path.join(output_location,grs.name)
            if(not os.path.isdir(curr_output_folder_name)): os.makedirs(curr_output_folder_name)

            ###################################################################################################
            # Gene Evidence Statistics
            ################################################################################################### 

            if(curr_input.gene_set):

                # Performing association of input region with gene_set
                grs = grs.gene_association(curr_input.gene_set, options.organism, options.promoter_length, options.maximum_association_length)

                # Writing gene-coordinate association
                output_file_name = os.path.join(curr_output_folder_name,output_association_name+".bed")
                output_file = open(output_file_name,"w")
                for gr in grs:
                    curr_gene_list = [e if e[0]!="." else e[1:] for gr.name.split(":")]
                    curr_prox_list = gr.data.split(":")
                    curr_name = ":".join([e[0]+"_"+e[1] for e in zip(curr_gene_list,curr_prox_list)])
                    output_file.write("\t".join([gr.chrom,gr.initial,gr.final,curr_name])+"\n")
                output_file.close()
                if(options.bigbed):
                    chrom_sizes_file = genome_data.get_chromosome_sizes()
                    bb_file_name = output_file_name+".bb"
                    try:
                        os.system(" ".join(["bedToBigBed", bed_file_name, chrom_sizes_file, bb_file_name, "-verbose=0"]))
                        os.remove(bed_file_name)
                    except Exception: pass # WARNING

                # Creating ev and nev sets
                curr_ev = GenomicRegionSet(grs.name+"_ev")
                curr_nev = GenomicRegionSet(grs.name+"_nev")

                # Populating ev and nev sets
                for gr in grs:
                    if(len([e for e in gr.name.split(":") if e[0]!="."]) > 0): curr_ev.add(gr)
                    else: curr_nev.add(gr)

                # Calculating statistics
                fisher_a_dict, fisher_b_dict = curr_ev.fishertable(mpbs_dict[original_name])
                fisher_c_dict, fisher_d_dict = curr_ev.fishertable(mpbs_dict[original_name])
                
                # Performing fisher test
                

                # Performing multiple test correction
                

                # Fetching ev and nev mpbs to write
                curr_ev_mpbs = # TODO
                curr_nev_mpbs = # TODO

                # Printing ev and nev mpbs
                

                # Printing statistics


            ###################################################################################################
            # Random Statistics
            ###################################################################################################

            # Evaluating random statistics

            # Iterating on grs



            ###################################################################################################
            # Writting
            ################################################################################################### 





    

    ###################################################################################################
    # Heatmap
    ###################################################################################################

    """

    """
    test_region = input_list[0].region_list[0]
    test_geneset = input_list[0].gene_set
    test_result = test_region.gene_association(gene_set=test_geneset)
    print "--------"
    for gr in test_result:
        print gr.chrom, gr.initial, gr.final
        print gr.name
        print gr.data
        print "--------"

    print "\n"
    all_genes, mapped_genes, all_proxs, mapped_proxs = test_region.filter_by_gene_association(test_geneset)
    print test_region.name
    print "--------"
    for gr in test_region:
        print gr.chrom, gr.initial, gr.final
        print gr.name
        print gr.data
        print "--------"
    for e in all_genes.genes: print e
    print "--------"
    for e in mapped_genes.genes: print e
    print "--------"
    for e in all_proxs: print e
    print "--------"
    for e in mapped_proxs: print e
    print "--------"
    """































