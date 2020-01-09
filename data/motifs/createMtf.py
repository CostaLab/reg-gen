###################################################################################################
##### Annotation File Standard (tab-separated):
# MATRIX_ID: The matrices' ID. It may change format for different repositories. (STRING)
# PWM_NAME: Name of the PWM inside the respective repository. (STRING - without the .pwm)
# SOURCE: The source repository of such matrix. (STRING)
# VERSION: The version of such matrix (1 for 'primary motif', 2 for 'secondary motif', etc). (INT)
# GENE_NAMES: Name of genes associated with such TF matrix. (LIST)
# GROUP: Name of factor "family" or "class" or "cluster", depending on repository. (STRING)
# UniProt: UniProt accession for the transcription factor. (STRING)
###################################################################################################
# * Mandatory fields: MATRIX_ID, SOURCE, VERSION, GENE_NAMES.
# * Fields with multiple entries should be separated by ';' (no spaces).
# * Fields with missing/non-existing/doesn't matter data should be filled with '.'
# * Co-binding should be represented by '+' (no spaces).
# * Group can be any string, and it may also contain whitespaces and punctuation (no tab).
###################################################################################################

# Import
import glob
import csv
import re
from optparse import OptionParser
from glob import glob
from os.path import basename

from MOODS import tools, parsers
import argparse

# Parameters
dataLocation = "./"
group = "."
fprList = [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
pseudocounts = 1.0

# Option Parser
usage_message = ""
parser = OptionParser(usage=usage_message)

parser.add_option("--all", dest="all", action="store_true", default=False, help="Create all mtf files")
parser.add_option("--hoc", dest="hoc", action="store_true", default=False, help="Create hocomoco.mtf")
parser.add_option("--jv", dest="jv", action="store_true", default=False, help="Create jaspar_vertebrates.mtf")
parser.add_option("--jp", dest="jp", action="store_true", default=False, help="Create jaspar_plants.mtf")
parser.add_option("--t", dest="t", action="store_true", default=False, help="Create transfac.mtf")
parser.add_option("--up", dest="up", action="store_true", default=False, help="Create uniprobe_primary.mtf")
parser.add_option("--us", dest="us", action="store_true", default=False, help="Create uniprobe_secondary.mtf")

options, arguments = parser.parse_args()

if options.all:
    options.hoc = True
    options.jv = True
    options.jp = True
    options.t = True
    options.up = True
    options.us = True

###################################################################################################
# HOCOMOCO
###################################################################################################
if options.hoc:
    # Fetching file names
    source = "hocomoco"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    hocomoco_anno = {}
    with open("hocomoco_anno.csv", "r") as f:
        csvf = csv.reader(f, dialect="unix")
        for l in csvf:
            hocomoco_anno[l[0]] = l[1:]
    for inputFileName in glob(inputLocation + "*.pwm"):
        ll = inputFileName.split("/")[-1].split(".")[0].split("_")
        matrix_id = ll[0]
        pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
        version = ".".join(pwm_name.split(".")[2:])
        gene_names = hocomoco_anno[pwm_name][0]
        group = hocomoco_anno[pwm_name][1]
        if not group:
            group = "."
        uniprot = hocomoco_anno[pwm_name][2]
        data_source = hocomoco_anno[pwm_name][3]
        taxGroup = "vertebrates"
        species = (pwm_name.split("_")[1]).split(".")[0]
        if species == "HUMAN":
            species = "Homo sapiens"
        elif species == "MOUSE":
            species = "Mus musculus"

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, pwm_name, version, gene_names, group, uniprot, data_source, taxGroup,
                             species, threshold])

    # Sorting results by ID
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()

###################################################################################################
# JASPAR VERTEBRATES
###################################################################################################
if options.jv:
    # Fetching file names
    source = "jaspar_vertebrates"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    jaspar_anno = {}
    with open("jaspar_anno.csv", "r") as f:
        csvf = csv.reader(f, dialect="unix")
        for l in csvf:
            if not l:
                continue
            jaspar_anno[l[0]] = l[1:]
    for inputFileName in glob(inputLocation + "*.pwm"):
        raw_name = inputFileName.split("/")[-1]
        full_name, pwm_name, _, _ = re.match("((.+?)(\(var\.(\d+)\))?).pwm", raw_name).groups()
        name_fields = pwm_name.split(".")
        matrix_id = name_fields[0]
        if "var" in full_name.split(".")[2]:
            version = ".".join([name_fields[1], (full_name.split(".")[3]).strip(")")])
        else:
            version = name_fields[1]
        gene_names = name_fields[2]
        group = jaspar_anno[full_name][1]
        if not group:
            group = "."
        uniprot = jaspar_anno[full_name][2]
        data_source = jaspar_anno[full_name][3]
        taxGroup = jaspar_anno[full_name][4]
        species = jaspar_anno[full_name][5]

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, full_name, version, gene_names, group, uniprot, data_source, taxGroup,
                             species, threshold])
    # Sorting results by ID
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()

###################################################################################################
# JASPAR PLANTS
###################################################################################################
if options.jp:
    # Fetching file names
    source = "jaspar_plants"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    jaspar_anno = {}
    with open("jaspar_anno.csv", "r") as f:
        csvf = csv.reader(f, dialect="unix")
        for l in csvf:
            if not l:
                continue
            jaspar_anno[l[0]] = l[1:]
    for inputFileName in glob(inputLocation + "*.pwm"):
        raw_name = inputFileName.split("/")[-1]
        full_name, pwm_name, _, _ = re.match("((.+?)(\(var\.(\d+)\))?).pwm", raw_name).groups()
        name_fields = pwm_name.split(".")
        matrix_id = name_fields[0]
        version = name_fields[1]
        gene_names = name_fields[2]
        group = jaspar_anno[full_name][1]
        if not group:
            group = "."
        uniprot = jaspar_anno[full_name][2]
        data_source = jaspar_anno[full_name][3]
        taxGroup = jaspar_anno[full_name][4]
        species = jaspar_anno[full_name][5]

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, pwm_name, version, gene_names, group, uniprot, data_source, taxGroup,
                             species, threshold])
    # Sorting results by ID
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()

###################################################################################################
# TRANSFAC PUBLIC
###################################################################################################
if options.t:
    # Fetching file names
    source = "transfac_public"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    for inputFileName in glob(inputLocation + "*.pwm"):
        ll = inputFileName.split("/")[-1].split(".")[0].split("_")
        matrix_id = ll[0]
        pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
        version = "1"
        gene_names = ll[1]

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, pwm_name, version, gene_names, ".", ".", ".", ".", ".", threshold])

    # Sorting results by ID
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()

###################################################################################################
# UNIPROBE PRIMARY
###################################################################################################
if options.up:
    # Fetching file names
    source = "uniprobe_primary"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    for inputFileName in glob(inputLocation + "*.pwm"):
        ll = inputFileName.split("/")[-1].split(".")[0].split("_")
        matrix_id = ll[0]
        pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
        version = ll[1]
        gene_names = ll[2]

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, pwm_name, version, gene_names, ".", ".", ".", ".", ".", threshold])

    # Sorting results by ID
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()

###################################################################################################
# UNIPROBE SECONDARY
###################################################################################################
if options.us:
    # Fetching file names
    source = "uniprobe_secondary"
    inputLocation = dataLocation + source + "/"
    resultMatrix = []
    for inputFileName in glob(inputLocation + "*.pwm"):
        ll = inputFileName.split("/")[-1].split(".")[0].split("_")
        matrix_id = ll[0]
        pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
        version = ll[1]
        gene_names = ll[2]

        # Creating PSSM
        pfm = parsers.pfm(inputFileName)
        bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
        pssm = tools.log_odds(pfm, bg, pseudocounts, 2)
        threshold_list = []

        # Evaluating thresholds
        for fpr in fprList:
            # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
            # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
            threshold = tools.threshold_from_p(pssm, bg, fpr)
            threshold_list.append(str(threshold))
        threshold = ",".join(threshold_list)

        resultMatrix.append([matrix_id, pwm_name, version, gene_names, ".", ".", ".", ".", ".", threshold])

    # Sorting results by ID and version
    resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

    # Writing to output file
    outputFileName = dataLocation + source + ".mtf"
    outputFile = open(outputFileName, "w")
    for resultVec in resultMatrix:
        outputFile.write("\t".join(resultVec) + "\n")
    outputFile.close()
