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

# Parameters
dataLocation = "./"
group = "."

###################################################################################################
# HOCOMOCO
###################################################################################################

# Fetching file names
source = "hocomoco"
inputLocation = dataLocation + source + "/"
resultMatrix = []
hocomoco_anno = {}
with open("hocomoco_anno.csv", "rb") as f:
    csvf = csv.reader(f)
    for l in csvf:
        hocomoco_anno[l[0]] = l[1:]
for inputFileName in glob.glob(inputLocation + "*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    gene_names = hocomoco_anno[pwm_name][0]
    group = hocomoco_anno[pwm_name][1]
    if not group:
        group = "."
    uniprot = hocomoco_anno[pwm_name][2]
    data_source = hocomoco_anno[pwm_name][3]
    resultMatrix.append([matrix_id, pwm_name, source, version, gene_names, group, uniprot, data_source])

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

# Fetching file names
source = "jaspar_vertebrates"
inputLocation = dataLocation + source + "/"
resultMatrix = []
jaspar_anno = {}
with open("jaspar_anno.csv", "rb") as f:
    csvf = csv.reader(f)
    for l in csvf:
        if not l:
            continue
        jaspar_anno[l[0]] = l[1:]
for inputFileName in glob.glob(inputLocation + "*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    if len(ll) > 4:
        version = ll[4]
    gene_names = jaspar_anno[pwm_name][0]
    group = jaspar_anno[pwm_name][1]
    if not group:
        group = "."
    uniprot = jaspar_anno[pwm_name][2]
    data_source = jaspar_anno[pwm_name][3]
    resultMatrix.append([matrix_id, pwm_name, source, version, gene_names, group, uniprot, data_source])

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

# Fetching file names
source = "transfac_public"
inputLocation = dataLocation + source + "/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation + "*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    gene_names = ll[1]
    resultMatrix.append([matrix_id, pwm_name, source, version, gene_names, ".", ".", ".", "."])

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

# Fetching file names
source = "uniprobe_primary"
inputLocation = dataLocation + source + "/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation + "*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = ll[1]
    gene_names = ll[2]
    resultMatrix.append([matrix_id, pwm_name, source, version, gene_names, ".", ".", ".", "."])

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

# Fetching file names
source = "uniprobe_secondary"
inputLocation = dataLocation + source + "/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation + "*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = ll[1]
    gene_names = ll[2]
    resultMatrix.append([matrix_id, pwm_name, source, version, gene_names, ".", ".", ".", "."])

# Sorting results by ID and version
resultMatrix = sorted(resultMatrix, key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation + source + ".mtf"
outputFile = open(outputFileName, "w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec) + "\n")
outputFile.close()
