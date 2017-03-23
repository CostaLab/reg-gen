
###################################################################################################
##### Annotation File Standard (tab-separated):
# MATRIX_ID: The matrices' ID. It may change format for different repositories. (STRING)
# PWM_NAME: Name of the PWM inside the respective repository. (STRING - without the .pwm)
# SOURCE: The source repository of such matrix. (STRING)
# VERSION: The version of such matrix (1 for 'primary motif', 2 for 'secondary motif', etc). (INT)
# GENE_NAMES: Name of genes associated with such TF matrix. (LIST)
# GROUP: <To be used in future - Will represent clusters of motifs>. (STRING)
###################################################################################################
# * Mandatory fields: MATRIX_ID, SOURCE, VERSION, GENE_NAMES.
# * Fields with multiple entries should be separated by ';' (no spaces).
# * Fields with missing/non-existing/doesn't matter data should be filled with '.'
# * Co-binding should be represented by '+' (no spaces).
###################################################################################################

# Import
import os
import sys
import glob

# Parameters
dataLocation = "./"
group = "."

###################################################################################################
# HOCOMOCO
###################################################################################################

# Fetching file names
# TODO: check if this still works for hocomoco v10
source = "hocomoco"
inputLocation = dataLocation+source+"/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation+"*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[1]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    gene_names = ll[0]
    resultMatrix.append([matrix_id,pwm_name,source,version,gene_names,group])

# Sorting results by ID
resultMatrix = sorted(resultMatrix ,key=lambda x: x[2])
resultMatrix = sorted(resultMatrix ,key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation+source+".mtf"
outputFile = open(outputFileName,"w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec)+"\n")
outputFile.close()

###################################################################################################
# JASPAR VERTEBRATES
###################################################################################################

# Fetching file names
source = "jaspar_vertebrates"
inputLocation = dataLocation+source+"/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation+"*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    if(len(ll) > 4): version = ll[4]
    gene_names = ll[2].replace("::","+")
    resultMatrix.append([matrix_id,pwm_name,source,version,gene_names,group])

# Sorting results by ID
resultMatrix = sorted(resultMatrix ,key=lambda x: x[2])
resultMatrix = sorted(resultMatrix ,key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation+source+".mtf"
outputFile = open(outputFileName,"w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec)+"\n")
outputFile.close()

###################################################################################################
# TRANSFAC PUBLIC
###################################################################################################

# Fetching file names
source = "transfac_public"
inputLocation = dataLocation+source+"/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation+"*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = "1"
    gene_names = ll[1]
    resultMatrix.append([matrix_id,pwm_name,source,version,gene_names,group])

# Sorting results by ID
resultMatrix = sorted(resultMatrix ,key=lambda x: x[2])
resultMatrix = sorted(resultMatrix ,key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation+source+".mtf"
outputFile = open(outputFileName,"w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec)+"\n")
outputFile.close()

###################################################################################################
# UNIPROBE PRIMARY
###################################################################################################

# Fetching file names
source = "uniprobe_primary"
inputLocation = dataLocation+source+"/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation+"*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = ll[1]
    gene_names = ll[2]
    resultMatrix.append([matrix_id,pwm_name,source,version,gene_names,group])

# Sorting results by ID
resultMatrix = sorted(resultMatrix ,key=lambda x: x[2])
resultMatrix = sorted(resultMatrix ,key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation+source+".mtf"
outputFile = open(outputFileName,"w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec)+"\n")
outputFile.close()

###################################################################################################
# UNIPROBE SECONDARY
###################################################################################################

# Fetching file names
source = "uniprobe_secondary"
inputLocation = dataLocation+source+"/"
resultMatrix = []
for inputFileName in glob.glob(inputLocation+"*.pwm"):
    ll = inputFileName.split("/")[-1].split(".")[0].split("_")
    matrix_id = ll[0]
    pwm_name = ".".join(inputFileName.split("/")[-1].split(".")[:-1])
    version = ll[1]
    gene_names = ll[2]
    resultMatrix.append([matrix_id,pwm_name,source,version,gene_names,group])

# Sorting results by ID and version
resultMatrix = sorted(resultMatrix ,key=lambda x: x[2])
resultMatrix = sorted(resultMatrix ,key=lambda x: x[0])

# Writing to output file
outputFileName = dataLocation+source+".mtf"
outputFile = open(outputFileName,"w")
for resultVec in resultMatrix:
    outputFile.write("\t".join(resultVec)+"\n")
outputFile.close()


