
# Import
import os
import sys
from glob import glob
from Bio import motifs
from os.path import basename

# Input
inFolder = sys.argv[1]
outFileName = sys.argv[2]

# Parameters
fprList = [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
pseudocounts = 0.1
background = {'A':0.25,'C':0.25,'G':0.25,'T':0.25}
precision = 10000

# Creating output file
outFile = open(outFileName,"w")
outFile.write("\t".join(["MOTIF"]+[str(e) for e in fprList])+"\n")

# Iterating on all PWMs
for pwmFileName in glob(inFolder+"*.pwm"):

  # Creating PSSM
  name = ".".join(basename(pwmFileName).split(".")[:-1])
  input_file = open(pwmFileName,"r")
  pfm = motifs.read(input_file, "pfm")
  pwm = pfm.counts.normalize(pseudocounts)
  input_file.close()
  pssm = pwm.log_odds(background)
  pssm_list = [pssm[e] for e in ["A","C","G","T"]]
  distribution = pssm.distribution(background=background, precision=precision)

  # Evaluating thresholds
  resVec = [name]
  for fpr in fprList:
    resVec.append(str(distribution.threshold_fpr(fpr)))
    
  # Writing results
  outFile.write("\t".join(resVec)+"\n")


