
# Import
import sys
from glob import glob
from Bio import motifs
from os.path import basename
import MOODS.tools
import MOODS.parsers

# Input
inFolder = sys.argv[1]
outFileName = sys.argv[2]

# Parameters
fprList = [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
pseudocounts = 1.0

outFile = open(outFileName, "w")
outFile.write("\t".join(["MOTIF"]+[str(e) for e in fprList])+"\n")

# Iterating on all PWMs
for pwmFileName in sorted(glob(inFolder+"*.pwm")):
    # Creating PSSM
    name = ".".join(basename(pwmFileName).split(".")[:-1])

    pfm = MOODS.parsers.pfm(pwmFileName)
    bg = MOODS.tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
    pssm = MOODS.tools.log_odds(pfm, bg, pseudocounts)

    # Evaluating thresholds
    resVec = [name]
    for fpr in fprList:
        resVec.append(str(MOODS.tools.threshold_from_p(pssm, bg, fpr)))

    # Writing results
    outFile.write("\t".join(resVec)+"\n")

outFile.close()
