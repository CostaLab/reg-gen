# Import
import sys
from glob import glob
from os.path import basename

from MOODS import tools, parsers

# Input
inFolder = sys.argv[1]
outFileName = sys.argv[2]

# Parameters
fprList = [0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001]
pseudocounts = 1.0

outFile = open(outFileName, "w")
outFile.write("\t".join(["MOTIF"] + [str(e) for e in fprList]) + "\n")

# Iterating on all PWMs
for pwmFileName in sorted(glob(inFolder + "*.pwm")):
    # Creating PSSM
    name = ".".join(basename(pwmFileName).split(".")[:-1])

    pfm = parsers.pfm(pwmFileName)
    bg = tools.flat_bg(len(pfm))  # total number of "points" to add, not per-row
    pssm = tools.log_odds(pfm, bg, pseudocounts, 2)

    # Evaluating thresholds
    resVec = [name]
    for fpr in fprList:
        # Note: this requires a modified version of MOODS. Only use it if you know what you are doing
        # resVec.append(str(tools.threshold_from_p(pssm, bg, fpr, 10000.0)))
        resVec.append(str(tools.threshold_from_p(pssm, bg, fpr)))

    # Writing results
    outFile.write("\t".join(resVec) + "\n")

outFile.close()
