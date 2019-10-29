#!/usr/bin/env python
# -*- coding: utf-8 -*-
#################################################################################################
params = []
params.append("")
params.append("Usage: protectionScore.py <FOOTPRINT_FILE> <MPBS_FILE> <DNASE_FILE> <GENOME_FILE> <OUTPUT_FILE>")
params.append("\nInput: ")
params.append("  FOOTPRINT_FILE: Footprint file (bed or bb). Output of rgt-hint program.")
params.append("  MPBS_FILE: Motif-predicted binding sites file name (bed or bb).")
params.append("             Output of rgt-motifanalysis --matching program.")
params.append("  DNASE_FILE: DNase-seq aligned reads file (bam).")
params.append("  GENOME_FILE: Genome fasta file.")
params.append("  OUTPUT_FILE: Name of output file where protection scores will be written.")
params.append("\nDescription: ")
params.append("  Evaluates the protection score for every MPBS in MPBS_FILE.")
params.append("  The protection score is evaluated only for MPBSs that")
params.append("  overlap with footprints. Also, the protection score is based")
params.append("  on the bias-corrected version of the DNase-seq signal in DNASE_FILE.")
params.append("  The correction will be performed based on k-mer frequencies estimated")
params.append("  on K562 (DNase-seq single-hit protocol). The output file consists of a plain")
params.append("  tab-separated text file in which each row has the MPBS name and its")
params.append("  protection score, respectively.")
params.append("")
#################################################################################################

# Import
import os
import sys
from math import log
from pysam import __version__ as ps_version
from pysam import Samfile
from pysam import Fastafile
from numpy import array
from rgt.HINT.pileupRegion import PileupRegion
if(len(sys.argv) != 6): 
    for e in params: print(e)
    sys.exit(0)

#################################################################################################
# INPUT
#################################################################################################

# Reading input
fpFileName = sys.argv[1]
mpbsFileName = sys.argv[2]
dnaseFileName = sys.argv[3]
genomeFileName = sys.argv[4]
outputFileName = sys.argv[5]

# Parameters
initial_clip = 1000
minValue = 100
halfWindow = 50
motifExt = 20
to_remove = []
scriptPath = os.path.dirname(os.path.realpath(__file__))
biasTableFFileName = os.path.join(scriptPath,"../data/fp_hmms/single_hit_bias_table_F.txt")
biasTableRFileName = os.path.join(scriptPath,"../data/fp_hmms/single_hit_bias_table_R.txt")

#################################################################################################
# BIAS CORRECTION FUNCTION
#################################################################################################

def bias_correction(bam, signal, fBiasDict, rBiasDict, genome_file_name, chrName, start, end):

  # Parameters
  window = 50
  defaultKmerValue = 1.0

  # Initialization
  fastaFile = Fastafile(genome_file_name)
  k_nb = len(list(fBiasDict.keys())[0])
  p1 = start; p2 = end
  p1_w = p1 - (window/2); p2_w = p2 + (window/2)
  p1_wk = p1_w - (k_nb/2); p2_wk = p2_w + (k_nb/2)

  # Raw counts
  nf = [0.0] * (p2_w-p1_w); nr = [0.0] * (p2_w-p1_w)
  for r in bam.fetch(chrName, p1_w, p2_w):
    if((not r.is_reverse) and (r.pos > p1_w)): nf[r.pos-p1_w] += 1.0
    if((r.is_reverse) and ((r.aend-1) < p2_w)): nr[r.aend-1-p1_w] += 1.0

  # Smoothed counts
  Nf = []; Nr = [];
  fSum = sum(nf[:window]); rSum = sum(nr[:window]);
  fLast = nf[0]; rLast = nr[0]
  for i in range((window/2),len(nf)-(window/2)):
    Nf.append(fSum)
    Nr.append(rSum)
    fSum -= fLast; fSum += nf[i+(window/2)]; fLast = nf[i-(window/2)+1]
    rSum -= rLast; rSum += nr[i+(window/2)]; rLast = nr[i-(window/2)+1]

  # Fetching sequence
  currStr = str(fastaFile.fetch(chrName, p1_wk-1, p2_wk-2)).upper()
  currRevComp = revcomp(str(fastaFile.fetch(chrName,p1_wk+2, p2_wk+1)).upper())

  # Iterating on sequence to create signal
  af = []; ar = []
  for i in range((k_nb/2),len(currStr)-(k_nb/2)+1):
    fseq = currStr[i-(k_nb/2):i+(k_nb/2)]
    rseq = currRevComp[len(currStr)-(k_nb/2)-i:len(currStr)+(k_nb/2)-i]
    try: af.append(fBiasDict[fseq])
    except Exception: af.append(defaultKmerValue)
    try: ar.append(rBiasDict[rseq])
    except Exception: ar.append(defaultKmerValue)

  # Calculating bias and writing to wig file
  fSum = sum(af[:window]); rSum = sum(ar[:window]);
  fLast = af[0]; rLast = ar[0]
  bias_corrected_signal = []
  for i in range((window/2),len(af)-(window/2)):
    nhatf = Nf[i-(window/2)]*(af[i]/fSum)
    nhatr = Nr[i-(window/2)]*(ar[i]/rSum)
    zf = log(nf[i]+1)-log(nhatf+1)
    zr = log(nr[i]+1)-log(nhatr+1)
    bias_corrected_signal.append(zf+zr)
    fSum -= fLast; fSum += af[i+(window/2)]; fLast = af[i-(window/2)+1]
    rSum -= rLast; rSum += ar[i+(window/2)]; rLast = ar[i-(window/2)+1]

  # Termination
  fastaFile.close()
  return bias_corrected_signal

def revcomp(s):
  revDict = dict([("A","T"),("T","A"),("C","G"),("G","C"),("N","N")])
  return "".join([revDict[e] for e in s[::-1]])

#################################################################################################
# COVERAGE OF BAM
#################################################################################################

covFileName = outputFileName+"_cov.txt"
to_remove.append(covFileName)
os.system("samtools view -c "+dnaseFileName+" > "+covFileName)
covFile = open(covFileName,"r")
rmp = 1000000.0 /float(covFile.readline().strip())
covFile.close()

#################################################################################################
# CONVERTING TO BED
#################################################################################################

# Footprint file
if(fpFileName.split(".")[-1] == "bb"):
  fpFileNameBed = ".".join(fpFileName.split(".")[:-1])+".bed"
  os.system("bigBedToBed "+fpFileName+" "+fpFileNameBed)
  to_remove.append(fpFileNameBed)
else: fpFileNameBed = fpFileName

# MPBS file
if(mpbsFileName.split(".")[-1] == "bb"):
  mpbsFileNameBed = ".".join(mpbsFileName.split(".")[:-1])+".bed"
  os.system("bigBedToBed "+mpbsFileName+" "+mpbsFileNameBed)
  to_remove.append(mpbsFileNameBed)
else: mpbsFileNameBed = mpbsFileName

#################################################################################################
# READING MPBS NAMES
#################################################################################################

mpbsList = []
tempFileName = outputFileName+"_tempmpbsname.txt"
to_remove.append(tempFileName)
os.system("cut -f 4 "+mpbsFileNameBed+" | sort | uniq > "+tempFileName)
tempFile = open(tempFileName,"r")
for line in tempFile:
  mpbsList.append(line.strip())
tempFile.close()
mpbsList = sorted(mpbsList)

#################################################################################################
# READING BIAS TABLE
#################################################################################################

biasTableF = dict()
biasTableFFile = open(biasTableFFileName,"r")
for line in biasTableFFile:
  ll = line.strip().split("\t")
  biasTableF[ll[0]] = float(ll[1])
biasTableFFile.close()

biasTableR = dict()
biasTableRFile = open(biasTableRFileName,"r")
for line in biasTableRFile:
  ll = line.strip().split("\t")
  biasTableR[ll[0]] = float(ll[1])
biasTableRFile.close()

#################################################################################################
# EVALUATING PROTECTION
#################################################################################################

# Initialization
protectDict = dict()
bam = Samfile(dnaseFileName,"rb")

# Iterating on MPBSs
for mpbs in mpbsList:

  # Fetching MPBSs
  grepFileName = outputFileName+"_grepmpbs.bed"
  to_remove.append(grepFileName)
  os.system("grep \"\t\""+mpbs+"\"\t\" "+mpbsFileNameBed+" | cut -f 1,2,3 | sort -k1,1 -k2,2n > "+grepFileName)

  # Intersect with footprints
  intFileName = outputFileName+"_fpint.bed"
  to_remove.append(intFileName)
  os.system("intersectBed -a "+grepFileName+" -b "+fpFileNameBed+" -wa -u > "+intFileName)

  # Iterating on MPBSs
  intFile = open(intFileName,"r")
  spr = 0.0
  counter = 0.0
  for line in intFile:

    # Fetching signal
    ll = line.strip().split("\t")
    mLen = int(ll[2]) - int(ll[1])
    mid = (int(ll[1])+int(ll[2]))/2
    p1 = max(mid - halfWindow,0)
    p2 = mid + halfWindow

    # Fetch raw signal
    pileup_region = PileupRegion(p1,p2,1)
    if(ps_version == "0.7.5"):
      bam.fetch(reference=ll[0], start=p1, end=p2, callback = pileup_region)
    else:
      iter = bam.fetch(reference=ll[0], start=p1, end=p2)
      for alignment in iter: pileup_region.__call__(alignment)
    raw_signal = array([min(e,initial_clip) for e in pileup_region.vector])
    
    # Std-based clipping
    mean = raw_signal.mean()
    std = raw_signal.std()
    clip_signal = [min(e, mean + (10 * std)) for e in raw_signal]

    # Bias Correction
    correctedSignal = bias_correction(bam, clip_signal, biasTableF, biasTableR, genomeFileName, ll[0], p1, p2)

    # Summing min value to signal
    stdzSignal = [e+minValue for e in correctedSignal]

    # Evaluating protection score
    signalHalfLen = len(stdzSignal)/2
    motifHalfLen = motifExt/2
    nc = sum(stdzSignal[signalHalfLen-motifHalfLen:signalHalfLen+motifHalfLen])
    nr = sum(stdzSignal[signalHalfLen+motifHalfLen:signalHalfLen+motifHalfLen+motifExt])
    nl = sum(stdzSignal[signalHalfLen-motifHalfLen-motifExt:signalHalfLen-motifHalfLen])
    spr += ((nr+nl)/2) - nc
    counter += 1.0

  intFile.close()

  # Updating protection
  protmean = "NA"
  if(counter > 0):
    protmean = spr/counter
  protectDict[mpbs] = str(protmean)

#################################################################################################
# WRITING RESULTS
#################################################################################################

# Writing protection score
outputFile = open(outputFileName,"w")
for m in mpbsList: outputFile.write("\t".join([m,protectDict[m]])+"\n")
outputFile.close()

# Termination
for e in to_remove: 
  if(os.path.isfile(e)): os.system("rm "+e)


