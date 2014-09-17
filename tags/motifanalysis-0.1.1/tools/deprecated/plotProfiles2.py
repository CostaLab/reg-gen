import sys
from os.path import basename
from rgt.CoverageSet import *



# A_P4      B_P4      C_P4     A_P13     B_P13     C_P13
#1.1982482 1.0613452 2.1236085 0.5630403 0.6752216 1.0295918

factors=[0.5630403,1.1982482,0.6752216,1.0613452,1.0295918,2.1236085]


bedFile= sys.argv[1]
bamFile= sys.argv[2]
bamFile2= sys.argv[3]
bamFile3= sys.argv[4]
bamFile4= sys.argv[5]
bamFile5=sys.argv[6]
bamFile6= sys.argv[7]

bedname=basename(bedFile)[:-4]


bed= SetGenomicRegions(bedname)
bed.readBed(bedFile)
bed.extend(19500,19500)

cov=CoverageSet(bed.name+"_"+basename(bamFile)[:-4],bed)
cov.coverageFromBam(bamFile,step=50,window=50)
#cov.normFactor(factors[0])
cov.normFPKM()
cov.plot(log=False)

cov2=CoverageSet(bed.name+"_"+basename(bamFile2)[:-4],bed)
cov2.coverageFromBam(bamFile2,step=50,window=50)
#cov2.normFactor(factors[1])
cov2.normFPKM()
cov2.plot(log=False)

cov.diff(cov2)
cov.plot(log=False,name=bedname+"_Diff_A_Late_Early")
#cov.abs()
#cov.plot(log=False,name=bedname+"_Abs_Diff_A_Late_Early")

cov=CoverageSet(bed.name+"_"+basename(bamFile3)[:-4],bed)
cov.coverageFromBam(bamFile3,step=50,window=50)
cov.normFPKM()
##cov.normFactor(factors[2])
cov.plot(log=False)
#
cov2=CoverageSet(bed.name+"_"+basename(bamFile4)[:-4],bed)
cov2.coverageFromBam(bamFile4,step=50,window=50)
#cov2.normFactor(factors[3])
cov2.normFPKM()
cov2.plot(log=False)
#
cov.diff(cov2)
cov.plot(log=False,name=bedname+"_Diff_B_Late_Early")
##cov.abs()
##cov.plot(log=True,name=bedname+"_Abs_Diff_B_Late_Early")
#
#
cov=CoverageSet(bed.name+"_"+basename(bamFile5)[:-4],bed)
cov.coverageFromBam(bamFile5,step=50,window=50)
cov.normFPKM()
#cov.normFactor(factors[4])
cov.plot(log=False)
#
cov2=CoverageSet(bed.name+"_"+basename(bamFile6)[:-4],bed)
cov2.coverageFromBam(bamFile6,step=50,window=50)
#cov2.normFactor(factors[5])
cov2.normFPKM()
cov2.plot(log=False)
#
cov.diff(cov2)
cov.plot(log=False,name=bedname+"_Diff_C_Late_Early")
##cov.abs()
##cov.plot(log=True,name=bedname+"_Abs_Diff_C_Late_Early")
