import sys
from rgt.SetGenomicRegions import *

bedFile = sys.argv[1]
geneList = sys.argv[2]
geneFile = sys.argv[3]
genomeFile = sys.argv[4]

bed = SetGenomicRegions(bedFile)
bed.filterByGeneAssociation(bedFile,geneList,geneFile,genomeFile)

bed.writeBed(bedFile+"_gene.bed")

