

import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
import scipy.stats

gprofilerline="http://biit.cs.ut.ee/gprofiler/index.cgi?organism=mmusculus&query="

back=False
designFile = sys.argv[1]
anotationPath = sys.argv[2]
genomeFile=anotationPath+"chrom.sizes"
geneFile=anotationPath+"association_file.bed"

exps=ExperimentalMatrix()
exps.read(designFile)

for region in exps.get_regionsets():

        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks, blad] = bed.filter_by_gene_association(region.fileName,None,geneFile,genomeFile,threshDist=50000)
        print mappedGenes
        print region.name+"\t"+("\t".join(bed.genes))
