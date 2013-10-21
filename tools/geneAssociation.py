import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from fisher import pvalue


designFile = sys.argv[1]
anotationPath = sys.argv[2]
genomeFile=anotationPath+"chrom.sizes"
geneFile=anotationPath+"association_file.bed"

exps=ExperimentalMatrix()
exps.read(designFile)

beds=[]
geneLists=[]

#this should be improved
bedGenes = GenomicRegionSet(geneFile)
bedGenes.read_bed(geneFile)
allgenes=[]
for r in bedGenes:
 allgenes.append(r.name)
allgenes=list(set(allgenes))



for g in exps.get_genesets():
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile)
        #print degenes,de_peak_genes, mappedGenes, totalPeaks
        prop_de=de_peak_genes/float(degenes)
        prop_back=mappedGenes/float(len(allgenes))
        a=de_peak_genes
        b=degenes-de_peak_genes
        c=mappedGenes-de_peak_genes
        d=len(allgenes)-b-c-a
        p= pvalue(a,b,c,len(allgenes)-b-c-a)
        print region.name,g.name,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail,p.left_tail


