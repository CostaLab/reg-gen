import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from fisher import pvalue

back=False
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
bedGenes.read(geneFile)
allgenes=[]
for r in bedGenes:
 allgenes.append(r.name)
allgenes=list(set(allgenes))

genesets=exps.get_genesets()

if len(sys.argv) > 3:
    back=True
    backGroundPeaks = sys.argv[3]
    backBed=GenomicRegionSet("BACK")
    backBed.read(backGroundPeaks)


backBed=GenomicRegionSet("BACK")    
backBed.read(backGroundPeaks)
backUP=GenomicRegionSet("BACKUP")
[back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks] = backUP.filter_by_gene_association(backGroundPeaks,genesets[0],geneFile,genomeFile)
prop_back=back_mappedGenes/float(len(allgenes))

for g in genesets:
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile)
        #print degenes
        #print bed.genes
        a=de_peak_genes
        b=degenes-de_peak_genes
        c=back_mappedGenes-de_peak_genes
        d=len(allgenes)-b-c-a
        prop_de=de_peak_genes/float(degenes)
        p= pvalue(a,b,c,d)
        print region.name,g.name,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail,p.left_tail


'''for g in genesets:
    if back:
        backBed=GenomicRegionSet("BACK")    
        backBed.read(backGroundPeaks)
        backUP=GenomicRegionSet("BACKUP")
        [back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks] = backUP.filter_by_gene_association(backGroundPeaks,g,geneFile,genomeFile)
        prop_back=len(backUP)/float(len(backBed))
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile)
        if back:
            a=len(bed)
            b=len(region)-a
            c=len(backUP)-a
            d=len(backBed)-a-b-c
            prop_de=len(bed)/float(len(region))            
        else:
            a=de_peak_genes
            b=degenes-de_peak_genes
            c=mappedGenes-de_peak_genes
            d=len(allgenes)-b-c-a
            prop_de=de_peak_genes/float(degenes)
            prop_back=mappedGenes/float(len(allgenes))
        p= pvalue(a,b,c,d)
        print region.name,g.name,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail,p.left_tail'''


