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
bedGenes.read_bed(geneFile)
allgenes=[]
for r in bedGenes:
 allgenes.append(r.name)
allgenes=list(set(allgenes))

genesets=exps.get_genesets()

if len(sys.argv) > 3:
    back=True
    backGroundPeaks = sys.argv[3]
    backBed=GenomicRegionSet("BACK")
    backBed.read_bed(backGroundPeaks)



for g in genesets:
    if back:
        backBed=GenomicRegionSet("BACK")    
        backBed.read_bed(backGroundPeaks)
        backUP=GenomicRegionSet("BACKUP")
        [back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks] = backUP.filter_by_gene_association(backGroundPeaks,g,geneFile,genomeFile)
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile)
        prop_de=de_peak_genes/float(degenes)
        prop_back=mappedGenes/float(len(allgenes))
        if back:
            a=len(bed)
            b=len(region)-a
            c=len(backUP)-a
            d=len(backBed)-a-b-c
            prop_de=len(bed)/float(len(region))
            prop_back=len(backUP)/float(len(backBed))
        else:
            a=de_peak_genes
            b=degenes-de_peak_genes
            c=mappedGenes-de_peak_genes
            d=len(allgenes)-b-c-a
        p= pvalue(a,b,c,d)
        print region.name,g.name,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail,p.left_tail



'''for g in genesets:
    if back:
        backBed=GenomicRegionSet("BACK")    
        backBed.read_bed(backGroundPeaks)
        [back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks] = backBed.filter_by_gene_association(backGroundPeaks,g,geneFile,genomeFile)
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile)
        prop_de=de_peak_genes/float(degenes)
        prop_back=mappedGenes/float(len(allgenes))
        if back:
            a=de_peak_genes
            b=mappedGenes-a
            c=back_de_peak_genes-a
            d=back_mappedGenes-a-b-c
            prop_de=de_peak_genes/float(mappedGenes)
            prop_back=back_de_peak_genes/float(back_mappedGenes)
        else:
            a=de_peak_genes
            b=degenes-de_peak_genes
            c=mappedGenes-de_peak_genes
            d=len(allgenes)-b-c-a
        p= pvalue(a,b,c,d)
        print region.name,g.name,a,b,c,d,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p.right_tail,p.left_tail'''

