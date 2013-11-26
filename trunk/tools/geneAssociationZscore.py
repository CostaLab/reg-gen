import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
#from fisher import pvalue
import scipy.stats

back=False
designFile = sys.argv[1]
anotationPath = sys.argv[2]
backGroundPeaks = sys.argv[3]
randomize = int(sys.argv[4])
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


for g in genesets:
    for region in exps.get_regionsets():
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks] = bed.filter_by_gene_association(region.fileName,g,geneFile,genomeFile,threshDist=500)
        randomRes=[]
        backBed=GenomicRegionSet("BACK")    
        backBed.read_bed(backGroundPeaks)
        for j,n in enumerate(range(randomize)):
            br=backBed.randomRegions(totalPeaks)
            br.write_bed(str(j)+"random.bed")
            backUP=GenomicRegionSet("BACKUP")
            [back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks] = backUP.filter_by_gene_association(str(j)+"random.bed",g,geneFile,genomeFile,threshDist=500)
            randomRes.append(back_de_peak_genes)
        randomRes=numpy.array(randomRes)
        a=de_peak_genes
        m=numpy.mean(randomRes)
        s=numpy.std(randomRes)
        z=(a-m)/s
        prop_de=de_peak_genes/float(degenes)
        prop_back=m/float(degenes)
        p= scipy.stats.norm.sf(z)
        print region.name,g.name,a,m,z,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p,degenes


