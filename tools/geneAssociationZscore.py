import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
#from fisher import pvalue
import scipy.stats


outdir=""

back=False
designFile = sys.argv[1]
anotationPath = sys.argv[2]
randomize = int(sys.argv[3])
backGroundPeaks=False
if len(sys.argv) > 4:
  backGroundPeaksName = sys.argv[4]
  backBed=GenomicRegionSet("BACK")
  backBed.read_bed(backGroundPeaksName)  
  backGroundPeaks=True
   
distance=50000
if len(sys.argv) > 5:
  distance=len(sys.argv[5])

if len(sys.argv) > 6:
  outdir=sys.argv[6]


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

if len(outdir)>0:
  outGene = open(outdir+"/mapped_genes.txt","w")
  
print genesets

for region in exps.get_regionsets():
    for j,n in enumerate(range(randomize)):
            if backGroundPeaks:
              print len(backBed), len(region)
              br=backBed.random_subregions(len(region))
            else:
              br=region.random_regions('hg19',total_size=len(region),overlap_result=True, overlap_input=True)

            br.write_bed(str(j)+"random.bed")
    for g in genesets:
        #print region,g
        bed = GenomicRegionSet("")
        [degenes,de_peak_genes, mappedGenes, totalPeaks,bla] = bed.filter_by_gene_association(region.fileName,g.genes,geneFile,genomeFile,threshDist=distance)
        randomRes=[]
        #backBed=GenomicRegionSet("BACK")    
        #backBed.read_bed(backGroundPeaks)
        for j,n in enumerate(range(randomize)):
            backUP=GenomicRegionSet("BACKUP")
            [back_de_genes,back_de_peak_genes, back_mappedGenes, back_totalPeaks,bla] = backUP.filter_by_gene_association(str(j)+"random.bed",g.genes,geneFile,genomeFile,threshDist=distance)
            randomRes.append(back_de_peak_genes)
            #print str(j)+"random.bed"
        randomRes=numpy.array(randomRes)
        #print randomRes
        a=de_peak_genes
        m=numpy.mean(randomRes)
        s=numpy.std(randomRes)
        z=(a-m)/s
        prop_de=de_peak_genes/float(degenes)
        prop_back=m/float(degenes)
        p= scipy.stats.norm.sf(z)
        print region.name,g.name,a,m,z,degenes,mappedGenes,len(allgenes),prop_de,prop_back,prop_de/prop_back,p,degenes

        if len(outdir)>0:
          outGene.write(region.name+"\t"+g.name+"\t"+("\t".join(bed.genes))+"\n")  
          bed.write_bed(outdir+"/"+g.name+"_"+region.name+".bed")  

