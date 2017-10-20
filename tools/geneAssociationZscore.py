import sys
import os.path
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
#from fisher import pvalue
import scipy.stats


outdir=""

back=False
designFile = sys.argv[1]
genomeName = sys.argv[2]
geneFile = sys.argv[3]
randomize = int(sys.argv[4])
backGroundPeaks=False
if len(sys.argv) > 5:
  backGroundPeaksName = sys.argv[6]
  backBed=GenomicRegionSet("BACK")
  backBed.read(backGroundPeaksName)  
  backGroundPeaks=True
   
distance=50000
if len(sys.argv) > 6:
  distance=len(sys.argv[6])

if len(sys.argv) > 7:
  outdir=sys.argv[7]


#genomeFile=anotationPath+"chrom.sizes"
#geneFile=anotationPath+"association_file.bed"

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

if len(outdir)>0:
  outGene = open(outdir+"/mapped_genes.txt","w")
  
print genesets

random_regions=[]

for region in exps.get_regionsets():
    for j,n in enumerate(range(randomize)):
            if backGroundPeaks:
              print len(backBed), len(region)
              br=backBed.random_subregions(len(region))
            else:
              br=region.random_regions('hg19',total_size=len(region),overlap_result=True, overlap_input=True)

            #br.write(str(j)+"random.bed")
            random_regions.append(br)
    for g in genesets:
        region_aux=deepcopy(region)
        [all_genes,mapped_genes,all_proxs,mapped_proxs] = region_aux.filter_by_gene_association(gene_set=g,organism=genomeName,threshDist=distance)

#            self in order to keep only the coordinates associated to genes which are in gene_set
#
#            - all_genes = GeneSet that contains all genes associated with the coordinates
#            - mapped_genes = GeneSet that contains the genes associated with the coordinates which are in gene_set
#            - all_proxs = List that contains all the proximity information of genes associated with the coordinates
#            - mapped_proxs = List that contains all the proximity information of genes associated with the coordinates which are in gene_set


        randomRes=[]
        #backBed=GenomicRegionSet("BACK")    
        #backBed.read(backGroundPeaks)
        for br in random_regions:
          random=deepcopy(br)
          [back_all_genes,back_mapped_genes,back_all_proxs,back_mapped_proxs] = random.filter_by_gene_association(gene_set=g,organism=genomeName,threshDist=distance)
          randomRes.append(len(back_mapped_genes))

        randomRes=numpy.array(randomRes)
        a=len(mapped_genes)
        m=numpy.mean(randomRes)
        s=numpy.std(randomRes)
        z=(a-m)/s
        #prop_de=a/float(degenes)
        #prop_back=m/float(degenes)
        p= scipy.stats.norm.sf(z)
        print region.name,g.name,a,m,z,len(all_genes),len(allgenes),p

        #if len(outdir)>0:
        #  outGene.write(region.name+"\t"+g.name+"\t"+("\t".join(bed.genes))+"\n")  
        #  bed.write(outdir+"/"+g.name+"_"+region.name+".bed")  

