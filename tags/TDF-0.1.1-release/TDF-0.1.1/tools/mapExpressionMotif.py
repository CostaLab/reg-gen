from rgt.MotifSet import MotifSet
from rgt.GeneSet import GeneSet
import sys

jaspar='/home/ivan/projects/reg-gen/data/motifs/jaspar_vertebrates.mtf'
uniprobe='/home/ivan/projects/reg-gen/data/motifs/uniprobe_primary.mtf'
internal='/home/ivan/projects/reg-gen/data/motifs/internal.mtf'

motif_set = MotifSet()
motif_set.read_file([jaspar,uniprobe,internal])

motifs=[(l.strip("\n")).split("\t")[0] for l in open(sys.argv[1])]

geneset_file=sys.argv[2]

search_mode=sys.argv[3]

genes=GeneSet("DC Genes")
genes.read_expression(geneset_file)
  
filtered=motif_set.filter_by_motifs(motifs)

[filtered_genes,g_m,m_g]=filtered.filter_by_genes(genes,search=search_mode)

genes_found=[]
not_found=[]
print "\t\t"+("\t".join(genes.cond))
for m in motifs:
  try:
    sel_genes=m_g[m]
    for g in sel_genes:
      print m+"\t"+g+"\t"+("\t".join([str(v) for v in genes.values[g]]))
      genes_found.append(g)
  except:
    not_found.append(m)

print not_found

import sets

print sets.Set(genes.genes).difference(genes_found)



#print filtered_genes

#print filtered.genes_map

#print filtered_genes.motifs_map

#print g_m

#print m_g
