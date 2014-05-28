from rgt.MotifSet import MotifSet
from rgt.GeneSet import GeneSet
import sys
import glob
import os.path

# motif databases

jaspar='/home/ivan/projetos/reg-gen/data/motifs/jaspar_vertebrates.mtf'
uniprobe='/home/ivan/projetos/reg-gen/data/motifs/uniprobe_primary.mtf'
internal='/home/ivan/projetos/reg-gen/data/motifs/internal.mtf'

   

enrichment_files=sys.argv[1]
geneset_file=sys.argv[2]
search_mode=sys.argv[3]
pvalue=float(sys.argv[4])
out=sys.argv[5]

# starting motif databases
motif_set = MotifSet()
motif_set.read_file([jaspar,uniprobe,internal])

# reading genes 
genes=GeneSet("genes")
genes.read(geneset_file)

# reading networks
for f in glob.glob(enrichment_files): 
  # use last dir name as name for condition
  condition=os.path.dirname(f)
  condition=condition.split("/")[-1]
  motif_set.read_motif_targets(f,condition,pvalue)

motif_set.write_cytoscape_network(genes,search_mode,out)

