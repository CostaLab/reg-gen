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
factor_file=sys.argv[2]
search_mode=sys.argv[3]
pvalue=float(sys.argv[4])
out=sys.argv[5]
filter_targets=[]
targets=None
if len(sys.argv) > 6:
  targets_file=sys.argv[6]
  # reading targets 
  targets=GeneSet("genes")
  targets.read(targets_file)


# starting motif databases
motif_set = MotifSet()
motif_set.read_file([jaspar,uniprobe,internal])

# reading genes 
factors=GeneSet("genes")
factors.read(factor_file)



# reading networks
for f in glob.glob(enrichment_files): 
  # use last dir name as name for condition
  condition=os.path.dirname(f)
  condition=condition.split("/")[-1]
  motif_set.read_motif_targets(f,condition,pvalue)

#print motif_set.motifs_map
motif_set.write_cytoscape_network(factors,search_mode,out,targets)

