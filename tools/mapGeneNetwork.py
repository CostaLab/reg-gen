import sys

from rgt.GeneSet import GeneSet
from rgt.MotifSet import MotifSet

# motif databases

# files with p-values
enrichment_files = sys.argv[1]
# tfs to include in the network
factor_file = sys.argv[2]
# search mode to map factors to motifs (exact or inexact)
search_mode = sys.argv[3]
# pvalue cutoff for definition of active factors
pvalue = float(sys.argv[4])
# output file
out = sys.argv[5]
# genes to be used as potential targets 
targets = None
if len(sys.argv) > 6:
    targets_file = sys.argv[6]
    # reading targets
    targets = GeneSet("genes")
    targets.read(targets_file)

# starting motif databases
if len(sys.argv) > 7:
    motif_set = MotifSet(preload_motifs=False)
    motif_set.read_mtf([sys.argv[7]])
else:
    motif_set = MotifSet(preload_motifs=True)

# reading genes 
factors = GeneSet("genes")
factors.read(factor_file)

# we only want a subset of the motif set
motif_set = motif_set.filter(factors.genes, key_type="gene_names", search=search_mode)

motif_set.read_enrichment(enrichment_files, pvalue)
motif_set.write_network(targets, out, pvalue)
