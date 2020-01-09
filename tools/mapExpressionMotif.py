

import sys

from rgt.GeneSet import GeneSet
from rgt.MotifSet import MotifSet

motifs = [(l.strip("\n")).split("\t")[0] for l in open(sys.argv[1])]
geneset_file = sys.argv[2]
search_mode = sys.argv[3]

# preload all available motifs from the repositories
motif_set = MotifSet(preload_motifs=True)

genes = GeneSet("DC Genes")
genes.read_expression(geneset_file)

# take only a subset of the motifs (using their exact names)
motif_set, _, _ = motif_set.filter(motifs, key_type="name")

# of these new motif set, take the subset of those matching these gene names
# (we only care about the motif2gene mapping)
_, m_g, _ = motif_set.filter(genes.genes, key_type="gene_names", search=search_mode)

genes_found = []
not_found = []
print("\t\t" + ("\t".join(genes.cond)))
for m in motifs:
    try:
        sel_genes = m_g[m]
        for g in sel_genes:
            print(m + "\t" + g + "\t" + ("\t".join([str(v) for v in genes.values[g]])))
            genes_found.append(g)
    except:
        not_found.append(m)

print(not_found)
print(set(genes.genes).difference(genes_found))
