from __future__ import print_function

import sys

from rgt.GeneSet import GeneSet
from rgt.MotifSet import MotifSet

jaspar = '/home/ivan/projects/reg-gen/data/motifs/jaspar_vertebrates.mtf'
uniprobe = '/home/ivan/projects/reg-gen/data/motifs/uniprobe_primary.mtf'
internal = '/home/ivan/projects/reg-gen/data/motifs/internal.mtf'

motif_set = MotifSet()
motif_set.read_mtf([jaspar, uniprobe, internal])

motifs = [(l.strip("\n")).split("\t")[0] for l in open(sys.argv[1])]

geneset_file = sys.argv[2]

search_mode = sys.argv[3]

genes = GeneSet("DC Genes")
genes.read_expression(geneset_file)

filtered = motif_set.filter(motifs, key_type="name")

[filtered_genes, m_g, g_m] = filtered.filter(genes.genes, key_type="gene_names", search=search_mode)

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

