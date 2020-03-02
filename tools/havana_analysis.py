
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
from rgt.AnnotationSet import AnnotationSet
import argparse 
import os        



##################################################################################
parser = argparse.ArgumentParser(description='Compare the number of genes among different filters on annotation', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', type=str, help="Gene list")
parser.add_argument('-organism', type=str, help="Define the organism")
args = parser.parse_args()


if os.path.isfile(args.i):
    print()

    for switch in [True, False]:
        print("HAVANA filter:\t\t"+str(switch))
        print("protein_coding filter:\t"+"True")
        print("known_only filter:\t"+"True")
        ann = AnnotationSet(args.organism, alias_source=args.organism,
                            filter_havana=switch, 
                            protein_coding=True, 
                            known_only=True)
        
        genes = GeneSet("genes")
        genes.read(args.i)
        print("\tInput gene number: \t"+str(len(genes)))
        print("\tFixing the names into Ensembl ID:")
        de_ensembl, unmap_gs, ensembl2symbol = ann.fix_gene_names(gene_set=genes, output_dict=True, mute_warn=True)
        print("\t\tMapped:\t\t"+str(len(de_ensembl)))
        print("\t\tUnmapped:\t"+str(len(unmap_gs)))

        genes.genes = de_ensembl
        
        de_prom, unmapped_gene_list = ann.get_promoters(promoter_length=1000,
                                                        gene_set=genes,
                                                        unmaplist=True)
        print("\tGetting promoters:")
        print("\t\tMapped:\t\t"+str(len(de_prom)))
        print("\t\tUnmapped:\t"+str(len(unmapped_gene_list)))

        de_prom.merge(namedistinct=True)
        print("\tMerging promoters by names:")
        print("\t\tMerged:\t\t"+str(len(de_prom)))


        print()
        nde_ensembl = [ g for g in list(ann.symbol_dict.keys()) if g not in de_ensembl ]
        print("\tBackground genes:\t"+str(len(nde_ensembl)))
        nde_gene = GeneSet("nde genes")
        nde_gene.genes = nde_ensembl

        nde_prom, unmapped_gene_list = ann.get_promoters(promoter_length=1000,
                                                         gene_set=nde_gene,
                                                         unmaplist=True)
        print("\tGetting background promoters:")
        print("\t\tMapped:\t\t"+str(len(nde_prom)))
        print("\t\tUnmapped:\t"+str(len(unmapped_gene_list)))

        nde_prom.merge(namedistinct=True)
        print("\tMerging background promoters by names:")
        print("\t\tMerged:\t\t"+str(len(nde_prom)))
