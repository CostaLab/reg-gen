import os
import sys
import pickle
import shutil
from rgt.GeneSet import GeneSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.AnnotationSet import AnnotationSet
from rgt.Util import GenomeData
from rgt.tdf.triplexTools import get_rna_region_str, rna_associated_gene, connect_rna

def load_dump(path, filename):
    file = open(os.path.join(path,filename),'r')
    object = pickle.load(file)
    file.close()
    print("\tLoading from file: "+filename)
    return object

def dump(object, path, filename):
    file = open(os.path.join(path,filename),'wb')
    pickle.dump(object,file)
    file.close()
    print("\tDump to file: "+filename)


def check_geneset_empty(geneset, filename, outdir):
    if len(geneset) == 0:
        print("Error: No genes are loaded from: " + filename)
        print("Please check the format.")
        try:
            shutil.rmtree(outdir)
        except:
            pass
        sys.exit(1)

class Input(object):
    """Process all the input files for TDF"""
    def __init__(self, pars):
        self.genome = GenomeData(organism=pars.organism)
        self.organism = pars.organism
        self.outdir = pars.o
        self.target_seq = os.path.join(pars.o, "dna_targets.fa")
        self.background_seq = os.path.join(pars.o, "dna_background.fa")
        self.pars = pars

    class DNA(object):
        def __init__(self):
            pass

        #     if format == "degenes":
        #         self.degenes(filename)
        #     elif format == "bed":
        #         self.bed_input()

        def degenes(self):
            dumpname = "_".join(["dump", Input.pars.de.rpartition("/")[-1].rpartition(".")[0],
                                 "".join([Input.pars.filter_havana,
                                          Input.pars.protein_coding ,
                                          Input.pars.known_only])])

            try:
                data = load_dump(path=Input.outfir, filename=dumpname)
                self.de_genes = data[0]
                self.nde_genes = data[1]
                self.de_promoters = data[2]
                self.nde_promoters = data[3]
                if Input.pars.score:
                    self.scores = data[4]
            except:

                #######################
                # Load DE genes
                self.de_genes = GeneSet("de genes")
                self.de_genes.read(Input.pars.de, score=Input.pars.score)
                check_geneset_empty(self.de_genes, Input.pars.de, Input.pars.o)


                ann = AnnotationSet(gene_source=Input.organism,
                                    alias_source=Input.organism,
                                    filter_havana=Input.pars.filter_havana,
                                    protein_coding=Input.pars.protein_coding,
                                    known_only=Input.pars.known_only)
                self.de_promoters, unmapped = ann.get_promoters(promoterLength=Input.pars.pl,
                                                      gene_set=self.de_genes,
                                                      regiondata=Input.pars.score,
                                                      unmaplist=True, variants=False)
                if Input.pars.score:
                    self.scores = self.de_promoters.get_score_dict()

                #######################
                # Non DE genes
                self.nde_genes = GeneSet("nde genes")
                self.nde_genes.genes = [g for g in ann.symbol_dict.values() if g not in self.de_genes]

                self.nde_promoters = ann.get_promoters(promoterLength=Input.pars.pl,
                                                       gene_set=self.nde_genes, unmaplist=False)
                # print2(summary, "   \t" + str(len(nde_prom)) + "\tmapped non-target promoters")
                self.nde_promoters.merge(namedistinct=True)

                #######################
                # Loading score
                data = [self.de_genes, self.nde_genes,
                        self.de_regide_promotersons, self.nde_promoters]
                if Input.pars.score:
                    data.append(self.scores)
                dump(object=data, path=Input.pars.o, filename=dumpname)

            pass

        def annotation(self):
            pass

        def bed_input(self):
            self.target_regions = GenomicRegionSet("targets")
            self.target_regions.read(Input.pars.bed)
            self.target_regions.remove_duplicates()
            if Input.pars.score:
                self.scores = self.target_regions.get_score_dict()
            self.nontarget_regions = GenomicRegionSet("background")
            self.nontarget_regions.read(Input.pars.bg)
            self.nontarget_regions.remove_duplicates()

    class RNA:
        def __init__(self):
            self.name = Input.pars.rn
            self.stat = {}

        def get_rna_info(self, filename, expfile=None):
            """Getting the rna region from the information header with the pattern:
                            REGION_chr3_51978050_51983935_-_
                        or  chr3:51978050-51983935 -    """
            self.regions = get_rna_region_str(filename)
            # print(self.rna_regions)
            if self.regions:
                r_genes = rna_associated_gene(rna_regions=self.regions,
                                              name=self.name, organism=Input.organism)
            if self.regions and len(self.regions[0]) == 5:
                self.expression = float(self.regions[0][-1])

            elif expfile:
                with open(expfile) as f:
                    for line in f:
                        l = line.strip().split()
                        if self.name == l[0].partition(".")[0]:
                            self.expression = l[1]
            else:
                self.expression = "n.a."

        def connect_rna(self):
            d = connect_rna(Input.pars.r, Input.pars.o, self.name)
            self.stat["exons"] = str(d[0])
            self.stat["seq_length"] = str(d[1])

