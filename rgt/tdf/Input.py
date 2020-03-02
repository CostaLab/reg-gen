# Python Libraries


import os
import sys
import pickle
import shutil
from ..GeneSet import GeneSet
from ..GenomicRegionSet import GenomicRegionSet
from ..AnnotationSet import AnnotationSet
from ..Util import GenomeData
from ..tdf.triplexTools import get_rna_region_str, rna_associated_gene, connect_rna


def load_dump(path, filename):
    f = open(os.path.join(path, filename), 'r')
    object_dump = pickle.load(f)
    f.close()
    print("\tLoading from file: " + filename)
    return object_dump


def dump(object, path, filename):
    f = open(os.path.join(path, filename), 'wb')
    pickle.dump(object, f)
    f.close()
    print("\tDump to file: " + filename)


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
        self.dna = self.DNA(pars=pars)

        self.rna = self.RNA(pars=pars)
        self.rna.get_rna_info()
        self.rna.connect_rna()

    class DNA(object):
        def __init__(self, pars):
            self.pars = pars
            self.scores = None

        #     if format == "degenes":
        #         self.degenes(filename)
        #     elif format == "bed":
        #         self.bed_input()

        def degenes(self):
            target_name = self.pars.de.rpartition("/")[-1].rpartition(".")[0]
            dumpname = ".".join(["tdf", "dump", self.pars.organism, target_name])
            try:
                data = load_dump(path=os.getcwd(), filename=dumpname)
                self.de_genes = data[0]
                self.nde_genes = data[1]
                self.target_regions = data[2]
                self.nontarget_regions = data[3]
                if self.pars.score:
                    self.scores = data[4]
            except:

                #######################
                # Load DE genes
                self.de_genes = GeneSet("de genes")
                self.de_genes.read(self.pars.de, score=self.pars.score)
                check_geneset_empty(self.de_genes, self.pars.de, self.pars.o)

                ann = AnnotationSet(gene_source=self.pars.organism,
                                    alias_source=self.pars.organism,
                                    filter_havana=self.pars.filter_havana,
                                    protein_coding=self.pars.protein_coding,
                                    known_only=self.pars.known_only)
                self.target_regions, unmapped = ann.get_promoters(promoter_length=self.pars.pl, tss=self.pars.tss,
                                                                  gene_set=self.de_genes,
                                                                  regiondata=self.pars.score,
                                                                  unmaplist=True, variants=False)
                if self.pars.score:
                    self.scores = self.target_regions.get_score_dict()

                #######################
                # Non DE genes
                self.nde_genes = GeneSet("nde genes")
                if self.pars.bg and self.pars.bg.endswith(".txt"):
                    self.nde_genes.read(self.pars.bg)
                else:
                    self.nde_genes.genes = [g for g in list(ann.symbol_dict.values()) if g not in self.de_genes]

                self.nontarget_regions = ann.get_promoters(promoter_length=self.pars.pl,
                                                           gene_set=self.nde_genes, unmaplist=False)
                # print2(summary, "   \t" + str(len(nde_prom)) + "\tmapped non-target promoters")
                self.nontarget_regions.merge(namedistinct=True)

                #######################
                # Loading score
                data = [self.de_genes, self.nde_genes, self.target_regions, self.nontarget_regions]
                if self.pars.score:
                    data.append(self.scores)
                dump(object=data, path=os.getcwd(), filename=dumpname)

            pass

        def annotation(self):
            pass

        def de_bed_input(self):
            self.target_regions = GenomicRegionSet("targets")
            self.target_regions.read(self.pars.bed)
            self.target_regions.remove_duplicates()
            if self.pars.score:
                self.scores = self.target_regions.get_score_dict()
            self.nontarget_regions = GenomicRegionSet("background")
            self.nontarget_regions.read(self.pars.bg)
            self.nontarget_regions.remove_duplicates()

        def bed_input(self, bed):
            self.target_regions = GenomicRegionSet("targets")
            self.target_regions.read(bed)
            self.target_regions.remove_duplicates()

            self.target_regions = self.target_regions.gene_association(organism=self.pars.organism)
            if self.pars.score:
                self.scores = self.target_regions.get_score_dict()

                # for rbs in self.rbss:
                #     tr = len(self.txp.merged_dict[rbs])
                #     self.counts_tr[rbs] = [tr, len(self.dna_region) - tr]
                #     self.counts_dbs[rbs] = len(self.txpf.merged_dict[rbs])

    class RNA:
        def __init__(self, pars):
            self.name = pars.rn
            self.pars = pars
            self.stat = {}

        def get_rna_info(self, expfile=None):
            """Get the RNA information such as loci from the FASTA header following the pattern as below:
                REGION_chr3_51978050_51983935_-_  or
                chr3:51978050-51983935 -
            """
            filename = self.pars.r
            self.regions = get_rna_region_str(filename)
            # print(self.rna_regions)
            if self.regions:
                r_genes = rna_associated_gene(rna_regions=self.regions,
                                              name=self.name, organism=self.pars.organism)
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
            d = connect_rna(self.pars.r, self.pars.o, self.name)
            self.num_exons = d[0]
            self.seq_length = d[1]
