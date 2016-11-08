# Python Libraries
from __future__ import print_function
import os
import time
import shutil
import datetime
from collections import *

# Local Libraries
import numpy
numpy.seterr(divide='ignore', invalid='ignore')

from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import scipy.cluster.hierarchy as sch
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter

# Distal Libraries
from rgt.GeneSet import GeneSet
from rgt.SequenceSet import SequenceSet
from rgt.AnnotationSet import AnnotationSet
from RNADNABindingSet import RNADNABindingSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.Util import SequenceType, Html, GenomeData, OverlapType
from triplexTools import dump, load_dump, print2, get_rna_region_str, connect_rna,\
    get_sequence, run_triplexator, dbd_regions, lineplot, value2str, rank_array,\
    split_gene_name, rna_associated_gene


# Color code for all analysis
target_color = "mediumblue"
nontarget_color = "darkgrey"
sig_color = "powderblue"


class PromoterTest:
    """Test the association between given triplex and differential expression genes"""
    def __init__(self, gene_list_file, bed, bg, organism, promoterLength, rna_name,
                 summary, temp, output, showdbs=None, score=False, scoreh=False,
                 filter_havana=True, protein_coding=False, known_only=True, gtf=None):
        """Initiation"""
        self.organism = organism
        if gtf: self.gtf = gtf
        else: self.gtf = organism
        genome = GenomeData(organism)
        self.rna_name = rna_name
        self.genome_path = genome.get_genome()
        self.showdbs = showdbs
        self.scores = None
        self.scoreh = scoreh
        self.motif = OrderedDict()
        self.sig_DBD = []
        self.stat = OrderedDict(name=rna_name, genome=self.gtf)

        # Input BED files
        if bed and bg:
            self.de_regions = GenomicRegionSet("de genes")
            self.nde_regions = GenomicRegionSet("nde genes")
            self.de_regions.read_bed(bed)
            self.de_regions.remove_duplicates()

            # score
            if score:
                self.scores = []
                for promoter in self.de_regions:
                    if promoter.data: self.scores.append(promoter.data.split("\t")[0])
                    else: self.scores = None

            try: self.de_regions = self.de_regions.gene_association(organism=self.organism)
            except: pass
            self.nde_regions.read_bed(bg)
            try: self.nde_regions = self.nde_regions.gene_association(organism=self.organism)
            except: pass
            self.nde_regions.remove_duplicates()

        # Input DE gene list
        else:
            if gtf:
                g = os.path.basename(gtf)
                dumpname = "_".join(["dump", gene_list_file.rpartition("/")[-1].rpartition(".")[0],
                                     g.rpartition(".")[0],
                                     filter_havana + protein_coding + known_only])
            else:
                dumpname = "_".join(["dump", gene_list_file.rpartition("/")[-1].rpartition(".")[0],
                                     filter_havana + protein_coding + known_only])
            try:
                data = load_dump(path=temp, filename=dumpname)
                self.de_gene = data[0]
                self.ensembl2symbol = data[1]
                self.nde_gene = data[2]
                self.de_regions = data[3]
                self.nde_regions = data[4]
                if score: self.scores = data[5]

            except:

                t1 = time.time()
                # Setting annotationSet
                if filter_havana== "T":
                    filter_havana = True
                else:
                    filter_havana = False
                if protein_coding == "T":
                    protein_coding = True
                else:
                    protein_coding = False
                if known_only == "T":
                    known_only = True
                else:
                    known_only = False

                ann = AnnotationSet(gene_source=self.gtf, alias_source=organism,
                                    filter_havana=filter_havana,
                                    protein_coding=protein_coding,
                                    known_only=known_only)

                t2 = time.time()
                print("\t" + str(datetime.timedelta(seconds=round(t2 - t1))))

                ######################################################
                # DE gene regions
                self.de_gene = GeneSet("de genes")

                if score:
                    self.de_gene.read_expression(geneListFile=gene_list_file, header=scoreh, valuestr=True)
                else:
                    self.de_gene.read(gene_list_file)
                # When there is no genes in the list
                if len(self.de_gene) == 0:
                    print("Error: No genes are loaded from: " + gene_list_file)
                    print("Please check the format.")
                    try:
                        shutil.rmtree(output)
                    except:
                        pass
                    sys.exit(1)

                print2(summary, "   \t" + str(len(self.de_gene)) + "\tinput genes ")
                # Generate a dict for ID transfer
                de_ensembl, unmap_gs, self.ensembl2symbol = ann.fix_gene_names(gene_set=self.de_gene, output_dict=True)
                print2(summary, "   \t" + str(len(de_ensembl)) + "\tgenes are mapped to Ensembl ID")
                if len(unmap_gs) > 0:
                    print2(summary, "   \t" + str(len(unmap_gs)) + "\tunmapped genes: " + ",".join(unmap_gs))

                # if "ENSG" in self.ensembl2symbol.keys()[0]: self.ensembl2symbol = None
                self.de_gene.genes = de_ensembl
                # print("After fixing")
                de_prom, unmapped_gene_list = ann.get_promoters(promoterLength=promoterLength,
                                                                gene_set=self.de_gene,
                                                                unmaplist=True)
                # print(len(unmapped_gene_list))
                print2(summary, "   \t" + str(len(de_prom)) + "\tmapped promoters")
                if len(unmapped_gene_list) > 0:
                    print2(summary, "   \t" + str(len(unmapped_gene_list)) + "\tunmapped promoters: " + ",".join(
                        unmapped_gene_list))

                if score: self.scores = []

                de_prom.merge(namedistinct=True, strand_specific=True)
                print2(summary, "   \t" + str(len(de_prom)) + "\tmerged promoters ")

                print2(summary, "   \t" + str(len(de_prom)) + "\tunique target promoters are loaded")
                for promoter in de_prom:
                    # if self.ensembl2symbol:
                    try:
                        gene_sym = self.ensembl2symbol[promoter.name]
                    except:
                        gene_sym = promoter.name

                    if score:
                        try:
                            try:
                                s = self.de_gene.values[promoter.name]
                            except:
                                s = self.de_gene.values[gene_sym.upper()]
                        except:
                            try:
                                print("Warning: " + promoter.name + "\tcannot be mapped to get its score.")
                            except:
                                print("Warning: " + gene_sym + "\tcannot be mapped to get its score.")
                            s = "0"
                        self.scores.append(s)

                self.de_regions = de_prom

                ######################################################
                # NonDE gene regions
                nde_ensembl = [g for g in ann.symbol_dict.keys() if g not in de_ensembl]
                print2(summary, "   \t" + str(len(nde_ensembl)) + "\tnon-target genes")
                self.nde_gene = GeneSet("nde genes")
                self.nde_gene.genes = nde_ensembl

                # Get promoters from nonDE gene
                nde_prom, unmapped_gene_list = ann.get_promoters(promoterLength=promoterLength,
                                                                 gene_set=self.nde_gene,
                                                                 unmaplist=True)
                print2(summary, "   \t" + str(len(nde_prom)) + "\tmapped non-target promoters")

                nde_prom.merge(namedistinct=True)
                print2(summary, "   \t" + str(len(nde_prom)) + "\tmerged non-target promoters")

                self.nde_regions = nde_prom

                print2(summary, "   \t" + str(len(nde_prom)) + "\tunique non-target promoters are loaded")
                # Loading score

                data = [self.de_gene, self.ensembl2symbol, self.nde_gene, self.de_regions, self.nde_regions]
                if score: data.append(self.scores)
                dump(object=data, path=temp, filename=dumpname)

            self.stat["target_regions"] = str(len(self.de_regions))
            self.stat["background_regions"] = str(len(self.nde_regions))

    def get_rna_region_str(self, rna):
        """Getting the rna region from the information header with the pattern:
                REGION_chr3_51978050_51983935_-_
            or  chr3:51978050-51983935 -    """
        self.rna_regions = get_rna_region_str(rna)

    def connect_rna(self, rna, temp):
        d = connect_rna(rna, temp, self.rna_name)
        self.stat["exons"] = str(d[0])
        self.stat["seq_length"] = str(d[1])

    def search_triplex(self, temp, l, e, c, fr, fm, of, mf, par, remove_temp=False, tp=False):
        print("    \tRunning Triplexator...")
        rna = os.path.join(temp, "rna_temp.fa")
        self.triplexator_p = [l, e, c, fr, fm, of, mf]

        # DE
        get_sequence(dir=temp, filename=os.path.join(temp, "de.fa"), regions=self.de_regions,
                     genome_path=self.genome_path)
        run_triplexator(ss=rna, ds=os.path.join(temp, "de.fa"),
                        output=os.path.join(temp, "de.txp"),
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, par=par, tp=tp)

        # non-DE
        get_sequence(dir=temp, filename=os.path.join(temp, "nde.fa"), regions=self.nde_regions,
                     genome_path=self.genome_path)
        run_triplexator(ss=rna, ds=os.path.join(temp, "nde.fa"),
                        output=os.path.join(temp, "nde.txp"),
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, par=par, tp=tp)

        # os.remove(os.path.join(args.o,"rna_temp.fa"))
        if remove_temp:
            os.remove(os.path.join(temp, "de.fa"))
            os.remove(os.path.join(temp, "nde.fa"))

    def count_frequency(self, temp, remove_temp, cutoff, l, obedp=False):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""

        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)
        # print([len_de, len_nde])

        self.frequency = {}
        self.frequency["promoters"] = {"de": OrderedDict(), "nde": OrderedDict()}

        ########################################################
        # Count the number of targeted promoters on each merged DNA Binding Domain
        print("\tCounting frequency of promoters on DBD...")
        self.txp_de = RNADNABindingSet("DE")

        self.txp_de.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=False)
        print("\t\t" + str(len(self.txp_de)) + "\tBinding de promoters")
        # txp_de.remove_duplicates()
        self.txp_de.merge_rbs(rm_duplicate=True, cutoff=cutoff,
                              region_set=self.de_regions, name_replace=self.de_regions)

        self.rbss = self.txp_de.merged_dict.keys()
        for rbs in self.rbss:
            # DE
            self.txp_de.merged_dict[rbs].remove_duplicates()
            l1 = len(self.txp_de.merged_dict[rbs])
            self.frequency["promoters"]["de"][rbs] = [l1, len_de - l1]

        self.txp_nde = RNADNABindingSet("non-DE")
        self.txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        print("\t\t" + str(len(self.txp_nde)) + "\tBinding nde promoters")
        # txp_nde.remove_duplicates()
        self.txp_nde.merge_rbs(rbss=self.rbss, region_set=self.nde_regions,
                               rm_duplicate=True)  # , asgene_organism=self.organism)

        for rbs in self.rbss:
            l2 = len(self.txp_nde.merged_dict[rbs])
            self.frequency["promoters"]["nde"][rbs] = [l2, len_nde - l2]

        self.stat["DBSs_target_all"] = str(len(self.txp_de))
        self.stat["DBSs_background_all"] = str(len(self.txp_nde))

        ########################################################
        # Count the number of hits on the promoters from each merged DBD

        print("\tCounting frequency of binding sites on DBD...")
        ####################
        # DE
        self.txp_def = RNADNABindingSet("DE")
        self.txp_def.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=True)
        self.txp_def.merge_rbs(rbss=self.rbss, rm_duplicate=True,
                               name_replace=self.de_regions)  # asgene_organism=self.organism
        print("\t\t" + str(len(self.txp_def)) + "\tBinding sites on de promoters")


        # Promoter profiling
        self.promoter = {"de": {},"nde": {}}

        self.promoter["de"]["rd"] = self.txp_def.sort_rd_by_regions(regionset=self.de_regions)
        # self.promoter["de"]["merged_dbs"] = {}
        self.promoter["de"]["dbs"] = {}
        self.promoter["de"]["dbs_coverage"] = {}

        for promoter in self.de_regions:
            dbs = self.promoter["de"]["rd"][promoter.toString()].get_dbs()
            # m_dbs = dbs.merge(w_return=True)
            self.promoter["de"]["dbs"][promoter.toString()] = len(dbs)
            # self.promoter["de"]["merged_dbs"][promoter.toString()] = len(m_dbs)
            self.promoter["de"]["dbs_coverage"][promoter.toString()] = float(dbs.total_coverage()) / len(promoter)

        ######################
        # nDE
        self.promoter["nde"]["dbs"] = {}
        self.promoter["nde"]["dbs_coverage"] = {}

        # ndef_dbs = GenomicRegionSet("nde_dbs")

        self.txp_nde = RNADNABindingSet("non-DE")
        self.txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        self.txp_nde.merge_rbs(rbss=self.rbss, rm_duplicate=True)  # , asgene_organism=self.organism)


        ndef_dbs = self.txp_nde.get_dbs()

        counts = self.nde_regions.counts_per_region(regionset=ndef_dbs)

        coverage = self.nde_regions.coverage_per_region(regionset=ndef_dbs)
        for i, p in enumerate(self.de_regions):
            self.promoter["nde"]["dbs"][p.toString()] = counts[i]
            # self.promoter["nde"]["merged_dbs"][p.toString()] = mcounts[i]
            self.promoter["nde"]["dbs_coverage"][p.toString()] = coverage[i]

        if self.showdbs:
            self.frequency["hits"] = {"de": OrderedDict(), "nde": OrderedDict()}
            numdbs_def = len(self.txp_def.get_dbs(rm_duplicate=True))
            numdbs_ndef = len(self.txp_ndef.get_dbs(rm_duplicate=True))
            for rbs in self.rbss:
                # DE
                l1 = len(self.txp_def.merged_dict[rbs])
                self.frequency["hits"]["de"][rbs] = [l1, numdbs_def - l1]
                # non-DE
                l2 = len(self.txp_ndef.merged_dict[rbs])
                self.frequency["hits"]["nde"][rbs] = [l2, numdbs_ndef - l2]

        if remove_temp:
            os.remove(os.path.join(temp, "de.txp"))

        os.remove(os.path.join(temp, "nde.txp"))

        if obedp:
            # try:
            output = self.de_regions.change_name_by_dict(convert_dict=self.ensembl2symbol)
            output.write_bed(filename=os.path.join(temp, obedp + "_target_promoters.bed"))

            self.txp_de.write_bed(filename=os.path.join(temp, obedp + "_target_promoters_dbs.bed"),
                                  remove_duplicates=False, convert_dict=self.ensembl2symbol)

            self.txp_def.write_bed(filename=os.path.join(temp, obedp + "_dbss.bed"),
                                   remove_duplicates=False, associated=self.organism)

    def fisher_exact(self, alpha):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        pvalues = []

        for rbs in self.frequency["promoters"]["de"]:
            table = numpy.array([self.frequency["promoters"]["de"][rbs], self.frequency["promoters"]["nde"][rbs]])
            # print(table)
            self.oddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
            pvalues.append(p)

        # correction
        if len(self.frequency["promoters"]["de"].keys()) > 1:
            reject, pvals_corrected = multiple_test_correction(pvalues, alpha=alpha, method='indep')
        else:
            pvals_corrected = pvalues
        for i, rbs in enumerate(self.frequency["promoters"]["de"].keys()):
            self.pvalue[rbs] = pvals_corrected[i]
            if pvals_corrected[i] < alpha:
                self.sig_DBD.append(rbs)

        if self.showdbs:
            self.hoddsratio = {}
            self.hpvalue = {}
            pvalues = []
            for rbs in self.frequency["hits"]["de"]:
                table = numpy.array([self.frequency["hits"]["de"][rbs], self.frequency["hits"]["nde"][rbs]])
                self.hoddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
                pvalues.append(p)
            # correction
            if len(self.frequency["hits"]["de"].keys()) > 1:
                reject, pvals_corrected = multiple_test_correction(pvalues, alpha=alpha, method='indep')
            else:
                pvals_corrected = pvalues
            for i, rbs in enumerate(self.frequency["hits"]["de"].keys()):
                self.hpvalue[rbs] = pvals_corrected[i]
                if pvals_corrected[i] < alpha:
                    self.sig_DBD.append(rbs)

    def dbd_regions(self, output):
        dbd_regions(exons=self.rna_regions, sig_region=self.sig_DBD,
                    rna_name=self.rna_name, output=output)
        self.stat["DBD_all"] = str(len(self.rbss))
        self.stat["DBD_sig"] = str(len(self.sig_DBD))

        sigDBD = GenomicRegionSet("DBD_sig")
        sigDBD.sequences = self.sig_DBD
        rbss = self.txp_de.get_rbs()
        overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
        self.stat["DBSs_target_DBD_sig"] = str(len(overlaps))
        rbss = self.txp_nde.get_rbs()
        overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
        self.stat["DBSs_background_DBD_sig"] = str(len(overlaps))

    def plot_lines(self, txp, rna, dirp, cut_off, log, ylabel, linelabel, filename, sig_region, ac=None, showpa=False):
        """Generate the plots for demonstration of RBS

        rna:     File path to the RNA sequence in FASTA format
        dir:     Directory where the figure is saved
        ac:      RNA accessibility data.
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(rna)
        if not self.rna_name: self.rna_name = rnas[0].name

        self.rna_len = rnas.total_len()

        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dirp=dirp, sig_region=sig_region,
                 cut_off=cut_off, log=log, ylabel=ylabel, linelabel=linelabel,
                 filename=filename, ac=ac, showpa=showpa, exons=self.rna_regions)

    def barplot(self, dirp, filename, sig_region, dbs=False):
        """Generate the barplot to show the difference between target promoters and non-target promoters"""

        def to_percent(y, position):
            # Ignore the passed in position. This has the effect of scaling the default
            # tick locations.
            s = str(100 * y)
            return s + '%'

        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6, 4))
        ind = range(len(self.rbss))
        width = 0.35

        if not dbs:
            propor_de = [float(b[0]) / len(self.de_regions) for b in self.frequency["promoters"]["de"].values()]
            propor_nde = [float(b[0]) / len(self.nde_regions) for b in self.frequency["promoters"]["nde"].values()]
        else:
            propor_de = [float(b[0]) / (b[0] + b[1]) for b in self.frequency["hits"]["de"].values()]
            propor_nde = [float(b[0]) / (b[0] + b[1]) for b in self.frequency["hits"]["nde"].values()]

        max_y = max([max(propor_de), max(propor_nde)]) * 1.2

        # Plotting
        rect = patches.Rectangle(xy=(1, 0), width=0.8, height=max_y, facecolor=sig_color,
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        for i, rbs in enumerate(self.rbss):
            if rbs in sig_region:
                rect = patches.Rectangle(xy=(i + 0.05, 0), width=0.9, height=max_y, facecolor=sig_color,
                                         edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
                ax.add_patch(rect)

        rects_de = ax.bar([i + 0.15 for i in ind], propor_de, width, color=target_color,
                          edgecolor="none", label="Target promoters")
        rects_nde = ax.bar([i + 0.15 + width for i in ind], propor_nde, width, color=nontarget_color,
                           edgecolor="none", label="Non-target promoters")

        # Legend
        tr_legend, = plt.plot([1, 1], color=target_color, linewidth=6, alpha=1)
        ntr_legend, = plt.plot([1, 1], color=nontarget_color, linewidth=6, alpha=1)
        ax.legend([tr_legend, ntr_legend, rect], ["Target promoters", "Non-target promoters", "Significant DBD"],
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  prop={'size': 9}, ncol=3)
        tr_legend.set_visible(False)
        ntr_legend.set_visible(False)

        tick_size = 8
        # Y axis
        ax.set_ylim([0, max_y])
        formatter = FuncFormatter(to_percent)
        # Set the formatter
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9)
        ax.set_ylabel("Proportion of promoters (%)", fontsize=9, rotation=90)

        # X axis
        ax.set_xlim([0, len(self.rbss)])
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax.set_xticks([i + 0.5 for i in range(len(self.rbss))])
        ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.rbss], rotation=35,
                           ha="right", fontsize=tick_size)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

        for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9)
        ax.set_xlabel(self.rna_name + " DNA Binding Domains", fontsize=9)

        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dirp, filename), facecolor='w', edgecolor='w',
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        pp = PdfPages(os.path.splitext(os.path.join(dirp, filename))[0] + '.pdf')
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def gen_html(self, directory, parameters, ccf, align=50, alpha=0.05):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = "Promoter Test: " + dir_name
        self.link_d = OrderedDict()
        self.link_d["RNA"] = "index.html"
        self.link_d["All promoters"] = "promoters.html"
        self.link_d["Sig promoters"] = "spromoters.html"
        self.link_d["Parameters"] = "parameters.html"

        if self.organism == "hg19":
            self.ani = "human"
        elif self.organism == "hg38":
            self.ani = "human"
        elif self.organism == "mm9":
            self.ani = "mouse"
        else:
            self.ani = None

        #############################################################
        # Index main page
        #############################################################
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_figure("plot_promoter.png", align="left", width="45%", more_images=["bar_promoter.png"])

        if self.showdbs:
            html.add_figure("plot_dbss.png", align="left", width="45%", more_images=["bar_dbss.png"])

        # Table of merged TBS on promoters
        if self.showdbs:
            header_list = [["", "", "Promoters", None, None, None, None, None, "DBSs", None, None, None, None, None],
                           ["#", "DBD",
                            "Target Promoter", None, "Non-target Promoter", None, "Statistics", None,
                            "Target Promoter", None, "Non-target Promoter", None, "Statistics", None],
                           [" ", " ",
                            "with DBS", "without DBS", "with DBS", "without DBS", "OR", "<i>p</i>-value",
                            "No. DBSs", "Other DBSs", "No. DBSs", "Other DBSs", "OR", "<i>p</i>-value"]]
            header_titles = [["", "", "Statistics on promoter level", None, None, None, None, None,
                              "Statistics on DBS level", None, None, None, None, None],
                             ["Rank of the talbe",
                              "DNA Binding Domain which is the functional region on RNA.",
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on promoters", None,
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on DBSs", None],
                             ["",
                              "",
                              "Number of target promoters which contain DBSs (DNA Binding Sites).",
                              "Number of target promoters which don't contain DBSs (DNA Binding Sites).",
                              "Number of non-target promoters which contain DBSs (DNA Binding Sites).",
                              "Number of non-target promoters which don't contain DBSs (DNA Binding Sites).",
                              "Odds Ratio", "P-value",
                              "Number of DBSs found in the target promoters.",
                              "Number of DBSs not found in the target promoters.",
                              "Number of DBSs found in the non-target promoters.",
                              "Number of DBSs not found in the non-target promoters.",
                              "Odds Ratio", "P-value"]
                             ]
            border_list = [" style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:2pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\""]
        else:
            header_list = [["#", "DBD", "Target Promoter", None, "Non-target Promoter", None, "Statistics", None],
                           [" ", " ", "with DBS", "without DBS", "with DBS", "without DBS", "OR", "<i>p</i>"]]
            header_titles = [["Rank of the talbe",
                              "DNA Binding Domain which is the functional region on RNA.",
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on promoters", None],
                             ["",
                              "",
                              "Number of target promoters which contain DBSs (DNA Binding Sites).",
                              "Number of target promoters which don't contain DBSs (DNA Binding Sites).",
                              "Number of non-target promoters which contain DBSs (DNA Binding Sites).",
                              "Number of non-target promoters which don't contain DBSs (DNA Binding Sites).",
                              "Odds Ratio", "P-value"]
                             ]
            border_list = ["style=\"border-right:1pt solid gray\"",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\""]

        type_list = 'ssssssssssssssssssss'
        col_size_list = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
        data_table = []
        rank = 0
        self.topDBD = ["-", 1]

        for rbs in self.frequency["promoters"]["de"].keys():
            if self.frequency["promoters"]["de"][rbs][0] < ccf: continue
            rank += 1
            if self.pvalue[rbs] < alpha:
                p_promoter = "<font color=\"red\">" + value2str(self.pvalue[rbs]) + "</font>"
            else:
                p_promoter = value2str(self.pvalue[rbs])

            if self.showdbs:
                if self.hpvalue[rbs] < alpha:
                    p_hit = "<font color=\"red\">" + value2str(self.hpvalue[rbs]) + "</font>"
                else:
                    p_hit = value2str(self.hpvalue[rbs])

            try:
                if self.pvalue[rbs] < self.topDBD[1]:
                    self.topDBD = [rbs.str_rna(pa=False), self.pvalue[rbs]]
            except:
                self.topDBD = [rbs.str_rna(pa=False), self.pvalue[rbs]]

            new_row = [str(rank),
                       rbs.str_rna(pa=False),
                       '<a href="' + "dbds_promoters.html#" + rbs.str_rna() +
                       '" style="text-align:left">' +
                       value2str(self.frequency["promoters"]["de"][rbs][0]) + '</a>',
                       value2str(self.frequency["promoters"]["de"][rbs][1]),
                       value2str(self.frequency["promoters"]["nde"][rbs][0]),
                       value2str(self.frequency["promoters"]["nde"][rbs][1]),
                       value2str(self.oddsratio[rbs]),
                       p_promoter]
            if self.showdbs:
                new_row += [value2str(self.frequency["hits"]["de"][rbs][0]),
                            value2str(self.frequency["hits"]["de"][rbs][1]),
                            value2str(self.frequency["hits"]["nde"][rbs][0]),
                            value2str(self.frequency["hits"]["nde"][rbs][1]),
                            value2str(self.hoddsratio[rbs]),
                            p_hit]

            data_table.append(new_row)

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titles, border_list=border_list, sortable=True)

        html.add_heading("Notes")
        html.add_list(["DBD stands for functional DNA Binding Domain on RNA.",
                       "RBS stands for RNA Binding Site on RNA.",
                       "DBS stands for DNA Binding Site on DNA."])

        #### Table for motif
        if self.motif:
            data_table = []
            header_list = ["DBD", "Motif"]
            for dbd, f in self.motif.iteritems():
                svg = '<object type="image/svg+xml" data="' + f + '">Your browser does not support SVG</object>'
                new_row = [dbd, svg]
                data_table.append(new_row)
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        ####
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "index.html"))

        #############################################################
        # RNA subpage: Profile of targeted promoters for each merged DNA Binding Domain
        #############################################################

        # dbds_promoters.html

        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if self.scores:
            if self.scoreh and self.de_gene:
                score_header = [self.de_gene.cond[0]]
            else:
                if "(" in self.scores[0]:
                    scs = self.scores[0].replace("(", "")
                    scs = scs.replace(")", "")
                    scs = scs.split(",")
                    # score_header = ["Score"] * len(scs)
                    score_header = ["Fold_change", "Filtered"]

                else:
                    score_header = ["Fold Change Score"]

            header_list = ["#", "Promoter", "Gene", "DBSs counts", "DBS coverage"]
            header_list += score_header
            header_list += ["Sum of Ranks"]
            header_titles = ["", "Target promoters", "Gene symbol",
                             "Number of DNA Binding sites locating within the promoter",
                             "The proportion of promoter covered by binding sites"]
            header_titles += [
                                 "Scores loaded by their absolute values from gene list or BED input. If there is annotation error for the gene names, it shows zero instead."] * len(
                score_header)
            header_titles += ["Sum up the ranks from left-hand side columns"]

        else:
            header_list = ["#", "Promoter", "Gene", "DBSs Count", "DBS coverage", "Sum of Ranks"]

            header_titles = ["", "Target promoters", "Gene symbol",
                             "Number of DNA Binding sites locating within the promoter",
                             "The proportion of promoter covered by binding sites",
                             "Sum up the ranks from left-hand side columns"]

        for rbsm in self.frequency["promoters"]["de"].keys():
            html.add_heading("DNA Binding Domain: " + rbsm.str_rna(), idtag=rbsm.str_rna())
            data_table = []
            # Calculate the ranking
            # rank_array
            rank_count = len(self.txp_de.merged_dict[rbsm]) - rank_array(
                [len(self.promoter["de"]["rd"][p.toString()]) for p in self.txp_de.merged_dict[rbsm]])
            rank_coverage = len(self.txp_de.merged_dict[rbsm]) - rank_array(
                [self.promoter["de"]["dbs_coverage"][p.toString()] for p in self.txp_de.merged_dict[rbsm]])

            if self.scores:
                multiple_scores = False
                if isinstance(self.scores[0], str):
                    if "(" in self.scores[0]:
                        def ranking(scores):
                            rank_score = len(self.txp_de.merged_dict[rbsm]) - rank_array(scores)
                            return rank_score

                        multiple_scores = True
                        scs = []
                        for i, promoter in enumerate(self.txp_de.merged_dict[rbsm]):
                            s = self.de_gene.values[promoter.name.upper()].replace("(", "")
                            s = s.replace(")", "")
                            s = s.split(",")
                            scs.append([float(ss) for ss in s])
                        ar = numpy.array(scs)
                        # ar.transpose()
                        # print(ar)
                        score_ar = ar.tolist()
                        rank_score = numpy.apply_along_axis(ranking, axis=0, arr=ar)
                        rank_score = rank_score.transpose()
                        rank_sum = numpy.sum(rank_score, axis=0).tolist()
                        rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_sum)]
                        rank_score = rank_score.tolist()

                    else:
                        try:
                            new_scores = []
                            for i, promoter in enumerate(self.txp_de.merged_dict[rbsm]):
                                try:
                                    try:
                                        s = self.de_gene.values[promoter.name]
                                    except:
                                        s = self.de_gene.values[self.ensembl2symbol[promoter.name].upper()]
                                except:
                                    s = new_scores.append(abs(float(s)))
                                if s == "Inf" or s == "inf":
                                    new_scores.append(float("inf"))
                                elif s == "-Inf" or s == "-inf":
                                    new_scores.append(-float("inf"))
                                else:
                                    new_scores.append(abs(float(s)))

                            scores = new_scores
                            rank_score = len(self.txp_de.merged_dict[rbsm]) - rank_array(scores)
                            rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]
                        except:
                            scores = [float(i) for i in self.scores]
                            rank_score = len(self.txp_de.merged_dict[rbsm]) - rank_array(scores)
                            rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

                else:
                    scores = [float(self.de_gene.values[p.name.upper()]) for p in self.txp_de.merged_dict[rbsm]]
                    rank_score = len(self.txp_de.merged_dict[rbsm]) - rank_array(scores)
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

            else:
                rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

            for i, promoter in enumerate(self.txp_de.merged_dict[rbsm]):
                # Add information
                if self.ani:
                    pr = '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism + "&position=" + promoter.chrom + "%3A" + str(
                        promoter.initial) + "-" + str(
                        promoter.final) + '" style="text-align:left">' + promoter.toString(space=True) + '</a>'

                elif self.organism == "tair10":
                    pr = "".join(['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=',
                                  promoter.chrom, "%3A", str(promoter.initial), "..", str(promoter.final),
                                  '" target="_blank">', promoter.toString(space=True), '</a>'])
                else:
                    pr = promoter.toString(space=True)

                try:
                    newline = [str(i + 1), pr,
                               split_gene_name(gene_name=self.ensembl2symbol[promoter.name], org=self.organism),
                               str(len(self.promoter["de"]["rd"][promoter.toString()])),
                               value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                               ]
                except:
                    newline = [str(i + 1), pr,
                               split_gene_name(gene_name=promoter.name, org=self.organism),
                               str(len(self.promoter["de"]["rd"][promoter.toString()])),
                               value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                               ]
                if self.scores:
                    if multiple_scores:

                        for j in range(len(score_ar[0])):

                            newline += [value2str(score_ar[i][j])]
                    else:
                        newline += [value2str(scores[i])]

                newline += ["<i>" + str(int(rank_sum[i])) + "</i>"]
                # print(newline)
                data_table.append(newline)
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_titles, sortable=True, border_list=None)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "dbds_promoters.html"))

        ################################################################
        ############# Parameters
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_heading("Parameters")
        header_list = ["Description", "Arguments", "Value"]

        if parameters.de:
            de = os.path.basename(parameters.de)
            bed = "False"
            bg = "False"
        else:
            de = "False"
            bed = os.path.basename(parameters.bed)
            bg = os.path.basename(parameters.bg)

        data_table = [["RNA sequence name", "-rn", parameters.rn],
                      ["Input RNA sequence file", "-r", os.path.basename(parameters.r)],
                      ["Input file for defferentially expression gene list", "-de", de],
                      ["Input BED file as promoters", "-bed", bed],
                      ["Input BED file as backgrounds", "-bg", bg],
                      ["Output directory", "-o", "/".join(parameters.o.partition("/")[-3:])],
                      ["Organism", "-organism", parameters.organism],
                      ["filter_havana", "-filter_havana", parameters.filter_havana],
                      ["protein_coding", "-protein_coding", parameters.protein_coding],
                      ["known_only", "-known_only", parameters.known_only],
                      ["Promoter length", "-pl", str(parameters.pl)],
                      ["Alpha level for rejection p value", "-a", str(parameters.a)],
                      ["Cut off value for filtering out the low counts of DBSs", "-ccf", str(parameters.ccf)],
                      ["Remove temporary files", "-rt", str(parameters.rt)],
                      ["Input file for RNA accecibility", "-ac", str(parameters.ac)],
                      ["Cut off value for RNA accecibility", "-accf", str(parameters.accf)],
                      ["Output the BED files for DNA binding sites.", "-obed", str(parameters.obed)],
                      ["Show parallel and antiparallel bindings in the plot separately.", "-showpa",
                       str(parameters.showpa)],
                      ["Minimum length", "-l", str(self.triplexator_p[0])],
                      ["Maximum error rate", "-e", str(self.triplexator_p[1])],
                      ["Tolerated number of consecutive errors", "-c", str(self.triplexator_p[2])],
                      ["Filtering repeats", "-fr", str(self.triplexator_p[3])],
                      ["Filtering mode", "-fm", str(self.triplexator_p[4])],
                      ["Output format", "-of", str(self.triplexator_p[5])],
                      ["Merge features", "-mf", str(self.triplexator_p[6])]]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True)
        html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory, "parameters.html"))

    def gen_html_genes(self, directory, align=50, alpha=0.05, nonDE=False):
        dir_name = os.path.basename(directory)
        html_header = "Promoter Test: " + dir_name
        self.ranktable = {}
        self.dbstable = {}
        type_list = 'sssssssssssssss'
        col_size_list = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

        #############################################################
        # Promoter centered
        #############################################################

        ##############################################################################################
        # promoters.html
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if self.scores:
            if self.scoreh and self.de_gene:
                score_header = [self.de_gene.cond[0]]
            else:
                if "(" in self.scores[0]:
                    scs = self.scores[0].replace("(", "")
                    scs = scs.replace(")", "")
                    scs = scs.split(",")
                    # score_header = ["Score"] * len(scs)
                    score_header = ["Fold_change", "Filtered"]

                else:
                    score_header = ["Fold Change Score"]
            header_listp = ["#", "Promoter", "Gene", "DBSs Count", "DBS coverage"]
            header_listp += score_header
            header_listp += ["Sum of Ranks"]

            header_titlesp = ["", "Target promoters", "Gene symbol",
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites"]
            header_titlesp += [
                                  "Scores loaded by their absolute values from gene list or BED input. If there is annotation error for the gene names, it shows zero instead."] * len(
                score_header)
            header_titlesp += ["Sum up the ranks from left-hand side columns"]

        else:
            header_listp = ["#", "Promoter", "Gene", "DBSs Count", "DBS coverage", "Sum of Ranks"]

            header_titlesp = ["", "Target promoters", "Gene symbol",
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites",
                              "Sum up the ranks from left-hand side columns"]

        html.add_heading("Target promoters")
        data_table = []

        if not self.de_regions.sorted: self.de_regions.sort()
        # Iterate by each gene promoter

        # Calculate the ranking
        rank_count = len(self.de_regions) - rank_array(
            [self.promoter["de"]["dbs"][p.toString()] for p in self.de_regions])
        rank_coverage = len(self.de_regions) - rank_array(
            [self.promoter["de"]["dbs_coverage"][p.toString()] for p in self.de_regions])

        if self.scores:
            multiple_scores = False
            if isinstance(self.scores[0], str):
                if "(" in self.scores[0]:
                    def ranking(scores):
                        rank_score = len(self.de_regions) - rank_array(scores)
                        return rank_score

                    multiple_scores = True
                    scs = []
                    for s in self.scores:
                        s = s.replace("(", "")
                        s = s.replace(")", "")
                        s = s.split(",")
                        scs.append([float(ss) for ss in s])
                    ar = numpy.array(scs)
                    # ar.transpose()
                    # print(ar)
                    score_ar = ar.tolist()
                    rank_score = numpy.apply_along_axis(ranking, axis=0, arr=ar)
                    rank_score = rank_score.transpose()
                    rank_sum = numpy.sum(rank_score, axis=0).tolist()
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_sum)]
                    rank_score = rank_score.tolist()

                else:
                    new_scores = []
                    for s in self.scores:
                        if s == "Inf" or s == "inf":
                            new_scores.append(float("inf"))
                        elif s == "-Inf" or s == "-inf":
                            new_scores.append(-float("inf"))
                        else:
                            new_scores.append(abs(float(s)))

                    scores = new_scores
                    rank_score = len(self.de_regions) - rank_array(scores)
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

            else:
                rank_score = len(self.de_regions) - rank_array(scores)
                rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

        for i, promoter in enumerate(self.de_regions):

            if self.promoter["de"]["dbs"][promoter.toString()] == 0:
                dbssount = str(0)
            else:
                dbssount = '<a href="promoters_dbds.html#' + promoter.toString() + '" style="text-align:left">' + str(
                    self.promoter["de"]["dbs"][promoter.toString()]) + '</a>'

            if self.ani:

                region_link = "".join(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=', self.organism,
                                       "&position=", promoter.chrom, "%3A", str(promoter.initial), "-",
                                       str(promoter.final), '" style="text-align:left" target="_blank">',
                                       promoter.toString(space=True), '</a>'])
            else:
                if self.organism == "tair10":
                    region_link = "".join(
                        ['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=',
                         promoter.chrom, "%3A", str(promoter.initial), "..", str(promoter.final),
                         '" target="_blank">',
                         promoter.toString(space=True), '</a>'])
                else:
                    region_link = promoter.toString(space=True)

            try:
                gn = self.ensembl2symbol[promoter.name]
                if not gn: gn = promoter.name
            except:
                gn = promoter.name

            self.ranktable[gn] = str(int(rank_sum[i]))
            self.dbstable[gn] = str(int(self.promoter["de"]["dbs"][promoter.toString()]))

            newline = [str(i + 1),
                       region_link,
                       split_gene_name(gene_name=gn, org=self.organism),
                       dbssount,
                       value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                       ]

            if self.scores:
                if multiple_scores:
                    # print(score_ar)
                    # print(len(score_ar))
                    for j in range(len(score_ar[0])):
                        # print(score_ar[j][i])
                        # print(score_ar[i][j])
                        newline += [value2str(score_ar[i][j])]
                else:
                    newline += [value2str(scores[i])]

            newline += ["<i>" + str(int(rank_sum[i])) + "</i>"]
            # print(newline)
            data_table.append(newline)

        # print(data_table)
        html.add_zebra_table(header_listp, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titlesp, border_list=None, sortable=True)
        html.add_heading("Notes")
        html.add_list(["DBS stands for DNA Binding Site on DNA.",
                       "DBS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "promoters.html"))

        ############################
        # Subpages for promoter centered page
        # promoters_dbds.html
        header_sub = ["#", "RBS", "DBS", "Strand", "Score", "Motif", "Orientation"]
        header_titles = ["", "RNA Binding Site", "DNA Binding Site", "Strand of DBS on DNA",
                         "Score of binding event", "Motif of binding by triple helix rule",
                         "Orientation of interaction between DNA and RNA. 'P'- Parallel; 'A'-Antiparallel"]
        header_list = header_sub
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for i, promoter in enumerate(self.de_regions):
            if self.promoter["de"]["dbs"][promoter.toString()] == 0:
                continue
            else:
                try:
                    gn = self.ensembl2symbol[promoter.name]
                except:
                    gn = promoter.name
                html.add_heading(split_gene_name(gene_name=gn, org=self.organism), idtag=promoter.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                                       "&position=" + promoter.chrom + "%3A" + str(promoter.initial) + "-" + str(
                    promoter.final) +
                                       '" style="margin-left:50">' +
                                       promoter.toString(space=True) + '</a>'])
                data_table = []

                for j, rd in enumerate(self.promoter["de"]["rd"][promoter.toString()]):
                    rbs = rd.rna.str_rna(pa=False)
                    for rbsm in self.sig_DBD:
                        # rbsm = rbsm.partition(":")[2].split("-")
                        if rd.rna.overlap(rbsm):
                            rbs = "<font color=\"red\">" + rd.rna.str_rna(pa=False) + "</font>"

                    data_table.append([str(j + 1),
                                       rbs,
                                       rd.dna.toString(space=True),
                                       rd.dna.orientation,
                                       rd.score,
                                       rd.motif, rd.orient])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     header_titles=header_titles, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "promoters_dbds.html"))

        ############################
        # Subpages for promoter centered page
        # spromoters_dbds.html
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        sig_promoter_count = {}
        for i, promoter in enumerate(self.de_regions):
            if self.promoter["de"]["dbs"][promoter.toString()] == 0:
                continue
            else:
                try:
                    gn = self.ensembl2symbol[promoter.name]
                except:
                    gn = promoter.name
                html.add_heading(split_gene_name(gene_name=gn, org=self.organism), idtag=promoter.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                                       "&position=" + promoter.chrom + "%3A" + str(promoter.initial) + "-" + str(
                    promoter.final) +
                                       '" style="margin-left:50">' +
                                       promoter.toString(space=True) + '</a>'])
                data_table = []

                for j, rd in enumerate(self.promoter["de"]["rd"][promoter.toString()]):
                    overlapping = False
                    for rbsm in self.sig_DBD:
                        if rd.rna.overlap(rbsm): overlapping = True
                    if overlapping:
                        data_table.append([str(j + 1),
                                           rd.rna.str_rna(pa=False),
                                           rd.dna.toString(space=True),
                                           rd.dna.orientation,
                                           rd.score,
                                           rd.motif, rd.orient])
                if overlapping:
                    sig_promoter_count[promoter] = len(data_table)
                    html.add_zebra_table(header_list, col_size_list, type_list, data_table,
                                         align=align, cell_align="left",
                                         header_titles=header_titles, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "spromoters_dbds.html"))

        ##############################################################################################
        # spromoters.html    for significant promoters
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        # Select promoters in sig DBD
        spromoters = sig_promoter_count.keys()
        if len(spromoters) == 0:
            html.add_heading("There is no significant DBD.")

        else:
            html.add_heading("Target promoters binded by significant DBD")
            # for rbsm in self.sig_region_promoter:
            #     spromoters = spromoters + [p for p in self.txp_de.merged_dict[rbsm]]
            # spromoters = list(set(spromoters))
            data_table = []

            # Iterate by each gene promoter

            # Calculate the ranking
            rank_count = len(spromoters) - rank_array([sig_promoter_count[p] for p in spromoters])
            rank_coverage = len(spromoters) - rank_array(
                [self.promoter["de"]["dbs_coverage"][p.toString()] for p in spromoters])

            if self.scores:
                multiple_scores = False
                sscores = []
                # de_genes_str = [g.name for g in self.de_gene.genes]
                for p in spromoters:

                    try:
                        gene_sym = self.ensembl2symbol[p.name].upper()
                    except:
                        gene_sym = p.name.upper()
                    try:
                        sscores.append(self.de_gene.values[gene_sym.upper()])
                    except:
                        sscores.append(0)

                if isinstance(sscores[0], str):
                    if "(" in sscores[0]:
                        def ranking(scores):
                            rank_score = len(spromoters) - rank_array(scores)
                            return rank_score

                        multiple_scores = True
                        scs = []

                        for s in sscores:
                            s = s.replace("(", "")
                            s = s.replace(")", "")
                            s = s.split(",")
                            scs.append([float(ss) for ss in s])
                        ar = numpy.array(scs)
                        # ar.transpose()
                        # print(ar)
                        score_ar = ar.tolist()
                        rank_score = numpy.apply_along_axis(ranking, axis=0, arr=ar)
                        rank_score = rank_score.transpose()
                        rank_sum = numpy.sum(rank_score, axis=0).tolist()
                        rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_sum)]
                        rank_score = rank_score.tolist()

                    else:
                        new_scores = []
                        for s in sscores:
                            if s == "Inf" or s == "inf":
                                new_scores.append(float("inf"))
                            elif s == "-Inf" or s == "-inf":
                                new_scores.append(-float("inf"))
                            elif isinstance(s, str):
                                new_scores.append(s)
                            else:
                                new_scores.append(abs(float(s)))

                        scores = new_scores
                        rank_score = len(spromoters) - rank_array(scores)
                        rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

                else:
                    rank_score = len(spromoters) - rank_array(sscores)
                    rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

            else:
                rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

            for i, promoter in enumerate(spromoters):

                dbssount = '<a href="spromoters_dbds.html#' + promoter.toString() + '" style="text-align:left">' + str(
                    sig_promoter_count[promoter]) + '</a>'

                if self.ani:
                    region_link = "".join(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=', self.organism,
                                           "&position=", promoter.chrom, "%3A", str(promoter.initial), "-",
                                           str(promoter.final), '" style="text-align:left" target="_blank">',
                                           promoter.toString(space=True), '</a>'])
                else:
                    if self.organism == "tair10":
                        region_link = "".join(
                            ['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=',
                             promoter.chrom, "%3A", str(promoter.initial), "..", str(promoter.final),
                             '" target="_blank">',
                             promoter.toString(space=True), '</a>'])
                    else:
                        region_link = promoter.toString(space=True)

                try:
                    gn = self.ensembl2symbol[promoter.name]
                    if not gn: gn = promoter.name
                except:
                    gn = promoter.name

                self.ranktable[gn] = str(int(rank_sum[i]))
                self.dbstable[gn] = str(sig_promoter_count[promoter])

                newline = [str(i + 1),
                           region_link,
                           split_gene_name(gene_name=gn, org=self.organism),
                           dbssount,
                           value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                           ]

                if self.scores:
                    if multiple_scores:
                        # print(score_ar)
                        # print(len(score_ar))
                        for j in range(len(score_ar[0])):
                            # print(score_ar[j][i])
                            # print(score_ar[i][j])
                            newline += [value2str(score_ar[i][j])]
                    else:
                        newline += [value2str(scores[i])]

                newline += ["<i>" + str(int(rank_sum[i])) + "</i>"]
                # print(newline)
                data_table.append(newline)

            # print(data_table)
            html.add_zebra_table(header_listp, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_titlesp, border_list=None, sortable=True)
            html.add_heading("Notes")
            html.add_list(["DBS stands for DNA Binding Site on DNA.",
                           "DBS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
            html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "spromoters.html"))

    def save_profile(self, output, bed, geneset):
        """Save statistics for comparison with other results"""
        pro_path = os.path.join(os.path.dirname(output), "profile.txt")
        exp = os.path.basename(output)
        # tag = os.path.basename(os.path.dirname(rnafile))
        if geneset:
            tar_reg = os.path.basename(geneset)
        else:
            tar_reg = os.path.basename(bed)
        # RNA name with region
        if self.rna_regions:
            trans = "Transcript:"
            for r in self.rna_regions:
                trans += " " + r[0] + ":" + str(r[1]) + "-" + str(r[2])
            rna = '<p title="' + trans + '">' + self.rna_name + "<p>"
        else:
            rna = self.rna_name
        # RNA associated genes
        r_genes = rna_associated_gene(rna_regions=self.rna_regions, name=self.rna_name, organism=self.organism)

        newlines = []
        # try:
        if os.path.isfile(pro_path):
            with open(pro_path, 'r') as f:
                new_exp = True
                for line in f:
                    line = line.strip()
                    line = line.split("\t")
                    if line[0] == exp:
                        newlines.append([exp, rna, output.split("_")[-1],
                                         self.organism, tar_reg, str(len(self.sig_DBD)),
                                         self.topDBD[0], value2str(self.topDBD[1]), r_genes])
                        new_exp = False
                    elif line[0] == "Experiment":
                        continue
                    else:
                        newlines.append(line)
                if new_exp:
                    newlines.append([exp, rna, output.split("_")[-1],
                                     self.organism, tar_reg, str(len(self.sig_DBD)),
                                     self.topDBD[0], value2str(self.topDBD[1]), r_genes])
        else:
            newlines.append([exp, rna, output.split("_")[-1],
                             self.organism, tar_reg, str(len(self.sig_DBD)),
                             self.topDBD[0], value2str(self.topDBD[1]), r_genes])

        newlines.sort(key=lambda x: float(x[7]))

        newlines = [["Experiment", "RNA_names", "Tag", "Organism", "Target_region", "No_sig_DBDs",
                     "Top_DBD", "p-value", "closest_genes"]] + newlines

        with open(pro_path, 'w') as f:
            for lines in newlines:
                print("\t".join(lines), file=f)

    def save_table(self, path, table, filename):
        """Save the summary rank into the table for heatmap"""
        "lncRNA_target_ranktable.txt"

        table_path = os.path.join(path, filename)
        # self.ranktable = {}
        rank_table = []

        if os.path.isfile(table_path):
            # Table exists
            f = open(table_path)
            for i, line in enumerate(f):
                line = line.strip().split()
                if i == 0:
                    if not line:
                        break
                        exist = False
                    if self.rna_name in line:
                        # lncRNA exists
                        exist = True
                        ind_rna = line.index(self.rna_name)
                        header = line
                    else:
                        exist = False
                        line.append(self.rna_name)
                        header = line

                else:
                    if exist and ind_rna:
                        try:
                            line[ind_rna + 1] = table[line[0]]
                            rank_table.append(line)
                        except:
                            rank_table.append(line)
                    else:
                        try:
                            line.append(table[line[0]])
                            rank_table.append(line)
                        except:
                            rank_table.append(line)
            f.close()

        else:
            # Table not exists
            header = ["gene", self.rna_name]
            for k, v in table.iteritems():
                rank_table.append([k, v])

        # Write into file
        g = open(table_path, "w")

        print("\t".join(header), file=g)
        for l in rank_table:
            # print(l)
            print("\t".join(l), file=g)
        g.close()

    def heatmap(self, table, temp):
        """Generate heatmap for comparisons among genes and lncRNAs"""

        # Generate random features and distance matrix.
        genes = []
        data = []
        if not os.path.isfile(os.path.join(temp, table)):
            print("*** No rank table is found.")
            return

        f = open(os.path.join(temp, table))

        for i, line in enumerate(f):
            line = line.strip().split()
            if i == 0:
                rnas = line[1:]
            else:
                genes.append(line[0])
                data.append([int(v) for v in line[1:]])

        f.close()
        data = numpy.array(data)
        # print(type(data))
        # print(data.shape)

        # Compute and plot first dendrogram.
        fig = plt.figure(figsize=(max(20, 4 + int(data.shape[1] * 1.6)),
                                  max(80, int(data.shape[0] * 0.2))))
        # fig.suptitle("Heatmap of summary ranks", fontsize=20, y=0.95)

        ax1 = fig.add_axes([0.09, 0.05, 0.2, 0.9])
        Y = sch.linkage(data, method='single')
        Z1 = sch.dendrogram(Y, orientation='right')
        ax1.set_xticks([])
        ax1.set_yticks([])

        # Compute and plot second dendrogram.
        ax2 = fig.add_axes([0.3, 1.1, 0.55, 0.04])
        # tdata = numpy.transpose(data)
        Y = sch.linkage(data.T, method='single')
        Z2 = sch.dendrogram(Y)
        ax2.set_xticks([])
        ax2.set_yticks([])

        # Plot distance matrix.
        axmatrix = fig.add_axes([0.3, 0.05, 0.55, 0.9])
        # axmatrix = fig.add_axes()
        idx1 = Z1['leaves']
        idx2 = Z2['leaves']
        data = data[idx1, :]
        data = data[:, idx2]
        im = axmatrix.matshow(data, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)

        axmatrix.set_xticks(range(data.shape[1]))
        axmatrix.set_xticklabels([rnas[i] for i in idx2], minor=False, ha="left")
        axmatrix.xaxis.set_label_position('top')
        axmatrix.xaxis.tick_top()
        plt.xticks(rotation=70, fontsize=12)

        axmatrix.set_yticks(range(data.shape[0]))
        axmatrix.set_yticklabels([genes[i] for i in idx1], minor=False)
        axmatrix.yaxis.set_label_position('right')
        axmatrix.yaxis.tick_right()
        plt.yticks(rotation=0, fontsize=12)

        # Plot colorbar.
        axcolor = fig.add_axes([0.1, 0.02, 0.8, 0.01])
        plt.colorbar(im, cax=axcolor, orientation='horizontal')
        axcolor.set_xlabel('summary rank')
        fig.savefig(os.path.join(temp, 'lncRNA_target_dendrogram.png'))
        fig.savefig(os.path.join(temp, 'lncRNA_target_dendrogram.pdf'), format="pdf")

