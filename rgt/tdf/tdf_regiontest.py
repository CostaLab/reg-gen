# Python Libraries
from __future__ import print_function
import os
import multiprocessing
from collections import *

# Local Libraries
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

# Distal Libraries
from rgt.SequenceSet import SequenceSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import SequenceType, Html, GenomeData
from triplexTools import get_rna_region_str, connect_rna,\
    dbd_regions, lineplot, value2str, rank_array,\
    split_gene_name, rna_associated_gene, find_triplex, random_each

# Color code for all analysis
target_color = "mediumblue"
nontarget_color = "darkgrey"
sig_color = "powderblue"

class RandomTest:
    def __init__(self, rna_fasta, rna_name, dna_region, organism, showdbs=False):
        self.organism = organism
        genome = GenomeData(organism)
        self.genome_path = genome.get_genome()
        # RNA: Path to the FASTA file
        self.rna_fasta = rna_fasta
        self.showdbs = showdbs

        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(self.rna_fasta)
        if rna_name:
            self.rna_name = rna_name
        else:
            self.rna_name = rnas[0].name

        # DNA: GenomicRegionSet
        self.dna_region = GenomicRegionSet(name="target")
        self.dna_region.read_bed(dna_region)
        self.dna_region = self.dna_region.gene_association(organism=self.organism)
        self.topDBD = []

    def get_rna_region_str(self, rna):
        """Getting the rna region from the information header with the pattern:
                REGION_chr3_51978050_51983935_-_"""

        self.rna_regions = get_rna_region_str(rna)

    def connect_rna(self, rna, temp):
        connect_rna(rna, temp, self.rna_name)

        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(os.path.join(temp, "rna_temp.fa"))
        if not self.rna_name: self.rna_name = rnas[0].name
        self.rna_len = rnas.total_len()

    def target_dna(self, temp, remove_temp, cutoff, l, e, c, fr, fm, of, mf, par, tp, obed=False):
        """Calculate the true counts of triplexes on the given dna regions"""
        self.triplexator_p = [ l, e, c, fr, fm, of, mf ]

        txp = find_triplex(rna_fasta=os.path.join(temp, "rna_temp.fa"), dna_region=self.dna_region,
                           temp=temp, organism=self.organism, remove_temp=remove_temp, tp=tp,
                           l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, par=par, genome_path=self.genome_path,
                           prefix="targeted_region", dna_fine_posi=False)
        txp.merge_rbs(rm_duplicate=True, region_set=self.dna_region, asgene_organism=self.organism, cutoff=cutoff)
        self.txp = txp
        txp.remove_duplicates()
        self.rbss = txp.merged_dict.keys()
        if len(self.rbss) == 0:
            print("ERROR: No potential binding event. Please change the parameters.")
            sys.exit(1)

        txpf = find_triplex(rna_fasta=os.path.join(temp, "rna_temp.fa"), dna_region=self.dna_region,
                            temp=temp, organism=self.organism, remove_temp=remove_temp, tp=tp,
                            l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, par=par, genome_path=self.genome_path,
                            prefix="dbs", dna_fine_posi=True)
        txpf.remove_duplicates()
        txpf.merge_rbs(rbss=self.rbss, rm_duplicate=True, asgene_organism=self.organism)
        self.txpf = txpf

        self.counts_tr = OrderedDict()
        self.counts_dbs = OrderedDict()

        for rbs in self.rbss:
            tr = len(self.txp.merged_dict[rbs])
            self.counts_tr[rbs] = [tr, len(self.dna_region) - tr]
            self.counts_dbs[rbs] = len(self.txpf.merged_dict[rbs])

        self.region_dbd = self.txpf.sort_rbs_by_regions(self.dna_region)

        self.region_dbs = self.txpf.sort_rd_by_regions(regionset=self.dna_region)
        self.region_dbsm = {}
        self.region_coverage = {}

        for region in self.dna_region:
            self.region_dbsm[region.toString()] = self.region_dbs[region.toString()].get_dbs().merge(w_return=True)
            self.region_coverage[region.toString()] = float(self.region_dbsm[region.toString()].total_coverage()) / len \
                (region)

        if obed:
            btr = self.txp.get_dbs()
            btr = btr.gene_association(organism=self.organism)
            btr.write_bed(os.path.join(temp, obed+ "_target_region_dbs.bed"))
            dbss = txpf.get_dbs()
            # dbss = dbss.gene_association(organism=self.organism)
            dbss.write_bed(os.path.join(temp, obed + "_dbss.bed"))

    def random_test(self, repeats, temp, remove_temp, l, e, c, fr, fm, of, mf, rm, par, tp, filter_bed, alpha):
        """Perform randomization for the given times"""
        self.repeats = repeats
        marks = numpy.round(numpy.linspace(0, repeats - 1, num=41)).tolist()
        print("random_test")
        print(par)
        # Prepare the input lists for multiprocessing
        mp_input = []
        for i in range(repeats):
            mp_input.append([str(i), os.path.join(temp, "rna_temp.fa"), self.dna_region,
                             temp, self.organism, self.rbss, str(marks.count(i)),
                             str(l), str(e), str(c), str(fr), str(fm), str(of), str(mf), str(rm),
                             filter_bed, self.genome_path, tp, par])
        # Multiprocessing
        print("\t\t|0%                  |                100%|")
        print("\t\t[", end="")
        pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
        mp_output = pool.map(random_each, mp_input)
        # print(mp_output)
        pool.close()
        pool.join()
        print("]")

        # Processing the result
        self.region_matrix = []
        self.dbss_matrix = []
        self.data = {"region": {"ave": [],
                                "sd": [],
                                "p": [],
                                "sig_region": [],
                                "sig_boolean": []},
                     "dbs": {"ave": [],
                             "sd": [],
                             "p": [],
                             "sig_region": [],
                             "sig_boolean": []}}

        region_counts = [v[0] for v in mp_output]
        dbss_counts = [v[1] for v in mp_output]

        for i, rbs in enumerate(self.rbss):

            counts_regions = [v[i] for v in region_counts]

            self.data["region"]["ave"].append(numpy.mean(counts_regions))
            self.data["region"]["sd"].append(numpy.std(counts_regions))
            num_sig = len([h for h in counts_regions if h > self.counts_tr[rbs][0]])
            p_region = float(num_sig) / repeats
            self.data["region"]["p"].append(p_region)
            self.region_matrix.append(counts_regions)

            if p_region < alpha:
                self.data["region"]["sig_region"].append(rbs)
                self.data["region"]["sig_boolean"].append(True)
            else:
                self.data["region"]["sig_boolean"].append(False)

            try:
                if p_region < self.topDBD[1]: self.topDBD = [rbs.str_rna(pa=False), p_region]
            except:
                self.topDBD = [rbs.str_rna(pa=False), p_region]

            # Analysis based on DBSs
            if self.showdbs:
                counts_dbss = [v[i] for v in dbss_counts]

                self.data["dbs"]["ave"].append(numpy.mean(counts_dbss))
                self.data["dbs"]["sd"].append(numpy.std(counts_dbss))
                num_sig = len([h for h in counts_dbss if h > self.counts_dbs[rbs]])
                p_dbs = float(num_sig) / repeats
                self.data["dbs"]["p"].append(p_dbs)
                self.dbss_matrix.append(counts_dbss)
                if p_dbs < alpha:
                    self.data["dbs"]["sig_region"].append(rbs)
                    self.data["dbs"]["sig_boolean"].append(True)
                else:
                    self.data["dbs"]["sig_boolean"].append(False)

        self.region_matrix = numpy.array(self.region_matrix)
        if self.showdbs: self.dbss_matrix = numpy.array(self.dbss_matrix)

    def dbd_regions(self, sig_region, output):
        """Generate the BED file of significant DBD regions and FASTA file of the sequences"""
        dbd_regions(exons=self.rna_regions, sig_region=sig_region, rna_name=self.rna_name, output=output)

    def lineplot(self, txp, dirp, ac, cut_off, log, ylabel, linelabel, showpa, sig_region, filename):
        """Generate lineplot for RNA"""

        lineplot(txp=txp, rnalen=self.rna_len, rnaname=self.rna_name, dirp=dirp, sig_region=sig_region,
                 cut_off=cut_off, log=log, ylabel=ylabel, linelabel=linelabel,
                 filename=filename, ac=ac, showpa=showpa)

    def boxplot(self, dir, matrix, sig_region, truecounts, sig_boolean, ylabel, filename):
        """Generate the visualized plot"""
        tick_size = 8
        label_size = 9

        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6, 4))
        max_y = int(max([matrix.max()] + truecounts) * 1.1) + 1
        min_y = max(int(matrix.min() * 0.9) - 1, 0)

        # Significant DBD
        rect = patches.Rectangle(xy=(1, 0), width=0.8, height=max_y, facecolor=sig_color,
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        for i, r in enumerate(sig_boolean):
            if r:
                rect = patches.Rectangle(xy=(i + 0.6, min_y), width=0.8, height=max_y, facecolor=sig_color,
                                         edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
                ax.add_patch(rect)

        # Plotting

        bp = ax.boxplot(matrix.transpose(), notch=False, sym='o', vert=True,
                        whis=1.5, positions=None, widths=None,
                        patch_artist=True, bootstrap=None)
        z = 10
        plt.setp(bp['boxes'], color=nontarget_color, alpha=1, edgecolor="none")
        plt.setp(bp['whiskers'], color='black', linestyle='-', linewidth=1, zorder=z, alpha=1)
        plt.setp(bp['fliers'], markerfacecolor='gray', color='white', alpha=0.3, markersize=1.8, zorder=z)
        plt.setp(bp['caps'], color='white', zorder=-1)
        plt.setp(bp['medians'], color='black', linewidth=1.5, zorder=z + 1)

        # Plot target regions
        plt.plot(range(1, len(self.rbss) + 1), truecounts, markerfacecolor=target_color,
                 marker='o', markersize=5, linestyle='None', markeredgecolor="white", zorder=z + 5)

        ax.set_xlabel(self.rna_name + " DNA Binding Domains", fontsize=label_size)
        ax.set_ylabel(ylabel, fontsize=label_size, rotation=90)

        ax.set_ylim([min_y, max_y])
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.rbss], rotation=35,
                           ha="right", fontsize=tick_size)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(tick_size)

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
        ax.tick_params(axis='y', which='both', left='on', right='off', labelbottom='off')

        # Legend
        dot_legend, = plt.plot([1, 1], color=target_color, marker='o', markersize=5, markeredgecolor="white",
                               linestyle='None')
        bp_legend, = plt.plot([1, 1], color=nontarget_color, linewidth=6, alpha=1)

        ax.legend([dot_legend, bp_legend, rect], ["Target Regions", "Non-target regions", "Significant DBD"],
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  prop={'size': 9}, ncol=3, numpoints=1)
        bp_legend.set_visible(False)
        dot_legend.set_visible(False)

        # f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(dir, filename + ".png"), facecolor='w', edgecolor='w',
                  bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)

        pp = PdfPages(os.path.join(dir, filename + '.pdf'))
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def gen_html(self, directory, parameters, obed, align=50, alpha=0.05, score=False):
        """Generate the HTML file"""
        dir_name = os.path.basename(directory)
        html_header = "Genomic Region Test: " + dir_name
        link_ds = OrderedDict()
        link_ds["RNA"] = "index.html"
        link_ds["Target regions"] = "target_regions.html"
        link_ds["Parameters"] = "parameters.html"

        if self.organism == "hg19":
            self.ani = "human"
        elif self.organism == "hg38":
            self.ani = "human"
        elif self.organism == "mm9":
            self.ani = "mouse"

        ##################################################
        # index.html

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        # Plots
        html.add_figure("lineplot_region.png", align="left", width="45%", more_images=["boxplot_regions.png"])
        if self.showdbs:
            html.add_figure("lineplot_dbs.png", align="left", width="45%", more_images=["boxplot_dbs.png"])

        if self.showdbs:
            header_list = [["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics",
                            "Target Regions", "Non-target Regions", None, "Statistics"],
                           ["", "", "with DBS", "without DBS", "with DBS (average)", "s.d.", "<i>p</i>-value",
                            "NO. DBSs", "NO. DBSs (average)", "s.d.", "<i>p</i>-value"]]
            header_titles = [["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                              "Regions from randomization", None, "Statistics based on target regions",
                              "Given target regions on DNA", "Regions from randomization", None,
                              "Statistics based on DNA Binding Sites"],
                             ["", "",
                              "Number of target regions with DBS binding",
                              "Number of target regions without DBS binding",
                              "Average number of regions from randomization with DBS binding",
                              "Standard deviation", "P value",
                              "Number of related DNA Binding Sites binding to target regions",
                              "Average number of DNA Binding Sites binding to random regions",
                              "Standard deviation", "P-value"]]
            border_list = [" style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:2pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\""]
        else:
            header_list = [["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics", None],
                           ["", "", "with DBS", "without DBS", "with DBS (average)", "s.d.", "<i>p</i>-value",
                            "z-score"]]
            header_titles = [["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                              "Regions from randomization", None, "Statistics based on target regions", None],
                             ["", "",
                              "Number of target regions with DBS binding",
                              "Number of target regions without DBS binding",
                              "Average number of regions from randomization with DBS binding",
                              "Standard deviation", "P value", "Z-score"]]
            border_list = [" style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", ""]

        type_list = 'ssssssssssssssss'
        col_size_list = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
        data_table = []

        for i, rbs in enumerate(self.rbss):
            if self.data["region"]["p"][i] < alpha:
                p_region = "<font color=\"red\">" + value2str(self.data["region"]["p"][i]) + "</font>"

            else:
                p_region = value2str(self.data["region"]["p"][i])
            zs = (self.counts_tr[rbs][0] - self.data["region"]["ave"][i]) / self.data["region"]["sd"][i]
            new_line = [str(i + 1),
                        rbs.str_rna(pa=False),
                        '<a href="dbd_region.html#' + rbs.str_rna() +
                        '" style="text-align:left">' + str(self.counts_tr[rbs][0]) + '</a>',
                        str(self.counts_tr[rbs][1]),
                        value2str(self.data["region"]["ave"][i]),
                        value2str(self.data["region"]["sd"][i]),
                        p_region,
                        value2str(zs)]
            if self.showdbs:
                if self.data["dbs"]["p"][i] < alpha:
                    p_dbs = "<font color=\"red\">" + value2str(self.data["dbs"]["p"][i]) + "</font>"
                else:
                    p_dbs = value2str(self.data["dbs"]["p"][i])

                new_line += [str(self.counts_dbs[rbs]),
                             value2str(self.data["dbs"]["ave"][i]),
                             value2str(self.data["dbs"]["sd"][i]),
                             p_dbs]
            data_table.append(new_line)

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, border_list=border_list, sortable=True)

        html.add_heading("Notes")
        html.add_list(["RNA name: " + self.rna_name,
                       "Randomization is performed for " + str(self.repeats) + " times.",
                       "DBD stands for DNA Binding Domain on RNA.",
                       "DBS stands for DNA Binding Site on DNA."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "index.html"))

        #############################################################
        # RNA subpage: Profile of targeted regions for each merged DNA Binding Domain
        #############################################################

        header_list = ["#", "Target Region",
                       "Associated Gene",
                       "No. of DBSs",
                       "DBS coverage"]
        header_titles = ["Rank", "Given target regions from BED files",
                         "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                         "Number of DNA Binding Sites locate within the region",
                         "The proportion of the region covered by DBS binding"]

        #########################################################
        # dbd_region.html
        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for rbsm in self.rbss:
            html.add_heading("DNA Binding Domain: " + rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for i, region in enumerate(self.txp.merged_dict[rbsm]):
                # Add information
                data_table.append([str(i + 1),
                                   '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                                   "&position=" + region.chrom + "%3A" + str(region.initial) + "-" + str(region.final) +
                                   '" style="text-align:left">' + region.toString(space=True) + '</a>',
                                   split_gene_name(gene_name=region.name, org=self.organism),
                                   str(len(self.region_dbs[region.toString()])),
                                   value2str(self.region_coverage[region.toString()])
                                   ])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 auto_width=True, header_titles=header_titles, sortable=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "dbd_region.html"))

        #############################################################
        # Targeted regions centered
        #############################################################

        ##############################################################################################
        # target_regions.html
        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if score:
            header_list = ["#", "Target region", "Associated Gene", "DBSs Count",
                           "DBS coverage", "Score", "Sum of ranks"]
            header_titles = ["Rank",
                             "Target regions loaded from the given BED file",
                             "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                             "Number of DNA Binding Sites within the region",
                             "The proportion of the region covered by DBS binding",
                             "Scores from BED file",
                             "Sum of all the left-hand-side ranks"]
        else:
            header_list = ["#", "Target region", "Associated Gene", "DBSs Count",
                           "DBS coverage", "Sum of ranks"]
            header_titles = ["Rank",
                             "Target regions loaded from the given BED file",
                             "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                             "Number of DNA Binding Sites within the region",
                             "The proportion of the region covered by DBS binding",
                             "Sum of all the left-hand-side ranks"]
        html.add_heading("Target Regions")
        data_table = []

        if not self.dna_region.sorted: self.dna_region.sort()

        # Calculate the ranking
        rank_count = len(self.dna_region) - rank_array([len(self.region_dbs[p.toString()]) for p in self.dna_region])
        rank_coverage = len(self.dna_region) - rank_array([self.region_coverage[p.toString()] for p in self.dna_region])

        if score:
            rank_score = len(self.dna_region) - rank_array([float(p.data.split("\t")[0]) for p in self.dna_region])
            rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]
            sum_rank = rank_array(rank_sum)  # method='min'
        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]
            sum_rank = rank_array(rank_sum)

        for i, region in enumerate(self.dna_region):
            dbs_counts = str(len(self.region_dbs[region.toString()]))
            dbs_cover = value2str(self.region_coverage[region.toString()])
            newline = [str(i + 1),
                       '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                       "&position=" + region.chrom + "%3A" + str(region.initial) + "-" + str(region.final) +
                       '" style="text-align:left">' + region.toString(space=True) + '</a>',
                       split_gene_name(gene_name=region.name, org=self.organism),
                       '<a href="region_dbs.html#' + region.toString() +
                       '" style="text-align:left">' + dbs_counts + '</a>',
                       dbs_cover]
            if score:
                dbs_score = str(region.data.split("\t")[0])
                newline += [dbs_score,
                            str(int(rank_sum[i] + 3))]
            else:
                ranking = str(int(rank_sum[i] + 2))
                newline += [ranking]

            data_table.append(newline)
            if score:
                region.data = "\t".join([dbs_counts, dbs_cover, dbs_score, ranking])
            else:
                region.data = "\t".join([dbs_counts, dbs_cover, ranking])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, sortable=True)
        html.add_heading("Notes")
        html.add_list(["All target regions without any bindings are ignored."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, "target_regions.html"))

        self.dna_region.sort_score()
        self.dna_region.write_bed(os.path.join(directory, obed + "_target_regions.bed"))

        ############################
        # Subpages for targeted region centered page
        # region_dbs.html
        header_list = ["RBS", "DBS", "Strand", "Score", "Motif", "Orientation"]

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for i, region in enumerate(self.dna_region):
            if len(self.region_dbs[region.toString()]) == 0:
                continue
            else:
                html.add_heading("Associated gene: " + split_gene_name(gene_name=region.name, org=self.organism),
                                 idtag=region.toString())
                html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                                       "&position=" + region.chrom + "%3A" + str(region.initial) + "-" + str(
                    region.final) +
                                       '" style="margin-left:50">' +
                                       region.toString(space=True) + '</a>'])
                data_table = []
                for rd in self.region_dbs[region.toString()]:
                    data_table.append([rd.rna.str_rna(pa=False),
                                       '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.organism +
                                       "&position=" + rd.dna.chrom + "%3A" + str(rd.dna.initial) + "-" + str(
                                           rd.dna.final) +
                                       '" style="text-align:left">' + rd.dna.toString(space=True) + '</a>',
                                       rd.dna.orientation,
                                       rd.score,
                                       rd.motif,
                                       rd.orient])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     auto_width=True)
        html.write(os.path.join(directory, "region_dbs.html"))

        ###############################################################################33
        ################ Parameters.html

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        html.add_heading("Parameters")
        header_list = ["Description", "Arguments", "Value"]

        data_table = [["RNA sequence name", "-rn", parameters.rn],
                      ["Input RNA sequence file", "-r", os.path.basename(parameters.r)],
                      ["Input BED file", "-bed", os.path.basename(parameters.bed)],
                      ["Output directory", "-o", os.path.basename(parameters.o)],
                      ["Organism", "-organism", parameters.organism],
                      ["Number of repitetion of andomization", "-n", str(parameters.n)],
                      ["Alpha level for rejection p value", "-a", str(parameters.a)],
                      ["Cut off value for filtering out the low counts of DBSs", "-ccf", str(parameters.ccf)],
                      ["Remove temporary files", "-rt", str(parameters.rt)],
                      ["Input BED file for masking in randomization", "-f", str(parameters.f)],
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

    def save_profile(self, output, bed):
        """Save some statistics for comparison with other results"""
        pro_path = os.path.join(os.path.dirname(output), "profile.txt")
        exp = os.path.basename(output)
        # tag = os.path.basename(os.path.dirname(rnafile))
        tar_reg = os.path.basename(bed)
        # RNA name with region
        if self.rna_regions:
            trans = "Transcript:"
            for r in self.rna_regions:
                trans += r[0] + ":" + str(r[1]) + "-" + str(r[2])
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
                                         self.organism, tar_reg, str(len(self.data["region"]["sig_region"])),
                                         self.topDBD[0], value2str(self.topDBD[1]), r_genes])
                        new_exp = False
                    else:
                        newlines.append(line)
                if new_exp:
                    newlines.append([exp, rna, output.split("_")[-1],
                                     self.organism, tar_reg, str(len(self.data["region"]["sig_region"])),
                                     self.topDBD[0], value2str(self.topDBD[1]), r_genes])
        else:
            newlines.append(["Experiment", "RNA_names", "Tag", "Organism", "Target_region", "No_sig_DBDs",
                             "Top_DBD", "p-value", "closest_genes"])
            newlines.append([exp, rna, output.split("_")[-1],
                             self.organism, tar_reg, str(len(self.data["region"]["sig_region"])),
                             self.topDBD[0], value2str(self.topDBD[1]), r_genes])

        with open(pro_path, 'w') as f:
            for lines in newlines:
                print("\t".join(lines), file=f)

