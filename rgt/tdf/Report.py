# Python Libraries


import os
import numpy
import natsort
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import MaxNLocator
from collections import OrderedDict
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter
from ..SequenceSet import SequenceSet
from ..GenomicRegionSet import GenomicRegionSet
from ..Util import Html, SequenceType

# Color code for all analysis
# target_color = "mediumblue"
target_color = "royalblue"
# nontarget_color = "darkgrey"
nontarget_color = "darkgrey"
sig_color = "lightgrey"
legend_fontsize=8

def uniq(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]

def read_ac(path, cut_off, rnalen):
    """Read the RNA accessibility file and output its positions and values

    The file should be a simple table with two columns:
    The first column is the position and the second one is the value
    '#' will be skipped

    """
    access = []
    with open(path) as f:
        i = 0
        while i < rnalen:
            for line in f:
                line = line.split()
                if not line:
                    continue
                elif line[0][0] == "#":
                    continue
                elif len(line) < 2:
                    continue
                else:
                    v = line[1]
                    if v == "NA":
                        access.append(0)
                    else:
                        try:
                            v = 2 ** (-float(v))
                        except:
                            continue
                        if v >= cut_off:
                            access.append(1)
                        else:
                            access.append(0)
                    i += 1
    return access

def value2str(value):
    if isinstance(value, str):
        try: value = float(value)
        except: return value
    if value == 0: return "0"
    if isinstance(value, int): return str(value)
    elif isinstance(value, float):
        if abs(value) >= 1000:
            try: r = "{}".format(int(value))
            except: r = "Inf"
        elif 1000 > abs(value) > 10: r = "{:.1f}".format(value)
        elif 10 > abs(value) >= 1: r = "{:.2f}".format(value)
        elif 1 > abs(value) >= 0.05: r = "{:.2f}".format(value)
        elif 0.05 > abs(value) > 0.0001: r = "{:.4f}".format(value)
        else: r = "{:.1e}".format(value)
        return r


def rank_array(a):
    try:
        a = numpy.array(a)
    except:
        a = numpy.array([float(b) for b in a])
    sa = numpy.searchsorted(numpy.sort(a), a)
    return sa

def region_link_internet(organism, region):
    ani = None
    if organism == "hg19":
        ani = "human"
    elif organism == "hg38":
        ani = "human"
    elif organism == "mm9":
        ani = "mouse"
    if ani:
        region_link = "".join(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=', organism,
                               "&position=", region.chrom, "%3A", str(region.initial), "-",
                               str(region.final), '" style="text-align:left" target="_blank">',
                               region.toString(space=False), '</a>'])
    else:
        if organism == "tair10":
            region_link = "".join(
                ['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=',
                 region.chrom, "%3A", str(region.initial), "..", str(region.final),
                 '" target="_blank">',
                 region.toString(space=False), '</a>'])
        else:
            region_link = region.toString(space=False)
    return region_link


def split_gene_name(gene_name, org):
    if gene_name == None:
        return ""
    if gene_name[0:2] == "chr":
        return gene_name

    if org == "hg19":
        ani = "Homo_sapiens"
    elif org == "hg38":
        ani = "Homo_sapiens"
    elif org == "mm9":
        ani = "Mus_musculus"
    else:
        ani = None

    if not ani:
        if org == "tair10":
            p1 = "".join(['<a href="https://www.arabidopsis.org/servlets/TairObject?name=', gene_name,
                          '&type=locus" target="_blank">', gene_name, '</a>'])
            return p1
        else:
            return gene_name
    else:
        p1 = '<a href="http://www.ensembl.org/' + ani + \
             "/Gene/Summary?g="
        p2 = '" target="_blank">'
        p3 = '</a>'

        if ":" in gene_name:
            genes = gene_name.split(":")
            genes = list(set(genes))
            result = []
            c = 0
            for i, g in enumerate(genes):
                if "(" in g:
                    d = g.partition('(')[2].partition(')')[0]
                    g = g.partition('(')[0]
                    if "-" in d:
                        result.insert(0, p1 + g + p2 + g + p3 + "(" + d + ")")
                    else:
                        result.append(p1 + g + p2 + g + p3 + "(" + d + ")")
                else:
                    c += 1
                    if c < 6:
                        result.append(p1 + g + p2 + g + p3)

            result = ",".join(result)

        elif gene_name == ".":
            result = "none"

        else:
            if "(" in gene_name:
                d = gene_name.partition('(')[2].partition(')')[0]
                g = gene_name.partition('(')[0]
                result = p1 + g + p2 + g + p3 + "(" + d + ")"
            else:
                result = p1 + gene_name + p2 + gene_name + p3

        return result



class Report(object):
    def __init__(self, pars, input, triplexes, stat):
        self.pars = pars
        self.input = input
        self.triplexes = triplexes
        self.stat = stat
        self.link_d = OrderedDict()

    def plot_lines(self, tpx, ylabel, linelabel, filename):
        """Generate the plots for demonstration of RBS
        """
        rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
        rnas.read_fasta(self.pars.r)

        rna_len = rnas.total_len()
        self.triplexes.autobinding.rna_track(rnalen=rna_len)

        self.lineplot(tpx=tpx, rnalen=rna_len, rnaname=self.pars.rn, dirp=self.pars.o,
                      sig_region=self.stat.sig_DBD, cut_off=self.pars.ccf, log=self.pars.log,
                      ylabel=ylabel, linelabel=linelabel, filename=filename,
                      ac=self.triplexes.autobinding.rna_track, showpa=self.pars.showpa,
                      exons=self.input.rna.regions)

    def lineplot(self, tpx, rnalen, rnaname, dirp, sig_region, cut_off, log, ylabel, linelabel,
                 filename, ac=None, showpa=False, exons=None):
        # Plotting
        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6, 4))

        # Extract data points
        x = list(range(rnalen))
        # print(rnalen)
        if log:
            all_y = [1] * rnalen
            p_y = [1] * rnalen
            a_y = [1] * rnalen
        else:
            all_y = [0] * rnalen
            p_y = [0] * rnalen
            a_y = [0] * rnalen

        tpx.remove_duplicates_by_dbs()
        for rd in tpx:
            # print(str(rd.rna.initial), str(rd.rna.final))
            if rd.rna.orientation == "P":
                for i in range(rd.rna.initial, rd.rna.final):
                    p_y[i] += 1
                    all_y[i] += 1
            if rd.rna.orientation == "A":
                for i in range(rd.rna.initial, rd.rna.final):
                    a_y[i] += 1
                    all_y[i] += 1
        # Log
        if log:
            all_y = numpy.log(all_y)
            p_y = numpy.log(p_y)
            a_y = numpy.log(a_y)
            max_y = max(all_y) + 0.5
            min_y = 1
            ylabel += "(log10)"
        else:
            max_y = float(max(all_y) * 1.1)
            min_y = 0

        if ac:
            # min_y = float(max_y * (-0.09))
            max_y = float(max_y * (1.05))

        # Plotting
        for rbs in sig_region:
            rect = patches.Rectangle(xy=(rbs.initial, 0), width=len(rbs), height=0.95*max_y, facecolor=sig_color,
                                     edgecolor="none", alpha=1, lw=None, label="Significant DBD", zorder=2)
            ax.add_patch(rect)

        # RNA accessbility
        if ac:
            if isinstance(ac, list):
                n_value = ac
            else:
                n_value = read_ac(ac, cut_off, rnalen=rnalen)
            drawing = False
            for i in x:
                if n_value[i] > 0:
                    if drawing:
                        continue
                    else:
                        last_i = i
                        drawing = True
                elif drawing:
                    # ax.add_patch(patches.Rectangle((last_i, min_y), i - last_i, -min_y,
                    #                                fill=True, color="silver", snap=False, linewidth=0,
                    #                                label="Autobinding"))

                    ax.add_patch(patches.Rectangle(xy=(last_i, 0.95*max_y), width=i - last_i, height=0.05*max_y,
                                                   fill=True, color="lightcoral", snap=False, linewidth=0,zorder=2,
                                                   label="Autobinding"))
                    drawing = False
                else:
                    continue

        lw = 1.5
        if showpa:
            ax.plot(x, all_y, color=target_color, alpha=1, lw=lw, label="Parallel + Anti-parallel")
            ax.plot(x, p_y, color="purple", alpha=1, lw=lw, label="Parallel")
            ax.plot(x, a_y, color="dimgrey", alpha=.8, lw=lw, label="Anti-parallel")
        else:
            ax.plot(x, all_y, color=target_color, alpha=0, lw=lw, markeredgecolor="none", zorder=10)
            ax.fill_between(x, 0, all_y, facecolor=target_color, alpha=1,  edgecolor="none", label=linelabel, zorder=10)
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.yaxis.set_ticks_position('left')
            ax.xaxis.set_ticks_position('bottom')



        # Legend
        handles, labels = ax.get_legend_handles_labels()
        legend_h = []
        legend_l = []
        for uniqlabel in uniq(labels):
            legend_h.append(handles[labels.index(uniqlabel)])
            legend_l.append(uniqlabel)
        ax.legend(legend_h, legend_l, fontsize=legend_fontsize,
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  ncol=3)

        # XY axis
        ax.set_xlim(left=0, right=rnalen)
        ax.set_ylim([min_y, max_y])
        for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9)
        ax.set_xlabel(rnaname + " sequence (nt)", fontsize=9)

        ax.set_ylabel(ylabel, fontsize=9, rotation=90)

        if None:
            if exons and len(exons) > 1:
                w = 0
                i = 0
                h = (max_y - min_y) * 0.02

                for exon in exons:
                    l = abs(exon[2] - exon[1])

                    # print([i,l,w])
                    # ax.axvline(x=w, color="gray", alpha=0.5, zorder=100)
                    if i % 2 == 0:
                        rect = matplotlib.patches.Rectangle((w, max_y - h), l, h, color="moccasin")
                    else:
                        rect = matplotlib.patches.Rectangle((w, max_y - h), l, h, color="gold")
                    ax.add_patch(rect)
                    i += 1
                    w += l
                ax.text(rnalen * 0.01, max_y - 2 * h, "exon boundaries", fontsize=5, color='black')

        # f.tight_layout(pad=1.08)

        f.savefig(os.path.join(dirp, filename), facecolor='w', edgecolor='w',
                  bbox_extra_artists=(plt.gci()),dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        ax.legend(legend_h, legend_l,fontsize=legend_fontsize,
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  prop={'size': 12}, ncol=3)
        pp = PdfPages(os.path.splitext(os.path.join(dirp,filename))[0] +'.pdf')
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def barplot(self, filename, dbs=False):
        """Generate the barplot to show the difference between target promoters and non-target promoters"""

        def to_percent(y, position):
            # Ignore the passed in position. This has the effect of scaling the default
            # tick locations.
            s = str(100 * y)
            return s + '%'

        f, ax = plt.subplots(1, 1, dpi=300, figsize=(6, 4))
        ind = list(range(len(self.stat.rbss)))
        width = 0.35

        if not dbs:
            propor_de = [float(b[0]) / len(self.input.dna.target_regions) for b in list(self.stat.frequency["promoters"]["de"].values())]
            propor_nde = [float(b[0]) / len(self.input.dna.nontarget_regions) for b in list(self.stat.frequency["promoters"]["nde"].values())]
        else:
            propor_de = [float(b[0]) / (b[0] + b[1]) for b in list(self.stat.frequency["hits"]["de"].values())]
            propor_nde = [float(b[0]) / (b[0] + b[1]) for b in list(self.stat.frequency["hits"]["nde"].values())]

        max_y = max([max(propor_de), max(propor_nde)]) * 1.2

        # Plotting
        rect = patches.Rectangle(xy=(1, 0), width=0.8, height=max_y, facecolor=sig_color,
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        for i, rbs in enumerate(self.stat.rbss):
            if rbs in self.stat.sig_DBD:
                rect = patches.Rectangle(xy=(i + 0.05, 0), width=0.9, height=max_y, facecolor=sig_color,
                                         edgecolor="none", lw=None, label="Significant DBD")
                ax.add_patch(rect)

        rects_de = ax.bar([i + 0.325 for i in ind], propor_de, width, color=target_color,
                          edgecolor="none", label="Target promoters")
        rects_nde = ax.bar([i + 0.325 + width for i in ind], propor_nde, width, color=nontarget_color,
                           edgecolor="none", label="Non-target promoters")

        # Legend
        # tr_legend, = plt.plot([1, 1], color=target_color, linewidth=6, alpha=1)
        # ntr_legend, = plt.plot([1, 1], color=nontarget_color, linewidth=6, alpha=1)
        ax.legend([rects_de, rects_nde, rect], ["Target promoters", "Non-target promoters", "Significant DBD"],
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,fontsize=legend_fontsize,
                  ncol=3)
        # tr_legend.set_visible(False)
        # ntr_legend.set_visible(False)

        tick_size = 8
        # Y axis
        ax.set_ylim([0, max_y])
        formatter = FuncFormatter(to_percent)
        # Set the formatter
        ax.yaxis.set_major_formatter(formatter)
        ax.tick_params(axis='y', which='both', left=True, right=False, labelbottom=False)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9)
        ax.set_ylabel("Proportion of promoters (%)", fontsize=9, rotation=90)

        # X axis
        ax.set_xlim([0, len(self.stat.rbss)])
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
        ax.set_xticks([i + 0.5 for i in range(len(self.stat.rbss))])
        # ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.stat.rbss], rotation=35,
        #                    ha="right", fontsize=tick_size)
        ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.stat.rbss], ha="center", fontsize=tick_size)
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)

        for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9)
        ax.set_xlabel(self.pars.rn + " DNA Binding Domains", fontsize=9)

        # f.tight_layout(pad=1.08)
        f.savefig(os.path.join(self.pars.o, filename), facecolor='w', edgecolor='w',
                  bbox_extra_artists=(plt.gci()), dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        ax.legend([rects_de, rects_nde, rect], ["Targets", "Non-targets", "Sig. DBD"],
                  fontsize=legend_fontsize,
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  prop={'size': 12}, ncol=3)
        pp = PdfPages(os.path.splitext(os.path.join(self.pars.o, filename))[0] + '.pdf')
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def boxplot(self, filename, matrix, sig_region, truecounts, sig_boolean, ylabel):
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
        plt.plot(list(range(1, len(self.stat.rbss) + 1)), truecounts, markerfacecolor=target_color,
                 marker='o', markersize=8, linestyle='None', markeredgecolor="white", zorder=z + 5, label="Target Regions")

        ax.set_xlabel(self.pars.rn + " DNA Binding Domains", fontsize=label_size)
        ax.set_ylabel(ylabel, fontsize=label_size, rotation=90)

        ax.set_ylim([min_y, max_y])
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))

        # ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.stat.rbss], rotation=35,
        #                    ha="right", fontsize=tick_size)
        ax.set_xticklabels([dbd.str_rna(pa=False) for dbd in self.stat.rbss], ha="center", fontsize=tick_size)
        for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(tick_size)

        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        ax.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
        ax.tick_params(axis='y', which='both', left=True, right=False, labelbottom=False)

        # Legend
        dot_legend, = plt.plot([1, 1], color=target_color, marker='o', markersize=5, markeredgecolor="white",
                               linestyle='None')
        # bp_legend, = plt.plot([1, 1], color=nontarget_color, linewidth=6, alpha=1)

        ax.legend([dot_legend, bp["boxes"][0], rect], ["Target Regions", "Non-target regions", "Significant DBD"],
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,fontsize=legend_fontsize,
                  ncol=3, numpoints=1)
        # bp_legend.set_visible(False)
        # dot_legend.set_visible(False)

        # f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.savefig(os.path.join(self.pars.o, filename), facecolor='w', edgecolor='w', bbox_extra_artists=(plt.gci()), dpi=300)
        # PDF
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
        ax.xaxis.label.set_size(14)
        ax.yaxis.label.set_size(14)
        ax.legend([dot_legend, bp["boxes"][0], rect], ["Targets", "Non-targets", "Sig. DBD"],
                  fontsize=legend_fontsize,
                  bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0.,
                  prop={'size': 12}, ncol=3)

        pp = PdfPages(os.path.splitext(os.path.join(self.pars.o, filename))[0] + '.pdf')
        pp.savefig(f, bbox_extra_artists=(plt.gci()), bbox_inches='tight')
        pp.close()

    def gen_html_promotertest(self):
        align = 50
        dir_name = os.path.basename(self.pars.o)
        # check_dir(directory)
        html_header = "Promoter Test: " + dir_name

        self.link_d["RNA"] = "index.html"
        self.link_d["Sig promoters"] = "spromoters.html"
        self.link_d["All promoters"] = "promoters.html"
        self.link_d["Autobinding"] = "autobinding.html"
        self.link_d["Parameters"] = "parameters.html"

        #############################################################
        # Index main page
        #############################################################
        html = Html(name=html_header, links_dict=self.link_d,
                    fig_dir=os.path.join(os.path.dirname(self.pars.o), "style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_figure(self.pars.rn+"_lineplot.png", align="left", width="45%",
                        more_images=[self.pars.rn+"_barplot.png"])

        if self.pars.showdbs:
            html.add_figure("plot_dbss.png", align="left", width="45%", more_images=["bar_dbss.png"])

        # Table of merged TBS on promoters
        if self.pars.showdbs:
            header_list = [["", "", "Promoters", None, None, None, None, None, "DBSs", None, None, None, None, None],
                           ["#", "DBD",
                            "Target Promoter", None, "Non-target Promoter", None, "Statistics", None,
                            "Target Promoter", None, "Non-target Promoter", None, "Statistics", None],
                           [" ", " ",
                            "with TTS", "without TTS", "with TTS", "without TTS", "OR", "<i>p</i>-value",
                            "No. TTSs", "Other TTSs", "No. TTSs", "Other TTSs", "OR", "<i>p</i>-value"]]
            header_titles = [["", "", "Statistics on promoter level", None, None, None, None, None,
                              "Statistics on TTS level", None, None, None, None, None],
                             ["Rank of the talbe",
                              "DNA Binding Domain which is the functional region on RNA.",
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on promoters", None,
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on TTSs", None],
                             ["",
                              "",
                              "Number of target promoters which contain TTSs.",
                              "Number of target promoters which don't contain TTSs.",
                              "Number of non-target promoters which contain TTSs.",
                              "Number of non-target promoters which don't contain TTSs.",
                              "Odds Ratio", "P-value",
                              "Number of TTSs found in the target promoters.",
                              "Number of TTSs not found in the target promoters.",
                              "Number of TTSs found in the non-target promoters.",
                              "Number of TTSs not found in the non-target promoters.",
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
            header_list = [["#", "DBD", "Target Promoter", None, "Non-target Promoter", None, "Statistics", None, "Autobinding"],
                           [" ", " ", "with TTS", "without TTS", "with TTS", "without TTS", "OR", "<i>p</i>", "Number"]]
            header_titles = [["Rank of the talbe",
                              "DNA Binding Domain which is the functional region on RNA.",
                              "Promoters of the differential expression genes.", None,
                              "Promoters of the non-differential expression genes.", None,
                              "Statistics based on promoters", None, "The RNA regions which bind to themselves"],
                             ["",
                              "",
                              "Number of target promoters which contain TTSs.",
                              "Number of target promoters which don't contain TTSs.",
                              "Number of non-target promoters which contain TTSs.",
                              "Number of non-target promoters which don't contain TTSs.",
                              "Odds Ratio", "P-value", "Number"]
                             ]
            border_list = ["style=\"border-right:1pt solid gray\"",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\"", "",
                           "style=\"border-right:1pt solid gray\"",
                           "style=\"border-right:1pt solid gray\""]

        type_list = 'ssssssssssssssssssss'
        col_size_list = [20] * 20
        data_table = []
        rank = 0
        self.topDBD = ["-", 1]

        for rbs in self.stat.frequency["promoters"]["de"]:
            if self.stat.frequency["promoters"]["de"][rbs][0] < self.pars.ccf: continue
            rank += 1
            if self.stat.pvalue[rbs] < self.pars.a:
                p_promoter = "<font color=\"red\">%s</font>" % (value2str(self.stat.pvalue[rbs]))
            else:
                p_promoter = value2str(self.stat.pvalue[rbs])

            if self.pars.showdbs:
                if self.stat.hpvalue[rbs] < self.pars.a:
                    p_hit = "<font color=\"red\">%s</font>" % (value2str(self.stat.hpvalue[rbs]))
                else:
                    p_hit = value2str(self.stat.hpvalue[rbs])

            try:
                if self.stat.pvalue[rbs] < self.topDBD[1]:
                    self.topDBD = [rbs.str_rna(pa=False), self.stat.pvalue[rbs]]
            except:
                self.topDBD = [rbs.str_rna(pa=False), self.stat.pvalue[rbs]]

            new_row = [str(rank),
                       rbs.str_rna(pa=False),
                       '<a href="dbds_promoters.html#%s" style="text-align:left">%s</a>' %
                       (rbs.str_rna(), value2str(self.stat.frequency["promoters"]["de"][rbs][0])),
                       value2str(self.stat.frequency["promoters"]["de"][rbs][1]),
                       value2str(self.stat.frequency["promoters"]["nde"][rbs][0]),
                       value2str(self.stat.frequency["promoters"]["nde"][rbs][1]),
                       value2str(self.stat.oddsratio[rbs]),
                       p_promoter,
                       '<a href="autobinding.html">' +
                       str(len(self.triplexes.autobinding.merged_dict[rbs]))+ '</a>']
            if self.pars.showdbs:
                new_row += [value2str(self.stat.frequency["hits"]["de"][rbs][0]),
                            value2str(self.stat.frequency["hits"]["de"][rbs][1]),
                            value2str(self.stat.frequency["hits"]["nde"][rbs][0]),
                            value2str(self.stat.frequency["hits"]["nde"][rbs][1]),
                            value2str(self.stat.hoddsratio[rbs]),
                            p_hit]

            data_table.append(new_row)
        # data_table = sorted(data_table, key=lambda x: x[-1])
        # data_table = sorted(data_table, key=lambda x: float(x[-1]))
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titles, border_list=border_list, sortable=True)
        html.add_heading("Notes")
        html.add_list(["DBD stands for functional DNA Binding Domain on RNA.",
                       "TFO stands for triplex forming oligonucleotide.",
                       "TTS stands for triplex target DNA site."])
        ####
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "index.html"))

        #############################################################
        # RNA subpage: Profile of targeted promoters for each merged DNA Binding Domain
        #############################################################

        # dbds_promoters.html

        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if self.input.dna.scores:
            if self.pars.scoreh and self.input.dna.de_gene:
                score_header = [self.input.dna.de_gene.cond[0]]
            else:
                score_header = ["Fold Change Score"]

            header_list = ["#", "Promoter", "Gene", "TTSs counts", "TTS coverage"]
            header_list += score_header
            header_list += ["Sum of Ranks"]
            header_titles = ["", "Target promoters", "Gene symbol",
                             "Number of DNA Binding sites locating within the promoter",
                             "The proportion of promoter covered by binding sites"]
            header_titles += ["Scores loaded by their absolute values from gene list or BED input. If there is annotation error for the gene names, it shows zero instead."] * len(
                score_header)
            header_titles += ["Sum up the ranks from left-hand side columns"]

        else:
            header_list = ["#", "Promoter", "Gene", "TTSs Count", "TTS coverage", "Sum of Ranks"]

            header_titles = ["", "Target promoters", "Gene symbol",
                             "Number of DNA Binding sites locating within the promoter",
                             "The proportion of promoter covered by binding sites",
                             "Sum up the ranks from left-hand side columns"]

        for rbsm in self.stat.frequency["promoters"]["de"]:
            html.add_heading("DNA Binding Domain: " + rbsm.str_rna(), idtag=rbsm.str_rna())
            data_table = []
            # Calculate the ranking
            # rank_array
            rank_count = len(self.stat.tpx_de.merged_dict[rbsm]) - rank_array(
                [len(self.stat.promoter["de"]["rd"][p.toString()]) for p in self.stat.tpx_de.merged_dict[rbsm]])
            rank_coverage = len(self.stat.tpx_de.merged_dict[rbsm]) - rank_array(
                [self.stat.promoter["de"]["dbs_coverage"][p.toString()] for p in self.stat.tpx_de.merged_dict[rbsm]])

            if self.input.dna.scores:
                rank_score = len(self.stat.tpx_de.merged_dict[rbsm]) - rank_array([self.input.dna.scores[p.toString()] for p in self.stat.tpx_de.merged_dict[rbsm] ])
                rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

            else:
                rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

            for i, promoter in enumerate(self.stat.tpx_de.merged_dict[rbsm]):
                # Add information
                region_link = region_link_internet(self.pars.organism, promoter)
                # try:
                #     newline = [str(i + 1), region_link,
                #                split_gene_name(gene_name=self.ensembl2symbol[promoter.name], org=self.pars.organism),
                #                str(len(self.promoter["de"]["rd"][promoter.toString()])),
                #                value2str(self.promoter["de"]["dbs_coverage"][promoter.toString()])
                #                ]
                # except:
                newline = [str(i + 1), region_link,
                           split_gene_name(gene_name=promoter.name, org=self.pars.organism),
                           str(len(self.stat.promoter["de"]["rd"][promoter.toString()])),
                           value2str(self.stat.promoter["de"]["dbs_coverage"][promoter.toString()])
                           ]
                if self.input.dna.scores:
                    # if multiple_scores:
                    #     for j in range(len(score_ar[0])):
                    #
                    #         newline += [value2str(score_ar[i][j])]
                    # else:
                    newline += [value2str(abs(self.input.dna.scores[promoter.toString()]))]

                newline += ["<i>" + str(int(rank_sum[i])) + "</i>"]
                # print(newline)
                data_table.append(newline)
            data_table = sorted(data_table, key=lambda x: x[-1])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_titles, sortable=True, border_list=None, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "dbds_promoters.html"))

        ################################################################
        ############# Autobinding
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_heading("Autobinding")
        header_list = ["#", "DBD", "RNA", "DNA", "Score", "Motif", "Orientation", "Sequence"]

        t = []
        for rbs in self.stat.rbss:
            for i, rd in enumerate(self.triplexes.autobinding):
                if rbs.overlap(rd.rna):
                    t.append([str(i), rbs.str_rna(pa=False),
                              str(rd.rna.initial) + "-" + str(rd.rna.final),
                              str(rd.dna.initial) + "-" + str(rd.dna.final),
                              rd.score, rd.motif, rd.orient, '<pre><font size="2">' + "\n".join(rd.match) + "</font></pre>"])
        if len(t) > 0:
            html.add_zebra_table(header_list, col_size_list, type_list, t, align=align, cell_align="left",
                                 sortable=True, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "autobinding.html"))

        ################################################################
        ############# Parameters
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_heading("Parameters")
        header_list = ["Description", "Arguments", "Value"]

        if self.pars.de:
            de = os.path.basename(self.pars.de)
            bed = "False"
            bg = "False"
        else:
            de = "False"
            bed = os.path.basename(self.pars.bed)
            bg = os.path.basename(self.pars.bg)

        data_table = [["RNA sequence name", "-rn", str(self.pars.rn)],
                      ["Input RNA sequence file", "-r", os.path.basename(self.pars.r)],
                      ["Input file for defferentially expression gene list", "-de", str(de)],
                      ["Input BED file as promoters", "-bed", str(bed)],
                      ["Input BED file as backgrounds", "-bg", str(bg)],
                      ["Output directory", "-o", "/".join(self.pars.o.partition("/")[-3:])],
                      ["Organism", "-organism", self.pars.organism],
                      ["filter_havana", "-filter_havana", str(self.pars.filter_havana)],
                      ["protein_coding", "-protein_coding", str(self.pars.protein_coding)],
                      ["known_only", "-known_only", str(self.pars.known_only)],
                      ["Promoter length", "-pl", str(self.pars.pl)],
                      ["Alpha level for rejection p value", "-a", str(self.pars.a)],
                      ["Cut off value for filtering out the DBD with low counts of triplexes", "-ccf", str(self.pars.ccf)],
                      ["Remove temporary files", "-rt", str(self.pars.rt)],
                      # ["Input file for RNA accecibility", "-ac", str(self.pars.ac)],
                      # ["Cut off value for RNA accecibility", "-accf", str(self.pars.accf)],
                      ["Output the BED files for DNA binding sites.", "-obed", str(self.pars.obed)],
                      ["Show parallel and antiparallel bindings in the plot separately.", "-showpa",
                       str(self.pars.showpa)],
                      ["Minimum length", "-l", str(self.pars.l)],
                      ["Maximum error rate", "-e", str(self.pars.e)],
                      ["Tolerated number of consecutive errors", "-c", str(self.pars.c)],
                      ["Filtering repeats", "-fr", str(self.pars.fr)],
                      ["Filtering mode", "-fm", str(self.pars.fm)],
                      ["Output format", "-of", str(self.pars.of)],
                      ["Merge features", "-mf", str(self.pars.mf)]]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, clean=True)
        # html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(self.pars.o, "parameters.html"))

    def gen_html_genes(self, align = 50, nonDE=False):

        dir_name = os.path.basename(self.pars.o)
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

        if self.input.dna.scores:
            if self.pars.scoreh and self.input.dna.de_gene:
                score_header = [self.input.dna.de_gene.cond[0]]
            else:
                # if "(" in self.scores.values()[0]:
                #     score_header = ["Fold_change", "Filtered"]
                # else:
                score_header = ["Fold Change Score"]
            header_listp = ["#", "Promoter", "Gene", "TTSs Count", "TTS coverage"]
            header_listp += score_header
            header_listp += ["Sum of Ranks"]

            header_titlesp = ["", "Target promoters", "Gene symbol",
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites"]
            header_titlesp += ["Scores loaded by their absolute values from gene list or BED input. If there is annotation error for the gene names, it shows zero instead."] * len(
                score_header)
            header_titlesp += ["Sum up the ranks from left-hand side columns"]

        else:
            header_listp = ["#", "Promoter", "Gene", "TTSs Count", "TTS coverage", "Sum of Ranks"]

            header_titlesp = ["", "Target promoters", "Gene symbol",
                              "Number of DNA Binding sites locating within the promoter",
                              "The proportion of promoter covered by binding sites",
                              "Sum up the ranks from left-hand side columns"]

        html.add_heading("Target promoters")
        data_table = []

        if not self.input.dna.target_regions.sorted: self.input.dna.target_regions.sort()
        # Iterate by each gene promoter

        # Calculate the ranking
        rank_count = len(self.input.dna.target_regions) - rank_array(
            [self.stat.promoter["de"]["dbs"][p.toString()] for p in self.input.dna.target_regions])
        rank_coverage = len(self.input.dna.target_regions) - rank_array(
            [self.stat.promoter["de"]["dbs_coverage"][p.toString()] for p in self.input.dna.target_regions])

        if self.input.dna.scores:

            rank_score = len(self.input.dna.target_regions) - rank_array([self.input.dna.scores[p.toString()] for p in self.input.dna.target_regions])
            rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

        for i, promoter in enumerate(self.input.dna.target_regions):

            if self.stat.promoter["de"]["dbs"][promoter.toString()] == 0:
                dbssount = str(0)
            else:
                dbssount = '<a href="promoters_dbds.html#' + promoter.toString() + '" style="text-align:left">' + \
                           str(self.stat.promoter["de"]["dbs"][promoter.toString()]) + '</a>'

            region_link = region_link_internet(self.pars.organism, promoter)

            # try:
            #     gn = self.ensembl2symbol[promoter.name]
            #     if not gn: gn = promoter.name
            # except:
            gn = promoter.name

            self.ranktable[gn] = str(int(rank_sum[i]))
            self.dbstable[gn] = str(int(self.stat.promoter["de"]["dbs"][promoter.toString()]))

            newline = [str(i + 1),
                       region_link,
                       split_gene_name(gene_name=gn, org=self.pars.organism),
                       dbssount,
                       value2str(self.stat.promoter["de"]["dbs_coverage"][promoter.toString()])
                       ]

            if self.input.dna.scores:
                # if multiple_scores:
                #     for j in range(len(score_ar[0])):
                #         newline += [value2str(abs(score_ar[i][j]))]
                # else:
                newline += [value2str(abs(self.input.dna.scores[promoter.toString()]))]

            newline += ["<i>" + str(int(rank_sum[i])) + "</i>"]
            # print(newline)
            data_table.append(newline)

        # print(data_table)
        data_table = natsort.natsorted(data_table, key=lambda x: x[-1])
        html.add_zebra_table(header_listp, col_size_list, type_list, data_table, align=align, cell_align="left",
                             header_titles=header_titlesp, border_list=None, sortable=True, clean=True)
        html.add_heading("Notes")
        html.add_list(["TTS stands for triplex target DNA site.",
                       "TTS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "promoters.html"))

        ############################
        # Subpages for promoter centered page
        # promoters_dbds.html
        header_sub = ["#", "TFO", "TTS", "Strand", "Score", "Motif", "Orientation", "Sequence"]
        header_titles = ["", "RNA Binding Site", "DNA Binding Site", "Strand of TTS on DNA",
                         "Score of binding event", "Motif of binding by triple helix rule",
                         "Orientation of interaction between DNA and RNA. 'P'- Parallel; 'A'-Antiparallel", "Binding Sequence between DNA and RNA"]
        header_list = header_sub
        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for i, promoter in enumerate(self.input.dna.target_regions):
            if self.stat.promoter["de"]["dbs"][promoter.toString()] == 0:
                continue
            else:
                gn = promoter.name
                html.add_heading(split_gene_name(gene_name=gn, org=self.pars.organism), idtag=promoter.toString())

                data_table = []

                for j, rd in enumerate(self.stat.promoter["de"]["rd"][promoter.toString()]):
                    rbs = rd.rna.str_rna(pa=False)
                    for rbsm in self.stat.sig_DBD:
                        # rbsm = rbsm.partition(":")[2].split("-")
                        if rd.rna.overlap(rbsm):
                            rbs = "<font color=\"red\">" + rbs + "</font>"
                    data_table.append([str(j + 1), rbs, rd.dna.toString(space=True),
                                       rd.dna.orientation, rd.score, rd.motif, rd.orient,
                                       '<pre><font size="2">' + "\n".join(rd.match) + "</font></pre>"])

                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     header_titles=header_titles, sortable=True, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "promoters_dbds.html"))

        ##############################################################################################
        # spromoters.html    for significant promoters
        spromoters = GenomicRegionSet("sig_promoters")
        self.stat.promoter["de"]["sig_dbs"] = {}
        self.stat.promoter["de"]["sig_dbs_coverage"] = {}
        for promoter in self.input.dna.target_regions:
            # for rd in self.promoter["de"]["rd"][promoter.toString()]:
            #     if rd.rna
            sig_bindings = self.stat.promoter["de"]["rd"][promoter.toString()].overlap_rbss(rbss=self.stat.sig_DBD)
            dbs = sig_bindings.get_dbs()
            if len(dbs) > 0:
                spromoters.add(promoter)
                m_dbs = dbs.merge(w_return=True)
                self.stat.promoter["de"]["sig_dbs"][promoter.toString()] = len(dbs)
                # self.promoter["de"]["merged_dbs"][promoter.toString()] = len(m_dbs)
                self.stat.promoter["de"]["sig_dbs_coverage"][promoter.toString()] = float(m_dbs.total_coverage()) / len(promoter)

        html = Html(name=html_header, links_dict=self.link_d,  # fig_dir=os.path.join(directory,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        # Select promoters in sig DBD
        # spromoters = self.promoter["de"]["sig_dbs"].keys()
        if len(spromoters) == 0:
            html.add_heading("There is no significant DBD.")
        else:
            html.add_heading("Target promoters bound by significant DBD")
            # for rbsm in self.sig_region_promoter:
            #     spromoters = spromoters + [p for p in self.tpx_de.merged_dict[rbsm]]
            # spromoters = list(set(spromoters))
            data_table = []

            # Iterate by each gene promoter

            # Calculate the ranking
            rank_count = len(spromoters) - rank_array([self.stat.promoter["de"]["sig_dbs"][p.toString()] for p in spromoters])
            rank_coverage = len(spromoters) - rank_array(
                [self.stat.promoter["de"]["sig_dbs_coverage"][p.toString()] for p in spromoters])

            if self.input.dna.scores:
                rank_score = len(spromoters) - rank_array([self.input.dna.scores[p.toString()] for p in spromoters ] )
                rank_sum = [x + y + z for x, y, z in zip(rank_count, rank_coverage, rank_score)]

            else:
                rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]

            for i, promoter in enumerate(spromoters):
                # try:
                #     gn = self.ensembl2symbol[promoter.name]
                # except:
                gn = promoter.name
                dbssount = '<a href="promoters_dbds.html#' + promoter.toString() + \
                           '" style="text-align:left">' + \
                           str(self.stat.promoter["de"]["sig_dbs"][promoter.toString()]) + '</a>'

                region_link = region_link_internet(self.pars.organism, promoter)

                self.ranktable[gn] = str(int(rank_sum[i]))
                self.dbstable[gn] = str(self.stat.promoter["de"]["sig_dbs"][promoter.toString()])

                newline = [str(i + 1),
                           region_link,
                           split_gene_name(gene_name=gn, org=self.pars.organism),
                           dbssount,
                           value2str(self.stat.promoter["de"]["sig_dbs_coverage"][promoter.toString()])
                           ]
                if self.input.dna.scores:
                    newline += [value2str(self.input.dna.scores[promoter.toString()])]

                newline += [str(int(rank_sum[i]))]
                # print(newline)
                data_table.append(newline)

            data_table = natsort.natsorted(data_table, key=lambda x: x[-1])
            html.add_zebra_table(header_listp, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_titlesp, border_list=None, sortable=True, clean=True)
            html.add_heading("Notes")
            html.add_list(["TTS stands for triplex target DNA site.",
                           "TTS coverage is the proportion of the promoter where has potential to form triple helices with the given RNA."])
            html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "spromoters.html"))

    def save_table(self, path, table, filename):
        """Save the summary rank into the table for heatmap"""

        table_path = os.path.join(path, filename)
        # self.ranktable = {}
        rank_table = []

        if os.path.isfile(table_path):
            # Table exists
            f = open(table_path)
            for i, line in enumerate(f):
                line = line.strip().split()
                if i == 0:
                    # if not line:
                    #     break
                    #     exist = False
                    if self.pars.rn in line:
                        # lncRNA exists
                        exist = True
                        ind_rna = line.index(self.pars.rn)
                        header = line
                    else:
                        exist = False
                        line.append(self.pars.rn)
                        header = line

                else:
                    if exist and ind_rna:
                        # try:
                        line[ind_rna + 1] = table[line[0]]
                        rank_table.append(line)
                        # except:
                        #     rank_table.append(line)
                    else:
                        # try:
                        line.append(table[line[0]])
                        rank_table.append(line)
                        # except:
                        #     rank_table.append(line)
            f.close()

        else:
            # Table not exists
            header = ["gene", self.pars.rn]
            for k, v in table.items():
                rank_table.append([k, v])

        # Write into file
        g = open(table_path, "w")

        print("\t".join(header), file=g)
        for l in rank_table:
            # print(l)
            print("\t".join(l), file=g)
        g.close()


    def gen_html_regiontest(self):
        """Generate the HTML file"""
        align = 50
        dir_name = os.path.basename(self.pars.o)
        html_header = "Genomic Region Test: " + dir_name
        link_ds = OrderedDict()
        link_ds["RNA"] = "index.html"
        link_ds["Sig Target Regions"] = "starget_regions.html"
        link_ds["Target Regions"] = "target_regions.html"
        link_ds["Autobinding"] = "autobinding.html"
        link_ds["Parameters"] = "parameters.html"

        ##################################################
        # index.html

        html = Html(name=html_header, links_dict=link_ds,
                    fig_dir=os.path.join(os.path.dirname(self.pars.o),"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        # Plots
        html.add_figure(self.pars.rn+"_lineplot.png", align="left", width="45%", more_images=[self.pars.rn+"_boxplot.png"])
        if self.pars.showdbs:
            html.add_figure("lineplot_dbs.png", align="left", width="45%", more_images=["boxplot_dbs.png"])

        if self.pars.showdbs:
            header_list = [["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics",
                            "Target Regions", "Non-target Regions", None, "Statistics"],
                           ["", "", "with TTS", "without TTS", "with TTS (average)", "s.d.", "<i>p</i>-value",
                            "NO. TTSs", "NO. TTSs (average)", "s.d.", "<i>p</i>-value"]]
            header_titles = [["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                              "Regions from randomization", None, "Statistics based on target regions",
                              "Given target regions on DNA", "Regions from randomization", None,
                              "Statistics based on DNA Binding Sites"],
                             ["", "",
                              "Number of target regions with triplex binding",
                              "Number of target regions without triplex binding",
                              "Average number of regions from randomization with triplex binding",
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
            header_list = [["#", "DBD", "Target Regions", None, "Non-target Regions", None, "Statistics", None, "Autobinding"],
                           ["", "", "with TTS", "without TTS", "with TTS (average)", "s.d.", "<i>p</i>-value",
                            "z-score", "Number"]]
            header_titles = [["Rank", "DNA Binding Domain", "Given target regions on DNA", None,
                              "Regions from randomization", None, "Statistics based on target regions", None,
                              "Regions bind to themselves"],
                             ["", "",
                              "Number of target regions with triplex binding",
                              "Number of target regions without triplex binding",
                              "Average number of regions from randomization with triplex binding",
                              "Standard deviation", "P value", "Z-score", ""]]
            border_list = [" style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"", "",
                           " style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"",
                           " style=\"border-right:1pt solid gray\"", ""]

        type_list = 'ssssssssssssssss'
        col_size_list = [50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50]
        data_table = []

        for i, rbs in enumerate(self.stat.rbss):
            if self.stat.data["region"]["p"][i] < self.pars.a:
                p_region = "<font color=\"red\">" + value2str(self.stat.data["region"]["p"][i]) + "</font>"

            else:
                p_region = value2str(self.stat.data["region"]["p"][i])
            zs = (self.stat.counts_tr[rbs][0] - self.stat.data["region"]["ave"][i]) / self.stat.data["region"]["sd"][i]
            new_line = [str(i + 1),
                        rbs.str_rna(pa=False),
                        '<a href="dbd_region.html#' + rbs.str_rna() +
                        '" style="text-align:left">' + str(self.stat.counts_tr[rbs][0]) + '</a>',
                        str(self.stat.counts_tr[rbs][1]),
                        value2str(self.stat.data["region"]["ave"][i]),
                        value2str(self.stat.data["region"]["sd"][i]),
                        p_region,
                        value2str(zs),
                        '<a href="autobinding.html">' +
                        str(len(self.triplexes.autobinding.merged_dict[rbs])) + '</a>'
                        ]
            if self.pars.showdbs:
                if self.stat.data["dbs"]["p"][i] < self.pars.a:
                    p_dbs = "<font color=\"red\">" + value2str(self.stat.data["dbs"]["p"][i]) + "</font>"
                else:
                    p_dbs = value2str(self.stat.data["dbs"]["p"][i])

                new_line += [str(self.stat.counts_dbs[rbs]),
                             value2str(self.stat.data["dbs"]["ave"][i]),
                             value2str(self.stat.data["dbs"]["sd"][i]),
                             p_dbs]
            data_table.append(new_line)

        data_table = natsort.natsorted(data_table, key=lambda x: x[6])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, border_list=border_list, sortable=True)

        html.add_heading("Notes")
        html.add_list(["RNA name: " + self.pars.rn,
                       "Randomization is performed for " + str(self.pars.n) + " times.",
                       "DBD stands for DNA Binding Domain on RNA.",
                       "TTS stands for triplex target site."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "index.html"))

        #############################################################
        # RNA subpage: Profile of targeted regions for each merged DNA Binding Domain
        #############################################################

        header_list = ["#", "Target Region",
                       "Associated Gene",
                       "No. of TTSs",
                       "TTS coverage"]
        header_titles = ["Rank", "Given target regions from BED files",
                         "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                         "Number of DNA Binding Sites locate within the region",
                         "The proportion of the region covered by triplex binding"]

        #########################################################
        # dbd_region.html
        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for rbsm in self.stat.rbss:
            html.add_heading("DNA Binding Domain: " + rbsm.str_rna(),
                             idtag=rbsm.str_rna())
            data_table = []
            for i, region in enumerate(self.stat.tpx.merged_dict[rbsm]):
                # Add information
                data_table.append([str(i + 1),
                                   '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.pars.organism +
                                   "&position=" + region.chrom + "%3A" + str(region.initial) + "-" + str(region.final) +
                                   '" style="text-align:left">' + region.toString(space=True) + '</a>',
                                   split_gene_name(gene_name=region.name, org=self.pars.organism),
                                   str(len(self.stat.region_dbs[region.toString()])),
                                   value2str(self.stat.region_coverage[region.toString()])
                                   ])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 auto_width=True, header_titles=header_titles, sortable=True, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "dbd_region.html"))

        #############################################################
        # Targeted regions centered
        #############################################################

        ##############################################################################################
        # target_regions.html
        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        if self.pars.score:
            header_list = ["#", "Target region", "Associated Gene", "TTSs Count", "Norm. TTSs",
                           "TTS coverage", "Score", "Sum of ranks"]
            header_titles = ["Rank",
                             "Target regions loaded from the given BED file",
                             "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                             "Number of DNA Binding Sites within the region",
                             "Normalized Number of DNA Binding Sites within the region (per 1000 bp)",
                             "The proportion of the region covered by TTS binding",
                             "Scores from BED file",
                             "Sum of all the left-hand-side ranks"]
        else:
            header_list = ["#", "Target region", "Associated Gene", "TTSs Count", "Norm. TTSs",
                           "TTS coverage", "Sum of ranks"]
            header_titles = ["Rank",
                             "Target regions loaded from the given BED file",
                             "Associated genes which is overlapping with the given region or close to it (less than 50000 bp)",
                             "Number of DNA Binding Sites within the region",
                             "Normalized Number of DNA Binding Sites within the region (per 1000 bp)",
                             "The proportion of the region covered by triplex binding",
                             "Sum of all the left-hand-side ranks"]
        html.add_heading("Target Regions")
        data_table = []

        if not self.input.dna.target_regions.sorted: self.input.dna.target_regions.sort()

        # Calculate the ranking
        rank_count = len(self.input.dna.target_regions) - rank_array([len(self.stat.region_dbs[p.toString()]) for p in self.input.dna.target_regions])
        rank_normcount = len(self.input.dna.target_regions) - rank_array([self.stat.region_normdbs[p.toString()] for p in self.input.dna.target_regions])
        rank_coverage = len(self.input.dna.target_regions) - rank_array([self.stat.region_coverage[p.toString()] for p in self.input.dna.target_regions])

        if self.pars.score:
            try:
                score_list = [float(p.data.split("\t")[0]) for p in self.input.dna.target_regions]
                rank_score = len(self.input.dna.target_regions) - rank_array([abs(s) for s in score_list])
                rank_sum = [x + y + z for x, y, z in zip(rank_normcount, rank_coverage, rank_score)]
                # sum_rank = rank_array(rank_sum)  # method='min'
            except ImportError:
                print("There is no score in BED file, please don't use '-score' argument.")
        else:
            rank_sum = [x + y for x, y in zip(rank_count, rank_coverage)]
            sum_rank = rank_array(rank_sum)

        for i, region in enumerate(self.input.dna.target_regions):
            dbs_counts = str(len(self.stat.region_dbs[region.toString()]))
            dbs_norm = "{0:.2f}".format(round(self.stat.region_normdbs[region.toString()],2))
            dbs_cover = value2str(self.stat.region_coverage[region.toString()])

            newline = [str(i + 1),
                       '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.pars.organism +
                       "&position=" + region.chrom + "%3A" + str(region.initial) + "-" + str(region.final) +
                       '" style="text-align:left">' + region.toString(space=True) + '</a>',
                       split_gene_name(gene_name=region.name, org=self.pars.organism),
                       '<a href="region_dbs.html#' + region.toString() +
                       '" style="text-align:left">' + dbs_counts + '</a>',
                       dbs_norm, dbs_cover]

            if self.pars.score:
                dbs_score = value2str(score_list[i])
                region.data = "\t".join([dbs_counts, dbs_norm, dbs_cover, dbs_score, str(rank_sum[i])])
                newline.append(dbs_score)
                newline.append(str(rank_sum[i]))
            else:
                region.data = "\t".join([dbs_counts, dbs_norm, dbs_cover, str(rank_sum[i])])
                newline.append(str(rank_sum[i]))
            data_table.append(newline)

        data_table = natsort.natsorted(data_table, key=lambda x: x[-1])
        # data_table = sorted(data_table, key=lambda x: x[-1])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, header_titles=header_titles, sortable=True, clean=True)
        html.add_heading("Notes")
        html.add_list(["All target regions without any bindings are ignored."])
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "target_regions.html"))

        self.input.dna.target_regions.sort_score()
        self.input.dna.target_regions.write(os.path.join(self.pars.o,  self.pars.rn+"_target_regions.bed"))



        ##############################################################################################
        # starget_regions.html    for significant target regions

        stargets = GenomicRegionSet("sig_targets")
        sig_dbs = {}
        sig_normdbs = {}
        sig_dbs_coverage = {}
        for i, r in enumerate(self.input.dna.target_regions):
            sig_bindings = self.stat.region_dbs[r.toString()].overlap_rbss(rbss=self.stat.data["region"]["sig_region"])
            dbs = sig_bindings.get_dbs()
            if len(dbs) > 0:
                stargets.add(r)
                m_dbs = dbs.merge(w_return=True)
                sig_dbs[r] = len(dbs)
                sig_normdbs[r] = len(dbs) * 1000 / len(r)
                # self.promoter["de"]["merged_dbs"][promoter.toString()] = len(m_dbs)
                sig_dbs_coverage[r] = float(m_dbs.total_coverage()) / len(r)

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        # Select promoters in sig DBD
        if len(self.stat.data["region"]["sig_region"]) == 0:
            html.add_heading("There is no significant DBD.")
        else:
            html.add_heading("Target regions bound by significant DBD")
            data_table = []
            # Calculate the ranking
            # rank_count = len(stargets) - rank_array([sig_dbs[p] for p in stargets])
            rank_normcount = len(stargets) - rank_array([sig_normdbs[p] for p in stargets])
            rank_coverage = len(stargets) - rank_array([sig_dbs_coverage[p] for p in stargets])
            if self.pars.score:
                score_list = [float(p.data.split("\t")[0]) for p in stargets]
                rank_score = len(stargets) - rank_array([abs(s) for s in score_list])
                rank_sum = [x + y + z for x, y, z in zip(rank_normcount, rank_coverage, rank_score)]
                sum_rank = rank_array(rank_sum)  # method='min'
            else:
                rank_sum = [x + y for x, y in zip(rank_normcount, rank_coverage)]
                sum_rank = rank_array(rank_sum)

            for i, region in enumerate(stargets):
                dbssount = '<a href="region_dbs.html#' + region.toString() + \
                           '" style="text-align:left">' + str(sig_dbs[region]) + '</a>'

                region_link = region_link_internet(self.pars.organism, region)

                newline = [str(i + 1), region_link,
                           split_gene_name(gene_name=region.name, org=self.pars.organism),
                           dbssount, "{0:.2f}".format(round(sig_normdbs[region],2)),
                           value2str(sig_dbs_coverage[region]) ]
                if self.pars.score:
                    dbs_score = value2str(score_list[i])
                    # region.data = "\t".join([dbs_counts, dbs_cover, dbs_score, str(sum_rank[i])])
                    newline.append(dbs_score)
                    newline.append(str(rank_sum[i]))
                    # print([dbs_score, str(sum_rank[i])])
                else:
                    # region.data = "\t".join([dbs_counts, dbs_cover, str(sum_rank[i])])
                    newline.append(str(rank_sum[i]))

                # newline += ["<i>" + str(rank_sum[i]) + "</i>"]
                # print(newline)
                data_table.append(newline)

            # print(data_table)
            # data_table = sorted(data_table, key=lambda x: x[-1])
            data_table = natsort.natsorted(data_table, key=lambda x: x[-1])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                 header_titles=header_titles, border_list=None, sortable=True, clean=True)
            html.add_heading("Notes")
            html.add_list(["TTS stands for triplex target DNA site.",
                           "TTS coverage is the proportion of the region where has potential to form triple helices with the given RNA."])
            html.add_fixed_rank_sortable()
            html.write(os.path.join(self.pars.o, "starget_regions.html"))

        ############################
        # Subpages for targeted region centered page
        # region_dbs.html
        header_list = ["TFO", "TTS", "Strand", "Score", "Motif", "Orientation", "Sequence"]
        header_titles = ["", "RNA Binding Site", "DNA Binding Site", "Strand of TTS on DNA",
                         "Score of binding event", "Motif of binding by triple helix rule",
                         "Orientation of interaction between DNA and RNA. 'P'- Parallel; 'A'-Antiparallel",
                         "Binding Sequence between DNA and RNA"]

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        for i, region in enumerate(self.input.dna.target_regions):
            if len(self.stat.region_dbs[region.toString()]) == 0:
                continue
            else:
                html.add_heading("Associated gene: " + split_gene_name(gene_name=region.name, org=self.pars.organism),
                                 idtag=region.toString())
                # html.add_free_content(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.pars.organism +
                #                        "&position=" + region.chrom + "%3A" + str(region.initial) +
                #                        "-" + str(region.final) + '" style="margin-left:50">' +
                #                        region.toString(space=True) + '</a>'])
                data_table = []
                for rd in self.stat.region_dbs[region.toString()]:
                    rbs = rd.rna.str_rna(pa=False)
                    for rbsm in self.stat.data["region"]["sig_region"]:
                        # rbsm = rbsm.partition(":")[2].split("-")
                        if rd.rna.overlap(rbsm):
                            rbs = "<font color=\"red\">" + rbs + "</font>"
                    data_table.append([rbs,
                                       '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=' + self.pars.organism +
                                       "&position=" + rd.dna.chrom + "%3A" + str(rd.dna.initial) + "-" + str(
                                           rd.dna.final) +
                                       '" style="text-align:left">' + rd.dna.toString(space=True) + '</a>',
                                       rd.dna.orientation, rd.score, rd.motif, rd.orient,
                                       '<pre><font size="2">' + "\n".join(rd.match) + "</font></pre>"])
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                                     header_titles=header_titles, auto_width=True, clean=True)
        html.write(os.path.join(self.pars.o, "region_dbs.html"))

        ################################################################
        ############# Autobinding
        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")

        html.add_heading("Autobinding")
        header_list = ["#", "DBD", "RNA", "DNA", "Score", "Motif", "Orientation", "Sequence"]
        t = []
        for rbs in self.stat.rbss:
            for i, rd in enumerate(self.triplexes.autobinding):
                if rbs.overlap(rd.rna):
                    t.append([str(i), rbs.str_rna(pa=False),
                              str(rd.rna.initial) + "-" + str(rd.rna.final),
                              str(rd.dna.initial) + "-" + str(rd.dna.final),
                              rd.score, rd.motif, rd.orient,
                              '<pre><font size="1">' + "\n".join(rd.match) + "</font></pre>"])
        if len(t) > 0:
            html.add_zebra_table(header_list, col_size_list, type_list, t, align=align, cell_align="left",
                                 sortable=True, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(self.pars.o, "autobinding.html"))

        ###############################################################################33
        ################ Parameters.html

        html = Html(name=html_header, links_dict=link_ds,  # fig_dir=os.path.join(self.pars.o,"style"),
                    fig_rpath="../style", RGT_header=False, other_logo="TDF", homepage="../index.html")
        html.add_heading("Parameters")
        header_list = ["Description", "Arguments", "Value"]

        data_table = [["RNA sequence name", "-rn", self.pars.rn],
                      ["Input RNA sequence file", "-r", os.path.basename(self.pars.r)],
                      ["Input BED file", "-bed", os.path.basename(self.pars.bed)],
                      ["Output directory", "-o", os.path.basename(self.pars.o)],
                      ["Organism", "-organism", self.pars.organism],
                      ["Number of repitetion of andomization", "-n", str(self.pars.n)],
                      ["Alpha level for rejection p value", "-a", str(self.pars.a)],
                      ["Cut off value for filtering out the DBD with low counts of TTSs", "-ccf", str(self.pars.ccf)],
                      ["Remove temporary files", "-rt", str(self.pars.rt)],
                      ["Input BED file for masking in randomization", "-f", str(self.pars.f)],
                      # ["Input file for RNA accecibility", "-ac", str(self.pars.ac)],
                      # ["Cut off value for RNA accecibility", "-accf", str(self.pars.accf)],
                      ["Output the BED files for DNA binding sites.", "-obed", str(self.pars.obed)],
                      ["Show parallel and antiparallel bindings in the plot separately.", "-showpa",
                       str(self.pars.showpa)],
                      ["Minimum length", "-l", str(self.pars.l)],
                      ["Maximum error rate", "-e", str(self.pars.e)],
                      ["Tolerated number of consecutive errors", "-c", str(self.pars.c)],
                      ["Filtering repeats", "-fr", str(self.pars.fr)],
                      ["Filtering mode", "-fm", str(self.pars.fm)],
                      ["Output format", "-of", str(self.pars.of)],
                      ["Merge features", "-mf", str(self.pars.mf)]]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left",
                             auto_width=True, clean=True)
        # html.add_free_content(['<a href="summary.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(self.pars.o, "parameters.html"))
