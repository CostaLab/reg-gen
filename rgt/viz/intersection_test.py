# Python Libraries


import itertools
import multiprocessing.pool
import matplotlib.ticker as mtick
from matplotlib_venn import venn3

# Local Libraries
# Distal Libraries
from .shared_function import *

# Local test

###########################################################################################
#                    Inersection test
###########################################################################################


def posi2set(regions, p):
    alln = list(range(len(regions)))
    inter_r = copy.deepcopy(regions[p[0]])

    for i in alln:
        # print("inter_r: "+inter_r.name)
        if i in p[1:]:
            inter_r = inter_r.intersect(regions[i], mode=OverlapType.OVERLAP)
        elif i == p[0]:
            pass
        else:
            inter_r = inter_r.subtract(regions[i], whole_region=False)
    # print("inter_r: "+inter_r.name)
    return inter_r


def posi2region(regions, p):
    # all = range(len(regions))
    new_r = GenomicRegionSet(name="")
    for r in p:
        new_r.combine(regions[r])
    return new_r


class Intersect:
    def __init__(self, reference_path, query_path, mode_count, organism):
        self.rEM, self.qEM = ExperimentalMatrix(), ExperimentalMatrix()
        self.rEM.read(reference_path)
        self.rEM.remove_empty_regionset()
        self.references = self.rEM.get_regionsets()
        self.referencenames = self.rEM.get_regionsnames()
        self.qEM.read(query_path)
        self.qEM.remove_empty_regionset()
        self.query = self.qEM.get_regionsets()
        self.querynames = self.qEM.get_regionsnames()
        self.mode_count = mode_count
        self.organism = organism
        self.sbar = None
        self.pbar = None
        self.test_d = None
        # self.background = None

    def background(self, path=None):
        """Given a bed file as the background for analysis"""
        # bgbed = GenomicRegionSet(name="Background")
        # if path:
        #     bgbed.read(path)
        #     # nq = []
        #     print("\tTrimming the queries by the given background: "+path)

        #     rlist = [ r.trim_by(background=bgbed) for r in self.references]
        #     self.references = rlist
        #     qlist = [ q.trim_by(background=bgbed) for q in self.query]
        #     self.query = qlist

        # self.background = bgbed

        if path:
            bg = GenomicRegionSet("background")
            bg.read(path)
            self.background = bg
            for ty in list(self.groupedreference.keys()):
                # self.background[ty] = bg
                rlist = [r.trim_by(background=bg) for r in self.groupedreference[ty]]
                self.groupedreference[ty] = rlist
                qlist = [q.trim_by(background=bg) for q in self.groupedquery[ty]]
                self.groupedquery[ty] = qlist

    def group_refque(self, groupby):
        self.groupedreference, self.groupedquery = group_refque(self.rEM, self.qEM,
                                                                groupby,
                                                                self.references,
                                                                self.query)
        remove_duplicates(self.groupedreference)
        remove_duplicates(self.groupedquery)

    def colors(self, colorby, definedinEM, ref_que="que"):
        """color_list is a Dict [query] : color """
        if ref_que == "que":
            self.color_list = color_groupded_region(self.qEM, self.groupedquery, colorby, definedinEM)
            if list(self.groupedquery.keys())[0] == "":
                self.color_tags = [n.name for n in self.groupedquery[""]]
            else:
                self.color_tags = gen_tags(self.qEM, colorby)
        elif ref_que == "ref":
            self.color_list = color_groupded_region(self.rEM, self.groupedreference, colorby, definedinEM)
            if list(self.groupedquery.keys())[0] == "":
                self.color_tags = [n.name for n in self.groupedquery[""]]
            else:
                self.color_tags = gen_tags(self.qEM, colorby)

    def colors_comb(self):
        """color_list is a list : color """

        if list(self.groupedquery.keys())[0] == "":
            self.color_tags = self.referencenames
        else:
            tags = []
            for t in [n.name for n in list(self.groupedreference.values())[0]]:
                nt = t.replace(list(self.groupedreference.keys())[0], "")
                nt = nt.replace("_", "")
                tags.append(nt)
            self.color_tags = tags
        self.color_list = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(self.color_tags))).tolist()

    def extend_ref(self, percentage):
        """percentage must be positive value"""
        for ty in self.groupedreference:
            for r in self.groupedreference[ty]:
                r.extend(left=percentage, right=percentage, percentage=False)
                r.merge()

    def count_intersect(self, threshold, frequency=True):

        self.counts = OrderedDict()
        self.rlen, self.qlen = {}, {}
        self.nalist = []

        if frequency:
            self.frequency = OrderedDict()

        # if self.mode_count == "bp":
        #    print2(self.parameter, "\n{0}\t{1}\t{2}\t{3}\t{4}".format("Reference","Length(bp)", "Query", "Length(bp)", "Length of Intersection(bp)"))
        # elif self.mode_count == "count":
        #    print2(self.parameter, "\n{0}\t{1}\t{2}\t{3}\t{4}".format("Reference","sequence_number", "Query", "sequence_number", "Number of Intersection"))

        for ty in list(self.groupedreference.keys()):
            self.counts[ty] = OrderedDict()
            self.rlen[ty], self.qlen[ty] = OrderedDict(), OrderedDict()
            if frequency:
                self.frequency[ty] = OrderedDict()

            for r in self.groupedreference[ty]:
                if r.total_coverage() == 0 and len(r) > 0:
                    self.nalist.append(r.name)
                    continue
                else:
                    self.counts[ty][r.name] = OrderedDict()
                    if self.mode_count == "bp":
                        rlen = r.total_coverage()
                    elif self.mode_count == "count":
                        rlen = len(r)
                    self.rlen[ty][r.name] = rlen

                    mp_input = []
                    for q in self.groupedquery[ty]:
                        if r.name == q.name:
                            continue
                        else:
                            mp_input.append([q, self.nalist, self.mode_count, self.qlen, threshold,
                                             self.counts, frequency, self.frequency, ty, r])
                    # q, nalist, mode_count, qlen_dict, threshold, counts, frequency, self_frequency, ty, r
                    pool = multiprocessing.Pool(processes=multiprocessing.cpu_count() - 1)
                    mp_output = pool.map(mp_count_intersect, mp_input)
                    pool.close()
                    pool.join()

                    # qname, nalist, qlen_dict[ty][q.name], counts[ty][r.name][q.name], self_frequency[ty][q.name].append(c[2])
                    for output in mp_output:
                        if output[1]:
                            self.nalist.append(output[1])
                        else:
                            self.qlen[ty][output[0]] = output[2]
                            self.counts[ty][r.name][output[0]] = output[3]
                            # print(r.name)
                            # print(output[0])
                            # print(output[3])
                            try:
                                self.frequency[ty][output[0]][r.name] = output[3][2]
                            except:
                                self.frequency[ty][output[0]] = {}
                                self.frequency[ty][output[0]][r.name] = output[3][2]
                                # print2(self.parameter, "{0}\t{1}\t{2}\t{3}\t{4}".format(r.name,rlen, q.name, qlen, c[2]))
                                # print(self.nalist)
                                # print(self.qlen)
                                # print(self.counts)
                                # print(self.frequency)
                                # sys.stdout.flush()

    def barplot(self, logt=False, percentage=False):
        f, axs = plt.subplots(len(list(self.counts.keys())), 1)
        f.subplots_adjust(left=0.3)
        self.xtickrotation, self.xtickalign = 0, "center"
        # if len(axs) == 1: axs = [axs]
        try:
            axs = axs.reshape(-1)
        except:
            axs = [axs]

        for ai, ax in enumerate(axs):
            if logt:
                ax.set_yscale('log')
                plus = 1
            else:
                ax.locator_params(axis='y', nbins=4)
                plus = 0
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.set_title(list(self.counts.keys())[ai], y=1)

            r_label = []
            for ind_r, r in enumerate(list(self.counts.values())[ai].keys()):
                for l in self.references:
                    if l.name == r:
                        lr = len(l)

                if len(axs) == 1:
                    r_label.append(r)
                else:
                    try:
                        r_label.append(self.rEM.get_type(r, "factor"))
                    except:
                        r_label.append(r)
                if len(r_label[-1]) > 15 or len(list(self.counts.values())[ai][r].keys()) * len(
                        list(self.counts.values())[ai].keys()) > 8:
                    self.xtickrotation, self.xtickalign = 50, "right"
                width = 0.8 / (len(list(self.counts.values())[ai][r].keys()) + 1)  # Plus one background
                for ind_q, q in enumerate(list(self.counts.values())[ai][r].keys()):
                    x = ind_r + ind_q * width + 0.1
                    if percentage:
                        # print(lr)
                        y = 100 * (list(self.counts.values())[ai][r][q][2] + plus) / lr
                    else:
                        y = list(self.counts.values())[ai][r][q][2] + plus  # intersect number

                    ax.bar(x, y, width=width, color=self.color_list[q], edgecolor="none", align='edge', log=logt)

            ax.yaxis.tick_left()
            ax.set_xticks([i + 0.5 - 0.5 * width for i in range(len(r_label))])
            ax.set_xticklabels(r_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=True)
            ax.set_xlim([0, len(list(self.counts.values())[ai].keys()) - 0.1])

            ax.legend(self.color_tags, loc='center left', handlelength=1, handletextpad=1,
                      columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            if percentage:
                fmt = '%.0f%%'  # Format you want the ticks, e.g. '40%'
                yticks = mtick.FormatStrFormatter(fmt)
                ax.yaxis.set_major_formatter(yticks)

                f.text(-0.025, 0.5, "Intersected percantage", rotation="vertical", va="center")
            else:
                f.text(-0.025, 0.5, "Intersected regions number", rotation="vertical", va="center")

        # f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        f.tight_layout()
        self.bar = f

    def stackedbar(self):
        f, axs = plt.subplots(len(list(self.counts.keys())), 1)
        f.subplots_adjust(left=0.3)
        # if len(axs) == 1: axs = [axs]
        try:
            axs = axs.reshape(-1)
        except:
            axs = [axs]

        for ai, ax in enumerate(axs):
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
            ax.locator_params(axis='y', nbins=2)
            ax.set_title(list(self.counts.keys())[ai], y=1)
            r_label = []

            for ind_r, r in enumerate(list(self.counts.values())[ai].keys()):
                if len(axs) == 1:
                    r_label.append(r)
                else:
                    try:
                        r_label.append(self.rEM.get_type(r, "factor"))
                    except:
                        r_label.append(r)
                width = 0.6
                bottom = 0
                reverse_q = list(list(self.counts.values())[ai][r].keys())[::-1]
                # ql = len(reverse_q) - 1
                rc = self.color_tags[::-1]
                for ind_q, q in enumerate(reverse_q):
                    x = ind_r
                    y = list(self.counts.values())[ai][r][q][2]  # intersect number
                    ax.bar(x, y, width=width, bottom=bottom, color=self.color_list[q],
                           edgecolor="none", align='center', label=rc[ind_q])
                    bottom = bottom + y
            ax.yaxis.tick_left()
            ax.set_xticks(list(range(len(r_label))))
            ax.set_xticklabels(r_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=True)
            ax.set_xlim([-0.5, ind_r + 0.5])

            # handles, labels = ax.get_legend_handles_labels()
            # nhandles = []
            # for t in self.color_tags:
            #    nhandles.append( handles[labels.index(t)] )
            # print(self.color_tags)
            handles, labels = ax.get_legend_handles_labels()
            by_label = OrderedDict(reversed(list(zip(labels, handles))))
            ax.legend(list(by_label.values()), list(by_label.keys()), loc='center left', handlelength=1, handletextpad=1,
                      columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            f.text(-0.025, 0.5, "Intersected regions number", rotation="vertical", va="center")

        # f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        f.tight_layout()
        self.sbar = f

    def percentagebar(self):
        self.color_list["No intersection"] = "0.7"
        f, axs = plt.subplots(len(list(self.counts.keys())), 1)
        f.subplots_adjust(left=0.3)
        # if len(axs) == 1: axs = [axs]
        try:
            axs = axs.reshape(-1)
        except:
            axs = [axs]
        self.percentage = []

        for ai, ax in enumerate(axs):
            # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.locator_params(axis='y', nbins=2)
            ax.set_title(list(self.counts.keys())[ai], y=1.05)
            r_label = []
            self.percentage.append({})
            for ind_r, r in enumerate(list(self.counts.values())[ai].keys()):
                if len(axs) == 1:
                    r_label.append(r)
                else:
                    try:
                        r_label.append(self.rEM.get_type(r, "factor"))
                    except:
                        r_label.append(r)
                width = 0.6
                bottom = 0
                if self.mode_count == "bp":
                    sumlength = self.rEM.objectsDict[r].total_coverage()
                elif self.mode_count == "count":
                    sumlength = len(self.rEM.objectsDict[r])
                self.percentage[ai][r] = {}

                if sumlength == 0:
                    for ind_q, q in enumerate(list(self.counts.values())[ai][r].keys()):
                        self.percentage[ai][r][q] = "ref is empty"
                else:
                    for ind_q, q in enumerate(list(self.counts.values())[ai][r].keys()):
                        x = ind_r
                        y = 100 * float(list(self.counts.values())[ai][r][q][2]) / sumlength  # percentage
                        ax.bar(x, y, width=width, bottom=bottom,
                               color=self.color_list[q], edgecolor="none", align='center')
                        bottom = bottom + y
                        self.percentage[ai][r][q] = y

                ax.bar(x, 100 - bottom, width=width, bottom=bottom,
                       color=self.color_list["No intersection"], align='center')

            ax.yaxis.tick_left()

            ax.set_xticks(list(range(len(r_label))))
            ax.set_xticklabels(r_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=True)
            ax.set_xlim([-0.5, ind_r + 0.5])
            ax.set_ylim([0, 100])

            legend_labels = self.color_tags + ["No intersection"]
            ax.legend(legend_labels, loc='center left', handlelength=1, handletextpad=1,
                      columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        f.text(-0.025, 0.5, "Proportion of intersected regions (%)", rotation="vertical", va="center")

        # f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        f.tight_layout()
        self.pbar = f

    def gen_html(self, directory, title, align, args):
        # fp = os.path.join(dir,outputname,title)
        # link_d = {title:"intersection.html"}
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = "Intersection Test: " + dir_name + "/" + title
        link_d = OrderedDict()
        link_d["Intersection test"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        html.add_figure("intersection_bar.png", align="center")
        html.add_figure("intersection_barp.png", align="center")
        if self.sbar:
            html.add_figure("intersection_stackedbar.png", align="center")

        header_list = ["#",
                       "Reference<br>name",
                       "Query<br>name",
                       "Reference<br>number",
                       "Query<br>number",
                       "Intersect.",
                       "Proportion<br>of Reference"]
        statistic_table = [["Reference_name", "Query_name", "Reference_number", "Query_number",
                            "Intersect.", "Proportion_of_Reference"]]
        if self.test_d:
            header_list += ["Average<br>intersect.", "Chi-square<br>statistic",
                            "Positive<br>Association<br>p-value", "Negative<br>Association<br>p-value"]
            statistic_table[0] += ["Average_intersect.", "Chi-square_statistic",
                                   "Positive_Association_p-value", "Negative_Association_p-value"]
        else:
            pass

        type_list = 'ssssssssssssssss'
        col_size_list = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]

        for ind_ty, ty in enumerate(self.counts.keys()):
            html.add_heading(ty, size=4, bold=False)
            data_table = []
            c = 0
            for ind_r, r in enumerate(self.counts[ty]):
                for ind_q, q in enumerate(self.counts[ty][r]):
                    if r == q:
                        continue
                    c += 1
                    pt = self.counts[ty][r][q][2] / self.rlen[ty][r]
                    internal_overlap = self.counts[ty][r][q][2]
                    if self.test_d:
                        aveinter = self.test_d[ty][r][q][0]
                        chisqua = value2str(self.test_d[ty][r][q][1])
                        pv = self.test_d[ty][r][q][2]
                        if isinstance(pv, str):
                            data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                               str(internal_overlap), "{:.2f}%".format(100 * pt),
                                               aveinter, chisqua, pv, "-"])
                        else:
                            npv = 1 - pv
                            if pv < 0.05:
                                if internal_overlap > aveinter:
                                    data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                                       str(internal_overlap), "{:.2f}%".format(100 * pt),
                                                       value2str(aveinter), chisqua,
                                                       "<font color=\"red\">" + value2str(pv) + "</font>",
                                                       value2str(npv)])
                                    statistic_table.append(
                                        [r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]), str(internal_overlap),
                                         "{:.2f}%".format(100 * pt),
                                         value2str(aveinter), chisqua, value2str(pv), value2str(npv)])
                                else:
                                    data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                                       str(internal_overlap), "{:.2f}%".format(100 * pt),
                                                       value2str(aveinter), chisqua, value2str(npv),
                                                       "<font color=\"red\">" + value2str(pv) + "</font>"])
                                    statistic_table.append(
                                        [r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]), str(internal_overlap),
                                         "{:.2f}%".format(100 * pt),
                                         value2str(aveinter), chisqua, value2str(npv), value2str(pv)])
                            elif self.test_d[ty][r][q][2] >= 0.05:
                                if internal_overlap > aveinter:
                                    data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                                       str(internal_overlap), "{:.2f}%".format(100 * pt),
                                                       value2str(aveinter), chisqua, value2str(pv), value2str(npv)])
                                    statistic_table.append(
                                        [r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]), str(internal_overlap),
                                         "{:.2f}%".format(100 * pt),
                                         value2str(aveinter), chisqua, value2str(pv), value2str(npv)])
                                else:
                                    data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                                       str(internal_overlap), "{:.2f}%".format(100 * pt),
                                                       value2str(aveinter), chisqua, value2str(npv), value2str(pv)])
                                    statistic_table.append(
                                        [r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]), str(internal_overlap),
                                         "{:.2f}%".format(100 * pt),
                                         value2str(aveinter), chisqua, value2str(npv), value2str(pv)])
                    else:
                        data_table.append([str(c), r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]),
                                           str(internal_overlap), "{:.2f}%".format(100 * pt)])
                        statistic_table.append([r, q, str(self.rlen[ty][r]), str(self.qlen[ty][q]), str(internal_overlap),
                                                "{:.2f}%".format(100 * pt)])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, sortable=True)
            output_array(statistic_table, directory=directory, folder=title, filename="statistics" + ty + ".txt")

        html.add_heading("Assumptions and hypothesis")
        list_ex = ['Positive association is defined by: True intersection number > Averaged random intersection.',
                   'Negative association is defined by: True intersection number < Averaged random intersection.']

        self.nalist = set(self.nalist)
        if len(self.nalist) > 0:
            list_ex.append(
                'The following region sets contain zero-length regions which cause error in intersection calculation, please check it:<br>' +
                '<font color=\"red\">' + ', '.join([s for s in self.nalist]) + '</font>')
        if self.test_d:
            list_ex.append('Randomly permutation for ' + str(self.test_time) + ' times.')
        html.add_list(list_ex)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))

        # Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        header_list = ["Description", "Argument", "Value"]
        data_table = [["Reference", "-r", args.r],
                      ["Query", "-q", args.q],
                      ["Output directory", "-o", os.path.basename(args.o)],
                      ["Experiment title", "-t", args.t],
                      # ["Grouping tag", "-g", args.g],
                      # ["Coloring tag", "-c", args.c],
                      # ["Background", "-bg", args.bg],
                      ["Organism", "-organism", args.organism]]

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_free_content([
            '<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(
            ['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory, title, "parameters.html"))

    def gen_html_comb(self, directory, title, align, args):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = "Combinatorial Test: " + dir_name + "/" + title
        link_d = OrderedDict()
        link_d["Combinatorial test"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d, fig_dir=os.path.join(directory, "style"),
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        # html.create_header()
        # html.add_heading(title)

        if self.sbar: html.add_figure("intersection_stackedbar.png", align="center")

        # html.add_free_content(['<p style=\"margin-left: '+str(align+150)+'">'+ '** </p>'])

        for ind_ty, ty in enumerate(self.groupedreference.keys()):
            html.add_heading(ty, size=4, bold=False)

            data_table = []
            header_list = ["Query<br>name", "Query<br>number", "Frequencies"]
            type_list = 'sssss'
            col_size_list = [10, 10, 30, 10]
            # print(self.frequency[ty].keys())
            sys.stdout.flush()

            for ind_q, q in enumerate(self.frequency[ty].keys()):
                html.add_figure("venn_" + ty + "_" + q + ".png", align="center", width="600")
                data_table.append([q, str(self.qlen[ty][q]),
                                   ",".join([str(v).rjust(7, " ") for v in list(self.frequency[ty][q].values())])])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align)

            header_list = ["Query 1", "Query 2", "Chi square", "p value"]
            data_table = []
            n = len(list(self.frequency[ty].keys()))
            for q in itertools.combinations(list(range(n)), 2):
                # print(q)
                q1 = list(self.frequency[ty].keys())[q[0]]
                q2 = list(self.frequency[ty].keys())[q[1]]

                sumqf = [x + y for x, y in zip(list(self.frequency[ty][q1].values()), list(self.frequency[ty][q2].values()))]
                nonzero = [i for i, e in enumerate(sumqf) if e != 0]
                qf1 = [list(self.frequency[ty][q1].values())[i] for i in nonzero]
                qf2 = [list(self.frequency[ty][q2].values())[i] for i in nonzero]

                chisq, p, dof, expected = stats.chi2_contingency([qf1, qf2])
                if p < 0.05:
                    data_table.append([q1, q2, value2str(chisq), "<font color=\"red\">" + value2str(p) + "</font>"])
                else:
                    data_table.append([q1, q2, value2str(chisq), value2str(p)])

            if len(data_table):
                html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align)
            else:
                html.add_free_content(["No overlapping regions found."])
        """
        header_list = ["Order"] + self.orig_refs
        type_list = 'sssss' * len(self.referencenames)
        col_size_list = [10,10,10,10] * len(self.referencenames)
        data_table = []
        print(self.referencenames)
        print(self.comb_ref_infor.keys())
        for i in range(len(self.referencenames)):
            data_table.append([i+1] + self.comb_ref_infor.values()[i])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align)
        """

        html.write(os.path.join(dir_name, title, "index.html"))

        # Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        header_list = ["Description", "Argument", "Value"]
        data_table = [["Reference", "-r", args.r],
                      ["Query", "-q", args.q],
                      ["Output directory", "-o", os.path.basename(args.o)],
                      ["Experiment title", "-t", args.t],
                      # ["Grouping tag", "-g", args.g],
                      # ["Coloring tag", "-c", args.c],
                      # ["Background", "-bg", args.bg],
                      ["Organism", "-organism", args.organism]]

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_free_content([
            '<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(
            ['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory, title, "parameters.html"))

    def combinatorial(self, background=None):
        def p2sign(plist, length):
            output = ["-"] * length
            for p in plist:
                output[p] = "+"
            return output

        new_refsp = OrderedDict()
        new_refs = OrderedDict()
        ref_names = []
        self.comb_ref_infor = {}

        for ty in list(self.groupedreference.keys()):
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []

            for i in range(1, n):
                new_refsp[ty].append(itertools.combinations(list(range(n)), i))
            for posi in new_refsp[ty]:
                # print(posi)
                posi = [list(i) for i in posi]

                for p in posi:
                    # print("   " + str(p))
                    pr = posi2set(self.groupedreference[ty], p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)
                    self.comb_ref_infor[pr.name] = p2sign(p, n)
            all_int = posi2set(self.groupedreference[ty], list(range(n)))
            new_refs[ty].append(all_int)
            ref_names.append(all_int.name)
            self.comb_ref_infor[all_int.name] = p2sign(list(range(n)), n)
            """
            # Background
            unions = GenomicRegionSet(name="")
            for r in self.groupedreference[ty]:
                unions.combine(r)
            unions.name = " + ".join([r.name for r in self.groupedreference[ty]])

            nonset = self.backgroung.subtract(unions)
            nonset.name = "!("+"+".join([r.name for r in self.groupedreference[ty]]) + ")"
            new_refs[ty].append(nonset)
            ref_names.append(nonset.name)
            """
        # self.comb_reference = new_refs
        self.groupedreference = copy.deepcopy(new_refs)
        self.orig_refs = copy.deepcopy(self.referencenames)
        self.referencenames = list(set(ref_names))

    def combine_regions(self, background=None):
        new_refsp = OrderedDict()
        new_refs = OrderedDict()
        ref_names = []

        for ty in list(self.groupedreference.keys()):
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []
            for i in range(1, n):
                new_refsp[ty].append(itertools.combinations(list(range(n)), i))
            for posi in new_refsp[ty]:
                posi = [list(i) for i in posi]
                for p in posi:
                    # print("   " + str(p))
                    pr = posi2region(self.groupedreference[ty], p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)

    def comb_stacked_plot(self):
        self.xtickrotation, self.xtickalign = 0, "center"
        f, axs = plt.subplots(1, len(list(self.frequency.keys())), sharey=True)
        f.subplots_adjust(left=0.3)
        # f.set_size_inches(18.5,15)
        # if len(axs) == 1: axs = [axs]
        try:
            axs = axs.reshape(-1)
        except:
            axs = [axs]

        for ai, ax in enumerate(axs):
            # ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ty = list(self.frequency.keys())[ai]
            ax.locator_params(axis='y', nbins=4)
            ax.set_title(list(self.frequency.keys())[ai], y=1)
            r_label = []
            q_label = []
            legends = []

            for ind_q, q in enumerate(self.frequency[ty].keys()):
                if len(axs) == 1:
                    q_label.append(q)
                else:
                    try:
                        q_label.append(self.qEM.get_type(q, "factor"))
                    except:
                        q_label.append(q)
                width = 0.6
                bottom = 0
                summ = sum(self.frequency[ty][q].values())
                for ind_r, r in enumerate(self.referencenames):
                    # for ind_r, rc in enumerate(self.frequency[ty][q]):
                    rc = self.frequency[ty][q][r]
                    if ind_q == 0:
                        # r = self.groupedreference[ty][ind_r].name
                        r_label.append(r)
                    x = ind_q
                    y = rc / summ  # intersect number
                    bar = ax.bar(x, y, width=width, bottom=bottom, color=self.color_list[ind_r], align='center')
                    bottom = bottom + y
                    if ind_q == 0:
                        legends.append(bar)
            ax.yaxis.tick_left()
            ax.set_xticks(list(range(len(q_label))))
            ax.set_xticklabels(q_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=True)
            ax.set_xlim([-0.5, len(q_label) + 0.5])
            ax.set_ylim([0, 1])

            legends.reverse()
            r_label.reverse()

            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)

        legend_name = reversed(self.color_tags)
        axs[-1].legend(legends, legend_name, loc='upper left', handlelength=1, handletextpad=1,
                       columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 1))
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            f.text(-0.025, 0.5, "Percentage", rotation="vertical", va="center")

        f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        self.sbar = f

    def comb_venn(self, directory):
        if len(self.references) == 2:
            print(2)
        elif len(self.references) == 3:
            for ind_ty, ty in enumerate(self.groupedreference.keys()):
                for q in self.query:
                    plt.figure(figsize=(6, 4))
                    plt.title("Venn Diagram: " + q.name)
                    freq = []
                    for r in self.groupedreference[ty]:
                        freq.append(self.frequency[ty][q.name][r.name])

                    # print([r.name for r in self.groupedreference[ty]])
                    self.venn = venn3(subsets=[freq[i] for i in [0, 1, 3, 2, 4, 5, 6]],
                                      set_labels=[n.name for n in self.references])
                    plt.annotate(str(len(q) - sum(freq)), xy=(0.1, 0.1), xytext=(-120, -120),
                                 ha='left', textcoords='offset points',
                                 bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1))
                    plt.savefig(os.path.join(directory, "venn_" + ty + "_" + q.name + ".png"))
                    plt.savefig(os.path.join(directory, "venn_" + ty + "_" + q.name + ".pdf"), format='pdf')
        else:
            print("*** For plotting Venn diagram, the number of references must be 2 or 3.")

    def stest(self, repeat, threshold, mp):

        print("\n\tIntersection random subsampling test:\n    Repeat " + str(repeat) + " times\n")
        self.test_time = repeat
        self.test_d = {}
        plist = OrderedDict()

        for ty in list(self.groupedreference.keys()):
            self.test_d[ty] = {}
            plist[ty] = OrderedDict()
            for r in self.groupedreference[ty]:
                if r.name in self.nalist:
                    continue
                print("\t" + r.name)
                self.test_d[ty][r.name] = {}
                plist[ty][r.name] = OrderedDict()
                print("\t.", end="")
                sys.stdout.flush()
                for q in self.groupedquery[ty]:
                    if r.name == q.name:
                        continue
                    else:
                        print(".", end="")
                        sys.stdout.flush()
                        if q.name in self.nalist:
                            continue
                        # True intersection
                        obs = self.counts[ty][r.name][q.name]
                        qn = q.name
                        if obs[2] == 0:
                            aveinter, chisq, p = "NA", "NA", "1"
                        else:
                            com = q.combine(r, change_name=False, output=True)
                            # Randomization
                            d = []

                            inp = [com, self.rlen[ty][r.name], self.mode_count, threshold]
                            mp_input = [inp for i in range(repeat)]

                            pool = multiprocessing.Pool(processes=mp)
                            mp_output = pool.map(mp_count_intersets, mp_input)
                            pool.close()
                            pool.join()

                            # for i in range(repeat):
                            #    random_r,random_q = com.random_split(size=self.rlen[ty][r.name])
                            #    d.append(random_r.intersect_count(random_q, mode_count=self.mode_count, threshold=threshold))
                            # d.append(count_intersect(random_r, random_q, mode_count=self.mode_count, threshold=threshold))
                            da = numpy.array(mp_output)

                            exp_m = numpy.mean(da, axis=0)
                            # print(exp_m)
                            # print(obs)
                            chisq, p, dof, expected = stats.chi2_contingency([exp_m, obs])
                            aveinter = exp_m[2]

                        plist[ty][r.name][qn] = p
                        self.test_d[ty][r.name][qn] = [aveinter, chisq, p]
                print()

            multiple_correction(plist)

            # c_p = 0
            for r in list(self.test_d[ty].keys()):
                if r in self.nalist:
                    continue
                for q in list(self.test_d[ty][r].keys()):
                    self.test_d[ty][r][q][2] = plist[ty][r][q]
