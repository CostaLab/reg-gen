# Python Libraries


import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import numpy

# Local Libraries
# Distal Libraries
from ..Util import Html
from ..CoverageSet import *
from ..ExperimentalMatrix import *
from .shared_function import output_array, group_refque, color_groupded_region, multiple_correction, value2str

# Local test
current_dir = os.getcwd()


###########################################################################################
#                    Projection test
###########################################################################################


class Projection:
    def __init__(self, reference_path, query_path, load_bed=True):
        # Reference
        self.rEM = ExperimentalMatrix()
        self.rEM.read(reference_path,load_bed=load_bed)
        # self.rEM.remove_empty_regionset()
        self.references = self.rEM.get_regionsets()
        self.referencenames = self.rEM.get_regionsnames()
        # Query
        self.qEM = ExperimentalMatrix()
        self.qEM.read(query_path)
        # self.qEM.remove_empty_regionset()
        self.query = self.qEM.get_regionsets()
        self.querynames = self.qEM.get_regionsnames()
        self.parameter = []
        self.background = None

    def group_refque(self, groupby=False):
        self.groupedreference, self.groupedquery = group_refque(self.rEM, self.qEM, groupby)

    def colors(self, colorby, definedinEM):
        ############# Color #####################################
        # self.color_list = colormap(self.qEM, colorby, definedinEM)
        self.color_list = color_groupded_region(self.qEM, self.groupedquery, colorby, definedinEM)
        # self.color_tags = gen_tags(self.qEM, colorby)
        # self.color_tags.append('Background')
        self.color_list['Background'] = '0.70'

    def ref_union(self):
        self.background = OrderedDict()
        for ty in list(self.groupedreference.keys()):
            self.background[ty] = GenomicRegionSet("union of references")
            for r in self.groupedreference[ty]:
                self.background[ty].combine(r)
            self.background[ty].merge()

        for ty in list(self.groupedreference.keys()):
            rlist = [r.trim_by(background=self.background[ty]) for r in self.groupedreference[ty]]
            self.groupedreference[ty] = rlist
            qlist = [q.trim_by(background=self.background[ty]) for q in self.groupedquery[ty]]
            self.groupedquery[ty] = qlist

    def set_background(self, bed_path):
        bg = GenomicRegionSet("background")
        bg.read(bed_path)
        self.background = OrderedDict()
        for ty in list(self.groupedreference.keys()):
            self.background[ty] = bg
            # rlist = [r.trim_by(background=bg) for r in self.groupedreference[ty]]
            # self.groupedreference[ty] = rlist
            #
            # qlist = [q.trim_by(background=bg) for q in self.groupedquery[ty]]
            # self.groupedquery[ty] = qlist

    def projection_test(self, organism):
        self.bglist = OrderedDict()
        self.qlist = OrderedDict()
        self.plist = OrderedDict()
        self.interq_list = OrderedDict()
        self.lenlist = {}
        # print2(self.parameter, "\nProjection test")
        # print2(self.parameter, "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}".format("Reference","Background", "Query", "Proportion", "p value"))

        all_p = {}
        for ty in list(self.groupedquery.keys()):
            # print(ty)
            self.bglist[ty] = OrderedDict()
            self.qlist[ty] = OrderedDict()
            self.plist[ty] = OrderedDict()
            self.interq_list[ty] = OrderedDict()
            if self.background:
                bgset = self.background[ty]
            else:
                bgset = None

            for i, r in enumerate(self.groupedreference[ty]):
                # print(r.name)
                self.bglist[ty][r.name] = OrderedDict()
                self.qlist[ty][r.name] = OrderedDict()
                self.plist[ty][r.name] = OrderedDict()
                self.interq_list[ty][r.name] = OrderedDict()
                self.lenlist[r.name] = len(r)
                for j, q in enumerate(self.groupedquery[ty]):
                    # print([ty, r.name, q.name])
                    if r.name == q.name:
                        continue
                    else:
                        rr = self.rEM.get_regionset(name=r.name)
                        # print([len(rr), len(q)])
                        print(".", end="")
                        bg, ratio, p, interq = rr.projection_test(q, organism, extra=True, background=bgset)
                        self.bglist[ty][r.name][q.name] = bg
                        self.qlist[ty][r.name][q.name] = ratio
                        self.plist[ty][r.name][q.name] = p
                        self.interq_list[ty][r.name][q.name] = interq
                        self.lenlist[q.name] = len(q)
                        # if r in self.backgrounds.keys(): pass
                        # else: self.backgrounds[r] = bg

        # multiple test correction
        multiple_correction(self.plist)

        for ty in list(self.groupedquery.keys()):
            for i, r in enumerate(self.groupedreference[ty]):
                for j, q in enumerate(self.groupedquery[ty]):
                    # if r.name == q.name: continue
                    # else:
                    # bg = self.bglist[ty][r.name][q.name]
                    # ratio = self.qlist[ty][r.name][q.name]
                    # p = self.plist[ty][r.name][q.name]
                    self.qlist[ty][r.name]['Background'] = list(self.bglist[ty][r.name].values())[0]

    def output_interq(self, directory):
        """Output the intersected query to the reference in BED format"""
        try:
            os.stat(os.path.dirname(directory))
        except:
            os.mkdir(os.path.dirname(directory))
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
        for ty in list(self.interq_list.keys()):
            if ty:
                g = ty + "_"
            else:
                g = ""
            for r in list(self.interq_list[ty].keys()):
                for q in list(self.interq_list[ty][r].keys()):
                    self.interq_list[ty][r][q].write_bed(os.path.join(directory, g + q + "_intersected_" + r + ".bed"))

    def plot(self, logt=None, pw=3, ph=3):

        tw = pw
        th = len(list(self.qlist.keys())) * ph
        f, ax = plt.subplots(len(list(self.qlist.keys())), 1, dpi=300, figsize=(tw, th))

        # f, ax = plt.subplots(len(self.qlist.keys()),1)
        try:
            ax = ax.reshape(-1)
        except:
            ax = [ax]
        # nm = len(self.groupedreference.keys()) * len(self.groupedreference.values()[0]) * len(self.groupedquery.values()[0])
        # if nm > 40:
        #     f.set_size_inches(nm * 0.2 +1 ,7)

        g_label = []
        for ind_ty, ty in enumerate(self.qlist.keys()):
            g_label.append(ty)
            r_label = []
            for ind_r, r in enumerate(self.qlist[ty].keys()):
                r_label.append(r)
                width = 0.8 / (len(list(self.qlist[ty][r].keys())) + 1)  # Plus one background
                for ind_q, q in enumerate(self.qlist[ty][r].keys()):
                    x = ind_r + ind_q * width + 0.1
                    y = self.qlist[ty][r][q]
                    if y == 0 and logt:
                        y = 0.000001
                    # print("    "+r+"     "+q+"     "+str(x)+"     "+str(y))
                    ax[ind_ty].bar(x, y, width=width, color=self.color_list[q], edgecolor="none",
                                   align='edge', log=logt, label=q)
            if logt:
                ax[ind_ty].set_yscale('log')
            else:
                ax[ind_ty].locator_params(axis='y', nbins=2)

            # ax[ind_ty].set_ylabel("Percentage of intersected regions",fontsize=12)
            ax[ind_ty].set_title(ty)
            ax[ind_ty].yaxis.tick_left()
            ax[ind_ty].set_ylabel('Percentage of intersected regions', fontsize=8)
            ax[ind_ty].set_xticks([i + 0.5 - 0.5 * width for i in range(len(r_label))])
            ax[ind_ty].set_xticklabels(r_label, rotation=30, ha="right", fontsize=8)
            ax[ind_ty].tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=True)

            handles, labels = ax[ind_ty].get_legend_handles_labels()
            # uniq_labels = unique(labels)
            uniq_labels = [q.name for q in self.groupedquery[ty]] + ["Background"]

            ax[ind_ty].legend([handles[labels.index(l)] for l in uniq_labels], uniq_labels,
                              loc='center left', handlelength=1, handletextpad=1,
                              columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax[ind_ty].spines[spine].set_visible(False)
        # f.text(-0.025, 0.5, "Percentage of intersected regions",fontsize=12, rotation="vertical", va="center")
        # f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        f.tight_layout()
        self.fig = f

    def heatmap(self):
        f, ax = plt.subplots(1, len(list(self.plist.keys())))
        try:
            ax = ax.reshape(-1)
        except:
            ax = [ax]

        g_label = []
        for ind_ty, ty in enumerate(self.plist.keys()):
            g_label.append(ty)
            r_label = []
            data = []
            for ind_r, r in enumerate(self.plist[ty].keys()):
                r_label.append(r)
                # data.append(self.plist[ty][r].values())
                for ind_q, q in enumerate(self.plist[ty][r].keys()):
                    pass
            da = numpy.array(data)
            da = da.transpose()
            # im = plt.imshow(da, cmap=ax[ind_r], vmin=, vmax, origin, extent, shape, filternorm, filterrad, imlim, resample, url, hold)

    def gen_html(self, directory, title, args, align=50):
        dir_name = os.path.basename(directory)
        statistic_table = []
        # check_dir(directory)
        html_header = "Projection Test: " + dir_name
        link_d = OrderedDict()
        link_d["Projection test"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        html.add_figure("projection_test.png", align="center")

        header_list = ["No.",
                       "Reference<br>name",
                       "Query<br>name",
                       "Reference<br>number",
                       "Query<br>number",
                       "Proportion",
                       "Background<br>proportion",
                       "Positive<br>association<br>p-value",
                       "Negative<br>association<br>p-value"]
        statistic_table.append(["Reference_name", "Query_name", "Reference_number",
                                "Query_number", "Proportion", "Background_proportion",
                                "Positive_association_p-value", "Negative_association_p-value"])
        type_list = 'ssssssssssssssss'
        col_size_list = [5, 10, 10, 10, 10, 10, 10, 15, 15]

        nalist = []
        for ind_ty, ty in enumerate(self.plist.keys()):
            html.add_heading(ty, size=4, bold=False)
            data_table = []
            for ind_r, r in enumerate(self.plist[ty].keys()):
                rlen = str(self.lenlist[r])
                for ind_q, q in enumerate(self.plist[ty][r].keys()):
                    qlen = str(self.lenlist[q])
                    backv = value2str(self.qlist[ty][r]['Background'])
                    propor = value2str(self.qlist[ty][r][q])
                    pv = self.plist[ty][r][q]
                    if pv == "na":
                        nalist.append(r)
                        continue
                    elif self.qlist[ty][r][q] < args.cfp:
                        continue
                    else:
                        pvn = 1 - pv

                        if self.plist[ty][r][q] < 0.05:
                            if self.qlist[ty][r]['Background'] < self.qlist[ty][r][q]:
                                data_table.append([str(ind_ty), r, q, rlen, qlen, propor, backv,
                                                   "<font color=\"red\">" + value2str(pv) + "</font>", value2str(pvn)])
                                statistic_table.append([r, q, rlen, qlen, propor, backv, value2str(pv), value2str(pvn)])
                            else:
                                data_table.append([str(ind_ty), r, q, rlen, qlen, propor, backv,
                                                   value2str(pvn), "<font color=\"red\">" + value2str(pv) + "</font>"])
                                statistic_table.append([r, q, rlen, qlen, propor, backv, value2str(pvn), value2str(pv)])
                        else:
                            data_table.append(
                                [str(ind_ty), r, q, rlen, qlen, propor, backv, value2str(pv), value2str(pvn)])
                            statistic_table.append([r, q, rlen, qlen, propor, backv, value2str(pv), value2str(pvn)])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, sortable=True)
            output_array(statistic_table, directory=directory, folder=title, filename="statistics" + ty + ".txt")

        header_list = ["Assumptions and hypothesis"]
        data_table = [['If the background proportion is too small, it may cause bias in p value.'],
                      [
                          'For projection test, the reference GenomicRegionSet should have non-zero length in order to calculate its background proportion.'],
                      ['P values are corrected by multiple test correction.'],
                      ['Positive association is defined by: Proportion > Background.'],
                      ['Negative association is defined by: Proportion < Background.']]

        nalist = set(nalist)
        if len(nalist) > 0:
            data_table.append([
                'The following references contain zero-length region which cause error in proportion calculation, please check it:<br>' +
                '     <font color=\"red\">' + ', '.join([s for s in nalist]) + '</font></p>'])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_fixed_rank_sortable()

        html.write(os.path.join(directory, os.path.join(title, "index.html")))

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
                      ["Organism", "-organism", args.organism],
                      ["Cutoff of proportion", "-cfp", str(args.cfp)]]

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")
        html.add_free_content([
            '<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(
            ['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory, os.path.join(title, "parameters.html")))

    def table(self, directory, folder):
        arr = numpy.array([["#reference", "query", "background", "proportion", "p-value"]])
        for ty in list(self.plist.keys()):
            for r in list(self.plist[ty].keys()):
                for q in list(self.plist[ty][r].keys()):
                    ar = numpy.array(
                        [[r, q, self.qlist[ty][r]['Background'], self.qlist[ty][r][q], self.plist[ty][r][q]]])
                    arr = numpy.vstack((arr, ar))
        output_array(arr, directory, folder, filename="output_table.txt")

    def distribution(self, organism):
        genome = GenomicRegionSet("genome")
        genome.get_genome_data(organism)
        all_cov = genome.total_coverage()
        self.chrom_list = []
        for ss in genome:
            self.chrom_list.append(ss.chrom)
        self.chrom_list.sort()

        # self.distriDict = OrderedDict()
        self.disperDict = OrderedDict()

        for ty in list(self.groupedreference.keys()):
            # self.distriDict[ty] = OrderedDict()
            self.disperDict[ty] = OrderedDict()
            # Reference
            for r in self.groupedreference[ty]:
                r.merge()
                len_r = r.total_coverage()
                # self.distriDict[ty][r.name] = []
                self.disperDict[ty][r.name] = []

                for ch in self.chrom_list:
                    rc = r.any_chrom(chrom=ch)
                    nr = sum([len(s) for s in rc])
                    # self.distriDict[ty][r.name].append(nr)
                    self.disperDict[ty][r.name].append(nr / len_r)

            # Query
            for q in self.groupedquery[ty]:
                q.merge()
                len_q = q.total_coverage()
                # self.distriDict[ty][q.name] = []
                self.disperDict[ty][q.name] = []

                for ch in self.chrom_list:
                    qc = q.any_chrom(chrom=chr)
                    nq = sum([len(s) for s in qc])
                    # self.distriDict[ty][q.name].append(nq)
                    self.disperDict[ty][q.name].append(nq / len_q)
            # Genome
            # self.distriDict[ty]["Genome"] = [len(genome.any_chrom(chrom=chr)) for chr in self.chrom_list]

            self.disperDict[ty]["Genome"] = [len(genome.any_chrom(chrom=chrom)[0]) / all_cov for chrom in self.chrom_list]

    def plot_distribution(self):
        def to_percentage(x, pos=0):
            return '{:.2f} %'.format(100 * x)

        self.fig = []

        for ty in list(self.disperDict.keys()):
            colors = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(list(self.disperDict[ty].keys())))).tolist()

            f, ax = plt.subplots()
            f.set_size_inches(10.5, 30)
            width = 0.9 / len(list(self.disperDict[ty].keys()))
            ind = np.arange(len(self.chrom_list))
            coverage = self.disperDict[ty]

            for ind_r, r in enumerate(self.disperDict[ty].keys()):
                ax.barh(ind + width * ind_r, self.disperDict[ty][r], width, color=colors[ind_r])

            plt.xlabel('Percentage')
            ax.xaxis.set_major_formatter(mtick.FuncFormatter(to_percentage))

            ax.minorticks_off()
            ax.set_yticks([x + 0.5 for x in range(len(self.chrom_list))])
            ax.set_yticklabels(self.chrom_list, rotation=0, ha="right")
            ax.tick_params(axis='y', which='both', top=False, bottom=False, labelbottom=True)

            ax.legend(list(self.disperDict[ty].keys()), loc='center left', handlelength=1, handletextpad=1,
                      columnspacing=2, borderaxespad=0., prop={'size': 10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right', 'left', 'bottom']:
                ax.spines[spine].set_visible(False)
            f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
            self.fig.append(f)

    def gen_html_distribution(self, outputname, title, align=50):
        fp = os.path.join(current_dir, outputname, title)
        link_d = {title: "distribution.html"}
        html = Html(name="Viz", links_dict=link_d, fig_dir=os.path.join(current_dir, outputname, "fig"),
                    other_logo="viz", homepage="../index.html")
        for i, f in enumerate(self.fig):
            html.add_figure("distribution_test_" + str(i) + ".png", align="center")

        html.add_free_content(['<p style=\"margin-left: ' + str(align + 150) + '">' +
                               '** </p>'])

        type_list = 'ssssssssssssssssssssssssssssssssssssssssssssss'
        col_size_list = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                         10, 10]
        data_table = []
        for ind_ty, ty in enumerate(self.disperDict.keys()):
            header_list = ["Chromosome"] + list(self.disperDict[ty].keys())
            html.add_heading(ty, size=4, bold=False)
            for i, ch in enumerate(self.chrom_list):
                # for ind_r,r in enumerate(self.disperDict[ty].keys()):

                data_table.append(
                    [ch] + ["{:.3f} %".format(100 * self.disperDict[ty][r][i]) for r in list(self.disperDict[ty].keys())])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align)

        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content([
            '<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(
            ['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.write(os.path.join(fp, "distribution.html"))
