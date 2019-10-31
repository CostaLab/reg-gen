# Python Libraries


# Local Libraries
# Distal Libraries
from .shared_function import *

# Local test


###########################################################################################
#                    Jaccard test
###########################################################################################

class Jaccard:
    def __init__(self, reference_path, query_path):
        self.rEM, self.qEM = ExperimentalMatrix(), ExperimentalMatrix()
        self.rEM.read(reference_path)
        self.qEM.read(query_path)
        self.references = self.rEM.get_regionsets()
        self.referencenames = self.rEM.get_regionsnames()
        self.query = self.qEM.get_regionsets()
        self.querynames = self.qEM.get_regionsnames()
        self.parameter = []

    def group_refque(self, groupby=False):
        self.groupedreference, self.groupedquery = group_refque(self.rEM, self.qEM, groupby)

    def colors(self, colorby, definedinEM):
        self.color_list = color_groupded_region(self.qEM, self.groupedquery, colorby, definedinEM)
        # self.color_list['Background'] = '0.70'

    def jaccard_test(self, runtime, organism):
        self.jlist = OrderedDict()
        self.realj = OrderedDict()
        self.plist = OrderedDict()
        self.rlen = {}
        self.qlen = {}
        self.rt = runtime
        self.nalist = []
        print2(self.parameter, "\nJaccard Test")
        print2(self.parameter,
               "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}".format("Reference", "Query", "Repeats", "True_Jaccard_index",
                                                                 "p-value", "Time"))
        for ty in list(self.groupedreference.keys()):
            self.jlist[ty] = OrderedDict()
            self.realj[ty] = OrderedDict()
            self.plist[ty] = OrderedDict()
            for i, r in enumerate(self.groupedreference[ty]):
                if r.total_coverage() == 0 and len(r) > 0:
                    self.nalist.append(r.name)
                    continue
                else:
                    if r.name not in list(self.rlen.keys()):
                        self.rlen[r.name] = len(r)
                    self.jlist[ty][r.name] = OrderedDict()
                    self.realj[ty][r.name] = OrderedDict()
                    self.plist[ty][r.name] = OrderedDict()
                    for j, q in enumerate(self.groupedquery[ty]):
                        ts = time.time()
                        # print(q.name + "      " + str(len(q.sepuences[0])))

                        # The real jaccard index from r and q
                        if q.total_coverage() == 0 and len(q) > 0:
                            self.nalist.append(q.name)
                            continue
                        else:
                            if q.name not in list(self.qlen.keys()):
                                self.qlen[q.name] = len(q)
                            self.jlist[ty][r.name][q.name] = []
                            self.realj[ty][r.name][q.name] = q.jaccard(r)
                            for k in range(runtime):
                                random = q.random_regions(organism=organism, multiply_factor=1, overlap_result=True,
                                                          overlap_input=True, chrom_M=False)
                                self.jlist[ty][r.name][q.name].append(r.jaccard(random))
                            # How many randomizations have higher jaccard index than the real index?
                            p = len([x for x in self.jlist[ty][r.name][q.name] if
                                     x > self.realj[ty][r.name][q.name]]) / runtime
                            self.plist[ty][r.name][q.name] = p
                            te = time.time()
                            print2(self.parameter, r.name + "\t" + q.name + "\tx" + str(runtime) + "\t" +
                                   value2str(self.realj[ty][r.name][q.name]) + "\t" + value2str(p) + "\t" +
                                   str(datetime.timedelta(seconds=round(te - ts))))

    def plot(self, logT=False, pw=3, ph=3):
        """ Return boxplot from the given tables.

        """
        self.fig = []
        self.xtickrotation, self.xtickalign = 0, "center"

        tw = pw
        # th = len(self.jlist.keys()) * ph
        # f, ax = plt.subplots(len(self.jlist[t].keys()), 1, dpi=300,
        #                       figsize=(tw, th) )
        legend_x = 1.05
        for it, t in enumerate(self.jlist.keys()):
            # f, axarr = plt.subplots(1, len(self.jlist[t].keys()), dpi=300, sharey=True)
            # legend_x = 1.05
            # nm = len(self.jlist.keys()) * len(self.jlist.values()[0]) * len(self.jlist.values()[0])
            # if nm > 30:
            #     f.set_size_inches(nm * 0.1 +1 ,nm * 0.1 +1)
            #     legend_x = 1.2
            #     self.xtickrotation, self.xtickalign = 70,"right"
            th = len(list(self.jlist[t].keys())) * ph
            f, axarr = plt.subplots(len(list(self.jlist[t].keys())), 1, dpi=300,
                                    figsize=(tw, th))
            try:
                axarr = axarr.reshape(-1)
            except:
                axarr = [axarr]
            plt.subplots_adjust(bottom=0.3)
            if logT:
                axarr[0].set_ylabel("Jaccard index (log)", fontsize=8)
            else:
                axarr[0].set_ylabel("Jaccard index (Intersect/Union)", fontsize=8)

            for i, r in enumerate(self.jlist[t].keys()):
                # axarr[i].set_title(r, y=0.94)
                if logT:
                    axarr[i].set_yscale('log')
                else:
                    axarr[i].locator_params(axis='y', nbins=4)
                axarr[i].tick_params(axis='y', direction='out')

                axarr[i].yaxis.tick_left()
                axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
                d = []  # Store data within group
                color_t = []  # Store tag for coloring boxes
                # x_ticklabels = []  # Store ticklabels
                axarr[i].set_xlabel(r, rotation=self.xtickrotation, ha=self.xtickalign)

                for j, q in enumerate(self.jlist[t][r].keys()):
                    d.append(self.jlist[t][r][q])
                    color_t.append(self.color_list[q])
                    # x_ticklabels.append(q)
                # Fine tuning boxplot
                bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None,
                                      patch_artist=True, bootstrap=None)
                z = 10  # zorder for bosplot
                plt.setp(bp['whiskers'], color='black', linestyle='-', linewidth=0.8, zorder=z)
                plt.setp(bp['fliers'], markerfacecolor='gray', color='white', alpha=0.3, markersize=1.8, zorder=z)
                plt.setp(bp['caps'], color='white', zorder=z)
                plt.setp(bp['medians'], color='black', linewidth=0.5, zorder=z + 1)
                legends = []
                for patch, color in zip(bp['boxes'], color_t):
                    patch.set_facecolor(color)  # When missing the data, the color patch will exceeds
                    patch.set_edgecolor("none")
                    patch.set_zorder(z)
                    legends.append(patch)
                axarr[i].scatter(x=list(range(1, 1 + len(list(self.jlist[t][r].keys())))),
                                 y=[y for y in list(self.realj[t][r].values())],
                                 s=10, c="red", edgecolors='none', zorder=2)
                # Fine tuning subplot
                # axarr[i].set_xticks(range(len(self.jlist[t][r].keys())))
                # plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
                # axarr[i].set_xticklabels(self.jlist[t][r].keys(), rotation=0, fontsize=10)

                # axarr[i].set_ylim(bottom=0.95)
                for spine in ['top', 'right', 'left', 'bottom']:
                    axarr[i].spines[spine].set_visible(False)
                axarr[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)

                if i > 0:
                    plt.setp(axarr[i].get_yticklabels(), visible=False)
                    # plt.setp(axarr[i].get_yticks(),visible=False)
                    axarr[i].minorticks_off()
                    axarr[i].tick_params(axis='y', which='both', left=False, right=False, labelbottom=False)

            plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
            axarr[-1].legend(legends[0:len(list(self.jlist[t][r].keys()))], list(self.jlist[t][r].keys()), loc='center left',
                             handlelength=1,
                             handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size': 10},
                             bbox_to_anchor=(legend_x, 0.5))
            # f.tight_layout(pad=2, h_pad=None, w_pad=None)
            self.fig.append(f)

    def gen_html(self, outputname, title, align=50):
        fp = os.path.join(current_dir, outputname, title)
        link_d = {title: "index.html"}
        html = Html(name="Viz", links_dict=link_d,
                    fig_dir=os.path.join(current_dir, outputname, "fig"), other_logo="viz",
                    homepage="../index.html")
        for i in range(len(self.fig)):
            html.add_figure("jaccard_test" + str(i + 1) + ".png", align="center")

        header_list = ["Reference<br>name",
                       "Query<br>name",
                       "Reference<br>number",
                       "Query<br>number",
                       "True<br>Jaccard<br>index",
                       "Average<br>random<br>Jaccard",
                       "Positive<br>Association<br>p-value",
                       "Negative<br>Association<br>p-value"]

        type_list = 'sssssssssssss'
        col_size_list = [10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10]
        data_table = []

        for ind_ty, ty in enumerate(self.jlist.keys()):
            # html.add_heading(ty, size = 4, bold = False)
            for ind_r, r in enumerate(self.jlist[ty].keys()):
                for ind_q, q in enumerate(self.jlist[ty][r].keys()):
                    rej = self.realj[ty][r][q]
                    rj = numpy.mean(self.jlist[ty][r][q])
                    p = self.plist[ty][r][q]
                    np = 1 - p
                    rl = str(self.rlen[r])
                    ql = str(self.qlen[q])

                    if self.plist[ty][r][q] < 0.05:
                        if self.realj[ty][r][q] > rj:
                            data_table.append([r, q, rl, ql, value2str(rej), value2str(rj),
                                               "<font color=\"red\">" + value2str(p) + "</font>",
                                               value2str(np)])
                        else:
                            data_table.append([r, q, rl, ql, value2str(rej), value2str(rj),
                                               value2str(np),
                                               "<font color=\"red\">" + value2str(p) + "</font>"])
                    else:
                        data_table.append([r, q, rl, ql, value2str(rej), value2str(rj), value2str(p), value2str(np)])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align)

        header_list = ["Assumptions and hypothesis"]
        data_table = [['Randomization was performed ' + str(self.rt) + ' times.'],
                      [
                          'For projection test, the reference and query should have non-zero length in order to calculate its Jaccard index.'],
                      ['Positive association is defined by: True Jaccard index > Averaged random Jaccard.'],
                      ['Negative association is defined by: True Jaccard index < Averaged random Jaccard.']]
        self.nalist = set(self.nalist)
        if len(self.nalist) > 0:
            data_table.append([
                'The following region sets contain zero-length regions which cause error in Jaccard index calculation, please check it:<br>' +
                '<font color=\"red\">' + ', '.join([s for s in self.nalist]) + '</font>'])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content([
            '<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(
            ['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.write(os.path.join(fp, "index.html"))

    def table(self, directory, folder):
        arr = numpy.array([["#reference", "query", "true_jaccard", "random_jaccard", "p-value"]])
        for ty in list(self.plist.keys()):
            for r in list(self.plist[ty].keys()):
                for q in list(self.plist[ty][r].keys()):
                    ar = numpy.array([[r, q, self.realj[ty][r][q], self.qlist[ty][r][q], self.plist[ty][r][q]]])
                    arr = numpy.vstack((arr, ar))
        output_array(arr, directory, folder, filename="output_table.txt")
