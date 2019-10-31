# Python Libraries



import os
from collections import OrderedDict

import numpy
from matplotlib import pyplot as plt
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from scipy.stats import mstats, mannwhitneyu

from .shared_function import output_array, gen_tags, tag_from_r, colormap, multiple_correction, value2str
from ..CoverageSet import CoverageSet
from ..ExperimentalMatrix import ExperimentalMatrix
from ..GenomicRegionSet import GenomicRegionSet
# Local Libraries
# Distal Libraries
from ..Util import Html


###########################################################################################
#                    Boxplot
###########################################################################################

class Boxplot:
    """
    input:
        exps: input experimental matrix
        title: Default = boxplot
        groupby: Group the data by the given factor in the header of experimental matrix

    output:
        parameters: list of records
        figs: a list of figure(s)
    """

    def __init__(self, EMpath, fields, title="boxplot", df=False):
        # Read the Experimental Matrix
        self.title = title
        self.exps = ExperimentalMatrix()
        self.exps.read(EMpath)
        for f in self.exps.fields:
            if f not in ['name', 'type', 'file', "reads", "regions", "factors"]:
                self.exps.match_ms_tags(f)
                self.exps.remove_name()
        self.beds = self.exps.get_regionsets()  # A list of GenomicRegionSets
        self.bednames = self.exps.get_regionsnames()
        self.reads = self.exps.get_readsfiles()
        self.readsnames = self.exps.get_readsnames()
        self.fieldsDict = self.exps.fieldsDict
        self.parameter = []
        self.df = df

    def combine_allregions(self):

        self.all_bed = GenomicRegionSet("All regions")
        for bed in self.beds:
            self.all_bed.combine(bed)
        self.all_bed.remove_duplicates()  # all_bed is sorted!!

    def bedCoverage(self):
        """ Return coverage matrix of multiple reads on one bed.
        bed --> GenomicRegionSet
        """
        c = []
        for rp in self.reads:
            print("    processing: ..." + rp[-45:])
            r = os.path.abspath(rp)  # Here change the relative path into absolute path
            cov = CoverageSet(r, self.all_bed)
            cov.coverage_from_genomicset(r)
            try:
                cov.normRPM()
            except:
                pass
            c.append(cov.coverage)
        self.all_table = numpy.transpose(c)

    def quantile_normalization(self):
        """ Return the np.array which contains the normalized values
        """
        rank_matrix = []
        for c in range(self.all_table.shape[1]):
            col = self.all_table[:, c]
            rank_col = mstats.rankdata(col)
            rank_matrix.append(rank_col)

        ranks = numpy.array(rank_matrix)
        trans_rank = numpy.transpose(ranks)

        # Calculate for means of ranks
        print("    Calculating for the mean of ranked data...")
        sort_matrix = numpy.sort(self.all_table, axis=0)
        means = []
        for r in range(self.all_table.shape[0]):
            row = [x for x in sort_matrix[r, :]]
            means.append(numpy.mean(row))

        # Replace the value by new means
        print("    Replacing the data value by normalized mean...")
        normalized_table = numpy.around(trans_rank)
        for i, v in enumerate(means):
            normalized_table[normalized_table == i + 1] = v
        # print(rounded_rank)
        self.norm_table = normalized_table

    def tables_for_plot(self):
        """ Return a Dict which stores all tables for each bed with file name as its key. """
        self.tableDict = OrderedDict()  # Storage all tables for each bed with bedname as the key
        conList = []  # Store containers of beds
        iterList = []

        for i, bed in enumerate(self.beds):
            self.tableDict[bed.name] = []
            bed.sort()
            conList.append(bed.__iter__())
            iterList.append(next(conList[-1]))

        for i, r in enumerate(self.all_bed.sequences):
            for j in range(len(self.beds)):
                while r > iterList[j]:
                    try:
                        iterList[j] = next(conList[j])
                    except:
                        break
                if r == iterList[j]:
                    self.tableDict[self.beds[j].name].append(self.norm_table[i])
                elif r < iterList[j]:
                    continue

    def print_plot_table(self, directory, folder):
        for i, bed in enumerate(self.tableDict.keys()):
            table = []
            header = ["loci"]
            for rp in self.reads:
                header.append(os.path.basename(rp))
            table.append(header)
            for j, re in enumerate(self.beds[i]):
                table.append(["_".join([re.chrom, str(re.initial), str(re.final)])] + self.tableDict[bed][j].tolist())
            # output_array(table, directory, folder, filename="table_" + bed + ".txt")
            output_array(table, directory, folder, filename="table_" + bed + ".txt")

    def group_tags(self, groupby, sortby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        self.tag_type = [groupby, sortby, colorby]

        if groupby == "None":
            self.group_tags = [""]
        else:
            self.group_tags = gen_tags(self.exps, groupby)
        if sortby == "None":
            self.sort_tags = [""]
        else:
            self.sort_tags = gen_tags(self.exps, sortby)
        if colorby == "None":
            self.color_tags = [""]
        else:
            self.color_tags = gen_tags(self.exps, colorby)

    def group_data(self, directory, folder, log=False):
        plotDict = OrderedDict()  # Extracting the data from different bed_bams file
        cuesbed = OrderedDict()  # Storing the cues for back tracking
        cuesbam = OrderedDict()
        for bedname in list(self.tableDict.keys()):
            plotDict[bedname] = OrderedDict()
            mt = numpy.array(self.tableDict[bedname])

            cuesbed[bedname] = set(tag_from_r(self.exps, self.tag_type, bedname))
            # cuesbed[bedname] = [tag for tag in self.exps.get_types(bedname) if tag in self.group_tags + self.sort_tags + self.color_tags]
            # print(cuesbed[bedname])

            for i, readname in enumerate(self.readsnames):
                plotDict[bedname][readname] = mt[:, i]
                # print(plotDict[bedname][readname])
                cuesbam[readname] = set(tag_from_r(self.exps, self.tag_type, readname))
                # print(cuesbam[readname])
                # cuesbam[readname] = [tag for tag in self.exps.get_types(readname) if tag in self.group_tags + self.sort_tags + self.color_tags]

        sortDict = OrderedDict()  # Storing the data by sorting tags
        for g in self.group_tags:
            # print("    "+g)
            sortDict[g] = OrderedDict()
            for a in self.sort_tags:
                # print("        "+a)
                sortDict[g][a] = OrderedDict()
                for c in self.color_tags:
                    # sortDict[g][a][c] = None
                    # print("            "+c)
                    for i, bed in enumerate(cuesbed.keys()):
                        # print([g, a, c])

                        if {g, a, c} >= cuesbed[bed]:
                            # print(cuesbed[bed])
                            sortDict[g][a][c] = []
                            for bam in list(cuesbam.keys()):

                                if {g, a, c} >= cuesbam[bam]:
                                    # print(cuesbam[bam])
                                    # print([bed, bam])

                                    if self.df:
                                        sortDict[g][a][c].append(plotDict[bed][bam])
                                        if len(sortDict[g][a][c]) == 2:
                                            if log:
                                                sortDict[g][a][c][0] = numpy.log(sortDict[g][a][c][0] + 1)
                                                sortDict[g][a][c][1] = numpy.log(sortDict[g][a][c][1] + 1)
                                                sortDict[g][a][c] = numpy.subtract(sortDict[g][a][c][0],
                                                                                   sortDict[g][a][c][1]).tolist()
                                            else:
                                                sortDict[g][a][c] = numpy.subtract(sortDict[g][a][c][0],
                                                                                   sortDict[g][a][c][1]).tolist()
                                    else:
                                        # print(plotDict[bed][bam])
                                        sortDict[g][a][c] = plotDict[bed][bam]
        self.sortDict = sortDict

    def color_map(self, colorby, definedinEM):
        self.colors = colormap(self.exps, colorby, definedinEM)

    def print_table(self, directory, folder):
        self.print_plot_table(directory,folder)
        # # self.printtable = OrderedDict()
        # table = []
        # maxn = 0
        # # table.append(["#group_tag", "sort_tag", "color_tag", "Signals"])
        # for i, g in enumerate(self.group_tags):
        #     for k, a in enumerate(self.sort_tags):
        #         for j, c in enumerate(self.color_tags):
        #             table.append([g, a, c] + [str(x) for x in self.sortDict[g][a][c]])
        #             c = len(self.sortDict[g][a][c]) + 3
        #             if c > maxn: maxn = c
        # for i, t in enumerate(table):
        #     if len(t) < maxn:
        #         table[i] = t + ["n.a."] * (maxn - len(t))

        # print

        # output_array(numpy.array([list(x) for x in zip(*table)]), directory, folder, filename="output_table.txt")

    def plot(self, title, scol, logT=False, ylim=False, pw=3, ph=4):
        """ Return boxplot from the given tables.

        """
        self.xtickrotation, self.xtickalign = 0, "center"
        if len(self.group_tags) < 2:
            ticklabelsize = pw * 1.5
        else:
            ticklabelsize = pw * 6
        tw = len(self.group_tags) * pw
        th = ph

        f, axarr = plt.subplots(1, len(self.group_tags), dpi=300, sharey=scol,
                                figsize=(tw, th))
        # f, axarr = plt.subplots(1, len(self.group_tags), dpi=300, sharey = scol)

        # nm = len(self.group_tags) * len(self.color_tags) * len(self.sort_tags)
        # if nm > 30:
        # f.set_size_inches(nm * 0.25 ,nm * 0.15)
        # legend_x = 1.2
        # self.xtickrotation, self.xtickalign = 70,"right"

        canvas = FigureCanvas(f)
        canvas.set_window_title(title)
        try:
            axarr = axarr.reshape(-1)
        except:
            axarr = [axarr]
        # plt.subplots_adjust(bottom=0.3)
        if logT:
            if self.df:
                axarr[0].set_ylabel("Read number difference (log)",
                                    fontsize=ticklabelsize + 1)
            else:
                axarr[0].set_ylabel("Read number (log)", fontsize=ticklabelsize + 1)
        else:
            if self.df:
                axarr[0].set_ylabel("Read number difference", fontsize=ticklabelsize + 1)
            else:
                axarr[0].set_ylabel("Read number", fontsize=ticklabelsize + 1)

        for i, g in enumerate(self.sortDict.keys()):
            # if self.df:
            #     axarr[i].set_title(g + "_df", y=1.02, fontsize=ticklabelsize + 2)
            # else:
            axarr[i].set_title(g, y=1.02, fontsize=ticklabelsize + 2)

            if logT and not self.df:
                axarr[i].set_yscale('log')
            else:
                axarr[i].locator_params(axis='y', nbins=4)

            axarr[i].tick_params(axis='y', direction='out')
            axarr[i].yaxis.tick_left()
            axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
            if ylim:
                axarr[i].set_ylim([-ylim, ylim])
            d = []  # Store data within group
            color_t = []  # Store tag for coloring boxes
            x_ticklabels = []  # Store ticklabels
            for j, a in enumerate(self.sortDict[g].keys()):
                # if len(a) > 10:
                # print(a)
                self.xtickrotation = 70
                self.xtickalign = "right"
                for k, c in enumerate(self.sortDict[g][a].keys()):
                    if not numpy.any(self.sortDict[g][a][c]):  # When there is no matching data, skip it
                        continue
                    else:
                        if self.df:
                            d.append(self.sortDict[g][a][c])
                        else:
                            d.append([x + 1 for x in self.sortDict[g][a][c]])
                        color_t.append(self.colors[k])
                        x_ticklabels.append(a)  # + "." + c

            # Fine tuning boxplot
            # print(d)
            bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None,
                                  widths=None, patch_artist=True, bootstrap=None)
            z = 10  # zorder for boxplot
            plt.setp(bp['whiskers'], color='black', linestyle='-', linewidth=0.8, zorder=z)
            plt.setp(bp['fliers'], markerfacecolor='gray', color='white', alpha=0.3, markersize=1.8, zorder=z)
            plt.setp(bp['caps'], color='white', zorder=z)
            plt.setp(bp['medians'], color='black', linewidth=1.5, zorder=z + 1)
            legends = []
            for patch, color in zip(bp['boxes'], color_t):
                patch.set_facecolor(color)  # When missing the data, the color patch will exceeds
                patch.set_edgecolor("none")
                patch.set_zorder(z)
                legends.append(patch)

            # Fine tuning subplot
            axarr[i].set_xticks([len(self.color_tags) * n + 1 + (len(self.color_tags) - 1) / 2 for n, s in
                                 enumerate(self.sortDict[g].keys())])
            # plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
            axarr[i].set_xticklabels(list(self.sortDict[g].keys()), rotation=self.xtickrotation,
                                     ha=self.xtickalign)
            # axarr[i].set_xticklabels(self.sortDict[g].keys(), rotation=70, ha=self.xtickalign, fontsize=10)

            # axarr[i].set_ylim(bottom=0.95)
            for spine in ['top', 'right', 'left', 'bottom']:
                axarr[i].spines[spine].set_visible(False)
            axarr[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=True)
            axarr[i].tick_params(labelsize=ticklabelsize + 1)
            if scol:
                # plt.setp(axarr[i].get_yticklabels(),visible=False)
                axarr[i].minorticks_off()

                # axarr[i].tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
            else:
                plt.setp(axarr[i].get_yticklabels(), visible=True)
                axarr[i].tick_params(axis='y', which='both', left=True, right=False, labelbottom=True)
                # plt.setp(axarr[i].get_yticks(),visible=False)

        axarr[-1].legend(legends[0:len(self.color_tags)], self.color_tags, loc='center left', handlelength=1,
                         handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size': ticklabelsize + 1},
                         bbox_to_anchor=(1.05, 0.5))
        # f.tight_layout(pad=2, h_pad=None, w_pad=None)
        # f.tight_layout()
        self.fig = f

    def gen_html(self, directory, title, align=50):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = title
        link_d = OrderedDict()
        link_d["Boxplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        # fp = os.path.join(dir,outputname,title)

        html.add_figure("boxplot.png", align="center")

        type_list = 'ssssssssssssssssssssssssssssssssssssssssssssss'

        #### Calculate p value ####
        plist = {}
        for g in list(self.sortDict.keys()):
            plist[g] = {}
            for s1 in list(self.sortDict[g].keys()):
                for c1 in list(self.sortDict[g][s1].keys()):
                    data1 = self.sortDict[g][s1][c1]
                    plist[g][s1 + c1] = {}
                    for s2 in list(self.sortDict[g].keys()):
                        for c2 in list(self.sortDict[g][s2].keys()):
                            if s2 == s1 and c2 == c1:
                                pass
                            else:
                                data2 = self.sortDict[g][s2][c2]
                                u, p_value = mannwhitneyu(data1, data2)
                                plist[g][s1 + c1][s2 + c2] = p_value

        print("Multiple test correction.")
        multiple_correction(plist)

        for g in list(self.sortDict.keys()):
            html.add_heading(g, size=4, bold=False)
            data_table = []
            col_size_list = [15]
            header_list = ["p-value"]
            for s in list(self.sortDict[g].keys()):
                for c in list(self.sortDict[g][s1].keys()):
                    header_list.append(s + "\n" + c)
                    col_size_list.append(15)

            for s1 in list(self.sortDict[g].keys()):
                for c1 in list(self.sortDict[g][s1].keys()):
                    row = [s1 + "\n" + c1]
                    for s2 in list(self.sortDict[g].keys()):
                        for c2 in list(self.sortDict[g][s2].keys()):
                            if s2 == s1 and c2 == c1:
                                row.append("-")
                            else:
                                p = plist[g][s1 + c1][s2 + c2]
                                if p > 0.05:
                                    row.append(value2str(p))
                                else:
                                    row.append("<font color=\"red\">" + value2str(p) + "</font>")
                    data_table.append(row)

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align + 50)

        # html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        header_list = ["Assumptions and hypothesis"]
        col_size_list = [50]
        data_table = [['All the regions among different BED files are normalized by quantile normalization.'],
                      [
                          'If there is any grouping problem, please check all the optional columns in input experimental matrix.']]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align, cell_align="left")

        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])
        html.write(os.path.join(directory, title, "parameters.html"))
