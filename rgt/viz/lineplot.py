# Python Libraries
from __future__ import print_function
from __future__ import division
import time
import numpy
import matplotlib.pyplot as plt

# Local Libraries
# Distal Libraries
from rgt.Util import Html
from rgt.CoverageSet import *
from rgt.ExperimentalMatrix import *
from shared_function import gen_tags, tag_from_r, print2, MyPool, compute_coverage, colormap, unique, output_array
# Local test
dir = os.getcwd()

###########################################################################################
#                    Lineplot
###########################################################################################


class Lineplot:
    def __init__(self, EMpath, title, annotation, organism, center, extend, rs, bs, ss, df, dft, fields, test, sense):

        # Read the Experimental Matrix
        self.title = title
        self.exps = ExperimentalMatrix()
        self.exps.read(EMpath, test=test)
        for f in self.exps.fields:
            if f not in ['name', 'type', 'file', "reads", "regions", "factors"]:
                self.exps.match_ms_tags(f, test=test)
                self.exps.remove_name()

        # if annotation:
        #     self.beds, self.bednames, self.annotation = annotation_dump(organism)

        # else:
        self.beds = self.exps.get_regionsets()  # A list of GenomicRegionSets
        self.bednames = self.exps.get_regionsnames()
        self.annotation = None

        self.reads = self.exps.get_readsfiles()
        self.readsnames = self.exps.get_readsnames()
        self.fieldsDict = self.exps.fieldsDict
        self.parameter = []
        self.center = center
        self.extend = extend
        self.rs = rs
        self.bs = bs
        self.ss = ss
        self.df = df
        self.dft = dft
        self.sense = sense

    def relocate_bed(self):
        self.processed_beds = []
        self.processed_bedsF = []  # Processed beds to be flapped

        for bed in self.beds:
            if self.center == 'bothends':
                newbed = bed.relocate_regions(center='leftend',
                                              left_length=self.extend + self.bs,
                                              right_length=self.extend + self.bs)
                self.processed_beds.append(newbed)
                newbedF = bed.relocate_regions(center='rightend',
                                               left_length=self.extend + self.bs,
                                               right_length=self.extend + self.bs)
                self.processed_bedsF.append(newbedF)
            elif self.center == 'upstream' or self.center == 'downstream':
                allbed = bed.relocate_regions(center=self.center,
                                              left_length=self.extend + self.bs,
                                              right_length=self.extend + self.bs)
                newbed = allbed.filter_strand(strand="+")
                self.processed_beds.append(newbed)
                newbedF = allbed.filter_strand(strand="-")
                self.processed_bedsF.append(newbedF)
            else:
                newbed = bed.relocate_regions(center=self.center,
                                              left_length=self.extend + int(0.5 * self.bs) + 2 * self.ss,
                                              right_length=self.extend + int(0.5 * self.bs) + 2 * self.ss)
                self.processed_beds.append(newbed)

    def group_tags(self, groupby, sortby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        self.tag_type = [sortby, groupby, colorby, self.dft]
        if "None" in self.tag_type: self.tag_type.remove("None")

        if groupby == "None":
            self.group_tags = [""]
        elif groupby == "regions" and self.annotation:
            self.group_tags = self.bednames
        else:
            self.group_tags = gen_tags(self.exps, groupby)

        if sortby == "None":
            self.sort_tags = [""]
        elif sortby == "regions" and self.annotation:
            self.sort_tags = self.bednames
        else:
            self.sort_tags = gen_tags(self.exps, sortby)

        if colorby == "None":
            self.color_tags = [""]
        elif colorby == "regions" and self.annotation:
            self.color_tags = self.bednames
        else:
            self.color_tags = gen_tags(self.exps, colorby)

        print("\tColumn labels:\t" + ",".join(self.group_tags))
        print("\tRow labels:\t" + ",".join(self.sort_tags))
        print("\tColor labels:\t" + ",".join(self.color_tags))

    def gen_cues(self):
        self.cuebed = OrderedDict()
        self.cuebam = OrderedDict()

        # if self.annotation:
        #     #all_tags = []
        #     #for dictt in self.exps.fieldsDict.values():
        #     #    for tag in dictt.keys():
        #     #        all_tags.append(tag)
        #     for bed in self.bednames:
        #     #    self.cuebed[bed] = set([bed]+all_tags)
        #         self.cuebed[bed] = set([bed])
        # else:
        for bed in self.bednames:
            self.cuebed[bed] = set(tag_from_r(self.exps, self.tag_type, bed))
            try:
                self.cuebed[bed].remove("None")
            except:
                pass
        for bam in self.readsnames:
            self.cuebam[bam] = set(tag_from_r(self.exps, self.tag_type, bam))

    def coverage(self, sortby, heatmap=False, logt=False, mp=0, log=False):

        def annot_ind(bednames, tags):
            """Find the index for annotation tag"""
            for ind, a in enumerate(bednames):
                if a in tags: return ind

        if mp>0: ts = time.time()
        normRPM = False
        # Calculate for coverage
        mp_input = []
        data = OrderedDict()

        bi = 0
        for s in self.sort_tags:
            data[s] = OrderedDict()
            for g in self.group_tags:
                data[s][g] = OrderedDict()
                for c in self.color_tags:
                    # if self.df: data[s][g][c] = []
                    data[s][g][c] = OrderedDict()
                    if not self.dft:
                        dfs = [c]
                    else:
                        dfs = self.exps.fieldsDict[self.dft].keys()
                    for d in dfs:
                        data[s][g][c][d] = defaultdict(list)
                        for bed in self.cuebed.keys():
                            # print(self.cuebed[bed])
                            # print(set([s,g,c,d]))
                            # print(self.cuebed[bed].issubset(set([s,g,c,d])))
                            if len(self.cuebed[bed].intersection(set([s, g, c, d]))) > 2 or self.cuebed[bed].issubset(
                                    set([s, g, c, d])):
                                # if self.cuebed[bed] <= set([s,g,c]):
                                for bam in self.cuebam.keys():

                                    # print(self.cuebam[bam])
                                    # print(set([s,g,c]))
                                    if self.cuebam[bam] <= set([s, g, c, d]):
                                        i = self.bednames.index(bed)
                                        j = self.readsnames.index(bam)
                                        # print(bed + "." + bam)

                                        # if len(self.processed_beds[i]) == 0:
                                        #     try:
                                        #         data[s][g][c][d].append(numpy.empty(1, dtype=object))
                                        #     except:
                                        #         data[s][g][c][d] = [numpy.empty(1, dtype=object)]
                                        #     continue
                                        #########################################################################
                                        if mp > 0:  # Multiple processing
                                            mp_input.append([self.processed_beds[i], self.reads[j],
                                                             self.rs, self.bs, self.ss, self.center, heatmap, logt,
                                                             s, g, c, d])
                                            data[s][g][c][d] = None

                                        #########################################################################
                                        else:  # Single thread
                                            ts = time.time()
                                            cov = CoverageSet(bed + "." + bam, self.processed_beds[i])

                                            # print(len(self.processed_beds[i]))
                                            if "Conservation" in [s,g,c,d]:
                                                cov.phastCons46way_score(stepsize=self.ss)

                                            elif ".bigwig" in self.reads[j].lower() or ".bw" in self.reads[j].lower():
                                                cov.coverage_from_bigwig(bigwig_file=self.reads[j], stepsize=self.ss)
                                            else:
                                                if not self.sense:
                                                    cov.coverage_from_bam(bam_file=self.reads[j],
                                                                          extension_size=self.rs, binsize=self.bs,
                                                                          stepsize=self.ss)
                                                    if normRPM: cov.normRPM()
                                                else:  # Sense specific
                                                    cov.coverage_from_bam(bam_file=self.reads[j],
                                                                          extension_size=self.rs, binsize=self.bs,
                                                                          stepsize=self.ss, get_sense_info=True,
                                                                          paired_reads=True)
                                                    cov.array_transpose()
                                                    if normRPM: cov.normRPM()

                                            # When bothends, consider the fliping end
                                            if self.center == 'bothends' or self.center == 'upstream' or self.center == 'downstream':
                                                if "Conservation" in [s,g,c,d]:
                                                    flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                    flap.phastCons46way_score(stepsize=self.ss)
                                                    ffcoverage = numpy.fliplr(flap.coverage)
                                                    cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
                                                elif ".bigwig" in self.reads[j].lower() or ".bw" in self.reads[j].lower():
                                                    flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                    flap.coverage_from_bigwig(bigwig_file=self.reads[j],
                                                                              stepsize=self.ss)
                                                    ffcoverage = numpy.fliplr(flap.coverage)
                                                    cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
                                                else:
                                                    flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                    if not self.sense:
                                                        flap.coverage_from_bam(self.reads[j], extension_size=self.rs,
                                                                               binsize=self.bs, stepsize=self.ss)
                                                        if normRPM: flap.normRPM()
                                                    else:  # Sense specific
                                                        flap.coverage_from_bam(bam_file=self.reads[j],
                                                                               extension_size=self.rs, binsize=self.bs,
                                                                               stepsize=self.ss, get_sense_info=True,
                                                                               paired_reads=True)
                                                        flap.array_transpose(flip=True)
                                                        if normRPM: flap.normRPM()
                                                    ffcoverage = numpy.fliplr(flap.coverage)
                                                    try: cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
                                                    except: pass

                                                    if self.sense:
                                                        cov.transpose_cov1 = numpy.concatenate((cov.transpose_cov1,
                                                                                                flap.transpose_cov1),axis=0)
                                                        cov.transpose_cov2 = numpy.concatenate((cov.transpose_cov2,
                                                                                                flap.transpose_cov2), axis=0)

                                            # Averaging the coverage of all regions of each bed file
                                            if heatmap:
                                                if logt:
                                                    data[s][g][c][d] = numpy.log10(numpy.vstack(
                                                        cov.coverage) + 1)  # Store the array into data list
                                                else:
                                                    data[s][g][c][d] = numpy.vstack(
                                                        cov.coverage)  # Store the array into data list
                                            else:
                                                if len(cov.coverage) == 0:
                                                    data[s][g][c][d] = None
                                                    print("** Warning: Cannot open " + self.reads[j])
                                                    continue
                                                else:
                                                    for i, car in enumerate(cov.coverage):
                                                        if i == 0: avearr = np.array(car, ndmin=2)
                                                        else:
                                                            # avearr = numpy.vstack((avearr, np.array(car, ndmin=2)))
                                                            try: avearr = numpy.vstack((avearr, np.array(car, ndmin=2)))
                                                            except: print(bed+"."+bam+"."+str(i))
                                                    if log:
                                                        avearr = numpy.log2(avearr+1)

                                                    avearr = numpy.average(avearr, axis=0)
                                                    if self.sense:
                                                        if log:
                                                            sense_1 = numpy.average(numpy.log2(cov.transpose_cov1+1), axis=0)
                                                            sense_2 = numpy.average(numpy.log2(cov.transpose_cov2+1), axis=0)
                                                        else:
                                                            sense_1 = numpy.average(cov.transpose_cov1,axis=0)
                                                            sense_2 = numpy.average(cov.transpose_cov2,axis=0)
                                                    cut_end = int(self.bs/self.ss)
                                                    avearr = avearr[cut_end:-cut_end]
                                                    data[s][g][c][d]["all"].append(avearr)

                                                    if self.sense:
                                                        sense_1 = sense_1[cut_end:-cut_end]
                                                        sense_2 = sense_2[cut_end:-cut_end]
                                                        data[s][g][c][d]["sense_1"].append(sense_1)
                                                        data[s][g][c][d]["sense_2"].append(sense_2)

                                            bi += 1
                                            te = time.time()
                                            print2(self.parameter,
                                                   "\t" + str(bi) + "\t" + "{0:30}\t--{1:<5.1f}s".format(
                                                       bed + "." + bam, ts - te))


        if mp > 0:
            pool = MyPool(mp)
            mp_output = pool.map(compute_coverage, mp_input)
            pool.close()
            pool.join()
            for s in data.keys():
                for g in data[s].keys():
                    for c in data[s][g].keys():
                        for d in data[s][g][c].keys():
                            for out in mp_output:
                                if out[0] == s and out[1] == g and out[2] == c and out[3] == d:
                                    if self.df:
                                        try:
                                            data[s][g][c][d][-1].append(out[4])
                                        except:
                                            data[s][g][c][d] = [[out[4]]]
                                    else:
                                        try:
                                            data[s][g][c][d].append(out[4])
                                        except:
                                            data[s][g][c][d] = [out[4]]
        if self.df:
            for s in data.keys():
                for g in data[s].keys():
                    for c in data[s][g].keys():
                        for d in data[s][g][c].keys():
                            if isinstance(data[s][g][c][d]["all"], list) and len(data[s][g][c][d]["all"]) > 1:
                                diff = numpy.subtract(data[s][g][c][d]["all"][0], data[s][g][c][d]["all"][1])
                                data[s][g][c][d]["df"].append(diff.tolist())
                            else:
                                print("Warning: There is no repetitive reads for calculating difference.\n"
                                      "         Please add one more entry in experimental matrix.")
        self.data = data

    def colormap(self, colorby, definedinEM):
        colors = colormap(self.exps, colorby, definedinEM, annotation=self.annotation)
        self.colors = {}
        for i, c in enumerate(self.color_tags):
            self.colors[c] = colors[i]

    def plot(self, groupby, colorby, output, printtable=False, scol=False, srow=False, w=2, h=2):

        rot = 50
        if len(self.data.values()[0].keys()) < 2:
            ticklabelsize = w * 1.5
        else:
            ticklabelsize = w * 3
        tw = len(self.data.values()[0].keys()) * w
        th = len(self.data.keys()) * (h * 0.8)

        f, axs = plt.subplots(len(self.data.keys()), len(self.data.values()[0].keys()), dpi=300,
                              figsize=(tw, th))

        yaxmax = [0] * len(self.data.values()[0])
        sx_ymax = [0] * len(self.data.keys())
        if self.df:
            yaxmin = [0] * len(self.data.values()[0])
            sx_ymin = [0] * len(self.data.keys())

        if printtable:
            bott = self.extend + int(0.5 * self.ss)
            pArr = [["Group_tag", "Sort_tag", "Color_tag", "Diff"] + [str(x) for x in range(-bott, bott + 10, self.ss)]]  # Header
        nit = len(self.data.keys())
        for it, s in enumerate(self.data.keys()):

            for i, g in enumerate(self.data[s].keys()):
                try:
                    ax = axs[it, i]
                except:
                    if len(self.data.keys()) == 1 and len(self.data[s].keys()) == 1:
                        ax = axs
                    elif len(self.data.keys()) == 1 and len(self.data[s].keys()) > 1:
                        ax = axs[i]
                    else:
                        ax = axs[it]

                if it == 0:
                    if self.df:
                        ax.set_title(g + "_df", fontsize=ticklabelsize + 2)
                    else:
                        ax.set_title(g, fontsize=ticklabelsize + 2)

                # Processing for future output
                for j, c in enumerate(self.data[s][g].keys()):

                    for k, d in enumerate(self.data[s][g][c].keys()):
                        if not self.data[s][g][c][d]:
                            continue
                        else:
                            if not self.sense:
                                if self.df: pt = self.data[s][g][c][d]["df"]
                                else: pt = self.data[s][g][c][d]["all"]

                                for l, y in enumerate(pt):
                                    # print(y)
                                    yaxmax[i] = max(numpy.amax(y), yaxmax[i])
                                    sx_ymax[it] = max(numpy.amax(y), sx_ymax[it])
                                    if self.df:
                                        yaxmin[i] = min(numpy.amin(y), yaxmin[i])
                                        sx_ymin[it] = min(numpy.amin(y), sx_ymin[it])

                                    x = numpy.linspace(-self.extend, self.extend, len(y))
                                    ax.plot(x, y, color=self.colors[c], lw=1, label=c)
                                    if it < nit - 1:
                                        ax.set_xticklabels([])
                                    # Processing for future output
                                    if printtable: pArr.append([g, s, c, d] + list(y))
                            else:
                                plt.text(0.5, 0.51, 'sense',transform=ax.transAxes,fontsize=ticklabelsize,
                                         horizontalalignment='center', verticalalignment='bottom')
                                plt.text(0.5, 0.49, 'anti-sense', transform=ax.transAxes,fontsize=ticklabelsize,
                                         horizontalalignment='center', verticalalignment='top')
                                plt.plot((-self.extend, self.extend), (0, 0), '0.1', linewidth=0.2)
                                # print(self.data[s][g][c][d])
                                for l, y in enumerate(self.data[s][g][c][d]["sense_1"]):
                                    # print(y)
                                    ymax1 = numpy.amax(y)
                                    yaxmax[i] = max(ymax1, yaxmax[i])
                                    sx_ymax[it] = max(ymax1, sx_ymax[it])
                                    x = numpy.linspace(-self.extend, self.extend, y.shape[0])
                                    ax.plot(x, y, color=self.colors[c], lw=1, label=c)
                                    if it < nit - 1: ax.set_xticklabels([])
                                    # Processing for future output
                                    if printtable: pArr.append([g, s, c, d, "+"] + list(y))

                                for l, y in enumerate(self.data[s][g][c][d]["sense_2"]):
                                    # print(y)
                                    ymax2 = numpy.amax(y)
                                    yaxmax[i] = max(ymax2, yaxmax[i])
                                    sx_ymax[it] = max(ymax2, sx_ymax[it])
                                    x = numpy.linspace(-self.extend, self.extend, y.shape[0])
                                    ax.plot(x, -y, color=self.colors[c], lw=1, label=c)
                                    if it < nit - 1: ax.set_xticklabels([])
                                    # Processing for future output
                                    if printtable: pArr.append([g, s, c, d, "-"] + list(y))
                                ym = 1.2 * max(max(yaxmax), max(sx_ymax))
                                ax.set_ylim([-ym, ym])

                ax.get_yaxis().set_label_coords(-0.1, 0.5)
                ax.set_xlim([-self.extend, self.extend])
                plt.setp(ax.get_xticklabels(), fontsize=ticklabelsize, rotation=rot)
                plt.setp(ax.get_yticklabels(), fontsize=ticklabelsize)


                ax.locator_params(axis='x', nbins=4)
                ax.locator_params(axis='y', nbins=3)
                # try:
                #
                # except:
                #     ax.locator_params(axis='y', nbins=2)
                #     pass
        if printtable:
            output_array(pArr, directory=output, folder=self.title, filename="plot_table.txt")

        for it, ty in enumerate(self.data.keys()):
            try:
                axs[it, 0].set_ylabel("{}".format(ty), fontsize=ticklabelsize + 1)
            except:
                try:
                    axs[it].set_ylabel("{}".format(ty), fontsize=ticklabelsize + 1)
                except:
                    axs.set_ylabel("{}".format(ty), fontsize=ticklabelsize + 1)

            for i, g in enumerate(self.data[ty].keys()):
                try: axx = axs[it, i]
                except:
                    try:
                        if len(self.data.keys()) == 1:
                            axx = axs[i]
                        else:
                            axx = axs[it]
                    except: axx = axs

                if self.df:
                    if scol and not srow:
                        ymin = yaxmin[i] - abs(yaxmin[i] * 0.2)
                        ymax = yaxmax[i] + abs(yaxmax[i] * 0.2)
                    elif srow and not scol:
                        ymin = sx_ymin[it] - abs(sx_ymin[it] * 0.2)
                        ymax = sx_ymax[it] + abs(sx_ymax[it] * 0.2)
                    elif scol and srow:
                        ymin = min(yaxmin[i], sx_ymin[it]) - abs(min(yaxmin[i], sx_ymin[it]) * 0.2)
                        ymax = max(yaxmax[i], sx_ymax[it]) + abs(max(yaxmax[i], sx_ymax[it]) * 0.2)

                else:
                    if scol and not srow: ymax = yaxmax[i] * 1.2
                    elif srow and not scol: ymax = sx_ymax[it] * 1.2
                    elif scol and srow: ymax = max(max(yaxmax), max(sx_ymax)) * 1.2
                    else:
                        ymax = axx.get_ylim()[1]
                    if self.sense: ymin = -ymax
                    else: ymin = 0

                try: axx.set_ylim([ymin, ymax])
                except: pass

        handles, labels = ax.get_legend_handles_labels()
        uniq_labels = unique(labels)

        plt.legend([handles[labels.index(l)] for l in uniq_labels], uniq_labels, loc='center left', handlelength=1,
                   handletextpad=1,
                   columnspacing=2, borderaxespad=0., prop={'size': ticklabelsize}, bbox_to_anchor=(1.05, 0.5))

        f.tight_layout()
        self.fig = f

    def gen_html(self, directory, title, align=50):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = dir_name + " / " + title
        link_d = OrderedDict()
        link_d["Lineplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        html.add_figure("lineplot.png", align="center", width="80%")

        html.write(os.path.join(directory, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        type_list = 'ssssssssss'
        col_size_list = [20, 20, 20, 20, 20, 20, 20, 20, 20]
        header_list = ["Assumptions and hypothesis"]
        data_table = []
        if self.annotation:
            data_table.append(
                ["Genomic annotation: TSS - Transcription Start Site; TTS - Transcription Termination Site."])
        data_table.append(["Directory:      " + directory.rpartition("/")[2]])
        data_table.append(["Title:          " + title])
        data_table.append(["Extend length:  " + str(self.extend)])
        data_table.append(["Read size:      " + str(self.rs)])
        data_table.append(["Bin size:       " + str(self.bs)])
        data_table.append(["Step size:      " + str(self.ss)])
        data_table.append(["Center mode:    " + self.center])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=align,
                             cell_align="left")

        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])

        html.write(os.path.join(directory, title, "parameters.html"))

    def hmsort(self, sort):
        if sort == None:
            pass
        elif sort == 0:
            for t in self.data.keys():
                for i, g in enumerate(self.data[t].keys()):
                    # print(numpy.sum(data[t][bed].values()[0], axis=1))
                    # print(len(numpy.sum(data[t][bed].values()[0], axis=1)))

                    sumarr = numpy.sum([numpy.sum(d, axis=1) for d in self.data[t][g].values()], axis=0)
                    # print(sumarr)
                    # sumarr = numpy.sum(sumarr, axis=1)
                    ind = stats.rankdata(sumarr, method='ordinal')  # The index for further sorting
                    # numpy.fliplr(ind)

                    for j, c in enumerate(self.data[t][g].keys()):
                        d = numpy.empty(shape=(self.data[t][g][c].shape))
                        for k, ranki in enumerate(ind):
                            d[-ranki, :] = self.data[t][g][c][k, :]
                        self.data[t][g][c] = d
        else:
            for t in self.data.keys():
                for i, g in enumerate(self.data[t].keys()):
                    sumarr = numpy.sum(self.data[t][g].values()[sort - 1], axis=1)
                    # print(sumarr)
                    # sumarr = numpy.sum(sumarr, axis=1)
                    ind = stats.rankdata(sumarr, method='ordinal')  # The index for further sorting
                    # list(ind)
                    # print(ind)
                    for j, c in enumerate(self.data[t][g].keys()):
                        d = numpy.empty(shape=(self.data[t][g][c].shape))
                        for k, ranki in enumerate(ind):
                            d[-ranki, :] = self.data[t][g][c][k, :]
                        self.data[t][g][c] = d
                        # print(data[t][bed].values()[0])

    def hmcmlist(self, colorby, definedinEM):
        # self.colors = colormaps(self.exps, colorby, definedinEM)
        self.colors = ["Reds", "Blues", "Oranges", "Greens", "Purples"]

    def heatmap(self, logt):
        tickfontsize = 6
        ratio = 10
        self.hmfiles = []
        self.figs = []
        for ti, t in enumerate(self.data.keys()):
            # fig.append(plt.figure())
            # rows = len(data[t].keys())
            columns = len(self.data[t].values()[0].keys())
            # fig, axs = plt.subplots(rows,columns, sharey=True, dpi=300)
            # matplotlib.pyplot.subplots_adjust(left=1, right=2, top=2, bottom=1)
            fig = plt.figure(t)
            plt.suptitle("Heatmap: " + t, y=1.05)
            rows = len(self.data[t].keys())

            # gs = gridspec.GridSpec(rows*ratio,columns)
            axs = numpy.empty(shape=(rows + 1, columns), dtype=object)

            for bi, g in enumerate(self.data[t].keys()):
                for bj, c in enumerate(self.data[t][g].keys()):
                    max_value = numpy.amax(self.data[t][g][c])
                    max_value = int(max_value)
                    axs[bi, bj] = plt.subplot2grid(shape=(rows * ratio + 1, columns), loc=(bi * ratio, bj),
                                                   rowspan=ratio)
                    if bi == 0: axs[bi, bj].set_title(c, fontsize=7)
                    # print(self.data[t][g][c])
                    # print(self.colors)
                    # print(bj)
                    # im = axs[bi, bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto',
                    #                        vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])

                    im = axs[bi, bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0, 1], aspect='auto',
                                            vmin=0, vmax=max_value, interpolation='nearest', cmap=plt.get_cmap("Blues"))

                    # for bi, g in enumerate(self.data[t].keys()):
                    #    for bj, c in enumerate(self.data[t][g].keys()):


                    # im = axs[bi, bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto',
                    #                        vmin=0, vmax=max_value, interpolation='nearest', cmap=cm.coolwarm)
                    axs[bi, bj].set_xlim([-self.extend, self.extend])
                    axs[bi, bj].set_xticks([-self.extend, 0, self.extend])
                    # axs[bi, bj].set_xticklabels([-args.e, 0, args.e]
                    plt.setp(axs[bi, bj].get_xticklabels(), fontsize=tickfontsize, rotation=0)
                    # plt.setp(axs[bi, bj].get_yticklabels(), fontsize=10)
                    # axs[bi, bj].locator_params(axis = 'x', nbins = 2)
                    # axs[bi, bj].locator_params(axis = 'y', nbins = 4)
                    for spine in ['top', 'right', 'left', 'bottom']:
                        axs[bi, bj].spines[spine].set_visible(False)
                    axs[bi, bj].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
                    axs[bi, bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')

                    # if bj > 0:
                    #    plt.setp(axs[bi, bj].get_yticklabels(),visible=False)
                    # plt.setp(axarr[i].get_yticks(),visible=False)
                    axs[bi, bj].minorticks_off()
                    if bj == 0:
                        # nregion = len(self.exps.objectsDict[g])
                        # axs[bi, bj].set_ylabel(self.exps.get_type(g,'factor')+" ("+str(nregion) + ")", fontsize=7)
                        axs[bi, bj].set_ylabel(g, fontsize=7)
                    if bi == rows - 1:
                        # divider = make_axes_locatable(axs[bi,bj])
                        # cax = divider.append_axes("bottom", size="5%", pad=0.5)
                        cbar_ax = plt.subplot2grid((rows * ratio + 4, columns), (rows * ratio + 3, bj))
                        # axs[rows,bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')


                        # cbar = grid.cbar_axes[i//2].colorbar(im)
                        # cbar = plt.colorbar(im, cax = axs[rows,bj], ticks=[0, max_value], orientation='horizontal')
                        # cbar = axs[rows,bj].imshow(range(int(max_value)), extent=[0, int(max_value),0,0], aspect=10, extent=[-self.extend, self.extend,0,0]
                        #                           vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])
                        # cbar = axs[rows,bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto',
                        #                    vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])
                        # cbar = axs[rows,bj].imshow([range(2*self.extend),range(2*self.extend),range(2*self.extend)],
                        #                           aspect='auto', vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj] )
                        # cbar.outline.set_linewidth(0.5)
                        # axs[rows,bj].set_ticks_position('none')
                        # axs[rows,bj].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
                        # axs[rows,bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')


                        # cbar.set_label('Amplitute of signal')
                        max_value = int(max_value)
                        # width = 0.4/rows
                        # cbar_ax = fig.add_axes([0.01 + bj/columns, 0, width, 0.01])
                        cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, max_value], orientation='horizontal')
                        cbar.ax.set_xticklabels([0, int(max_value)])
                        if logt:
                            cbar.ax.set_xticklabels(['0', '{:1.1f}'.format(max_value)],
                                                    fontsize=tickfontsize)  # horizontal colorbar
                            cbar.set_label('log10', fontsize=tickfontsize)
                            # else:
                            # cbar.ax.set_xticklabels(['0', int(max_value)], fontsize=tickfontsize)# horizontal colorbar
                            # pass
                            # cbar.outline.set_linewidth(0.1)

            # fig.tight_layout()
            # fig.tight_layout(pad=1.08, h_pad=None, w_pad=None)
            # fig.tight_layout(pad=1, h_pad=1, w_pad=1)
            self.figs.append(fig)
            self.hmfiles.append("heatmap" + "_" + t)

    def gen_htmlhm(self, outputname, title, align=50):
        dir_name = os.path.basename(outputname)
        # check_dir(directory)
        html_header = title
        link_d = OrderedDict()
        link_d["Lineplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        # Each row is a plot with its data
        for name in self.hmfiles:
            html.add_figure(name + ".png", align="center")
        html.write(os.path.join(outputname, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])
        html.write(os.path.join(outputname, title, "parameters.html"))

