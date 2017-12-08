# Python Libraries
from __future__ import print_function
from __future__ import division
import time
import numpy
import matplotlib.pyplot as plt
import os
from collections import OrderedDict, defaultdict
from scipy.stats import stats

# Local Libraries
# Distal Libraries
from rgt.Util import Html
from rgt.CoverageSet import CoverageSet
from rgt.ExperimentalMatrix import ExperimentalMatrix
from shared_function import gen_tags, tag_from_r, print2, MyPool, compute_coverage, colormap, unique, output_array


# Local test
# dir_path = os.getcwd()


###########################################################################################
#                    Lineplot
###########################################################################################


class Lineplot:
    def __init__(self, em_path, title, annotation, organism, center, extend, rs, bs, ss,
                 df, dft, fields, test, sense, strand, flipnegative):

        # Read the Experimental Matrix
        self.title = title
        self.exps = ExperimentalMatrix()
        self.exps.read(em_path, test=test, add_region_len=True)
        for f in self.exps.fields:
            if f not in ['name', 'type', 'file', "reads", "regions", "factors"]:
                self.exps.match_ms_tags(f, test=test)
                self.exps.remove_name()
        # print(self.exps.fieldsDict)
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
        self.strand = strand
        self.flipnegative = flipnegative

        self.processed_beds = []
        self.processed_bedsF = []  # Processed beds to be flapped

    def relocate_bed(self):

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
                allbed = bed.relocate_regions(center=self.center,
                                              left_length=self.extend + int(0.5 * self.bs) + 2 * self.ss,
                                              right_length=self.extend + int(0.5 * self.bs) + 2 * self.ss)
                self.processed_beds.append(allbed)

    def group_tags(self, groupby, rowby, columnby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        self.tag_type = [groupby, rowby, columnby, colorby, self.dft]
        if "None" in self.tag_type:
            self.tag_type.remove("None")

        if groupby == "None":
            self.group_tags = [""]
        elif groupby == "regions" and self.annotation:
            self.group_tags = self.bednames
        else:
            self.group_tags = gen_tags(self.exps, groupby, region_len=False)

        if rowby == "None":
            self.row_tags = [""]
        elif rowby == "regions" and self.annotation:
            self.row_tags = self.bednames
        else:
            self.row_tags = gen_tags(self.exps, rowby, region_len=False)

        if columnby == "None":
            self.column_tags = [""]
        elif columnby == "regions" and self.annotation:
            self.column_tags = self.bednames
        else:
            self.column_tags = gen_tags(self.exps, columnby, region_len=False)

        if colorby == "None":
            self.color_tags = [""]
        elif colorby == "regions" and self.annotation:
            self.color_tags = self.bednames
        else:
            self.color_tags = gen_tags(self.exps, colorby, region_len=False)

        print("\tGroup labels:\t" + ",".join(self.group_tags))
        print("\tRow labels:\t" + ",".join(self.row_tags))
        print("\tColumn labels:\t" + ",".join(self.column_tags))
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
            # print(self.cuebed[bed])
            try:
                self.cuebed[bed].remove("None")
            except:
                pass
        for bam in self.readsnames:
            self.cuebam[bam] = set(tag_from_r(self.exps, self.tag_type, bam))

    def coverage(self, sortby, heatmap=False, logt=False, mp=0, log=False, average=False):

        def annot_ind(bednames, tags):
            """Find the index for annotation tag"""
            for ind, a_region in enumerate(bednames):
                if a_region in tags: return ind

        if mp > 0:
            ts = time.time()
        normRPM = False
        # Calculate for coverage
        mp_input = []
        data = OrderedDict()

        bi = 0
        for g in self.group_tags:
            data[g] = OrderedDict()
            for r in self.row_tags:
                data[g][r] = OrderedDict()
                for c in self.column_tags:
                    data[g][r][c] = OrderedDict()
                    for cc in self.color_tags:
                        # if self.df: data[s][g][c] = []
                        data[g][r][c][cc] = OrderedDict()
                        if not self.dft:
                            dfs = [cc]
                        else:
                            if self.dft == "regions":
                                dfs = self.exps.get_regionsnames()
                            elif self.dft == "reads":
                                dfs = self.exps.get_readsnames()
                            else:
                                dfs = self.exps.fieldsDict[self.dft].keys()

                        for d in dfs:
                            data[g][r][c][cc][d] = defaultdict(list)
                            for bed in self.cuebed.keys():
                                # print(self.cuebed[bed])
                                # print(set([s,g,c,d]))
                                # print(self.cuebed[bed].issubset(set([s,g,c,d])))
                                if len(self.cuebed[bed].intersection({g, r, c, cc, d})) > 2 or self.cuebed[
                                    bed].issubset(
                                    {g, r, c, cc, d}):
                                    # if self.cuebed[bed] <= set([s,g,c]):
                                    for bam in self.cuebam.keys():

                                        # print(self.cuebam[bam])
                                        # print(set([s,g,c]))
                                        if self.cuebam[bam] <= {g, r, c, cc, d}:
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
                                                                 g, r, c, cc, d])
                                                data[g][r][c][cc][d] = None

                                            #########################################################################
                                            else:  # Single thread
                                                ts = time.time()
                                                cov = CoverageSet(bed + "." + bam, self.processed_beds[i])

                                                # print(len(self.processed_beds[i]))
                                                if "Conservation" in [g, r, c, cc, d]:
                                                    cov.phastCons46way_score(stepsize=self.ss)

                                                elif ".bigwig" in self.reads[j].lower() or ".bw" in self.reads[
                                                    j].lower():
                                                    cov.coverage_from_bigwig(bigwig_file=self.reads[j],
                                                                             stepsize=self.ss)
                                                else:
                                                    if not self.sense and not self.strand:
                                                        cov.coverage_from_bam(bam_file=self.reads[j],
                                                                              extension_size=self.rs, binsize=self.bs,
                                                                              stepsize=self.ss)
                                                        if normRPM: cov.normRPM()
                                                    else:  # Sense specific
                                                        cov.coverage_from_bam(bam_file=self.reads[j],
                                                                              extension_size=self.rs, binsize=self.bs,
                                                                              stepsize=self.ss,
                                                                              get_sense_info=self.sense,
                                                                              get_strand_info=self.strand,
                                                                              paired_reads=True)
                                                        cov.array_transpose()
                                                        if normRPM: cov.normRPM()

                                                if self.center == "midpoint" and self.flipnegative:
                                                    for k, re in enumerate(self.processed_beds[i]):
                                                        if re.orientation == "-":
                                                            # print(k)
                                                            # print(cov.coverage[k])
                                                            cov.coverage[k] = cov.coverage[k][::-1]

                                                # When bothends, consider the fliping end
                                                if self.center == 'bothends' or self.center == 'upstream' or self.center == 'downstream':
                                                    if "Conservation" in [g, r, c, cc, d]:
                                                        flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                        flap.phastCons46way_score(stepsize=self.ss)
                                                        ffcoverage = numpy.fliplr(flap.coverage)
                                                        cov.coverage = numpy.concatenate((cov.coverage, ffcoverage),
                                                                                         axis=0)
                                                    elif ".bigwig" in self.reads[j].lower() or ".bw" in self.reads[
                                                        j].lower():
                                                        flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                        flap.coverage_from_bigwig(bigwig_file=self.reads[j],
                                                                                  stepsize=self.ss)
                                                        ffcoverage = numpy.fliplr(flap.coverage)
                                                        cov.coverage = numpy.concatenate((cov.coverage, ffcoverage),
                                                                                         axis=0)
                                                    else:
                                                        flap = CoverageSet("for flap", self.processed_bedsF[i])
                                                        if not self.sense:
                                                            flap.coverage_from_bam(self.reads[j],
                                                                                   extension_size=self.rs,
                                                                                   binsize=self.bs, stepsize=self.ss)
                                                            if normRPM: flap.normRPM()
                                                        else:  # Sense specific
                                                            flap.coverage_from_bam(bam_file=self.reads[j],
                                                                                   extension_size=self.rs,
                                                                                   binsize=self.bs,
                                                                                   stepsize=self.ss,
                                                                                   get_sense_info=True,
                                                                                   paired_reads=True)
                                                            flap.array_transpose(flip=True)
                                                            if normRPM: flap.normRPM()
                                                        ffcoverage = numpy.fliplr(flap.coverage)
                                                        try:
                                                            cov.coverage = numpy.concatenate((cov.coverage, ffcoverage),
                                                                                             axis=0)
                                                        except:
                                                            pass

                                                        if self.sense:
                                                            cov.transpose_cov1 = numpy.concatenate((cov.transpose_cov1,
                                                                                                    flap.transpose_cov1),
                                                                                                   axis=0)
                                                            cov.transpose_cov2 = numpy.concatenate((cov.transpose_cov2,
                                                                                                    flap.transpose_cov2),
                                                                                                   axis=0)

                                                # Averaging the coverage of all regions of each bed file
                                                if heatmap:
                                                    if logt:
                                                        data[g][r][c][cc][d] = numpy.log10(numpy.vstack(
                                                            cov.coverage) + 1)  # Store the array into data list
                                                    else:
                                                        data[g][r][c][cc][d] = numpy.vstack(
                                                            cov.coverage)  # Store the array into data list
                                                else:
                                                    if len(cov.coverage) == 0:
                                                        data[g][r][c][cc][d] = None
                                                        print("** Warning: Cannot open " + self.reads[j])
                                                        continue
                                                    else:
                                                        for i, car in enumerate(cov.coverage):
                                                            if i == 0:
                                                                avearr = numpy.array(car, ndmin=2)
                                                            else:
                                                                # avearr = numpy.vstack((avearr, np.array(car, ndmin=2)))
                                                                try:
                                                                    avearr = numpy.vstack(
                                                                        (avearr, numpy.array(car, ndmin=2)))
                                                                except:
                                                                    print(bed + "." + bam + "." + str(i))
                                                        if log:
                                                            avearr = numpy.log10(avearr + 1)

                                                        avearr = numpy.average(avearr, axis=0)
                                                        if self.sense or self.strand:
                                                            if log:
                                                                sense_1 = numpy.average(
                                                                    numpy.log2(cov.transpose_cov1 + 1), axis=0)
                                                                sense_2 = numpy.average(
                                                                    numpy.log2(cov.transpose_cov2 + 1), axis=0)
                                                            else:
                                                                sense_1 = numpy.average(cov.transpose_cov1, axis=0)
                                                                sense_2 = numpy.average(cov.transpose_cov2, axis=0)
                                                        cut_end = int(self.bs / self.ss)
                                                        avearr = avearr[cut_end:-cut_end]
                                                        data[g][r][c][cc][d]["all"].append(avearr)

                                                        if self.sense or self.strand:
                                                            sense_1 = sense_1[cut_end:-cut_end]
                                                            sense_2 = sense_2[cut_end:-cut_end]
                                                            data[g][r][c][cc][d]["sense_1"].append(sense_1)
                                                            data[g][r][c][cc][d]["sense_2"].append(sense_2)

                                                bi += 1
                                                te = time.time()
                                                print2(self.parameter,
                                                       "\t" + str(bi) + "\t" + "{0:30}\t--{1:<5.1f}\tsec".format(
                                                           bed + "." + bam, ts - te))

        if mp > 0:
            pool = MyPool(mp)
            mp_output = pool.map(compute_coverage, mp_input)
            pool.close()
            pool.join()
            for g in data.keys():
                for r in data[g].keys():
                    for c in data[g][r].keys():
                        for cc in data[g][r][c].keys():
                            for d in data[g][r][c][cc].keys():
                                for out in mp_output:
                                    if out[0] == g and out[1] == r and out[2] == c and out[3] == cc and out[3] == d:
                                        if self.df:
                                            try:
                                                data[g][r][c][cc][d][-1].append(out[4])
                                            except:
                                                data[g][r][c][cc][d] = [[out[4]]]
                                        else:
                                            try:
                                                data[g][r][c][cc][d].append(out[4])
                                            except:
                                                data[g][r][c][cc][d] = [out[4]]
        if average:
            for g in data.keys():
                for r in data[g].keys():
                    for c in data[g][r].keys():
                        for cc in data[g][r][c].keys():
                            for d in data[g][r][c][cc].keys():
                                if isinstance(data[g][r][c][cc][d]["all"], list) and len(
                                        data[g][r][c][cc][d]["all"]) > 1:
                                    a = numpy.array(data[g][r][c][cc][d]["all"])
                                    averaged_array = numpy.array(numpy.average(a, axis=0))
                                    # print(averaged_array)
                                    # sys.exit(1)
                                    data[g][r][c][cc][d]["all"] = [averaged_array]
                                    # print(len(data[s][g][c][d]["all"]))
        if self.df:
            for g in data.keys():
                for r in data[g].keys():
                    for c in data[g][r].keys():
                        for cc in data[g][r][c].keys():
                            for d in data[g][r][c][cc].keys():
                                if isinstance(data[g][r][c][cc][d]["all"], list) and len(
                                        data[g][r][c][cc][d]["all"]) > 1:
                                    diff = numpy.subtract(data[g][r][c][cc][d]["all"][0],
                                                          data[g][r][c][cc][d]["all"][1])
                                    data[g][r][c][cc][d]["df"].append(diff.tolist())
                                else:
                                    print("Warning: There is no repetitive reads for calculating difference.\n"
                                          "         Please add one more entry in experimental matrix.")

        self.data = data

    def colormap(self, colorby, definedinEM):
        colors = colormap(self.exps, colorby, definedinEM, annotation=self.annotation)
        self.colors = {}
        for i, c in enumerate(self.color_tags):
            self.colors[c] = colors[i]

    def plot(self, output, printtable=False, scol=False, srow=False, w=2, h=2, ylog=False):
        linewidth = 1
        self.fig = []

        rot = 30
        if len(self.data.values()[0].keys()) < 2:
            ticklabelsize = w * 1.5
        else:
            ticklabelsize = w * 3
        tw = len(self.data.values()[0].keys()) * w
        th = len(self.data.keys()) * (h * 0.8)

        for g in self.group_tags:

            f, axs = plt.subplots(len(self.data[g].keys()), len(self.data[g].values()[0].keys()), dpi=300,
                                  figsize=(tw, th))

            yaxmax = [0] * len(self.data[g].values()[0])
            sx_ymax = [0] * len(self.data[g].keys())
            if self.df:
                yaxmin = [0] * len(self.data[g].values()[0])
                sx_ymin = [0] * len(self.data[g].keys())

            if printtable:
                bott = self.extend + int(0.5 * self.ss)
                pArr = [["Group_tag", "Sort_tag", "Color_tag", "Diff"] + [str(x) for x in
                                                                          range(-bott, bott + 10, self.ss)]]  # Header
            nit = len(self.data[g].keys())
            for ir, r in enumerate(self.data[g].keys()):

                for ic, c in enumerate(self.data[g][r].keys()):
                    try:
                        ax = axs[ir, ic]
                    except:
                        if len(self.data[g].keys()) == 1 and len(self.data[g][r].keys()) == 1:
                            ax = axs
                        elif len(self.data[g].keys()) == 1 and len(self.data[g][r].keys()) > 1:
                            ax = axs[ic]
                        else:
                            ax = axs[ir]

                    if ir == 0:
                        if self.df:
                            ax.set_title(c + "_df", fontsize=ticklabelsize + 2)
                        else:
                            ax.set_title(c, fontsize=ticklabelsize + 2)

                    # Processing for future output
                    for j, cc in enumerate(self.data[g][r][c].keys()):

                        for k, d in enumerate(self.data[g][r][c][cc].keys()):
                            if not self.data[g][r][c][cc][d]:
                                continue
                            else:
                                if not self.sense and not self.strand:
                                    if self.df:
                                        pt = self.data[g][r][c][cc][d]["df"]
                                    else:
                                        pt = self.data[g][r][c][cc][d]["all"]

                                    for l, y in enumerate(pt):
                                        # print(y)
                                        yaxmax[ic] = max(numpy.amax(y), yaxmax[ic])
                                        sx_ymax[ir] = max(numpy.amax(y), sx_ymax[ir])
                                        if self.df:
                                            yaxmin[ic] = min(numpy.amin(y), yaxmin[ic])
                                            sx_ymin[ir] = min(numpy.amin(y), sx_ymin[ir])

                                        x = numpy.linspace(-self.extend, self.extend, len(y))
                                        ax.plot(x, y, color=self.colors[cc], lw=linewidth, label=cc)
                                        if ir < nit - 1:
                                            ax.set_xticklabels([])
                                        # Processing for future output
                                        if printtable:
                                            pArr.append([g, r, c, cc, d] + list(y))
                                else:
                                    if self.sense:
                                        plt.text(0.5, 0.51, 'sense', transform=ax.transAxes, fontsize=ticklabelsize,
                                                 horizontalalignment='center', verticalalignment='bottom')
                                        plt.text(0.5, 0.49, 'anti-sense', transform=ax.transAxes,
                                                 fontsize=ticklabelsize,
                                                 horizontalalignment='center', verticalalignment='top')
                                    elif self.strand:
                                        plt.text(0.5, 0.51, 'Forward strand', transform=ax.transAxes,
                                                 fontsize=ticklabelsize,
                                                 horizontalalignment='center', verticalalignment='bottom')
                                        plt.text(0.5, 0.49, 'Reverse strand', transform=ax.transAxes,
                                                 fontsize=ticklabelsize,
                                                 horizontalalignment='center', verticalalignment='top')
                                    plt.plot((-self.extend, self.extend), (0, 0), '0.1', linewidth=0.2)
                                    # print(self.data[s][g][c][d])
                                    for l, y in enumerate(self.data[g][r][c][cc][d]["sense_1"]):
                                        # print(y)
                                        ymax1 = numpy.amax(y)
                                        yaxmax[ic] = max(ymax1, yaxmax[ic])
                                        sx_ymax[ir] = max(ymax1, sx_ymax[ir])
                                        x = numpy.linspace(-self.extend, self.extend, y.shape[0])
                                        ax.plot(x, y, color=self.colors[c], lw=linewidth, label=c)
                                        if ir < nit - 1: ax.set_xticklabels([])
                                        # Processing for future output
                                        if printtable: pArr.append([g, r, c, cc, d, "+"] + list(y))

                                    for l, y in enumerate(self.data[g][r][c][cc][d]["sense_2"]):
                                        # print(y)
                                        ymax2 = numpy.amax(y)
                                        yaxmax[ic] = max(ymax2, yaxmax[ic])
                                        sx_ymax[ir] = max(ymax2, sx_ymax[ir])
                                        x = numpy.linspace(-self.extend, self.extend, y.shape[0])
                                        ax.plot(x, -y, color=self.colors[c], lw=linewidth, label=c)
                                        if ir < nit - 1: ax.set_xticklabels([])
                                        # Processing for future output
                                        if printtable: pArr.append([g, r, c, cc, d, "-"] + list(y))
                                    ym = 1.2 * max(max(yaxmax), max(sx_ymax))
                                    ax.set_ylim([-ym, ym])

                    # ax.get_yaxis().set_label_coords(-0.1, 0.5)
                    ax.set_xlim([-self.extend, self.extend])
                    plt.setp(ax.get_xticklabels(), fontsize=ticklabelsize, rotation=rot, ha='right')
                    plt.setp(ax.get_yticklabels(), fontsize=ticklabelsize)
                    ax.locator_params(axis='x', nbins=4)
                    ax.locator_params(axis='y', nbins=3)

            if printtable:
                output_array(pArr, directory=output, folder=self.title, filename="plot_table_" + g + ".txt")

            handles = []
            labels = []
            for ir, r in enumerate(self.data[g].keys()):
                if ylog:
                    nr = r + " (log10)"
                else:
                    nr = r
                try:
                    axs[ir, 0].set_ylabel("{}".format(nr), fontsize=ticklabelsize + 1)
                    axs[ir, 0].get_yaxis().set_label_coords(-0.1, 0.5)
                except:
                    try:
                        axs[ir].set_ylabel("{}".format(nr), fontsize=ticklabelsize + 1)
                        axs[ir].get_yaxis().set_label_coords(-0.1, 0.5)
                    except:
                        axs.set_ylabel("{}".format(nr), fontsize=ticklabelsize + 1)
                        axs.get_yaxis().set_label_coords(-0.1, 0.5)

                for ic, c in enumerate(self.data[g][r].keys()):
                    try:
                        axx = axs[ir, ic]
                    except:
                        try:
                            if len(self.data[g].keys()) == 1:
                                axx = axs[ic]
                            else:
                                axx = axs[ir]
                        except:
                            axx = axs

                    if self.df:
                        if scol and not srow:
                            ymin = yaxmin[ic] - abs(yaxmin[ic] * 0.2)
                            ymax = yaxmax[ic] + abs(yaxmax[ic] * 0.2)
                        elif srow and not scol:
                            ymin = sx_ymin[ir] - abs(sx_ymin[ir] * 0.2)
                            ymax = sx_ymax[ir] + abs(sx_ymax[ir] * 0.2)
                        elif scol and srow:
                            ymin = min(yaxmin[ic], sx_ymin[ir]) - abs(min(yaxmin[ic], sx_ymin[ir]) * 0.2)
                            ymax = max(yaxmax[ic], sx_ymax[ir]) + abs(max(yaxmax[ic], sx_ymax[ir]) * 0.2)

                    else:
                        if scol and not srow:
                            ymax = yaxmax[ic] * 1.2
                        elif srow and not scol:
                            ymax = sx_ymax[ir] * 1.2
                        elif scol and srow:
                            ymax = max(max(yaxmax), max(sx_ymax)) * 1.2
                        else:
                            ymax = axx.get_ylim()[1]
                        if self.sense or self.strand:
                            ymin = -ymax
                        else:
                            ymin = 0

                    try:
                        axx.set_ylim([ymin, ymax])
                    except:
                        pass
                    hand, lab = axx.get_legend_handles_labels()
                    handles += hand
                    labels += lab
            # handles, labels = ax.get_legend_handles_labels()
            uniq_labels = unique(labels)

            plt.legend([handles[labels.index(l)] for l in uniq_labels], uniq_labels, loc='center left', handlelength=1,
                       handletextpad=1,
                       columnspacing=2, borderaxespad=0., prop={'size': ticklabelsize}, bbox_to_anchor=(1.05, 0.5))

            f.tight_layout()
            self.fig.append(f)

    def gen_html(self, directory, title, align=50):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = dir_name + " / " + title
        link_d = OrderedDict()
        link_d["Lineplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        for g in self.group_tags:
            html.add_heading(heading=g)
            html.add_figure("lineplot_" + g + ".png", align="center", width="80%")

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
        if not sort:
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
                        d = numpy.empty(shape=self.data[t][g][c].shape)
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
                        d = numpy.empty(shape=self.data[t][g][c].shape)
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
