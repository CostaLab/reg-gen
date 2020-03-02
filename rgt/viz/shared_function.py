# Python Libraries


import datetime
import fnmatch
import multiprocessing.pool
import pickle
import re
import time
from collections import defaultdict

from shutil import copyfile

import matplotlib
import matplotlib.cm as cmx
import matplotlib.colors as colormat
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages

from ..AnnotationSet import AnnotationSet
from ..CoverageSet import *
from ..ExperimentalMatrix import *
# Local Libraries
# Distal Libraries
from ..Util import Html
from ..motifanalysis.Statistics import multiple_test_correction

if sys.version_info >= (3, 1):
    from urllib.request import urlopen
else:
    import urllib.request, urllib.error, urllib.parse

# Local test
current_dir = os.getcwd()


###########################################################################################
#                    Universal functions
###########################################################################################


def print2(parameter, string):
    """ Show the message on the console and also save in a list for future backup. """
    print(string)
    parameter.append(string)


def unique(a):
    try:
        b = [matplotlib.colors.to_hex(x) for x in a]
        seen = set()
        return [seen.add(x) or x for x in b if x not in seen]
    except:
        seen = set()
        return [seen.add(x) or x for x in a if x not in seen]


def purge(directory, pattern):
    for f in os.listdir(directory):
        if re.search(pattern, f):
            os.remove(os.path.join(directory, f))


def gen_tags(exps, tag, region_len=False):
    """Generate the unique tags from the EM according to the given tag. """
    if "factor" not in exps.fields:
        exps.add_factor_col()
    if tag == "reads":
        try:
            l = [exps.get_type(i, "factor") for i in exps.get_readsnames()]
        except:
            # print("You must define 'factor' column in experimental matrix for grouping.")
            sys.exit(1)
    elif tag == "regions":
        # try:
        if not region_len:
            l = [exps.get_type(i, "factor") for i in exps.get_regionsnames()]
        else:
            l = [exps.get_type(i, "factor") + "(" + str(len(exps.get_regionsets()[i])) + ")" for i, n in
                 enumerate(exps.get_regionsnames())]

        # except:
        #     # print("You must define 'factor' column in experimental matrix for grouping.")
        #     sys.exit(1)
    else:
        l = exps.fieldsDict[tag]
        try:
            l = exps.fieldsDict[tag]
        except:
            print('Cannot find the column "' + tag + '"')
            sys.exit(1)
    try:
        l.remove("ALL")
    except:
        pass
    return unique(l)


def tag_from_r(exps, tag_type, name):
    tags = []
    for t in tag_type:
        if t == "reads" or t == "regions": t = "factor"
        try:
            tags.append(exps.get_type(name, t))
        except:
            pass

    return [x for x in tags if x is not None]


def colormap(exps, colorby, definedinEM, annotation=None):
    """Generate the self.colors in the format which compatible with matplotlib"""
    if definedinEM:
        if colorby == "reads":
            color_res = []
            for i in exps.get_readsnames():
                c = exps.get_type(i, "color")
                if c[0] == "(":
                    rgb = [eval(j) for j in c.strip('()').split(',')]
                    color_res.append([v / 255 for v in rgb])
                else:
                    color_res.append(c)
        elif colorby == "regions":
            color_res = []
            for i in exps.get_regionsnames():
                c = exps.get_type(i, "color")
                if c[0] == "(":
                    rgb = [eval(j) for j in c.strip('()').split(',')]
                    color_res.append([v / 255 for v in rgb])
                else:
                    color_res.append(c)
        else:
            color_res = []
            for i in list(exps.fieldsDict[colorby].values()):
                c = exps.get_type(i[0], "color")
                if c[0] == "(":
                    rgb = [float(j) for j in c.strip('()').split(',')]
                    color_res.append([v / 255 for v in rgb])
                else:
                    color_res.append(c)

    else:
        if annotation:
            color_res = plt.cm.Set1(numpy.linspace(0, 1, len(annotation))).tolist()
        else:
            # colors = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan', 'orange']
            # colors = plt.cm.jet(numpy.linspace(0.1, 0.9, len(gen_tags(exps, colorby)))).tolist()
            if colorby == "reads":
                ks = []
                for gr in exps.get_readsnames():
                    ks.append(exps.get_type(name=gr, field="factor"))
                n = len(set(ks))
            elif colorby == "regions":
                ks = []
                for gr in exps.get_regionsnames():
                    ks.append(exps.get_type(name=gr, field="factor"))
                n = len(set(ks))
            else:
                n = len(list(exps.fieldsDict[colorby].keys()))
            # print(n)

            cmap = plt.cm.get_cmap('jet')
            print(cmap)
            sys.exit()
            color_res = cmap(numpy.arange(n))

            # if n < 8:
            #     indn = np.linspace(0, 32, 256)
            #     color_res = [plt.cm.Set1(indn[i]) for i in range(n)]
            # else:
            #     set1 = plt.get_cmap('Set1')
            #     cNorm = colormat.Normalize(vmin=0, vmax=n)
            #     scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=set1)
            #     color_res = [scalarMap.to_rgba(d) for d in range(n)]

            # color_res = plt.cm.Set1(numpy.linspace(0.1, 0.9, n)).tolist()
            # print(len(plt.cm.Set1().tolist()))
            #
            # np.linspace(0, 1, 9)

    color_res = unique(color_res)
    return color_res


def colormaps(exps, colorby, definedinEM):
    """Generate a list of colormaps in the format which compatible with matplotlib"""
    if definedinEM:
        if colorby == "reads":
            colors = []
            for i in exps.get_readsnames():
                c = exps.get_type(i, "color")
                colors.append(exps.get_type(i, "color"))
        elif colorby == "regions":
            colors = []
            for i in exps.get_regionsnames():
                c = exps.get_type(i, "color")
                colors.append(exps.get_type(i, "color"))
        else:
            colors = [exps.get_type(i, "color") for i in exps.fieldsDict[colorby]]
    else:
        if colorby == "reads" or colorby == "regions":
            n = len(list(exps.fieldsDict["factor"].keys()))
        else:
            n = len(list(exps.fieldsDict[colorby].keys()))
        colors = plt.cm.Set1(numpy.linspace(0, n - 1, n)).tolist()
        # if len(exps.get_regionsnames()) < 20:
        #    colors = ['Blues', 'Oranges', 'Greens', 'Reds',  'Purples', 'Greys', 'YlGnBu', 'gist_yarg', 'GnBu',
        #              'OrRd', 'PuBu', 'PuRd', 'RdPu', 'YlGn', 'BuGn', 'YlOrBr', 'BuPu','YlOrRd','PuBuGn','binary']
        # else:
        # colors = plt.cm.Set2(numpy.linspace(0.1, 0.9, len(exps.get_regionsnames()))).tolist()
        #    colors = plt.cm.Spectral(numpy.linspace(0.1, 0.9, len(exps.get_regionsnames()))).tolist()
    return colors


def color_groupded_region(EM, grouped_region, colorby, definedinEM):
    """Generate the self.colors in the format which compatible with matplotlib"""
    if definedinEM:
        colors = OrderedDict()
        for ty in list(grouped_region.keys()):
            for q in list(grouped_region[ty].keys()):
                c = EM.get_type(q.name, "color")
                if c[0] == "(":
                    rgb = [eval(j) for j in c.strip('()').split(',')]
                    colors[q.name] = rgb
                else:
                    colors[q.name] = c
    else:
        colors = OrderedDict()
        qs = []
        if colorby == "regions":
            for ty in list(grouped_region.keys()):
                for q in grouped_region[ty]:
                    qs.append(q.name)
            qs = list(set(qs))
            # Accent Spectral hsv Set1
            colormap = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(qs))).tolist()

            for i, q in enumerate(qs):
                colors[q] = colormap[i]
        else:
            types = list(EM.fieldsDict[colorby].keys())

            colormap = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(types))).tolist()

            for ty in list(grouped_region.keys()):
                for q in grouped_region[ty]:
                    i = types.index(EM.get_type(q.name, colorby))
                    colors[q.name] = colormap[i]
    return colors


def output_array(array, directory, folder, filename):
    """ Write a txt file from the given array. """
    pd = os.path.join(current_dir, directory, folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)

    f = open(os.path.join(pd, filename), "w")
    for i, line in enumerate(array):
        f.write(("\t".join(str(j) for j in line)) + "\n")
    f.close()


def remove_duplicates(grouped_dict):
    for ty in grouped_dict:
        for r in grouped_dict[ty]:
            r.remove_duplicates()


def group_refque(rEM, qEM, groupby, rRegion=None, qRegion=None):
    """ Group regionsets of rEM and qEM according to groupby """
    groupedreference = OrderedDict()  # Store all bed names according to their types
    groupedquery = OrderedDict()  # Store all bed names according to their types
    if rRegion:
        rregs = rRegion
    else:
        rregs = rEM.get_regionsets()
    if qRegion:
        qregs = qRegion
    else:
        qregs = qEM.get_regionsets()

    if groupby:
        for r in rregs:
            ty = rEM.get_type(r.name, groupby)
            try:
                groupedreference[ty].append(r)
            except:
                groupedreference[ty] = [r]

        for q in qregs:
            ty = qEM.get_type(q.name, groupby)
            try:
                groupedquery[ty].append(q)
            except:
                groupedquery[ty] = [q]
    else:
        groupedreference[""] = rregs
        groupedquery[""] = qregs
    return groupedreference, groupedquery


def count_intersect(reference, query, mode_count="count", threshold=False):
    bed1 = copy.deepcopy(reference)
    bed2 = copy.deepcopy(query)
    if mode_count == "count":
        if threshold:
            if bed1.total_coverage() == 0:
                print(
                    "\n ** Warning : " + bed1.name + " has no length (only points) for finding intersection with given threshold.")
                sys.exit(1)
            if bed2.total_coverage() == 0:
                print(
                    "\n ** Warning : " + bed2.name + " has no length (only points) for finding intersection with given threshold.")
                sys.exit(1)
            if 50 >= threshold > 0:
                bed1.extend(-threshold, -threshold, percentage=True)
            elif threshold > 50 or threshold < 0:
                print("\n **** Threshold should be the percentage between 0 and 50. ****\n")
                sys.exit(1)
        ##if bed1.total_coverage() == 0: bed1.extend(0,1)
        # if bed2.total_coverage() == 0: bed2.extend(0,1)
        intersect_r = bed1.intersect(bed2, mode=OverlapType.OVERLAP)
        c_inter = len(intersect_r)
        intersect_1 = bed1.intersect(intersect_r, mode=OverlapType.ORIGINAL)
        c_12 = len(bed1) - len(intersect_1)
        intersect_2 = bed2.intersect(intersect_r, mode=OverlapType.ORIGINAL)
        c_21 = len(bed2) - len(intersect_2)
        # print(c_12, c_21, c_inter)
        return c_12, c_21, c_inter

    elif mode_count == "bp":
        intersect_r = bed1.intersect(bed2, mode=OverlapType.OVERLAP)
        len_inter = intersect_r.total_coverage()
        allbed1 = bed1.total_coverage()
        allbed2 = bed2.total_coverage()
        len_12 = allbed1 - len_inter
        len_21 = allbed2 - len_inter
        return len_12, len_21, len_inter


def value2str(value):
    if isinstance(value, str): return value
    if value == 0: return "0"
    if isinstance(value, int):
        return str(value)
    elif isinstance(value, float):
        if value >= 1000:
            r = "{}".format(int(value))
        elif 1000 > value > 10:
            r = "{:.1f}".format(value)
        elif 10 > value >= 1:
            r = "{:.2f}".format(value)
        elif 0.9999 > value > 0.0001:
            r = "{:.4f}".format(value)
        elif 1 > value > 0.9999:
            r = "1"
        else:
            r = "{:.1e}".format(value)
        return r


def multiple_correction(dic):
    """
    dic[ty][r][q] = p
    """
    for ty in list(dic.keys()):
        all_p = []
        rn = len(list(dic[ty].keys()))
        qn = len(list(dic[ty].values())[0].keys())
        cue = {}
        i = 0
        if rn == 1 and qn == 1: return
        # get all p values from the dictionary
        for r in list(dic[ty].keys()):
            for q in list(dic[ty][r].keys()):

                if isinstance(dic[ty][r][q], str):
                    pass
                else:
                    all_p.append(dic[ty][r][q])
                    cue[ty + r + q] = i
                    i = i + 1
        # correction
        reject, pvals_corrected = multiple_test_correction(all_p, alpha=0.05, method='indep')
        # modify all p values
        for ir, r in enumerate(dic[ty].keys()):
            for iq, q in enumerate(dic[ty][r].keys()):
                try:
                    dic[ty][r][q] = pvals_corrected[cue[ty + r + q]]
                except:
                    pass


def compute_coverage(inputs):
    """
    [1]  bedname, bamname, processed_beds, processed_bedsF,
    [5]  g, r, c, cc, d, read_file,
    [11] rs, bs, ss, center, heatmap, logt,
    [17] sense, strand, flipnegative, center, outside, extend]
    """
    ts = time.time()
    normRPM = True

    [ bedname, bamname, processed_beds, processed_bedsF, g, r, c, cc, d, read_file,
      rs, bs, ss, center_end, heatmap, logt, sense, strand, flipnegative, outside, extend ] = inputs

    res = defaultdict(list)

    if len(processed_beds) == 0:
        res["all"].append(numpy.zeros(5))
        return [g, r, c, cc, d, res]
    else:
        cov = CoverageSet(bedname + "." + bamname, processed_beds)

        if "Conservation" in [g, r, c, cc, d]:
            cov.phastCons46way_score(stepsize=ss)
        elif ".bigwig" in read_file.lower() or ".bw" in read_file.lower():
            cov.coverage_from_bigwig(bigwig_file=read_file, stepsize=ss)
        else:
            if not sense and not strand:
                cov.coverage_from_bam(bam_file=read_file, extension_size=rs, binsize=bs, stepsize=ss)
                if normRPM: cov.normRPM()
            else:  # Sense specific
                cov.coverage_from_bam(bam_file=read_file, extension_size=rs, binsize=bs, stepsize=ss,
                                      get_sense_info=sense, get_strand_info=strand, paired_reads=True)
                cov.array_transpose()
                if normRPM: cov.normRPM()

        if center_end == "midpoint" and flipnegative:
            for k, re in enumerate(processed_beds):
                if re.orientation == "-":
                    cov.coverage[k] = cov.coverage[k][::-1]

        # When bothends, consider the fliping end
        # if center_end == 'bothends' or center_end == 'upstream' or center_end == 'downstream':
        if center_end == 'bothends':
            if "Conservation" in [g, r, c, cc, d]:
                flap = CoverageSet("for flap", processed_bedsF)
                flap.phastCons46way_score(stepsize=ss)
                ffcoverage = numpy.fliplr(flap.coverage)
                cov.coverage = numpy.concatenate((cov.coverage, ffcoverage),
                                                 axis=0)
            elif ".bigwig" in read_file.lower() or ".bw" in read_file.lower():
                flap = CoverageSet("for flap", processed_bedsF)
                flap.coverage_from_bigwig(bigwig_file=read_file, stepsize=ss)
                ffcoverage = numpy.fliplr(flap.coverage)
                cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
            else:
                flap = CoverageSet("for flap", processed_bedsF)
                if not sense:
                    flap.coverage_from_bam(read_file, extension_size=rs, binsize=bs, stepsize=ss)
                    if normRPM: flap.normRPM()
                else:  # Sense specific
                    flap.coverage_from_bam(bam_file=read_file, extension_size=rs, binsize=bs,
                                           stepsize=ss, get_sense_info=True, paired_reads=True)
                    flap.array_transpose(flip=True)
                    if normRPM: flap.normRPM()
                ffcoverage = numpy.fliplr(flap.coverage)
                try:
                    cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
                except:
                    pass

                if sense:
                    cov.transpose_cov1 = numpy.concatenate((cov.transpose_cov1, flap.transpose_cov1), axis=0)
                    cov.transpose_cov2 = numpy.concatenate((cov.transpose_cov2, flap.transpose_cov2), axis=0)
        # Extend outside
        if outside:
            new_arrays = []
            for ar in cov.coverage:
                ss_side = int(extend / ss) - 1
                # if len(ar) < 2*ss_side
                left_ar = ar[0:ss_side]
                right_ar = ar[-ss_side:]
                rest = ar[ss_side:-ss_side]

                # print([len[ar], len(left_ar), len(rest), len(right_ar)])
                try:
                    xp = numpy.linspace(0, ss_side, len(rest))
                    rest = numpy.interp(list(range(ss_side)), xp=xp, fp=rest)
                    nar = numpy.concatenate((left_ar, rest))
                    nar = numpy.concatenate((nar, right_ar))
                    new_arrays.append(nar)
                except:
                    print([ss_side, extend, ss, len(ar), rest])
            cov.coverage = new_arrays

        # Averaging the coverage of all regions of each bed file
        if heatmap:
            if logt:
                res = numpy.log10(numpy.vstack(cov.coverage) + 1)  # Store the array into data list
            else:
                res = numpy.vstack(cov.coverage)  # Store the array into data list
        else:
            if len(cov.coverage) == 0:
                res = None
                print("** Warning: Cannot open " + read_file)
            else:
                for i, car in enumerate(cov.coverage):

                    if i == 0:
                        avearr = numpy.array(car, ndmin=2)
                    else:

                        if avearr.shape[1] == len(car):
                            avearr = numpy.vstack((avearr, numpy.array(car, ndmin=2)))
                if logt:
                    avearr = numpy.log10(avearr + 1)

                avearr = numpy.average(avearr, axis=0)
                if sense or strand:
                    if logt:
                        sense_1 = numpy.average(numpy.log2(cov.transpose_cov1 + 1), axis=0)
                        sense_2 = numpy.average(numpy.log2(cov.transpose_cov2 + 1), axis=0)
                    else:
                        sense_1 = numpy.average(cov.transpose_cov1, axis=0)
                        sense_2 = numpy.average(cov.transpose_cov2, axis=0)
                cut_end = int(bs / ss)
                avearr = avearr[cut_end:-cut_end]
                res["all"].append(avearr)

                if sense or strand:
                    sense_1 = sense_1[cut_end:-cut_end]
                    sense_2 = sense_2[cut_end:-cut_end]
                    res["sense_1"].append(sense_1)
                    res["sense_2"].append(sense_2)

        result = [g, r, c, cc, d, res]  # Store the array into data list
        te = time.time()
        print("\tComputing " + bedname + " . " + bamname + "\t\t" +
              str(datetime.timedelta(seconds=round(te - ts))))
        return result


def mp_count_intersets(inps):
    # [ com, self.rlen[ty][r.name], self.mode_count, threshold ]
    random_r, random_q = inps[0].random_split(size=inps[1])
    d = random_r.intersect_count(random_q, mode_count=inps[2], threshold=inps[3])
    return d


def mp_count_intersect(inputs):
    # q, nalist, mode_count, qlen_dict, threshold, counts, frequency, self_frequency, ty, r
    q = inputs[0]
    nalist = inputs[1]
    mode_count = inputs[2]
    # qlen_dict = inputs[3]
    threshold = inputs[4]
    # counts = inputs[5]
    frequency = inputs[6]
    # self_frequency = inputs[7]
    ty = inputs[8]
    r = inputs[9]

    output = [q.name]

    if q.total_coverage() == 0 and len(q) > 0:
        output.append(q.name)

    else:
        output.append(None)

        if mode_count == "bp":
            qlen = q.total_coverage()
        elif mode_count == "count":
            qlen = len(q)

        output.append(qlen)
        # Define different mode of intersection and count here
        print(".", end="")
        sys.stdout.flush()
        c = r.intersect_count(q, mode_count=mode_count, threshold=threshold)

        output.append(c)
        # counts[ty][r.name][q.name] = c
        if frequency:
            output.append(c[2])
            # try: self_frequency[ty][q.name].append(c[2])
            # except:
            #    self_frequency[ty][q.name] = [c[2]]

    return output


def get_url(url, filename):
    if sys.version_info >= (3, 1):
        u = urlopen(url)
        print("** Loading " + url, end="")
        f = open(filename, 'wb')
        print(u.read(), file=f)
        f.close()
        print("\t... Done")
    else:
        u = urllib.request.urlopen(url)
        print("** Loading " + url, end="")
        f = open(filename, 'wb')
        print(u.read(), file=f)
        f.close()
        print("\t... Done")


def load_dump(path, filename):
    print("\tLoading from file: " + filename)
    f = open(os.path.join(path, filename), 'r')
    ob = pickle.load(f)
    f.close()
    return ob


def dump(object, path, filename):
    print("\tDump to file: " + filename)
    f = open(os.path.join(path, filename), 'wb')
    pickle.dump(object, f)
    f.close()


def read_gtf(gd, organism):
    try:
        anno = load_dump(gd.get_annotation_dump_dir(), "gtf.dump")
        print("ee")
    except:
        anno = AnnotationSet(organism)
        dump(anno, gd.get_annotation_dump_dir(), "gtf.dump")
        print("qq")
    return anno


def annotation_dump(organism):
    print("\nLoading genetic annotation data...\n")
    gd = GenomeData(organism)
    beds = []
    # TSS
    try:
        tss = load_dump(gd.get_annotation_dump_dir(), "tss.dump")
    except:
        anno = read_gtf(gd, organism)
        tss = anno.get_tss()
        dump(tss, gd.get_annotation_dump_dir(), "tss.dump")
    beds.append(tss)
    print("\tTSS: " + str(len(tss)))
    # TTS
    try:
        tts = load_dump(gd.get_annotation_dump_dir(), "tts.dump")
    except:
        anno = read_gtf(gd, organism)
        tts = anno.get_tts()
        dump(tts, gd.get_annotation_dump_dir(), "tts.dump")
    beds.append(tts)
    print("\tTTS: " + str(len(tts)))

    # exon
    try:
        exon = load_dump(gd.get_annotation_dump_dir(), "exon.dump")
    except:
        anno = read_gtf(gd, organism)
        exon = anno.get_exons(start_site=False, end_site=False)
        dump(exon, gd.get_annotation_dump_dir(), "exon.dump")
    print("\texon: " + str(len(exon)))
    # exons
    try:
        exons = load_dump(gd.get_annotation_dump_dir(), "exons.dump")
    except:
        anno = read_gtf(gd, organism)
        exons = anno.get_exons(start_site=True, end_site=False)
        dump(exons, gd.get_annotation_dump_dir(), "exons.dump")
    beds.append(exons)
    print("\texons: " + str(len(exons)))
    # exone
    try:
        exone = load_dump(gd.get_annotation_dump_dir(), "exone.dump")
    except:
        anno = read_gtf(gd, organism)
        exone = anno.get_exons(start_site=False, end_site=True)
        dump(exone, gd.get_annotation_dump_dir(), "exone.dump")
    beds.append(exone)
    print("\texone: " + str(len(exone)))

    # bednames = ["TSS", "TTS", "Exon start site", "Exon end site", "Intron start site", "Intron end site"]
    bednames = ["TSS", "TTS", "Exon start site", "Exon end site"]
    annotation = bednames
    return beds, bednames, annotation

def output(f, directory, folder, filename, extra=None, pdf=False, show=None):
    """Output the file in the defined folder """
    pd = os.path.normpath(os.path.join(current_dir, directory, folder))
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)

        # Saving
    if not extra:
        f.savefig(os.path.join(pd, filename), facecolor='w', edgecolor='w',
                  bbox_inches='tight', dpi=400)
    else:
        f.savefig(os.path.join(pd, filename), facecolor='w', edgecolor='w',
                  bbox_extra_artists=extra, bbox_inches='tight', dpi=400)

    if pdf:
        try:
            pp = PdfPages(os.path.join(pd, filename) + '.pdf')
            pp.savefig(f, bbox_extra_artists=extra, bbox_inches='tight')
            pp.close()
        except:
            print("ERROR: Problem in PDF conversion. Skipped.")
    if show:
        plt.show()


def output_parameters(parameter, directory, folder, filename):
    pd = os.path.join(current_dir, directory, folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)
    if parameter:
        with open(os.path.join(pd, "parameters.txt"), 'w') as f:
            for s in parameter:
                print(s, file=f)


def copy_em(em, directory, folder, filename="experimental_matrix.txt"):
    copyfile(em, os.path.join(current_dir, directory, folder, filename))


def list_all_index(path):
    """Creat an 'index.html' in the defined directory """
    dirname = os.path.basename(path)
    parentdir = os.path.basename(os.path.dirname(path))

    # link_d = {"List":"index.html"}
    link_d = {}
    ####
    for root, dirnames, filenames in os.walk(os.path.dirname(path)):
        for filename in fnmatch.filter(filenames, 'index.html'):
            if root.split('/')[-2] == parentdir:
                link_d[root.split('/')[-1]] = "../" + root.split('/')[-1] + "/index.html"
    link_d = OrderedDict(sorted(list(link_d.items()), key=lambda key_value: key_value[0]))

    ###

    html = Html(name="Directory: " + dirname, links_dict=link_d,
                fig_dir=os.path.join(path, "style"), fig_rpath="./style", RGT_header=False, other_logo="viz")
    header_list = ["No.", "Experiments"]
    html.add_heading("All experiments in: " + dirname + "/")
    data_table = []
    type_list = 'ssss'
    col_size_list = [10, 10, 10]
    c = 0
    for root, dirnames, filenames in os.walk(path):
        # roots = root.split('/')
        for filename in fnmatch.filter(filenames, '*.html'):
            if filename == 'index.html' and root.split('/')[-1] != dirname:
                # print(root)
                c += 1
                data_table.append([str(c),
                                   '<a href="' + os.path.join(root.split('/')[-1], filename) + '"><font size="4">' +
                                   root.split('/')[-1] + "</a>"])
                # print(link_d[roots[-1]])
    html.add_zebra_table(header_list, col_size_list, type_list, data_table, align=50, cell_align="left", sortable=True)
    html.add_fixed_rank_sortable()
    html.write(os.path.join(path, "index.html"))


def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try:
        os.stat(path)
    except:
        os.mkdir(path)


def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero
    
    http://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(multiprocessing.pool.Pool):
    def __reduce__(self):
        pass

    Process = NoDaemonProcess


def walklevel(some_dir, level=1):
    some_dir = some_dir.rstrip(os.path.sep)
    assert os.path.isdir(some_dir)
    num_sep = some_dir.count(os.path.sep)
    for root, dirs, files in os.walk(some_dir):
        yield root, dirs, files
        num_sep_this = root.count(os.path.sep)
        if num_sep + level <= num_sep_this:
            del dirs[:]
