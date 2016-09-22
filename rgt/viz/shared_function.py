# Python Libraries
from __future__ import print_function
from __future__ import division
import time, datetime
import multiprocessing.pool
import urllib2
import re
import matplotlib.pyplot as plt

# Local Libraries
# Distal Libraries
from rgt.ExperimentalMatrix import *
from rgt.CoverageSet import *
from rgt.motifanalysis.Statistics import multiple_test_correction

# Local test
dir = os.getcwd()

###########################################################################################
#                    Universal functions
###########################################################################################


def print2(parameter, string):
    """ Show the message on the console and also save in a list for future backup. """
    print(string)
    parameter.append(string)


def unique(a):
    seen = set()
    return [seen.add(x) or x for x in a if x not in seen]


def purge(dir, pattern):
    for f in os.listdir(dir):
        if re.search(pattern, f):
            os.remove(os.path.join(dir, f))


def gen_tags(exps, tag):
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
        try:
            l = [exps.get_type(i, "factor") for i in exps.get_regionsnames()]
        except:
            # print("You must define 'factor' column in experimental matrix for grouping.")
            sys.exit(1)
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
            colors = []
            for i in exps.get_readsnames():
                c = exps.get_type(i, "color")
                if c[0] == "(":
                    rgb = [eval(j) for j in c.strip('()').split(',')]
                    colors.append(rgb)
                else:
                    colors.append(c)
        elif colorby == "regions":
            colors = []
            for i in exps.get_regionsnames():
                c = exps.get_type(i, "color")
                if c[0] == "(":
                    rgb = [eval(j) for j in c.strip('()').split(',')]
                    colors.append([v / 255 for v in rgb])
                else:
                    colors.append(c)
        else:
            colors = []
            for i in exps.fieldsDict[colorby].values():
                c = exps.get_type(i[0], "color")
                if c[0] == "(":
                    rgb = [float(j) for j in c.strip('()').split(',')]
                    colors.append([v / 255 for v in rgb])
                else:
                    colors.append(c)

    else:
        if annotation:
            colors = plt.cm.Set1(numpy.linspace(0, 1, len(annotation))).tolist()
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
                n = len(exps.fieldsDict[colorby].keys())
            # print(n)
            # colors = plt.cm.Spectral(numpy.linspace(0.1, 0.9, n)).tolist()
            colors = plt.cm.Set1(numpy.linspace(0, 1, n)).tolist()
    return colors


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
            n = len(exps.fieldsDict["factor"].keys())
        else:
            n = len(exps.fieldsDict[colorby].keys())
        colors = plt.cm.Set1(numpy.linspace(0, 1, n)).tolist()


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
        for ty in grouped_region.keys():
            for q in grouped_region[ty].keys():
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
            for ty in grouped_region.keys():
                for q in grouped_region[ty]:
                    qs.append(q.name)
            qs = list(set(qs))
            # Accent Spectral hsv Set1
            colormap = plt.cm.Spectral(numpy.linspace(0.1, 0.9, len(qs))).tolist()

            for i, q in enumerate(qs):
                colors[q] = colormap[i]
        else:
            types = EM.fieldsDict[colorby].keys()

            colormap = plt.cm.Spectral(numpy.linspace(0.1, 0.9, len(types))).tolist()

            for ty in grouped_region.keys():
                for q in grouped_region[ty]:
                    i = types.index(EM.get_type(q.name, colorby))
                    colors[q.name] = colormap[i]
    return colors


def output_array(array, directory, folder, filename):
    """ Write a txt file from the given array. """
    pd = os.path.join(dir, directory, folder)
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
    if (isinstance(value, str)): return value
    if value == 0: return "0"
    if (isinstance(value, int)):
        return str(value)
    elif (isinstance(value, float)):
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
    for ty in dic.keys():
        all_p = []
        rn = len(dic[ty].keys())
        qn = len(dic[ty].values()[0].keys())
        cue = {}
        i = 0
        if rn == 1 and qn == 1: return
        # get all p values from the dictionary
        for r in dic[ty].keys():
            for q in dic[ty][r].keys():

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


def compute_coverage(input):
    """
    bed, bam, rs, bs, ss, center, heatmap, logt, s, g, c, d
    """
    ts = time.time()
    cov = CoverageSet(input[0].name + ".", input[0])
    if ".bigWig" in input[1] or ".bw" in input[1]:
        cov.coverage_from_bigwig(bigwig_file=input[1], stepsize=input[4])
    else:
        cov.coverage_from_bam(bam_file=input[1], extension_size=input[2], binsize=input[3], stepsize=input[4])
        cov.normRPM()
    # When bothends, consider the fliping end
    if input[5] == 'bothends':
        flap = CoverageSet("for flap", input[0])
        flap.coverage_from_bam(input[1], extension_size=input[2], binsize=input[3], stepsize=input[4])
        ffcoverage = numpy.fliplr(flap.coverage)
        cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
    # Averaging the coverage of all regions of each bed file
    if input[6]:
        if input[7]:
            result = numpy.log10(numpy.vstack(cov.coverage))  # Store the array into data list
        else:
            result = numpy.vstack(cov.coverage)  # Store the array into data list
    else:
        # print(cov.coverage)
        for i, car in enumerate(cov.coverage):
            car = numpy.delete(car, [0, 1])
            if i == 0:
                avearr = np.array(car)
                lenr = car.shape[0]
            elif car.shape[0] == lenr:
                avearr = numpy.vstack((avearr, car))
            else:
                pass
        # avearr = numpy.array(cov.coverage)
        # print(avearr)
        # print(avearr.shape)
        avearr = numpy.average(avearr, axis=0)
        # numpy.transpose(avearr)
        result = [input[8], input[9], input[10], input[11], avearr]  # Store the array into data list
    te = time.time()
    print("\tComputing " + os.path.basename(input[1]) + " . " + input[0].name + "\t\t" + str(
        datetime.timedelta(seconds=round(te - ts))))
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

    # qname, nalist, qlen_dict[ty][q.name], counts[ty][r.name][q.name], self_frequency[ty][q.name].append(c[2])
    return output


def get_url(url, filename):
    # file_name = url.split('/')[-1]
    u = urllib2.urlopen(url)
    print("** Loading " + url, end="")
    f = open(filename, 'wb')
    print(u.read(), file=f)
    f.close()
    print("\t... Done")


def load_dump(path, filename):
    print("\tLoading from file: " + filename)
    file = open(os.path.join(path, filename), 'r')
    object = pickle.load(file)
    file.close()
    return object


def dump(object, path, filename):
    print("\tDump to file: " + filename)
    file = open(os.path.join(path, filename), 'wb')
    pickle.dump(object, file)
    file.close()


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
    Process = NoDaemonProcess