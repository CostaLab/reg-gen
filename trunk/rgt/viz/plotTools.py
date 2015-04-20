# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os
import numpy
from scipy.stats import mstats, wilcoxon, mannwhitneyu, rankdata
import time, datetime, argparse
from collections import *
import copy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FuncFormatter 
from matplotlib import cm
import itertools
import pickle
import multiprocessing

# Local Libraries
# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.ExperimentalMatrix import *
from rgt.AnnotationSet import *
from rgt.Util import GenomeData, OverlapType, Html
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
        
def gen_tags(exps, tag):
    """Generate the unique tags from the EM according to the given tag. """
    if tag == "reads":
        try: l = [exps.get_type(i,"factor") for i in exps.get_readsnames()]
        except: 
            print("You must define 'factor' column in experimental matrix for grouping.")
            sys.exit(1)
    elif tag == "regions":
        try: l = [exps.get_type(i,"factor") for i in exps.get_regionsnames()]
        except: 
            print("You must define 'factor' column in experimental matrix for grouping.")
            sys.exit(1)
    else:
        try: l = exps.fieldsDict[tag]
        except: 
            print('Cannot find the column "' + tag +'"')
            sys.exit(1)
    return unique(l)

def tag_from_r(exps, tag_type, name):
    tags = []
    for type in tag_type:
        if type == "reads" or type == "regions": type = "factor" 
        try:
            tags.append(exps.get_type(name,type))
        except: pass
    return tags
        
def colormap(exps, colorby, definedinEM, annotation=None):
    """Generate the self.colors in the format which compatible with matplotlib"""
    if definedinEM:
        if colorby == "reads":
            colors = []
            for i in exps.get_readsnames():
                c = exps.get_type(i,"color")
                if c[0] == "(":
                    rgb = [ eval(j) for j in c.strip('()').split(',')]
                    colors.append(rgb)
                else:
                    colors.append(c)
        elif colorby == "regions":
            colors = []
            for i in exps.get_regionsnames():
                c = exps.get_type(i,"color")
                if c[0] == "(":
                    rgb = [ eval(j) for j in c.strip('()').split(',')]
                    colors.append([v/255 for v in rgb])
                else:
                    colors.append(c)
        else:
            colors = []
            for i in exps.fieldsDict[colorby].values():
                c = exps.get_type(i[0],"color")
                if c[0] == "(":
                    rgb = [ float(j) for j in c.strip('()').split(',')]
                    colors.append([v/255 for v in rgb])
                else:
                    colors.append(c)
            
    else:
        if annotation:
            colors = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(annotation))).tolist()
        else:
            #colors = [ 'lightgreen', 'pink', 'cyan', 'lightblue', 'tan', 'orange']
            #colors = plt.cm.jet(numpy.linspace(0.1, 0.9, len(gen_tags(exps, colorby)))).tolist()
            colors = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(gen_tags(exps, colorby)))).tolist()
    return colors

def colormaps(exps, colorby, definedinEM):
    """Generate a list of colormaps in the format which compatible with matplotlib"""
    if definedinEM:
        if colorby == "reads":
            colors = []
            for i in exps.get_readsnames():
                c = exps.get_type(i,"color")
                colors.append(exps.get_type(i,"color"))
        elif colorby == "regions":
            colors = []
            for i in exps.get_regionsnames():
                c = exps.get_type(i,"color")
                colors.append(exps.get_type(i,"color"))
        else:
            colors = [exps.get_type(i,"color") for i in exps.fieldsDict[colorby]]
    else:
        
        if len(exps.get_regionsnames()) < 20:
            colors = ['Blues', 'Oranges', 'Greens', 'Reds',  'Purples', 'Greys', 'YlGnBu', 'gist_yarg', 'GnBu', 
                      'OrRd', 'PuBu', 'PuRd', 'RdPu', 'YlGn', 'BuGn', 'YlOrBr', 'BuPu','YlOrRd','PuBuGn','binary']
        else:
            colors = plt.cm.Set2(numpy.linspace(0.1, 0.9, len(exps.get_regionsnames()))).tolist()
    return colors

def color_groupded_region(EM, grouped_region, colorby, definedinEM):
    """Generate the self.colors in the format which compatible with matplotlib"""
    if definedinEM:
        colors = OrderedDict()
        for ty in grouped_region.keys():
            for q in grouped_region[ty].keys():
                c = EM.get_type(q.name,"color")
                if c[0] == "(":
                    rgb = [ eval(j) for j in c.strip('()').split(',')]
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
            colormap = plt.cm.Spectral(numpy.linspace(0, 1, len(qs))).tolist()
            for i, q in enumerate(qs):
                colors[q] = colormap[i]
        else:
            types = EM.fieldsDict[colorby].keys()
            colormap = plt.cm.Spectral(numpy.linspace(0, 1, len(types))).tolist()
            for ty in grouped_region.keys():
                for q in grouped_region[ty]: 
                    i = types.index(EM.get_type(q.name, colorby))
                    colors[q.name] = colormap[i]
    return colors 

def output_array(array, directory, folder, filename):
    """ Write a txt file from the given array. """
    pd = os.path.join(dir,directory,folder)
    try:
        os.stat(os.path.dirname(pd))
    except:
        os.mkdir(os.path.dirname(pd))
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)   
             
    f = open(os.path.join(pd,filename),"w")
    for i,line in enumerate(array):
        f.write(("\t".join(str(j) for j in line))+"\n")
    f.close()

def remove_duplicates(grouped_dict):
    for ty in grouped_dict:
        for r in grouped_dict[ty]:
            r.remove_duplicates()

def group_refque(rEM, qEM, groupby):
    """ Group regionsets of rEM and qEM according to groupby """
    groupedreference = OrderedDict()  # Store all bed names according to their types
    groupedquery = OrderedDict()  # Store all bed names according to their types
    if groupby:
        for r in rEM.get_regionsets():
            ty = rEM.get_type(r.name,groupby)
            try: groupedreference[ty].append(r)
            except: groupedreference[ty] =[r]
        
        for q in qEM.get_regionsets():
            ty = qEM.get_type(q.name,groupby)
            try: groupedquery[ty].append(q)
            except: groupedquery[ty] =[q]
    else:
        groupedreference["All region sets without grouping"] = rEM.get_regionsets()
        groupedquery["All region sets without grouping"] = qEM.get_regionsets()
    return groupedreference, groupedquery

def count_intersect(reference, query, mode_count="count", threshold=False):
    bed1 = copy.deepcopy(reference)
    bed2 = copy.deepcopy(query)
    if mode_count=="count":
        if threshold:
            if bed1.total_coverage() == 0:
                print("\n ** Warning : "+ bed1.name +" has no length (only points) for finding intersection with given threshold.")
                sys.exit(1)
            if bed2.total_coverage() == 0:
                print("\n ** Warning : "+ bed2.name +" has no length (only points) for finding intersection with given threshold.")
                sys.exit(1)
            if 50 >= threshold > 0:
                bed1.extend(-threshold,-threshold,percentage=True)
            elif threshold > 50 or threshold < 0:
                print("\n **** Threshold should be the percentage between 0 and 50. ****\n")
                sys.exit(1)
        ##if bed1.total_coverage() == 0: bed1.extend(0,1)
        #if bed2.total_coverage() == 0: bed2.extend(0,1)
        intersect_r = bed1.intersect(bed2, mode=OverlapType.ORIGINAL)
        #intersect_r.remove_duplicates()
        c_inter = len(intersect_r)
        c_12 = len(bed1) - c_inter
        c_21 = len(bed2) - c_inter
        #print(c_12, c_21, c_inter)
        return c_12, c_21, c_inter
        
    elif mode_count=="bp":
        intersect_r = bed1.intersect(bed2, mode=OverlapType.OVERLAP)
        len_inter = intersect_r.total_coverage()
        allbed1 = bed1.total_coverage()
        allbed2 = bed2.total_coverage()
        len_12 = allbed1 - len_inter
        len_21 = allbed2 - len_inter
        return len_12, len_21, len_inter

def count_intersect3(bedA, bedB, bedC, m="OVERLAP"):
    mode = eval("OverlapType."+m)
    #ABC
    AB = bedA.intersect(bedB, mode)
    BC = bedB.intersect(bedC, mode)
    AC = bedA.intersect(bedC, mode)
    ABC = AB.intersect(BC, mode)
    Abc = bedA.subtract(AB.combine(AC))
    aBc = bedB.subtract(AB.combine(BC))
    ABc = AB.subtract(ABC)
    abC = bedC.subtract(AC.combine(BC))
    AbC = AC.subtract(ABC)
    aBC = BC.subtract(ABC)
    return len(Abc), len(aBc), len(ABc), len(abC), len(AbC), len(aBC), len(ABC)

def value2str(value):
    if (isinstance(value,str)): return value
    if value == 0: return "0"
    if(isinstance(value,int)): return str(value)
    elif(isinstance(value,float)):
        if value >= 1000: r = "{}".format(int(value))
        elif 1000 > value > 10: r = "{:.1f}".format(value)
        elif 10 > value >= 1: r = "{:.2f}".format(value)
        elif 1 > value > 0.0001: r = "{:.4f}".format(value)
        else: r = "{:.1e}".format(value)
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
                    cue[ty+r+q] = i
                    i = i + 1
        # correction
        reject, pvals_corrected = multiple_test_correction(all_p, alpha=0.05, method='indep')
        # modify all p values
        for ir, r in enumerate(dic[ty].keys()):
            for iq, q in enumerate(dic[ty][r].keys()):
                try: dic[ty][r][q] = pvals_corrected[cue[ty+r+q]]
                except: 
                    pass

def compute_coverage(input):
    """
    bed, bam, rs, bs, ss, center, heatmap, logt, s, g, c
    """
    
    ts = time.time()
    cov = CoverageSet(input[0].name+".", input[0])
    cov.coverage_from_bam(bam_file=input[1], read_size = input[2], binsize = input[3], stepsize = input[4])
    cov.normRPM()
    # When bothends, consider the fliping end
    if input[5] == 'bothends':
        flap = CoverageSet("for flap", input[0])
        flap.coverage_from_bam(input[1], read_size = input[2], binsize = input[3], stepsize = input[4])
        ffcoverage = numpy.fliplr(flap.coverage)
        cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
    # Averaging the coverage of all regions of each bed file
    if input[6]:
        if input[7]:
            result = numpy.log10(numpy.vstack(cov.coverage)) # Store the array into data list
        else:
            result = numpy.vstack(cov.coverage) # Store the array into data list
    else:
        #print(cov.coverage)
        for i, car in enumerate(cov.coverage):
            car = numpy.delete(car, [0,1])
            if i == 0:
                avearr = np.array(car)
                lenr = car.shape[0]
            elif car.shape[0] == lenr:
                avearr = numpy.vstack((avearr, car))
            else:
                pass
        #avearr = numpy.array(cov.coverage)
        #print(avearr)
        #print(avearr.shape)
        avearr = numpy.average(avearr, axis=0)
        #numpy.transpose(avearr)
        result = [input[8], input[9], input[10], avearr] # Store the array into data list
    te = time.time()
    print("\tComputing "+os.path.basename(input[1])+" . "+input[0].name + "\t\t"+str(datetime.timedelta(seconds=round(te-ts))))
    return result    
###########################################################################################
#                    Projection test
###########################################################################################

class Projection:
    def __init__(self, reference_path, query_path):
        # Reference
        self.rEM = ExperimentalMatrix()
        self.rEM.read(reference_path)
        self.references = self.rEM.get_regionsets()
        self.referencenames = self.rEM.get_regionsnames()
        # Query
        self.qEM = ExperimentalMatrix()
        self.qEM.read(query_path)
        self.query = self.qEM.get_regionsets()
        self.querynames = self.qEM.get_regionsnames()
        self.parameter = []
        
    def group_refque(self, groupby=False):
        self.groupedreference, self.groupedquery = group_refque(self.rEM, self.qEM, groupby)
    
    def colors(self, colorby, definedinEM):
        ############# Color #####################################
        #self.color_list = colormap(self.qEM, colorby, definedinEM)
        self.color_list = color_groupded_region(self.qEM, self.groupedquery, colorby, definedinEM)
        #self.color_tags = gen_tags(self.qEM, colorby)
        #self.color_tags.append('Background')
        self.color_list['Background'] = '0.70'
    
    def ref_union(self):
        self.background = OrderedDict()
        for ty in self.groupedreference.keys():
            self.background[ty] = GenomicRegionSet("union of references")
            for r in self.groupedreference[ty]:
                self.background[ty].combine(r)
            self.background[ty].merge()
            
    def background(self, bed_path):
        bg = GenomicRegionSet("background")
        bg.read_bed(bed_path)

        self.background = OrderedDict()
        for ty in self.groupedreference.keys():
            self.background[ty] = bg
        
    def projection_test(self, organism):
        self.bglist = OrderedDict()
        self.qlist = OrderedDict()
        self.plist = OrderedDict()
        self.lenlist = {}
        #print2(self.parameter, "\nProjection test")
        #print2(self.parameter, "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}".format("Reference","Background", "Query", "Proportion", "p value"))
        
        all_p = {}
        for ty in self.groupedquery.keys():
            self.bglist[ty] = OrderedDict()
            self.qlist[ty] = OrderedDict()
            self.plist[ty] = OrderedDict()
            try:
                if self.background: bgset = self.background[ty]
            except: bgset = None
            
            for i, r in enumerate(self.groupedreference[ty]):
                self.bglist[ty][r.name] = OrderedDict()
                self.qlist[ty][r.name] = OrderedDict()
                self.plist[ty][r.name] = OrderedDict()
                self.lenlist[r.name] = len(r)
                for j, q in enumerate(self.groupedquery[ty]):
                    #print(r.name, q.name, sep="\t")
                    bg, ratio, p = r.projection_test(q, organism, extra=True, background=bgset)
                    self.bglist[ty][r.name][q.name] = bg
                    self.qlist[ty][r.name][q.name] = ratio
                    self.plist[ty][r.name][q.name] = p
                    self.lenlist[q.name] = len(q)
                    #if r in self.backgrounds.keys(): pass
                    #else: self.backgrounds[r] = bg
         
        # multiple test correction       
        multiple_correction(self.plist)
        
        for ty in self.groupedquery.keys():
            for i, r in enumerate(self.groupedreference[ty]):
                for j, q in enumerate(self.groupedquery[ty]):
                    bg = self.bglist[ty][r.name][q.name]
                    ratio = self.qlist[ty][r.name][q.name]
                    p = self.plist[ty][r.name][q.name]
                    #print(p)
                    #if len(q) == 0:
                    #    note = "Empty query!"
                    #elif p < 0.05 and bg > ratio: 
                    #    note = "Negatively unassociated!"
                    #elif p < 0.05 and bg < ratio:
                    #    note = "Positively associated!"
                    #else:
                    #    note = ""
                    #print2(self.parameter, r.name+"\t"+value2str(bg)+"\t"+q.name+"\t"+value2str(ratio)+"\t"+value2str(p)+"\t"+note)
                    
                    self.qlist[ty][r.name]['Background'] = self.bglist[ty][r.name][q.name]

    def plot(self, logt=None):
        f, ax = plt.subplots(len(self.qlist.keys()),1)
        try: ax = ax.reshape(-1)
        except: ax = [ax]
        nm = len(self.groupedreference.keys()) * len(self.groupedreference.values()[0]) * len(self.groupedquery.values()[0])
        if nm > 40:
            f.set_size_inches(nm * 0.2 +1 ,7)
            
        g_label = []
        for ind_ty, ty in enumerate(self.qlist.keys()):
            g_label.append(ty)
            r_label = []   
            for ind_r,r in enumerate(self.qlist[ty].keys()):
                r_label.append(r)
                width = 0.8/(len(self.qlist[ty][r].keys())+1) # Plus one background
                for ind_q, q in enumerate(self.qlist[ty][r].keys()):
                    x = ind_r + ind_q*width + 0.1
                    y = self.qlist[ty][r][q]
                    if y == 0 and logt: y = 0.000001
                    #print("    "+r+"     "+q+"     "+str(x)+"     "+str(y))
                    ax[ind_ty].bar(x, y, width=width, color=self.color_list[q], edgecolor="none", 
                                   align='edge', log=logt)
            if logt:
                ax[ind_ty].set_yscale('log')
            else:
                ax[ind_ty].locator_params(axis = 'y', nbins = 2)
                
            #ax[ind_ty].set_ylabel("Percentage of intersected regions",fontsize=12)
            ax[ind_ty].set_title(ty)
            ax[ind_ty].yaxis.tick_left()
            ax[ind_ty].set_xticks([i + 0.5 - 0.5*width for i in range(len(r_label))])
            ax[ind_ty].set_xticklabels(r_label,rotation=40, ha="right")
            ax[ind_ty].tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
            ax[ind_ty].legend(self.qlist[ty][r].keys(), loc='center left', handlelength=1, handletextpad=1, 
                      columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax[ind_ty].spines[spine].set_visible(False)
        f.text(-0.025, 0.5, "Percentage of intersected regions",fontsize=12, rotation="vertical", va="center")
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        self.fig = f

    def heatmap(self):
        f, ax = plt.subplots(1, len(self.plist.keys()))
        try: ax = ax.reshape(-1)
        except: ax = [ax]
        
        g_label = []
        for ind_ty, ty in enumerate(self.plist.keys()):
            g_label.append(ty)
            r_label = []
            data = []
            for ind_r,r in enumerate(self.plist[ty].keys()):
                r_label.append(r)
                #data.append(self.plist[ty][r].values())
                for ind_q, q in enumerate(self.plist[ty][r].keys()):
                    pass
            da = numpy.array(data)
            da = da.transpose()
            #im = plt.imshow(da, cmap=ax[ind_r], vmin=, vmax, origin, extent, shape, filternorm, filterrad, imlim, resample, url, hold)


    def gen_html(self, directory, title, args, align=50):
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = "Projection Test: "+dir_name
        link_d = OrderedDict()
        link_d["Projection test"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        html.add_figure("projection_test.png", align="center")
        
        header_list = ["Reference<br>name",
                       "Query<br>name", 
                       "Reference<br>number",
                       "Query<br>number", 
                       "Proportion",
                       "Background<br>proportion",
                       "Positive<br>association<br>p-value",
                       "Negative<br>association<br>p-value"]
        
        type_list = 'sssssssssss'
        col_size_list = [10,10,10,10,10,10,15,15]
        
        nalist = []
        for ind_ty, ty in enumerate(self.plist.keys()):
            html.add_heading(ty, size = 4, bold = False)
            data_table = []
            for ind_r,r in enumerate(self.plist[ty].keys()):
                rlen = str(self.lenlist[r])
                for ind_q, q in enumerate(self.plist[ty][r].keys()):
                    qlen = str(self.lenlist[q])
                    backv = value2str(self.qlist[ty][r]['Background'])
                    propor = value2str(self.qlist[ty][r][q])
                    pv = self.plist[ty][r][q]
                    if pv == "na": 
                        nalist.append(r)
                        continue
                    else:
                        pvn = 1-pv
                    
                        if self.plist[ty][r][q] < 0.05:
                            if self.qlist[ty][r]['Background'] <  self.qlist[ty][r][q]:
                                data_table.append([r,q,rlen,qlen,propor,backv,
                                                   "<font color=\"red\">"+value2str(pv)+"</font>", value2str(pvn)])
                            else:
                                data_table.append([r,q,rlen,qlen,propor,backv,
                                                   value2str(pvn), "<font color=\"red\">"+value2str(pv)+"</font>"])
                        else:
                            data_table.append([r,q,rlen,qlen,propor,backv,value2str(pv),value2str(pvn)])

            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, sortable=True)

        
        header_list=["Assumptions and hypothesis"]
        data_table = [['If the background proportion is too small, it may cause bias in p value.'],
                      ['For projection test, the reference GenomicRegionSet should have non-zero length in order to calculate its background proportion.'],
                      ['P values are corrected by multiple test correction.'],
                      ['Positive association is defined by: Proportion > Background.'],
                      ['Negative association is defined by: Proportion < Background.']]
        
        nalist = set(nalist)
        if len(nalist) > 0:
            data_table.append(['The following references contain zero-length region which cause error in proportion calculation, please check it:<br>'+
                               '     <font color=\"red\">'+', '.join([s for s in nalist])+'</font></p>'])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, cell_align="left")
        html.add_fixed_rank_sortable()
        
        html.write(os.path.join(directory,os.path.join(title,"index.html")))

        # Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        header_list = ["Description", "Argument", "Value"]
        data_table = [["Reference", "-r", args.r ],
                      ["Query", "-q", args.q],
                      ["Output directory", "-o", os.path.basename(args.o)],
                      ["Experiment title", "-t", args.t],
                      #["Grouping tag", "-g", args.g],
                      #["Coloring tag", "-c", args.c],
                      #["Background", "-bg", args.bg],
                      ["Organism", "-organism", args.organism]]

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, cell_align="left")
        html.add_free_content(['<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory,os.path.join(title,"parameters.html")))

    def table(self, directory, folder):
        arr = numpy.array([["#reference", "query", "background", "proportion", "p-value"]])
        for ty in self.plist.keys():
            for r in self.plist[ty].keys():
                for q in self.plist[ty][r].keys():
                    ar = numpy.array([[r, q, self.qlist[ty][r]['Background'], self.qlist[ty][r][q],self.plist[ty][r][q]]])
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
        
        #self.distriDict = OrderedDict()
        self.disperDict = OrderedDict()
        
        for ty in self.groupedreference.keys():
            #self.distriDict[ty] = OrderedDict()
            self.disperDict[ty] = OrderedDict()
            # Reference
            for r in self.groupedreference[ty]:
                r.merge()
                len_r = r.total_coverage()
                #self.distriDict[ty][r.name] = []
                self.disperDict[ty][r.name] = []
                
                for ch in self.chrom_list:
                    rc = r.any_chrom(chrom=ch)
                    nr = sum([len(s) for s in rc])
                    #self.distriDict[ty][r.name].append(nr)
                    self.disperDict[ty][r.name].append(nr/len_r)
                    
            # Query
            for q in self.groupedquery[ty]:
                q.merge()
                len_q = q.total_coverage()
                #self.distriDict[ty][q.name] = []
                self.disperDict[ty][q.name] = []
                
                for ch in self.chrom_list:
                    qc = q.any_chrom(chrom=chr)
                    nq = sum([len(s) for s in qc])
                    #self.distriDict[ty][q.name].append(nq)
                    self.disperDict[ty][q.name].append(nq/len_q)
            # Genome
            #self.distriDict[ty]["Genome"] = [len(genome.any_chrom(chrom=chr)) for chr in self.chrom_list]
            
            self.disperDict[ty]["Genome"] = [len(genome.any_chrom(chrom=chr)[0])/all_cov for chr in self.chrom_list]
        
    def plot_distribution(self):
        def to_percentage(x, pos=0): 
            return '{:.2f} %'.format(100*x)
         
        self.fig = []
        
        for ty in self.disperDict.keys():
            colors = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(self.disperDict[ty].keys()))).tolist()
            
            f, ax = plt.subplots()
            f.set_size_inches(10.5, 30)
            width = 0.9/len(self.disperDict[ty].keys())
            ind = np.arange(len(self.chrom_list))
            coverage = self.disperDict[ty]
            
            for ind_r, r in enumerate(self.disperDict[ty].keys()):
            
                ax.barh(ind + width*ind_r, self.disperDict[ty][r], width, color=colors[ind_r])
            
            plt.xlabel('Percentage')
            ax.xaxis.set_major_formatter(FuncFormatter(to_percentage)) 
            
            ax.minorticks_off()
            ax.set_yticks([ x + 0.5 for x in range(len(self.chrom_list))])
            ax.set_yticklabels(self.chrom_list,rotation=0, ha="right")
            ax.tick_params(axis='y', which='both', top='off', bottom='off', labelbottom='on')
            
            ax.legend(self.disperDict[ty].keys(), loc='center left', handlelength=1, handletextpad=1, 
                      columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right', 'left', 'bottom']:
                ax.spines[spine].set_visible(False)
            f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
            self.fig.append(f)

    def gen_html_distribution(self, outputname, title, align=50):
        fp = os.path.join(dir,outputname,title)
        link_d = {title:"distribution.html"}
        html = Html(name="Viz", links_dict=link_d, fig_dir=os.path.join(dir,outputname,"fig"), other_logo="viz")
        for i, f in enumerate(self.fig):
            html.add_figure("distribution_test_"+str(i)+".png", align="center")
        
        html.add_free_content(['<p style=\"margin-left: '+str(align+150)+'">'+
                               '** </p>'])
        
        type_list = 'ssssssssssssssssssssssssssssssssssssssssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]
        data_table = []
        for ind_ty, ty in enumerate(self.disperDict.keys()):
            header_list = ["Chromosome"] + self.disperDict[ty].keys()
            html.add_heading(ty, size = 4, bold = False)
            for i, ch in enumerate(self.chrom_list):
            #for ind_r,r in enumerate(self.disperDict[ty].keys()):
            
                data_table.append([ch]+["{:.3f} %".format(100 * self.disperDict[ty][r][i]) for r in self.disperDict[ty].keys()])
                  
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align)
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.write(os.path.join(fp,"distribution.html"))
        
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
        #self.color_list['Background'] = '0.70'
    
    def jaccard_test(self, runtime, organism):
        self.jlist = OrderedDict()
        self.realj = OrderedDict()
        self.plist = OrderedDict()
        self.rlen = {}
        self.qlen = {}
        self.rt = runtime
        self.nalist = []
        print2(self.parameter, "\nJaccard Test")
        print2(self.parameter, "{0:s}\t{1:s}\t{2:s}\t{3:s}\t{4:s}\t{5:s}".format("Reference","Query","Repeats", "True_Jaccard_index", "p-value", "Time"))
        for ty in self.groupedreference.keys():
            self.jlist[ty] = OrderedDict()
            self.realj[ty] = OrderedDict()
            self.plist[ty] = OrderedDict()
            for i, r in enumerate(self.groupedreference[ty]):
                if r.total_coverage() == 0 and len(r)>0:
                    self.nalist.append(r.name)
                    continue
                else:
                    if r.name not in self.rlen.keys(): self.rlen[r.name] = len(r)
                    self.jlist[ty][r.name] = OrderedDict()
                    self.realj[ty][r.name] = OrderedDict()
                    self.plist[ty][r.name] = OrderedDict()
                    for j, q in enumerate(self.groupedquery[ty]):
                        ts = time.time()
                        #print(q.name + "      " + str(len(q.sepuences[0])))
                        
                        # The real jaccard index from r and q
                        if q.total_coverage() == 0 and len(q)>0:
                            self.nalist.append(q.name)
                            continue
                        else:
                            if q.name not in self.qlen.keys(): self.qlen[q.name] = len(q)
                            self.jlist[ty][r.name][q.name] = []
                            self.realj[ty][r.name][q.name] = q.jaccard(r)
                            for k in range(runtime):
                                random = q.random_regions(organism=organism, multiply_factor=1, overlap_result=True, overlap_input=True, chrom_M=False)
                                self.jlist[ty][r.name][q.name].append(r.jaccard(random))
                            # How many randomizations have higher jaccard index than the real index?
                            p = len([x for x in self.jlist[ty][r.name][q.name] if x > self.realj[ty][r.name][q.name]])/runtime
                            self.plist[ty][r.name][q.name] = p
                            te = time.time()
                            print2(self.parameter, r.name +"\t"+ q.name +"\tx"+str(runtime)+"\t"+ 
                                   value2str(self.realj[ty][r.name][q.name]) +"\t"+ value2str(p) +"\t"+ 
                                   str(datetime.timedelta(seconds=round(te-ts))))    
    def plot(self, logT=False):
        """ Return boxplot from the given tables.
        
        """
        self.fig = []
        self.xtickrotation, self.xtickalign = 0,"center"
        
        for it, t in enumerate(self.jlist.keys()):
            f, axarr = plt.subplots(1, len(self.jlist[t].keys()), dpi=300, sharey=True)
            legend_x = 1.05
            nm = len(self.jlist.keys()) * len(self.jlist.values()[0]) * len(self.jlist.values()[0])
            if nm > 30:
                f.set_size_inches(nm * 0.1 +1 ,nm * 0.1 +1)
                legend_x = 1.2
                self.xtickrotation, self.xtickalign = 70,"right"
            try: axarr = axarr.reshape(-1)
            except: axarr = [axarr]
            plt.subplots_adjust(bottom=0.3)
            if logT:
                axarr[0].set_ylabel("Jaccard index (log)")
            else:
                axarr[0].set_ylabel("Jaccard index (Intersect/Union)")
            
            for i, r in enumerate(self.jlist[t].keys()):
                #axarr[i].set_title(r, y=0.94)
                if logT:
                    axarr[i].set_yscale('log')
                axarr[i].tick_params(axis='y', direction='out')
                axarr[i].yaxis.tick_left()
                axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
                d = []  # Store data within group
                color_t = []  # Store tag for coloring boxes
                #x_ticklabels = []  # Store ticklabels
                axarr[i].set_xlabel(r, rotation=self.xtickrotation, ha=self.xtickalign)
                
                
                
                for j,q in enumerate(self.jlist[t][r].keys()):
                    d.append(self.jlist[t][r][q])
                    color_t.append(self.color_list[q])
                    #x_ticklabels.append(q)
                # Fine tuning boxplot
                axarr[i].scatter(x=range(len(self.jlist[t][r].keys())), y=[y for y in self.realj[t][r].values()], s=2, c="red", edgecolors='none')
                bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, widths=None, 
                                      patch_artist=True, bootstrap=None)
                z = 10 # zorder for bosplot
                plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8,zorder=z)
                plt.setp(bp['fliers'], markerfacecolor='gray',color='white', alpha=0.3,markersize=1.8,zorder=z)
                plt.setp(bp['caps'],color='white', zorder=z)
                plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
                legends = []
                for patch, color in zip(bp['boxes'], color_t):
                    patch.set_facecolor(color) # When missing the data, the color patch will exceeds
                    patch.set_zorder(z)
                    legends.append(patch)
                
                # Fine tuning subplot
                #axarr[i].set_xticks(range(len(self.jlist[t][r].keys())))
                #plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
                #axarr[i].set_xticklabels(self.jlist[t][r].keys(), rotation=0, fontsize=10)
                
                #axarr[i].set_ylim(bottom=0.95)
                for spine in ['top', 'right', 'left', 'bottom']:
                    axarr[i].spines[spine].set_visible(False)
                axarr[i].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
                
                if i > 0:
                    plt.setp(axarr[i].get_yticklabels(),visible=False)
                    #plt.setp(axarr[i].get_yticks(),visible=False)
                    axarr[i].minorticks_off()
                    axarr[i].tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
            
            plt.setp([a.get_yticklabels() for a in axarr[1:]], visible=False)
            axarr[-1].legend(legends[0:len(self.jlist[t][r].keys())], self.jlist[t][r].keys(), loc='center left', handlelength=1, 
                     handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10},
                     bbox_to_anchor=(legend_x, 0.5))
            f.tight_layout(pad=2, h_pad=None, w_pad=None)
            self.fig.append(f)
  
    def gen_html(self, outputname, title, align=50):
        fp = os.path.join(dir,outputname,title)
        link_d = {title:"jaccard.html"}
        html = Html(name="Viz", links_dict=link_d, fig_dir=os.path.join(dir,outputname,"fig"), other_logo="viz")
        for i in range(len(self.fig)):
            html.add_figure("jaccard_test"+str(i+1)+".png", align="center")
        
        header_list = ["Reference<br>name",
                       "Query<br>name", 
                       "Reference<br>number", 
                       "Query<br>number", 
                       "True<br>Jaccard<br>index",
                       "Average<br>random<br>Jaccard",
                       "Positive<br>Association<br>p-value",
                       "Negative<br>Association<br>p-value"]
        
        type_list = 'sssssssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10,10,10,10,10]
        data_table = []
        
        for ind_ty, ty in enumerate(self.jlist.keys()):
            html.add_heading(ty, size = 4, bold = False)
            for ind_r,r in enumerate(self.jlist[ty].keys()):
                for ind_q, q in enumerate(self.jlist[ty][r].keys()):
                    rej = self.realj[ty][r][q]
                    rj = numpy.mean(self.jlist[ty][r][q])
                    p = self.plist[ty][r][q]
                    np = 1-p
                    rl = str(self.rlen[r])
                    ql = str(self.qlen[q])
                    
                    if self.plist[ty][r][q] < 0.05:
                        if self.realj[ty][r][q] > rj:
                            data_table.append([r,q,rl,ql,value2str(rej),value2str(rj),
                                               "<font color=\"red\">"+value2str(p)+"</font>",
                                               value2str(np)])
                        else:
                            data_table.append([r,q,rl,ql,value2str(rej),value2str(rj),
                                               value2str(np),
                                               "<font color=\"red\">"+value2str(p)+"</font>"])
                    else:
                        data_table.append([r,q,rl,ql,value2str(rej),value2str(rj), value2str(p),value2str(np)])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align)
        
        
        header_list=["Assumptions and hypothesis"]
        data_table = [['Randomization was performed '+ str(self.rt)+ ' times.'],
                      ['For projection test, the reference and query should have non-zero length in order to calculate its Jaccard index.'],
                      ['Positive association is defined by: True Jaccard index > Averaged random Jaccard.'],
                      ['Negative association is defined by: True Jaccard index < Averaged random Jaccard.']]
        self.nalist = set(self.nalist)
        if len(self.nalist)>0:
            data_table.append(['The following region sets contain zero-length regions which cause error in Jaccard index calculation, please check it:<br>'+
                               '<font color=\"red\">'+', '.join([s for s in self.nalist])+'</font>'])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, cell_align="left")
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.write(os.path.join(fp,"jaccard.html"))
        
    def table(self, directory, folder):
        arr = numpy.array([["#reference", "query", "true_jaccard", "random_jaccard", "p-value"]])
        for ty in self.plist.keys():
            for r in self.plist[ty].keys():
                for q in self.plist[ty][r].keys():
                    ar = numpy.array([[r, q, self.realj[ty][r][q], self.qlist[ty][r][q],self.plist[ty][r][q]]])
                    arr = numpy.vstack((arr, ar))
        output_array(arr, directory, folder, filename="output_table.txt")
        
###########################################################################################
#                    Inersection test
###########################################################################################

class Intersect:
    def __init__(self, reference_path, query_path, mode_count, organism):
        self.rEM, self.qEM = ExperimentalMatrix(), ExperimentalMatrix()
        self.rEM.read(reference_path)
        self.references = self.rEM.get_regionsets()
        self.referencenames = self.rEM.get_regionsnames()
        self.qEM.read(query_path)
        self.query = self.qEM.get_regionsets()
        self.querynames = self.qEM.get_regionsnames()
        self.mode_count = mode_count
        self.organism = organism
        self.sbar = None
        self.pbar = None
        self.test_d = None

    def background(self,path=None):
        """Given a bed file as the background for analysis"""
        bgbed = GenomicRegionSet(name="Background")
        if path:
            bgbed.read_bed(path)
        else:
            bgbed.get_genome_data(organism=self.organism)
        self.backgroung = bgbed

    def group_refque(self, groupby):
        self.groupedreference, self.groupedquery = group_refque(self.rEM, self.qEM, groupby)
        remove_duplicates(self.groupedreference)
        remove_duplicates(self.groupedquery)
        
    def colors(self, colorby, definedinEM, ref_que = "que"):
        """color_list is a Dict [query] : color """
        if ref_que == "que":
            self.color_list = color_groupded_region(self.qEM, self.groupedquery, colorby, definedinEM)
            if self.groupedquery.keys()[0] == "All region sets without grouping":
                self.color_tags = [n.name for n in self.groupedquery["All region sets without grouping"]]
            else:
                self.color_tags = gen_tags(self.qEM, colorby)
        elif ref_que == "ref":
            self.color_list = color_groupded_region(self.rEM, self.groupedreference, colorby, definedinEM)
            if self.groupedquery.keys()[0] == "All region sets without grouping":
                self.color_tags = [n.name for n in self.groupedquery["All region sets without grouping"]]
            else:
                self.color_tags = gen_tags(self.qEM, colorby)
            
    def colors_comb(self):
        """color_list is a list : color """
        
        if self.groupedquery.keys()[0] == "All region sets without grouping":
            self.color_tags = self.referencenames
        else:
            tags = []
            for t in [n.name for n in self.groupedreference.values()[0]]:
                nt = t.replace(self.groupedreference.keys()[0],"")
                nt = nt.replace("_","")
                tags.append(nt)
            self.color_tags = tags
        self.color_list = plt.cm.Set1(numpy.linspace(0.1, 0.9, len(self.color_tags))).tolist()
    
    def extend_ref(self,percentage):
        """percentage must be positive value"""
        for ty in self.groupedreference:
            for r in self.groupedreference[ty]:
                r.extend(left=percentage,right=percentage,percentage=True)
        
    def count_intersect(self, threshold, frequency=False):
        self.counts = OrderedDict()
        self.rlen, self.qlen = {}, {}
        self.nalist = []
        if frequency: self.frequency = OrderedDict()
        
        #if self.mode_count == "bp":
        #    print2(self.parameter, "\n{0}\t{1}\t{2}\t{3}\t{4}".format("Reference","Length(bp)", "Query", "Length(bp)", "Length of Intersection(bp)"))
        #elif self.mode_count == "count":
        #    print2(self.parameter, "\n{0}\t{1}\t{2}\t{3}\t{4}".format("Reference","sequence_number", "Query", "sequence_number", "Number of Intersection"))
        
        for ty in self.groupedreference.keys():
            self.counts[ty] = OrderedDict()
            self.rlen[ty], self.qlen[ty] = OrderedDict(), OrderedDict()
            if frequency: self.frequency[ty] = OrderedDict()
            
            for r in self.groupedreference[ty]:
                if r.total_coverage() == 0 and len(r)>0:
                    self.nalist.append(r.name)
                    continue
                else:
                    self.counts[ty][r.name] = OrderedDict()
                    if self.mode_count == "bp": rlen = r.total_coverage()
                    elif self.mode_count == "count": rlen = len(r)
                    self.rlen[ty][r.name] = rlen
                    
                    
                    for q in self.groupedquery[ty]:
                        if q.total_coverage() == 0 and len(q)>0:
                            self.nalist.append(q.name)
                            continue
                        else:
                            if self.mode_count == "bp": qlen = q.total_coverage()
                            elif self.mode_count == "count": qlen = len(q)
                            self.qlen[ty][q.name] = qlen
                            # Define different mode of intersection and count here
                            c = count_intersect(r,q, mode_count=self.mode_count, threshold=threshold)
                            self.counts[ty][r.name][q.name] = c
                            if frequency: 
                                try: self.frequency[ty][q.name].append(c[2])
                                except: 
                                    self.frequency[ty][q.name] = [c[2]]
                            
                            #print2(self.parameter, "{0}\t{1}\t{2}\t{3}\t{4}".format(r.name,rlen, q.name, qlen, c[2]))
            
    def barplot(self, logt=False, percentage=False):
        f, axs = plt.subplots(len(self.counts.keys()),1)
        f.subplots_adjust(left=0.3)
        self.xtickrotation, self.xtickalign = 0,"center"
        #if len(axs) == 1: axs = [axs]
        try: axs = axs.reshape(-1)
        except: axs = [axs]
        
        for ai, ax in enumerate(axs):
            if logt:
                ax.set_yscale('log')
                plus = 1
            else:
                ax.locator_params(axis = 'y', nbins = 4)
                plus = 0
                ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.set_title(self.counts.keys()[ai], y=1)
            r_label = []   
            for ind_r,r in enumerate(self.counts.values()[ai].keys()):
                for l in self.references:
                    if l.name == r: lr = len(l)
                    
                if len(axs) == 1: 
                    r_label.append(r)
                else: 
                    try: r_label.append(self.rEM.get_type(r,"factor"))
                    except: r_label.append(r)
                if len(r_label[-1]) > 15 or len(self.counts.values()[ai][r].keys())*len(self.counts.values()[ai].keys()) > 8: 
                    self.xtickrotation, self.xtickalign = 50,"right"
                width = 0.8/(len(self.counts.values()[ai][r].keys())+1) # Plus one background
                for ind_q, q in enumerate(self.counts.values()[ai][r].keys()):
                    x = ind_r + ind_q*width + 0.1
                    if percentage:
                        y = (self.counts.values()[ai][r][q][2] + plus)/lr
                    else:
                        y = self.counts.values()[ai][r][q][2] + plus # intersect number
                    
                    ax.bar(x, y, width=width, color=self.color_list[q], edgecolor="none", align='edge', log=logt)
                    
            ax.yaxis.tick_left()
            ax.set_xticks([i + 0.5 - 0.5*width for i in range(len(r_label))])
            ax.set_xticklabels(r_label,fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
            ax.set_xlim([0, len(self.counts.values()[ai].keys())-0.1])
            
            ax.legend(self.color_tags, loc='center left', handlelength=1, handletextpad=1, 
                      columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            if percentage:
                f.text(-0.025, 0.5, "Intersected percantage", rotation="vertical", va="center")
            else:
                f.text(-0.025, 0.5, "Intersected regions number", rotation="vertical", va="center")
            
        f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        self.bar = f

    def stackedbar(self):
        f, axs = plt.subplots(len(self.counts.keys()),1)
        f.subplots_adjust(left=0.3)
        #if len(axs) == 1: axs = [axs]
        try: axs = axs.reshape(-1)
        except: axs = [axs]
        
        for ai, ax in enumerate(axs):
            ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.locator_params(axis = 'y', nbins = 2)
            ax.set_title(self.counts.keys()[ai], y=1)
            r_label = []   
            for ind_r,r in enumerate(self.counts.values()[ai].keys()):
                if len(axs) == 1: r_label.append(r)
                else: 
                    try: r_label.append(self.rEM.get_type(r,"factor"))
                    except: r_label.append(r)
                width = 0.6
                bottom = 0
                for ind_q, q in enumerate(self.counts.values()[ai][r].keys()):
                    x = ind_r
                    y = self.counts.values()[ai][r][q][2] # intersect number
                    ax.bar(x, y, width=width, bottom=bottom, color=self.color_list[q], edgecolor="none", align='center')
                    bottom = bottom + y
            ax.yaxis.tick_left()
            ax.set_xticks(range(len(r_label)))
            ax.set_xticklabels(r_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
            ax.set_xlim([-0.5, ind_r+0.5])
            
            ax.legend(self.color_tags, loc='center left', handlelength=1, handletextpad=1, 
                      columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            f.text(-0.025, 0.5, "Intersected regions number", rotation="vertical", va="center")
    
        f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        self.sbar = f

    def percentagebar(self):
        self.color_list["No intersection"] = "0.7"
        f, axs = plt.subplots(len(self.counts.keys()),1)
        f.subplots_adjust(left=0.3)
        #if len(axs) == 1: axs = [axs]
        try: axs = axs.reshape(-1)
        except: axs = [axs]
        self.percentage = []
        
        for ai, ax in enumerate(axs):
            #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ax.locator_params(axis = 'y', nbins = 2)
            ax.set_title(self.counts.keys()[ai], y=1.05)
            r_label = []
            self.percentage.append({})
            for ind_r,r in enumerate(self.counts.values()[ai].keys()):
                if len(axs) == 1: r_label.append(r)
                else: 
                    try: r_label.append(self.rEM.get_type(r,"factor"))
                    except: r_label.append(r)
                width = 0.6
                bottom = 0
                if self.mode_count == "bp":
                    sumlength = self.rEM.objectsDict[r].total_coverage()
                elif self.mode_count == "count":
                    sumlength = len(self.rEM.objectsDict[r])
                self.percentage[ai][r] = {}
                
                if sumlength == 0:
                    for ind_q, q in enumerate(self.counts.values()[ai][r].keys()):
                        self.percentage[ai][r][q] = "ref is empty"
                else:
                    for ind_q, q in enumerate(self.counts.values()[ai][r].keys()):
                        x = ind_r
                        y = int(self.counts.values()[ai][r][q][2])/sumlength # percentage
                        ax.bar(x, y, width=width, bottom=bottom, color=self.color_list[q], edgecolor="none", align='center')
                        bottom = bottom + y
                        self.percentage[ai][r][q] = y
                ax.bar(x,1-bottom,width=width, bottom=bottom, color=self.color_list["No intersection"], align='center')
            ax.yaxis.tick_left()
            ax.set_xticks(range(len(r_label)))
            ax.set_xticklabels(r_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
            ax.set_xlim([-0.5, ind_r+0.5])
            ax.set_ylim([0,1])
            legend_labels = self.color_tags + ["No intersection"]
            ax.legend(legend_labels, loc='center left', handlelength=1, handletextpad=1, 
                      columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
        f.text(-0.025, 0.5, "Proportion of intersected regions (%)", rotation="vertical", va="center")
    
        f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        self.pbar = f

    def gen_html(self, directory, title, align, args):
        #fp = os.path.join(dir,outputname,title)
        #link_d = {title:"intersection.html"}
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = "Intersection Test: "+dir_name
        link_d = OrderedDict()
        link_d["Intersection test"] = "index.html"
        link_d["Parameters"] = "parameters.html"

        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        
        html.add_figure("intersection_bar.png", align="center")
        if self.sbar: html.add_figure("intersection_stackedbar.png", align="center")
        html.add_figure("intersection_barp.png", align="center")
        
        header_list = ["Reference<br>name",
                       "Query<br>name", 
                       "Reference<br>number", 
                       "Query<br>number", 
                       "Intersect.",
                       "Proportion <br>of Reference"]
       
        if self.test_d: 
            header_list += ["Average<br>intersect.", "Chi-square<br>statistic", 
                            "Positive<br>Association<br>p-value", "Negative<br>Association<br>p-value"]
        else: pass
        
        type_list = 'ssssssssssssssss'
        col_size_list = [10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10]
        
        for ind_ty, ty in enumerate(self.counts.keys()):
            html.add_heading(ty, size = 4, bold = False)
            data_table = []
            for ind_r,r in enumerate(self.counts[ty]):
                for ind_q, q in enumerate(self.counts[ty][r]):
                    if r == q: continue
                    pt = self.counts[ty][r][q][2]/self.rlen[ty][r]
                    intern = self.counts[ty][r][q][2]
                    if self.test_d:
                        aveinter = self.test_d[ty][r][q][0]
                        chisqua = value2str(self.test_d[ty][r][q][1])
                        pv = self.test_d[ty][r][q][2]
                        if isinstance(pv, str):
                            data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                              str(intern), "{:.2f}%".format(100*pt),
                                              aveinter, chisqua, pv,"-"])
                        else: 
                            npv = 1 - pv
                            if pv < 0.05:
                                if intern > aveinter:
                                    data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                                       str(intern), "{:.2f}%".format(100*pt),
                                                       value2str(aveinter), chisqua, "<font color=\"red\">"+value2str(pv)+"</font>", value2str(npv)])
                                else:
                                    data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                                       str(intern), "{:.2f}%".format(100*pt),
                                                       value2str(aveinter), chisqua, value2str(npv), "<font color=\"red\">"+value2str(pv)+"</font>"])
                            elif self.test_d[ty][r][q][2] >= 0.05:
                                if intern > aveinter:
                                    data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                                       str(intern), "{:.2f}%".format(100*pt),
                                                       value2str(aveinter), chisqua, value2str(pv),value2str(npv)])
                                else:
                                    data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                                       str(intern), "{:.2f}%".format(100*pt),
                                                       value2str(aveinter), chisqua, value2str(npv),value2str(pv)])
                    else:
                        data_table.append([r, q,str(self.rlen[ty][r]), str(self.qlen[ty][q]), 
                                           str(intern), "{:.2f}%".format(100*pt)])
        
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, sortable=True)
        
        html.add_heading("Assumptions and hypothesis")
        list_ex = ['Positive association is defined by: True intersection number > Averaged random intersection.',
                   'Negative association is defined by: True intersection number < Averaged random intersection.']
        
        self.nalist = set(self.nalist)
        if len(self.nalist)>0:
            list_ex.append('The following region sets contain zero-length regions which cause error in intersection calculation, please check it:<br>'+
                           '<font color=\"red\">'+', '.join([s for s in self.nalist])+'</font>')
        if self.test_d:
            list_ex.append('Randomly permutation for '+ str(self.test_time)+ ' times.')
        html.add_list(list_ex)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))

        # Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        header_list = ["Description", "Argument", "Value"]
        data_table = [["Reference", "-r", args.r ],
                      ["Query", "-q", args.q],
                      ["Output directory", "-o", os.path.basename(args.o)],
                      ["Experiment title", "-t", args.t],
                      #["Grouping tag", "-g", args.g],
                      #["Coloring tag", "-c", args.c],
                      #["Background", "-bg", args.bg],
                      ["Organism", "-organism", args.organism]]

        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, cell_align="left")
        html.add_free_content(['<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See details</a>'])
        html.write(os.path.join(directory, title,"parameters.html"))
    
    def gen_html_comb(self, outputname, title, align):
        fp = os.path.join(dir,outputname,title)
        link_d = {title:"combinatorial.html"}
        html = Html(name="Viz", links_dict=link_d, fig_dir=os.path.join(dir,outputname,"fig"))
        #html.create_header()
        #html.add_heading(title)
        
        if self.sbar: html.add_figure("intersection_stackedbar.png", align="center")
        
        html.add_free_content(['<p style=\"margin-left: '+str(align+150)+'">'+
                               '** </p>'])
        
        for ind_ty, ty in enumerate(self.groupedreference.keys()):
            html.add_heading(ty, size = 4, bold = False)
            
            data_table = []
            header_list = ["Query<br>name", "Query<br>number", "Frequencies"]
            type_list = 'sssss'
            col_size_list = [10,10,30,10] 
        
            for ind_q, q in enumerate(self.frequency[ty].keys()):
                data_table.append([q,  str(self.qlen[ty][q]), ",".join(["%05d" % v for v in self.frequency[ty][q]])])
                
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align)
            
            header_list = ["Query 1", "Query 2", "Chi square", "p value"]
            data_table = []  
            for ind_q1, q1 in enumerate(self.frequency[ty].keys()[:-1]):
                for ind_q2, q2 in enumerate(self.frequency[ty].keys()[ind_q1+1:]):
                    if q1 == q2: continue
                    else:
                        chisq, p, dof, expected = stats.chi2_contingency([self.frequency[ty][q1],self.frequency[ty][q2]])
                        if p < 0.05:
                            data_table.append([q1,q2,value2str(chisq), "<font color=\"red\">"+value2str(p)+"</font>"])
                        else:
                            data_table.append([q1,q2,value2str(chisq), value2str(p)])
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align)
            
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
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="reference_experimental_matrix.txt" style="margin-left:100">See reference experimental matrix</a>'])
        html.add_free_content(['<a href="query_experimental_matrix.txt" style="margin-left:100">See query experimental matrix</a>'])
        html.write(os.path.join(fp,"combinatorial.html"))
    
    def posi2region(self, regions, p):
        all = range(len(regions))
        new_r = GenomicRegionSet(name="")
        for r in p:
            new_r.combine(regions[r])
        return new_r
    
    def posi2set(self,regions, p):
        all = range(len(regions))
        inter_r = copy.deepcopy(regions[p[0]])

        for i in all:
            #print("inter_r: "+inter_r.name)
            if i in p[1:]:
                inter_r = inter_r.intersect(regions[i],mode=OverlapType.OVERLAP)
            elif i == p[0]: pass
            else:
                inter_r = inter_r.subtract(regions[i], whole_region=False)
        #print("inter_r: "+inter_r.name)
        return inter_r
    
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
        
        for ty in self.groupedreference.keys():
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []
            
            for i in range(1,n):
                new_refsp[ty].append(itertools.combinations(range(n),i))
            for posi in new_refsp[ty]:
                posi = [list(i) for i in posi]
                
                for p in posi:
                    #print("   " + str(p))
                    pr = self.posi2set(self.groupedreference[ty],p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)
                    self.comb_ref_infor[pr.name] = p2sign(p,n)
            all_int = self.posi2set(self.groupedreference[ty],range(n))
            new_refs[ty].append(all_int)
            ref_names.append(all_int.name)
            self.comb_ref_infor[all_int.name] = p2sign(range(n),n)
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
        #self.comb_reference = new_refs
        self.groupedreference = copy.deepcopy(new_refs)
        self.orig_refs = copy.deepcopy(self.referencenames)
        self.referencenames = list(set(ref_names))
    
    def combine_regions(self, background=None):
        new_refsp = OrderedDict()
        new_refs = OrderedDict()
        ref_names = []
        
        for ty in self.groupedreference.keys():
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []
            for i in range(1,n):
                new_refsp[ty].append(itertools.combinations(range(n),i))
            for posi in new_refsp[ty]:
                posi = [list(i) for i in posi]
                for p in posi:
                    print("   " + str(p))
                    pr = self.posi2region(self.groupedreference[ty],p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)
    
    def comb_stacked_plot(self):
        self.xtickrotation, self.xtickalign = 0,"center"
        f, axs = plt.subplots(1, len(self.frequency.keys()), sharey = True)
        f.subplots_adjust(left=0.3)
        #f.set_size_inches(18.5,15)
        #if len(axs) == 1: axs = [axs]
        try: axs = axs.reshape(-1)
        except: axs = [axs]
        
        for ai, ax in enumerate(axs):
            #ax.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
            ty = self.frequency.keys()[ai]
            ax.locator_params(axis = 'y', nbins = 4)
            ax.set_title(self.frequency.keys()[ai], y=1)
            r_label = []
            q_label = []   
            legends = []
            
            for ind_q,q in enumerate(self.frequency[ty].keys()):
                if len(axs) == 1: 
                    q_label.append(q)
                else: 
                    try: q_label.append(self.qEM.get_type(q,"factor"))
                    except: q_label.append(q)
                width = 0.6
                bottom = 0
                summ = sum(self.frequency[ty][q])
                for ind_r, rc in enumerate(self.frequency[ty][q]):
                    if ind_q == 0: 
                        r = self.groupedreference[ty][ind_r].name
                        r_label.append(r)
                    x = ind_q
                    y = rc/summ # intersect number
                    bar = ax.bar(x, y, width=width, bottom=bottom, color=self.color_list[ind_r], align='center')
                    bottom = bottom + y
                    if ind_q == 0: legends.append(bar) 
            ax.yaxis.tick_left()
            ax.set_xticks(range(len(q_label)))
            ax.set_xticklabels(q_label, fontsize=9, rotation=self.xtickrotation, ha=self.xtickalign)
            ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
            ax.set_xlim([-0.5, ind_q+0.5])
            
            legends.reverse()
            r_label.reverse()
            
            for spine in ['top', 'right']:  # 'left', 'bottom'
                ax.spines[spine].set_visible(False)
                
        legend_name = reversed(self.color_tags)
        axs[-1].legend(legends, legend_name, loc='upper left', handlelength=1, handletextpad=1, 
                  columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 1))
        if self.mode_count == "bp":
            f.text(-0.025, 0.5, "Intersected regions (bp)", rotation="vertical", va="center")
        elif self.mode_count == "count":
            f.text(-0.025, 0.5, "Percentage", rotation="vertical", va="center")
    
        f.tight_layout(pad=0.5, h_pad=None, w_pad=0.5)
        self.sbar = f
        
        
    def stest(self,repeat,threshold):
        print("\nIntersection random subsampling test:\n    Repeat "+str(repeat)+" times\n")
        self.test_time = repeat
        self.test_d = {}
        plist = OrderedDict()
        
        for ty in self.groupedreference.keys():
            self.test_d[ty] = {}
            plist[ty] = OrderedDict()
            for r in self.groupedreference[ty]:
                if r.name in self.nalist: continue
                self.test_d[ty][r.name] = {}
                plist[ty][r.name] = OrderedDict()
                for q in self.groupedquery[ty]:
                    if q.name in self.nalist: continue
                    # True intersection
                    obs = self.counts[ty][r.name][q.name]
                    qn = q.name
                    if obs[2] == 0:
                        aveinter, chisq, p = "NA", "NA", "NA"
                    else:
                        com = q.combine(r, change_name=False, output=True)
                        # Randomization
                        d = []
                        for i in range(repeat):
                            random_r,random_q = com.random_split(size=self.rlen[ty][r.name])                           
                            d.append(count_intersect(random_r, random_q, mode_count=self.mode_count, threshold=threshold))
                        da = numpy.array(d)
                        
                        exp_m = numpy.mean(da, axis=0)
                        chisq, p, dof, expected = stats.chi2_contingency([exp_m,obs])
                        aveinter = exp_m[2]
                    plist[ty][r.name][qn] = p
                    self.test_d[ty][r.name][qn] = [aveinter, chisq, p]
                    
            multiple_correction(plist)
            
            #c_p = 0
            for r in self.test_d[ty].keys():
                if r in self.nalist: continue
                for q in self.test_d[ty][r].keys():
                    self.test_d[ty][r][q][2] = plist[ty][r][q]
                    
###########################################################################################
#                    Proportion plot
###########################################################################################

class Proportion:

    def __init__(self, queEM, refEM, organism):
        """Initiate the proportion plot"""
        self.exp_que = ExperimentalMatrix()
        self.exp_que.read(queEM)
        self.exp_ref = ExperimentalMatrix()
        self.exp_ref.read(refEM)
        
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
    def __init__(self,EMpath, title="boxplot", df=False):
        # Read the Experimental Matrix
        self.title = title
        self.exps = ExperimentalMatrix()
        self.exps.read(EMpath)
        self.beds = self.exps.get_regionsets() # A list of GenomicRegionSets
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
        self.all_bed.remove_duplicates() #all_bed is sorted!!
        
    
    def bedCoverage(self):
        """ Return coverage matrix of multiple reads on one bed. 
        bed --> GenomicRegionSet
        """
        c=[]
        for rp in self.reads:
            print("    processing: "+rp)
            r = os.path.abspath(rp)   # Here change the relative path into absolute path
            cov = CoverageSet(r,self.all_bed)
            cov.coverage_from_genomicset(r)
            cov.normRPM()
            c.append(cov.coverage)
        self.all_table = numpy.transpose(c)
      
    def quantile_normalization(self):
        """ Return the np.array which contains the normalized values
        """
        rank_matrix = []
        for c in range(self.all_table.shape[1]):
            col = self.all_table[:,c]
            rank_col = mstats.rankdata(col)
            rank_matrix.append(rank_col)
    
        ranks = numpy.array(rank_matrix)
        trans_rank = numpy.transpose(ranks)
        
        # Calculate for means of ranks
        print("    Calculating for the mean of ranked data...")
        sort_matrix = numpy.sort(self.all_table,axis=0)
        means = []
        for r in range(self.all_table.shape[0]):
            row = [x for x in sort_matrix[r,:]]
            means.append(numpy.mean(row))
    
        # Replace the value by new means
        print("    Replacing the data value by normalized mean...")
        normalized_table = numpy.around(trans_rank)
        for i, v  in enumerate(means):
            normalized_table[normalized_table == i+1] = v
        #print(rounded_rank)
        self.norm_table = normalized_table

    def tables_for_plot(self):
        """ Return a Dict which stores all tables for each bed with file name as its key. """
        self.tableDict = OrderedDict() # Storage all tables for each bed with bedname as the key
        conList = []   # Store containers of beds
        iterList = []
        
        for i,bed in enumerate(self.beds):
            self.tableDict[bed.name] = []
            bed.sort()
            conList.append(bed.__iter__())
            iterList.append(conList[-1].next())
            
        for i, r in enumerate(self.all_bed.sequences):
            for j in range(len(self.beds)):
                while r > iterList[j]:
                    try:
                        iterList[j] = conList[j].next()
                    except:
                        break
                if r == iterList[j]:
                    self.tableDict[self.beds[j].name].append(self.norm_table[i])
                elif r < iterList[j]:
                    continue

    def print_plot_table(self, directory, folder):
        for i,bed in enumerate(self.tableDict.keys()):
            table = []
            header = ["chrom", "initial", "final"]
            for rp in self.reads:
                header.append(os.path.basename(rp))
            table.append(header)
            for j, re in enumerate(self.beds[i]):
                table.append([re.chrom, re.initial, re.final] + self.tableDict[bed][j].tolist())
            output_array(table, directory, folder, filename="table_"+bed+".txt")  
        
        
    def group_tags(self, groupby, sortby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        self.tag_type = [groupby, sortby, colorby]
        
        if groupby == "None":
            self.group_tags = ["All region sets without grouping"]
        else:
            self.group_tags = gen_tags(self.exps, groupby)
        if sortby == "None":
            self.sort_tags = ["All region sets without grouping"]
        else:
            self.sort_tags = gen_tags(self.exps, sortby)
        if colorby == "None":
            self.color_tags = ["All region sets without grouping"]
        else:
            self.color_tags = gen_tags(self.exps, colorby)
    
    
        
    def group_data(self, directory, folder, log=False):  
        plotDict = OrderedDict()  # Extracting the data from different bed_bams file
        cuesbed = OrderedDict()   # Storing the cues for back tracking
        cuesbam = OrderedDict()
        for bedname in self.tableDict.keys():
            plotDict[bedname] = OrderedDict()
            mt = numpy.array(self.tableDict[bedname])
            
            cuesbed[bedname] = set(tag_from_r(self.exps, self.tag_type, bedname))
            #cuesbed[bedname] = [tag for tag in self.exps.get_types(bedname) if tag in self.group_tags + self.sort_tags + self.color_tags]
            
            for i,readname in enumerate(self.readsnames):
                plotDict[bedname][readname] = mt[:,i]
                #print(plotDict[bedname][readname])
                cuesbam[readname] = set(tag_from_r(self.exps, self.tag_type, readname))
                #cuesbam[readname] = [tag for tag in self.exps.get_types(readname) if tag in self.group_tags + self.sort_tags + self.color_tags]
            
        sortDict = OrderedDict()  # Storing the data by sorting tags
        for g in self.group_tags:
            #print("    "+g)
            sortDict[g] = OrderedDict()
            for a in self.sort_tags:
                #print("        "+a)
                sortDict[g][a] = OrderedDict()
                for c in self.color_tags:
                    #sortDict[g][a][c] = None
                    #print("            "+c)
                    for i, bed in enumerate(cuesbed.keys()):
                        if set([g,a,c]) >= cuesbed[bed]:
                            sortDict[g][a][c] = []
                            for bam in cuesbam.keys():
                                if set([g,a,c]) >= cuesbam[bam]:
                                    if self.df:
                                        sortDict[g][a][c].append(plotDict[bed][bam])
                                        if len(sortDict[g][a][c]) == 2: 
                                            bam2 = bam
                                            if log:
                                                sortDict[g][a][c][0] = numpy.log(sortDict[g][a][c][0])
                                                sortDict[g][a][c][1] = numpy.log(sortDict[g][a][c][1])
                                                sortDict[g][a][c] = numpy.subtract(sortDict[g][a][c][0],sortDict[g][a][c][1]).tolist()
                                            else:
                                                sortDict[g][a][c] = numpy.subtract(sortDict[g][a][c][0],sortDict[g][a][c][1]).tolist()
                                        else: 
                                            bam1 = bam
                                    else:
                                        sortDict[g][a][c] = plotDict[bed][bam]
        self.sortDict = sortDict

    def color_map(self, colorby, definedinEM):
        self.colors = colormap(self.exps, colorby, definedinEM)
        
    
    def print_table(self, directory, folder):
        self.printtable = OrderedDict()
        table = []
        table.append(["#group_tag", "sort_tag", "color_tag", "Signals"])
        for i, g in enumerate(self.group_tags):
            for k, a in enumerate(self.sort_tags):
                for j, c in enumerate(self.color_tags):
                    table.append([g, a, c] + [str(x) for x in self.sortDict[g][a][c]])
        #print(table)
        output_array(table, directory, folder, filename="output_table.txt")  
                    
    def plot(self, title, sy, logT=False, ylim=False):
        """ Return boxplot from the given tables.
        
        """
        f, axarr = plt.subplots(1, len(self.group_tags), dpi=300, sharey = sy)
        self.xtickrotation, self.xtickalign = 0,"center"

        nm = len(self.group_tags) * len(self.color_tags) * len(self.sort_tags)
        if nm > 30:
            f.set_size_inches(nm * 0.25 ,nm * 0.15)
            #legend_x = 1.2
            self.xtickrotation, self.xtickalign = 70,"right"
        
        canvas = FigureCanvas(f)
        canvas.set_window_title(title)
        try: axarr = axarr.reshape(-1)
        except: axarr = [axarr]
        #plt.subplots_adjust(bottom=0.3)
        if logT:
            if self.df: axarr[0].set_ylabel("Count number difference (log)")
            else: axarr[0].set_ylabel("Count number (log)")
        else:
            
            if self.df: axarr[0].set_ylabel("Count number difference")
            else: axarr[0].set_ylabel("Count number")
            
            
        for i, g in enumerate(self.sortDict.keys()):
            if self.df: axarr[i].set_title(g+"_df", y=1.02)
            else: axarr[i].set_title(g, y=1.02)
            
            
            if logT and not self.df: 
                axarr[i].set_yscale('log')
            else:
                axarr[i].locator_params(axis = 'y', nbins = 4)

            axarr[i].tick_params(axis='y', direction='out')
            axarr[i].yaxis.tick_left()
            axarr[i].yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.7, zorder=1)
            if ylim:
                axarr[i].set_ylim([-ylim, ylim])
            d = []  # Store data within group
            color_t = []  # Store tag for coloring boxes
            x_ticklabels = []  # Store ticklabels
            for j, a in enumerate(self.sortDict[g].keys()):
                for k, c in enumerate(self.sortDict[g][a].keys()):
                    if self.sortDict[g][a][c] == None:  # When there is no matching data, skip it
                        continue
                    else:
                        if self.df: 
                            d.append(self.sortDict[g][a][c])
                        else:
                            d.append([x+1 for x in self.sortDict[g][a][c]])
                        color_t.append(self.colors[k])
                        x_ticklabels.append(a)  #  + "." + c
            # Fine tuning boxplot
            #print(d)
            bp = axarr[i].boxplot(d, notch=False, sym='o', vert=True, whis=1.5, positions=None, 
                                  widths=None, patch_artist=True, bootstrap=None)
            z = 10 # zorder for bosplot
            plt.setp(bp['whiskers'], color='black',linestyle='-',linewidth=0.8,zorder=z)
            plt.setp(bp['fliers'], markerfacecolor='gray',color='white',alpha=0.3,markersize=1.8,zorder=z)
            plt.setp(bp['caps'],color='white',zorder=z)
            plt.setp(bp['medians'], color='black', linewidth=1.5,zorder=z+1)
            legends = []
            for patch, color in zip(bp['boxes'], color_t):
                patch.set_facecolor(color) # When missing the data, the color patch will exceeds
                patch.set_edgecolor("none") 
                patch.set_zorder(z)
                legends.append(patch)
                
            # Fine tuning subplot
            axarr[i].set_xticks([len(self.color_tags)*n + 1 + (len(self.color_tags)-1)/2 for n,s in enumerate(self.sortDict[g].keys())])
            #plt.xticks(xlocations, sort_tags, rotation=90, fontsize=10)
            axarr[i].set_xticklabels(self.sortDict[g].keys(), self.xtickrotation, ha=self.xtickalign, fontsize=10)
            
            #axarr[i].set_ylim(bottom=0.95)
            for spine in ['top', 'right', 'left', 'bottom']:
                axarr[i].spines[spine].set_visible(False)
            axarr[i].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
            
            if sy: 
                #plt.setp(axarr[i].get_yticklabels(),visible=False)
                axarr[i].minorticks_off()
                
                #axarr[i].tick_params(axis='y', which='both', left='off', right='off', labelbottom='off')
            else: 
                plt.setp(axarr[i].get_yticklabels(),visible=True)
                axarr[i].tick_params(axis='y', which='both', left='on', right='off', labelbottom='on')
            #plt.setp(axarr[i].get_yticks(),visible=False)
                
                
        axarr[-1].legend(legends[0:len(self.color_tags)], self.color_tags, loc='center left', handlelength=1, 
                 handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10},
                 bbox_to_anchor=(1.05, 0.5))
        f.tight_layout(pad=2, h_pad=None, w_pad=None)
        self.fig = f
    
    def gen_html(self, directory, title, align=50):
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = title
        link_d = OrderedDict()
        link_d["Boxplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"


        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        #fp = os.path.join(dir,outputname,title)
        
        html.add_figure("boxplot.png", align="center")
        
        
        type_list = 'ssssssssssssssssssssssssssssssssssssssssssssss'
        
        #### Calculate p value ####
        plist = {}
        for g in self.sortDict.keys():
            plist[g] = {}
            for s1 in self.sortDict[g].keys():
                for c1 in self.sortDict[g][s1].keys():
                    data1 = self.sortDict[g][s1][c1]
                    plist[g][s1+c1] = {}
                    for s2 in self.sortDict[g].keys():
                        for c2 in self.sortDict[g][s2].keys():
                            if s2 == s1 and c2 == c1: pass
                            else:
                                data2 = self.sortDict[g][s2][c2]
                                u, p_value = mannwhitneyu(data1, data2)
                                plist[g][s1+c1][s2+c2] = p_value
        
        print("Multiple test correction.")
        multiple_correction(plist)
        
        for g in self.sortDict.keys():
            html.add_heading(g, size = 4, bold = False)
            data_table = []
            col_size_list = [15]
            header_list = ["p-value"]
            for s in self.sortDict[g].keys():
                for c in self.sortDict[g][s1].keys():
                    header_list.append(s+"\n"+c)        
                    col_size_list.append(15)
            
            for s1 in self.sortDict[g].keys():
                for c1 in self.sortDict[g][s1].keys():
                    row = [s1+"\n"+c1]
                    for s2 in self.sortDict[g].keys():
                        for c2 in self.sortDict[g][s2].keys():
                            if s2 == s1 and c2 == c1: 
                                row.append("-")
                            else:
                                p = plist[g][s1+c1][s2+c2]
                                if p > 0.05:
                                    row.append(value2str(p))
                                else:
                                    row.append("<font color=\"red\">"+value2str(p)+"</font>")
                    data_table.append(row)
        
            html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align+50)
        
        #html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        header_list=["Assumptions and hypothesis"]
        col_size_list = [50]
        data_table = [['All the regions among different BED files are normalized by quantile normalization.'],
                      ['If there is any grouping problem, please check all the optional columns in input experimental matrix.']]
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, cell_align="left")
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])
        html.write(os.path.join(directory, title, "parameters.html"))

        
###########################################################################################
#                    Lineplot 
###########################################################################################

def load_dump(path, filename):
    print("\tLoading from file: "+filename)
    file = open(os.path.join(path,filename),'r')
    object = pickle.load(file)
    file.close()
    return object

def dump(object, path, filename):
    print("\tDump to file: "+filename)
    file = open(os.path.join(path,filename),'wb')
    pickle.dump(object,file)
    file.close()

def read_gtf(gd, organism):
    try:
        anno = load_dump(gd.get_annotation_dump_dir(),"gtf.dump")
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
    try: tss = load_dump(gd.get_annotation_dump_dir(),"tss.dump")
    except:
        anno = read_gtf(gd, organism)
        tss = anno.get_tss()
        dump(tss, gd.get_annotation_dump_dir(),"tss.dump")
    beds.append(tss)
    print("\tTSS: "+str(len(tss)))
    # TTS
    try: tts = load_dump(gd.get_annotation_dump_dir(),"tts.dump")
    except:
        anno = read_gtf(gd, organism)
        tts = anno.get_tts()
        dump(tts, gd.get_annotation_dump_dir(),"tts.dump")
    beds.append(tts)
    print("\tTTS: "+str(len(tts)))
    
    
    # exon
    try: exon = load_dump(gd.get_annotation_dump_dir(),"exon.dump")
    except:
        anno = read_gtf(gd, organism)
        exon = anno.get_exons(start_site=False, end_site=False)
        dump(exon, gd.get_annotation_dump_dir(),"exon.dump")
    print("\texon: "+str(len(exon)))
    # exons
    try: exons = load_dump(gd.get_annotation_dump_dir(),"exons.dump")
    except:
        anno = read_gtf(gd, organism)
        exons = anno.get_exons(start_site=True, end_site=False)
        dump(exons, gd.get_annotation_dump_dir(),"exons.dump")
    beds.append(exons)
    print("\texons: "+str(len(exons)))
    # exone
    try: exone = load_dump(gd.get_annotation_dump_dir(),"exone.dump")
    except:
        anno = read_gtf(gd, organism)
        exone = anno.get_exons(start_site=False, end_site=True)
        dump(exone, gd.get_annotation_dump_dir(),"exone.dump")
    beds.append(exone)
    print("\texone: "+str(len(exone)))
    
    
    #bednames = ["TSS", "TTS", "Exon start site", "Exon end site", "Intron start site", "Intron end site"]
    bednames = ["TSS", "TTS", "Exon start site", "Exon end site"]
    
    annotation = bednames

    return beds, bednames, annotation


class Lineplot:
    def __init__(self, EMpath, title, annotation, organism, center, extend, rs, bs, ss, df):
        # Read the Experimental Matrix
        self.title = title
        self.exps = ExperimentalMatrix()
        self.exps.read(EMpath)
        if annotation:
            self.beds, self.bednames, self.annotation = annotation_dump(organism)

        else:
            self.beds = self.exps.get_regionsets() # A list of GenomicRegionSets
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
    
    def relocate_bed(self):
        processed_beds = []
        processed_bedsF = [] # Processed beds to be flapped
        
        for bed in self.beds:
            if self.center == 'bothends':
                newbed = bed.relocate_regions(center='leftend', left_length=self.extend+int(0.5*self.bs), right_length=self.extend+int(0.5*self.bs))
                processed_beds.append(newbed)
                newbedF = bed.relocate_regions(center='rightend', left_length=self.extend+int(0.5*self.bs), right_length=self.extend+int(0.5*self.bs))
                processed_bedsF.append(newbedF)
            else:
                newbed = bed.relocate_regions(center=self.center, left_length=self.extend+int(0.5*self.bs), right_length=self.extend+int(0.5*self.bs))
                processed_beds.append(newbed)
        self.processed_beds = processed_beds
        self.processed_bedsF = processed_bedsF
        
    def group_tags(self, groupby, sortby, colorby):
        """Generate the tags for the grouping of plot
        Parameters:
            groupby = 'reads','regions','cell',or 'factor'
            colorby = 'reads','regions','cell',or 'factor'
            sortby = 'reads','regions','cell',or 'factor'
        """
        self.tag_type = [sortby, groupby, colorby]
        if groupby == "None":
            self.group_tags = ["All region sets without grouping"]
        elif groupby == "regions" and self.annotation:
            self.group_tags = self.bednames
        else:
            self.group_tags = gen_tags(self.exps, groupby)
            
        if sortby == "None":
            self.sort_tags = ["All region sets without grouping"]
        elif sortby == "regions" and self.annotation:
            self.sort_tags = self.bednames
        else:
            self.sort_tags = gen_tags(self.exps, sortby)
            
        if colorby == "None":
            self.color_tags = ["All region sets without grouping"]
        elif colorby == "regions" and self.annotation:
            self.color_tags = self.bednames
        else:
            self.color_tags = gen_tags(self.exps, colorby)
        
    def gen_cues(self):
        self.cuebed = OrderedDict()
        self.cuebam = OrderedDict()
        
        if self.annotation:
            #all_tags = []
            #for dictt in self.exps.fieldsDict.values():
            #    for tag in dictt.keys():
            #        all_tags.append(tag) 
            for bed in self.bednames:
            #    self.cuebed[bed] = set([bed]+all_tags)
                self.cuebed[bed] = set([bed])
        else:
            for bed in self.bednames:
                self.cuebed[bed] = set(tag_from_r(self.exps, self.tag_type,bed))
        for bam in self.readsnames:
            self.cuebam[bam] = set(tag_from_r(self.exps, self.tag_type,bam))
        
        
    def coverage(self, sortby, heatmap=False, logt=False, mp=False):
        
        def annot_ind(bednames, tags):
            """Find the index for annotation tag"""
            for ind, a in enumerate(bednames):
                if a in tags: return ind

        if mp: ts = time.time()
        # Calculate for coverage
        mp_input = []
        data = OrderedDict()
        totn = len(self.sort_tags) * len(self.group_tags) * len(self.color_tags)
        if self.df: totn = totn * 2
        bi = 0
        for s in self.sort_tags:
            data[s] = OrderedDict()
            for g in self.group_tags:
                data[s][g] = OrderedDict()
                for c in self.color_tags:
                    #if self.df: data[s][g][c] = []
                    for bed in self.cuebed.keys():
                        if self.cuebed[bed] <= set([s,g,c]):
                            for bam in self.cuebam.keys():
                                if self.cuebam[bam] <= set([s,g,c]):
                                    if mp: 
                                        if self.annotation:
                                            i = annot_ind(self.bednames, [s,g,c])
                                        else: i = self.bednames.index(bed)
                                        j = self.readsnames.index(bam)
                                        mp_input.append([ self.processed_beds[i], self.reads[j], 
                                                          self.rs, self.bs, self.ss, self.center, heatmap, logt,
                                                          s, g, c])
                                        if self.df: data[s][g][c] = []
                                        else: data[s][g][c] = 0
                                    else:
                                        
                                        ts = time.time()
                                        if self.annotation:
                                            i = annot_ind(self.bednames, [s,g,c])
                                        else: i = self.bednames.index(bed)
                                        j = self.readsnames.index(bam)
                                        
                                        cov = CoverageSet(bed+"."+bam, self.processed_beds[i])
                                        cov.coverage_from_bam(bam_file=self.reads[j], read_size = self.rs, binsize = self.bs, stepsize = self.ss)
                                        cov.normRPM()
                                        # When bothends, consider the fliping end
                                        if self.center == 'bothends':
                                            flap = CoverageSet("for flap", self.processed_bedsF[i])
                                            flap.coverage_from_bam(self.reads[j], read_size = self.rs, binsize = self.bs, stepsize = self.ss)
                                            ffcoverage = numpy.fliplr(flap.coverage)
                                            cov.coverage = numpy.concatenate((cov.coverage, ffcoverage), axis=0)
                                        # Averaging the coverage of all regions of each bed file
                                        if heatmap:
                                            if logt:
                                                data[s][g][c] = numpy.log10(numpy.vstack(cov.coverage)) # Store the array into data list
                                            else:
                                                data[s][g][c] = numpy.vstack(cov.coverage) # Store the array into data list
                                        else:
                                            #print(cov.coverage)
                                            for i, car in enumerate(cov.coverage):
                                                car = numpy.delete(car, [0,1])
                                                if i == 0:
                                                    avearr = np.array(car)
                                                    lenr = car.shape[0]
                                                elif car.shape[0] == lenr:
                                                    avearr = numpy.vstack((avearr, car))
                                                else:
                                                    pass
                                            
                                            avearr = numpy.average(avearr, axis=0)
                                            
                                            if self.df: 
                                                try: data[s][g][c].append(avearr)
                                                except: data[s][g][c] = [avearr]
                                            else: data[s][g][c] = avearr # Store the array into data list
                                        bi += 1
                                        te = time.time()
                                        print2(self.parameter, "     ("+str(bi)+"/"+str(totn)+") Computing\t" + "{0:30}   --{1:<6.1f}secs".format(bed+"."+bam, ts-te))
        if mp: 
            pool = multiprocessing.Pool(processes=multiprocessing.cpu_count())
            mp_output = pool.map(compute_coverage, mp_input)
            pool.close()
            pool.join()
        
            for s in data.keys():
                for g in data[s].keys():
                    for c in data[s][g].keys():
                        for out in mp_output:
                            if out[0] == s and out[1] == g and out[2] == c:
                                if self.df:
                                    data[s][g][c].append(out[3])
                                else:
                                    data[s][g][c] = out[3]     
            te = time.time()
            
        if self.df:
            for s in data.keys():
                for g in data[s].keys():
                    for c in data[s][g].keys():
                        #print(s+"  "+ g+"  "+ c)
                        #print(len(data[s][g][c]))
                        try: 
                            d = numpy.subtract(data[s][g][c][0],data[s][g][c][1])
                            data[s][g][c] = d
                        except: pass

        self.data = data
        
    def colormap(self, colorby, definedinEM):
        self.colors = colormap(self.exps, colorby, definedinEM, annotation=self.annotation)
        
    def plot(self, groupby, colorby, output, printtable=False, sy=False, sx=False):
        
        rot = 50
        ticklabelsize = 7
        
        f, axs = plt.subplots(len(self.data.keys()),len(self.data.values()[0].keys()), figsize=(8.27, 11.69), dpi=300) # 
        #if len(self.data.keys()) == 1 and len(self.data.values()[0]) == 1: 
        #    axs=numpy.array([[axs,None],[None,None]])
        
        yaxmax = [0]*len(self.data.values()[0])
        sx_ymax = [0]*len(self.data.keys())
        if self.df: 
            yaxmin = [0]*len(self.data.values()[0])
            sx_ymin = [0]*len(self.data.keys())

        for it, s in enumerate(self.data.keys()):
            for i,g in enumerate(self.data[s].keys()):
                
                try: ax = axs[it,i]
                except: 
                    if len(self.data.keys()) == 1 and len(self.data[s].keys()) == 1:
                        ax = axs
                    elif len(self.data.keys()) == 1 and len(self.data[s].keys()) > 1:
                        ax = axs[i]
                    else:
                        ax = axs[it]
                              
                if it == 0:
                    if self.df:
                        ax.set_title(g+"_df",fontsize=11)
                    else:
                        ax.set_title(g,fontsize=11)
                    
                # Processing for future output
                if printtable:
                    pArr = numpy.array(["Name","X","Y"]) # Header
                    
                for j, c in enumerate(self.data[s][g].keys()):
                    y = self.data[s][g][c]
                    
                    try: 
                        yaxmax[i] = max(numpy.amax(y), yaxmax[i])
                        sx_ymax[i] = max(numpy.amax(y), sx_ymax[i])
                        if self.df: 
                            yaxmin[i] = min(numpy.amin(y), yaxmin[i])
                            sx_ymin[i] = min(numpy.amin(y), sx_ymin[i])

                    except: continue
                    x = numpy.linspace(-self.extend, self.extend, len(y))
                    ax.plot(x,y, color=self.colors[j], lw=1)
                    # Processing for future output
                    if printtable:
                        [bed] = [bed for bed in self.bednames if [g,c,s] in self.cuebed[bed]]
                        name = numpy.array(*len(x))
                        xvalue = numpy.array(x)
                        yvalue = numpy.array(y)
                        conArr = numpy.vstack([name,xvalue,yvalue])
                        conArr = numpy.transpose(conArr)
                        pArr = numpy.vstack([pArr, conArr])
                if printtable:
                    [bam] = [bam for bam in self.readsnames if [g,c,s] in self.cuebam[bam]]
                    output_array(pArr, directory = output, folder ="lineplot_tables",filename=s+"_"+bam)
                
                ax.set_xlim([-self.extend, self.extend])
                plt.setp(ax.get_xticklabels(), fontsize=ticklabelsize, rotation=rot)
                plt.setp(ax.get_yticklabels(), fontsize=ticklabelsize)
                ax.locator_params(axis = 'x', nbins = 4)
                ax.locator_params(axis = 'y', nbins = 4)
                try: axs[0,-1].legend(self.color_tags, loc='center left', handlelength=1, handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
                except: axs[-1].legend(self.color_tags, loc='center left', handlelength=1, handletextpad=1, columnspacing=2, borderaxespad=0., prop={'size':10}, bbox_to_anchor=(1.05, 0.5))
                
        for it,ty in enumerate(self.data.keys()):
            try: 
                axs[it,0].set_ylabel("{}".format(ty),fontsize=12)
            except:
                if len(self.data.keys()) == 1:
                    axs[0].set_ylabel("{}".format(ty),fontsize=12)
                else:
                    axs[it].set_ylabel("{}".format(ty),fontsize=12)
                    
            if sy:
                for i,g in enumerate(self.data[ty].keys()):
                    if self.df:
                        try: axs[it,i].set_ylim([yaxmin[i] - abs(yaxmin[i]*0.2), yaxmax[i] + abs(yaxmax[i]*0.2)])
                        except:
                            if len(self.data.keys()) == 1:
                                axs[i].set_ylim([yaxmin[i] - abs(yaxmin[i]*0.2), yaxmax[i] + abs(yaxmax[i]*0.2)])
                            else:
                                axs[it].set_ylim([yaxmin[i] - abs(yaxmin[i]*0.2), yaxmax[i] + abs(yaxmax[i]*0.2)])
                    else:
                        try: axs[it,i].set_ylim([0, yaxmax[i]*1.2])
                        except:
                            if len(self.data.keys()) == 1:
                                axs[i].set_ylim([0, yaxmax[i]*1.2])
                            else:
                                axs[it].set_ylim([0, yaxmax[i]*1.2])
            elif sx:
                for i,g in enumerate(self.data[ty].keys()):
                    if self.df:
                        try: axs[it,i].set_ylim([sx_ymin[it] - abs(sx_ymin[it]*0.2), sx_ymax[it] + abs(sx_ymax[it]*0.2)])
                        except:
                            if len(self.data.keys()) == 1:
                                axs[i].set_ylim([sx_ymin[it] - abs(sx_ymin[it]*0.2), sx_ymax[it] + abs(sx_ymax[it]*0.2)])
                            else:
                                axs[it].set_ylim([sx_ymin[it] - abs(sx_ymin[it]*0.2), sx_ymax[it] + abs(sx_ymax[it]*0.2)])
                    else:
                        try: axs[it,i].set_ylim([0, sx_ymax[it]*1.2])
                        except:
                            if len(self.data.keys()) == 1:
                                axs[i].set_ylim([0, sx_ymax[it]*1.2])
                            else:
                                axs[it].set_ylim([0, sx_ymax[it]*1.2])

                
        f.tight_layout(pad=1.08, h_pad=None, w_pad=None)
        self.fig = f

    def gen_html(self, directory, title, align=50):
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = title
        link_d = OrderedDict()
        link_d["Lineplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"


        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        html.add_figure("lineplot.png", align="center")
        
        
        html.write(os.path.join(directory, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        type_list = 'ssssssssss'
        col_size_list = [20,20,20,20,20,20,20,20,20]
        header_list=["Assumptions and hypothesis"]
        data_table = [[]]
        if self.annotation:
            data_table.append("Genomic annotation: TSS - Transcription Start Site; TTS - Transcription Termination Site.")
            
        html.add_zebra_table(header_list, col_size_list, type_list, data_table, align = align, 
                             cell_align="left")
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])

        html.write(os.path.join(directory, title, "parameters.html"))

        
    def hmsort(self,sort):
        if sort == None:
            pass
        elif sort == 0:
            for t in self.data.keys():
                for i, g in enumerate(self.data[t].keys()):
                    #print(numpy.sum(data[t][bed].values()[0], axis=1))
                    #print(len(numpy.sum(data[t][bed].values()[0], axis=1)))
                    
                    sumarr = numpy.sum([numpy.sum(d, axis=1) for d in self.data[t][g].values()], axis=0)
                    #print(sumarr)
                    #sumarr = numpy.sum(sumarr, axis=1)
                    ind = stats.rankdata(sumarr,method='ordinal') # The index for further sorting
                    #numpy.fliplr(ind)
                    
                    for j, c in enumerate(self.data[t][g].keys()):
                        d = numpy.empty(shape=(self.data[t][g][c].shape))
                        for k, ranki in enumerate(ind):
                            d[-ranki,:] = self.data[t][g][c][k,:]
                        self.data[t][g][c] = d
        else:
            for t in self.data.keys():
                for i, g in enumerate(self.data[t].keys()):
                    sumarr = numpy.sum(self.data[t][g].values()[sort - 1], axis=1)
                    #print(sumarr)
                    #sumarr = numpy.sum(sumarr, axis=1)
                    ind = stats.rankdata(sumarr,method='ordinal') # The index for further sorting
                    #list(ind)
                    #print(ind)
                    for j, c in enumerate(self.data[t][g].keys()):
                        d = numpy.empty(shape=(self.data[t][g][c].shape))
                        for k, ranki in enumerate(ind):
                            d[-ranki,:] = self.data[t][g][c][k,:]
                        self.data[t][g][c] = d
                    #print(data[t][bed].values()[0])
 
    def hmcmlist(self, colorby, definedinEM):
        self.colors = colormaps(self.exps, colorby, definedinEM)
    
    def heatmap(self, logt):
        tickfontsize = 6
        ratio = 10
        self.hmfiles = []
        self.figs = []
        for ti, t in enumerate(self.data.keys()):
            #fig.append(plt.figure())
            #rows = len(data[t].keys())
            columns = len(self.data[t].values()[0].keys())
            #fig, axs = plt.subplots(rows,columns, sharey=True, dpi=300)
            #matplotlib.pyplot.subplots_adjust(left=1, right=2, top=2, bottom=1)
            fig = plt.figure(t)
            plt.suptitle("Heatmap: "+t, y=1.05)
            rows = len(self.data[t].keys())
            
            #gs = gridspec.GridSpec(rows*ratio,columns)
            axs = numpy.empty(shape=(rows+1,columns), dtype=object)
    
            for bi, g in enumerate(self.data[t].keys()):
                for bj, c in enumerate(self.data[t][g].keys()):
                    max_value = numpy.amax(self.data[t][g][c])
                    max_value = int(max_value)
                    axs[bi, bj] = plt.subplot2grid(shape=(rows*ratio+1, columns), loc=(bi*ratio, bj), rowspan=ratio)
                    if bi == 0: axs[bi, bj].set_title(c, fontsize=7)
                    #print(self.data[t][g][c])
                    im = axs[bi, bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto', 
                                            vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])
                    
            for bi, g in enumerate(self.data[t].keys()):
                for bj, c in enumerate(self.data[t][g].keys()):
                    
                    
                    #im = axs[bi, bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto', 
                    #                        vmin=0, vmax=max_value, interpolation='nearest', cmap=cm.coolwarm)
                    axs[bi, bj].set_xlim([-self.extend, self.extend])
                    axs[bi, bj].set_xticks([-self.extend, 0, self.extend])
                    #axs[bi, bj].set_xticklabels([-args.e, 0, args.e]
                    plt.setp(axs[bi, bj].get_xticklabels(), fontsize=tickfontsize, rotation=0)
                    #plt.setp(axs[bi, bj].get_yticklabels(), fontsize=10)
                    #axs[bi, bj].locator_params(axis = 'x', nbins = 2)
                    #axs[bi, bj].locator_params(axis = 'y', nbins = 4)
                    for spine in ['top', 'right', 'left', 'bottom']:
                        axs[bi, bj].spines[spine].set_visible(False)
                    axs[bi, bj].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='on')
                    axs[bi, bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
                    
                    #if bj > 0:
                    #    plt.setp(axs[bi, bj].get_yticklabels(),visible=False)
                    #plt.setp(axarr[i].get_yticks(),visible=False)
                    axs[bi, bj].minorticks_off()
                    if bj == 0:
                        #nregion = len(self.exps.objectsDict[g])
                        #axs[bi, bj].set_ylabel(self.exps.get_type(g,'factor')+" ("+str(nregion) + ")", fontsize=7)
                        axs[bi, bj].set_ylabel(g, fontsize=7)
                    if bi == rows-1:
                        #divider = make_axes_locatable(axs[bi,bj])
                        #cax = divider.append_axes("bottom", size="5%", pad=0.5)
                        cbar_ax = plt.subplot2grid((rows*ratio+4, columns), (rows*ratio+3, bj))
                        #axs[rows,bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
                        
                        
                        #cbar = grid.cbar_axes[i//2].colorbar(im)     
                        #cbar = plt.colorbar(im, cax = axs[rows,bj], ticks=[0, max_value], orientation='horizontal')
                        #cbar = axs[rows,bj].imshow(range(int(max_value)), extent=[0, int(max_value),0,0], aspect=10, extent=[-self.extend, self.extend,0,0]
                        #                           vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])
                        #cbar = axs[rows,bj].imshow(self.data[t][g][c], extent=[-self.extend, self.extend, 0,1], aspect='auto', 
                        #                    vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj])
                        #cbar = axs[rows,bj].imshow([range(2*self.extend),range(2*self.extend),range(2*self.extend)], 
                        #                           aspect='auto', vmin=0, vmax=max_value, interpolation='nearest', cmap=self.colors[bj] )
                        #cbar.outline.set_linewidth(0.5)
                        #axs[rows,bj].set_ticks_position('none')
                        #axs[rows,bj].tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
                        #axs[rows,bj].tick_params(axis='y', which='both', left='off', right='off', labelleft='off')
                        if logt:
                            cbar.ax.set_xticklabels(['0', '{:1.1f}'.format(max_value)], fontsize=tickfontsize)# horizontal colorbar
                            cbar.set_label('log10', fontsize=tickfontsize)
                        else:
                            #cbar.ax.set_xticklabels(['0', int(max_value)], fontsize=tickfontsize)# horizontal colorbar
                            pass
                        #cbar.set_label('Amplitute of signal')
                        max_value = int(max_value)
                        #width = 0.4/rows
                        #cbar_ax = fig.add_axes([0.01 + bj/columns, 0, width, 0.01])
                        cbar = plt.colorbar(im, cax=cbar_ax, ticks=[0, max_value], orientation='horizontal')
                        cbar.ax.set_xticklabels([0, int(max_value)])
                        cbar.outline.set_linewidth(0.1)
                        
            fig.tight_layout()
            #fig.tight_layout(pad=1.08, h_pad=None, w_pad=None)
            #fig.tight_layout(pad=1, h_pad=1, w_pad=1)
            self.figs.append(fig)
            self.hmfiles.append("heatmap"+ "_" + t)

    def gen_htmlhm(self, outputname, title, align=50):
        dir_name = os.path.basename(directory)
        #check_dir(directory)
        html_header = title
        link_d = OrderedDict()
        link_d["Lineplot"] = "index.html"
        link_d["Parameters"] = "parameters.html"


        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        
        # Each row is a plot with its data
        for name in self.hmfiles:
            html.add_figure(name+".png", align="center")
        html.write(os.path.join(directory, title, "index.html"))

        ## Parameters
        html = Html(name=html_header, links_dict=link_d,
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")
        
        html.add_free_content(['<a href="parameters.txt" style="margin-left:100">See parameters</a>'])
        html.add_free_content(['<a href="experimental_matrix.txt" style="margin-left:100">See experimental matrix</a>'])
        html.write(os.path.join(directory, title, "parameters.html"))





