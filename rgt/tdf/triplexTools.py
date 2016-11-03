# Python Libraries
from __future__ import print_function
import os
import pysam
import shutil
import pickle
from ctypes import *
from collections import *

# Local Libraries
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator


# Distal Libraries
from rgt.SequenceSet import SequenceSet
from rgt.viz.plotTools import output_array
from rgt.GenomicRegion import GenomicRegion
from RNADNABindingSet import RNADNABindingSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.Util import SequenceType, Html, ConfigurationFile, GenomeData, Library_path

# Color code for all analysis
target_color = "mediumblue"
nontarget_color = "darkgrey"
sig_color = "powderblue"

####################################################################################
####################################################################################

def print2(summary, string):
    """ Show the message on the console and also save in summary. """
    print(string)
    summary.append(string)


def output_summary(summary, directory, filename):
    """Save the summary log file into the defined directory"""
    pd = os.path.join(dir, directory)
    try:
        os.stat(pd)
    except:
        os.mkdir(pd)
    if summary:
        with open(os.path.join(pd, filename), 'w') as f:
            print("********* RGT Triplex: Summary information *********",
                  file=f)
            for s in summary:
                print(s, file=f)


def check_dir(path):
    """Check the availability of the given directory and creat it"""
    try: os.stat(path)
    except:
        try: os.mkdir(path)
        except: pass


def try_int(s):
    "Convert to integer if possible."
    try:
        return int(s)
    except:
        return s


def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))


def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))


def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())


def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)


def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp


def list_all_index(path, link_d=None, show_RNA_ass_gene=False):
    """Creat an 'index.html' in the defined directory """

    dirname = os.path.basename(path)

    if link_d:
        pass
    else:
        link_d = {"List": "index.html"}

    html = Html(name="Directory: " + dirname, links_dict=link_d,
                fig_rpath="./style", fig_dir=os.path.join(path, "style"),
                RGT_header=False, other_logo="TDF", homepage="../index.html")

    html.add_heading("All experiments in: " + dirname + "/")

    data_table = []
    type_list = 'sssssssssssss'
    col_size_list = [20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20]
    c = 0
    if show_RNA_ass_gene:
        header_list = ["No.", "Experiments", "RNA", "Closest genes",
                       "No sig. DBD",
                       "Top DBD", "p-value", "Organism", "Target region"]
    else:
        header_list = ["No.", "Experiments", "RNA", "No sig. DBD",
                       "Top DBD", "p-value", "Organism",  # "Condition",
                       "Target region"]

    profile_f = open(os.path.join(path, "profile.txt"), 'r')
    profile = {}
    for line in profile_f:
        line = line.strip()
        line = line.split("\t")
        profile[line[0]] = line[1:]

    for i, exp in enumerate(profile.keys()):
        # print(exp)
        c += 1

        try:
            if profile[exp][5] == "-":
                new_line = [str(c), exp, profile[exp][0]]
            else:
                new_line = [str(c),
                            '<a href="' + os.path.join(exp, "index.html") + \
                            '">' + exp + "</a>", profile[exp][0]]

            if show_RNA_ass_gene:
                new_line.append(
                    split_gene_name(gene_name=profile[exp][7],
                                    org=profile[exp][2])
                )

            if profile[exp][6] == "-":
                new_line += [profile[exp][4],
                             profile[exp][5], profile[exp][6],
                             profile[exp][2], profile[exp][3]]

            elif float(profile[exp][6]) < 0.05:
                new_line += [profile[exp][4],
                             profile[exp][5],
                             "<font color=\"red\">" + \
                             profile[exp][6] + "</font>",
                             profile[exp][2], profile[exp][3]]
            else:
                new_line += [profile[exp][4],
                             profile[exp][5], profile[exp][6],
                             profile[exp][2], profile[exp][3]]
            data_table.append(new_line)
        except:
            if exp != "Experiment":
                print("Error in loading profile: " + exp)
            continue

    html.add_zebra_table(header_list, col_size_list, type_list, data_table,
                         align=10, cell_align="left", sortable=True)

    html.add_fixed_rank_sortable()
    html.write(os.path.join(path, "index.html"))


def revise_index(root, show_RNA_ass_gene=False):
    "Revise other index.html in the same project"

    dir_list = {}
    plist = {}
    for item in os.listdir(root):
        h = os.path.join(root, item, "index.html")
        pro = os.path.join(root, item, "profile.txt")
        if os.path.isfile(pro):
            dir_list[os.path.basename(item)] = "../" + item + "/index.html"
            plist[os.path.basename(item)] = h
    dir_list = OrderedDict(sorted(dir_list.items()))
    # print(dir_list)
    for d, p in plist.iteritems():
        list_all_index(path=os.path.dirname(p),
                       link_d=dir_list, show_RNA_ass_gene=show_RNA_ass_gene)


def gen_heatmap(path):
    """Generate the heatmap to show the sig RNA in among conditions"""

    def fmt(x, pos):
        # a, b = '{:.0e}'.format(x).split('e')
        # b = int(b)
        # return r'${} \times 10^{{{}}}$'.format(a, b)
        a = -numpy.log10(x)
        return '{:.0f}'.format(a)

    matrix = OrderedDict()
    rnas = []
    for item in os.listdir(path):
        # print(item)
        # print(os.path.isdir(os.path.join(path,item)))
        if not os.path.isdir(os.path.join(path, item)): continue
        if item == "style": continue
        # if item == "index.html": continue
        matrix[item] = {}
        pro = os.path.join(path, item, "profile.txt")
        with open(pro) as f:
            for line in f:
                line = line.strip().split("\t")
                if line[0] == "Experiment": continue
                if line[6] == "-": continue
                matrix[item][line[0]] = float(line[7])
                rnas.append(line[0])
    rnas = list(set(rnas))
    # rnas.sort()

    # Convert into array
    ar = []
    exps = natsorted(matrix.keys())
    rnas = natsorted(rnas)
    # print(exps)
    for exp in exps:
        row = []
        for rna in rnas:
            try:
                row.append(matrix[exp][rna])
            except:
                row.append(1)
        ar.append(row)
    ar = numpy.array(ar)
    ar = numpy.transpose(ar)

    # print(ar.shape)
    data = ar[~numpy.all(ar == 1, axis=1)]
    # print(data.shape)


    fig = plt.figure(figsize=(len(matrix.keys()) * 1.5, len(rnas) * 2.5))
    # fig = plt.figure()
    # ax1 = fig.add_axes([0.09,0.2,0.2,0.6])
    # Y = sch.linkage(data, method='single')
    # Z1 = sch.dendrogram(Y, orientation='right')
    # ax1.set_xticks([])
    # ax1.set_yticks([])

    # Compute and plot second dendrogram.
    # ax2 = fig.add_axes([0.3, 1.1, 0.55,0.04])
    # Y = sch.linkage(data.T, method='single')
    # Z2 = sch.dendrogram(Y)
    # ax2.set_xticks([])
    # ax2.set_yticks([])
    bounds = []
    for n in range(-8, 0):
        # print(10**n)
        bounds.append(10 ** n)
    norm = colors.BoundaryNorm(bounds, plt.cm.YlOrRd_r.N)

    # Plot distance matrix.
    axmatrix = fig.add_axes([0.1, 0.2, 0.8, 0.6])
    # axmatrix = fig.add_axes()
    # idx1 = Z1['leaves']
    # idx2 = Z2['leaves']
    # data = data[idx1,:]
    # data = data[:,idx2]
    im = axmatrix.matshow(data, aspect='auto', origin='lower', cmap=plt.cm.YlOrRd_r, norm=norm)

    axmatrix.set_xticks(range(data.shape[1]))
    axmatrix.set_xticklabels(exps, minor=False, ha="left")
    axmatrix.xaxis.set_label_position('top')
    axmatrix.xaxis.tick_top()
    plt.xticks(rotation=40, fontsize=10)

    axmatrix.set_yticks(range(data.shape[0]))
    axmatrix.set_yticklabels(rnas, minor=False)
    # axmatrix.set_yticklabels( [ rnas[i] for i in idx1 ], minor=False)
    axmatrix.yaxis.set_label_position('right')
    axmatrix.yaxis.tick_right()
    plt.yticks(rotation=0, fontsize=10)
    # axmatrix.tight_layout()

    # Plot colorbar.
    axcolor = fig.add_axes([0.1, 0.1, 0.8, 0.02])

    plt.colorbar(im, cax=axcolor, orientation='horizontal', norm=norm,
                 boundaries=bounds, ticks=bounds, format=matplotlib.ticker.FuncFormatter(fmt))
    axcolor.set_xlabel('p value (-log10)')
    lmats = ar.tolist()
    for i, r in enumerate(rnas):
        lmats[i] = [r] + [str(x) for x in lmats[i]]
    lmats = [["p-value"] + exps] + lmats

    output_array(array=lmats, directory=path, folder="", filename='matrix_p.txt')
    # os.path.join(path,'matrix_p.txt'), lmats, delimiter='\t')
    fig.savefig(os.path.join(path, 'condition_lncRNA_dendrogram.png'))
    fig.savefig(os.path.join(path, 'condition_lncRNA_dendrogram.pdf'), format="pdf")


def generate_rna_exp_pv_table(root, multi_corr=True):
    "Generate p value table for Experiments vs RNA in the same project"

    nested_dict = lambda: defaultdict(nested_dict)
    # nested_dict = lambda: defaultdict(lambda: 'n.a.')

    data = nested_dict()
    rnas = []

    for item in os.listdir(root):
        pro = os.path.join(root, item, "profile.txt")
        if os.path.isfile(pro):
            with open(pro) as f:
                for line in f:
                    if line.startswith("Experiment"):
                        continue
                    else:
                        line = line.strip().split("\t")
                        data[item][line[0]] = float(line[7])
                        rnas.append(line[0])

    exp_list = sorted(data.keys())
    rnas = sorted(list(set(rnas)))

    pvs = []
    for rna in rnas:
        for exp in exp_list:
            if data[exp][rna]: pvs.append(data[exp][rna])
    if multi_corr:
        reject, pvals_corrected = multiple_test_correction(pvs, alpha=0.05, method='indep')
    else:
        pvals_corrected = pvs

    with open(os.path.join(root, "table_exp_rna_pv.txt"), "w") as t:
        print("\t".join(["RNA_ID"] + exp_list), file=t)
        i = 0
        for rna in rnas:
            newline = [rna]
            for exp in exp_list:
                if data[exp][rna]:
                    newline.append(str(pvals_corrected[i]))
                    i += 1
                else:
                    newline.append("n.a.")
            print("\t".join(newline), file=t)


def value2str(value):
    if (isinstance(value,str)): return value
    if value == 0: return "0"
    if(isinstance(value,int)): return str(value)
    elif(isinstance(value,float)):
        if abs(value) >= 1000: 
            try: r = "{}".format(int(value))
            except: r = "Inf"
        elif 1000 > abs(value) > 10: r = "{:.1f}".format(value)
        elif 10 > abs(value) >= 1: r = "{:.2f}".format(value)
        elif 1 > abs(value) >= 0.05: r = "{:.2f}".format(value)
        elif 0.05 > abs(value) > 0.0001: r = "{:.4f}".format(value)
        else: r = "{:.1e}".format(value)
        return r


def uniq(seq):
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if not (x in seen or seen_add(x))]


def random_each(input):
    """Return the counts of DNA Binding sites with randomization
    For multiprocessing. 
    Input contains:
    0       1               2                3     4              5          6                    
    str(i), self.rna_fasta, self.dna_region, temp, self.organism, self.rbss, str(marks.count(i)),
    number, rna,            region,          temp, organism,      rbss,      number of mark

    7  8  9  10  11  12  13  14  15          16                17   18
    l, e, c, fr, fm, of, mf, rm, filter_bed, self.genome_path, tp   par
    """
    import sys
    # Filter BED file
    if input[15]:
        random = input[2].random_regions(organism=input[4], multiply_factor=1,
                                         overlap_result=True, overlap_input=True,
                                         chrom_X=True, chrom_M=False, filter_path=input[15])
    else:
        random = input[2].random_regions(organism=input[4], multiply_factor=1,
                                         overlap_result=True, overlap_input=True,
                                         chrom_X=True, chrom_M=False)
    
    txp = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3],
                       organism=input[4], prefix=str(input[0]), remove_temp=True,
                       l=int(input[7]), e=int(input[8]), c=input[9], fr=input[10],
                       fm=input[11], of=input[12], mf=input[13], rm=input[14],
                       par=input[18], genome_path=input[16],
                       dna_fine_posi=False, tp=input[17])

    txp.merge_rbs(rbss=input[5], rm_duplicate=True)

    txpf = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[4], prefix=str(input[0]), remove_temp=True, 
                       l=int(input[7]), e=int(input[8]),  c=input[9], fr=input[10], 
                       fm=input[11], of=input[12], mf=input[13], rm=input[14], 
                       par=input[18], genome_path=input[16],
                       dna_fine_posi=True, tp=input[17])

    txpf.merge_rbs(rbss=input[5], rm_duplicate=True)
    sys.stdout.flush()
    print("".join(["="]*int(input[6])), end="")

    return [ [len(tr) for tr in txp.merged_dict.values() ], [len(dbss) for dbss in txpf.merged_dict.values()] ]


def get_sequence(dir, filename, regions, genome_path):
    """
    Fetch sequence into FASTA file according to the given BED file
    """
    genome = pysam.Fastafile(genome_path)
    with open(os.path.join(dir, filename), 'w') as output:
        for region in regions:
            print(">"+ region.toString(), file=output)
            # if region.orientation == "-":
            #     seq = Seq(genome.fetch(region.chrom, max(0, region.initial), region.final), IUPAC.unambiguous_dna)
            #     seq = seq.reverse_complement()
            #     print(seq, file=output)
            # else:
            print(genome.fetch(region.chrom, max(0, region.initial), region.final), file=output)


def find_triplex(rna_fasta, dna_region, temp, organism, l, e, dna_fine_posi, genome_path, prefix="", remove_temp=False, 
                 c=None, fr=None, fm=None, of=None, mf=None, rm=None, par="", tp=False):
    """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
    
    # Generate FASTA 
    get_sequence(dir=temp, filename="dna_"+prefix+".fa", regions=dna_region, genome_path=genome_path)

    # Triplexator
    run_triplexator(ss=rna_fasta, ds=os.path.join(temp,"dna_"+prefix+".fa"), 
                    output=os.path.join(temp, "dna_"+prefix+".txp"), 
                    l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, rm=rm, par=par, tp=tp)
    # Read txp
    txp = RNADNABindingSet("dna")
    txp.read_txp(os.path.join(temp, "dna_"+prefix+".txp"), dna_fine_posi=dna_fine_posi)
    txp.remove_duplicates()

    if remove_temp:
        os.remove(os.path.join(temp,"dna_"+prefix+".fa"))
        os.remove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp

def run_triplexator(ss, ds, output, l=None, e=None, c=None, fr=None, fm=None, of=None, mf=None, rm=None, par="", tp=False):
    """Perform Triplexator"""
    #triplexator_path = check_triplexator_path()
    # triplexator -ss -ds -l 15 -e 20 -c 2 -fr off -fm 0 -of 1 -rm
    triclass = Library_path()
    triplex_lib_path = triclass.get_triplexator()
    triplex_lib  = cdll.LoadLibrary(triplex_lib_path)

    arguments = ""
    if ss: arguments += "-ss "+ss+" "
    if ds: arguments += "-ds "+ds+" "
    if l: arguments += "-l "+str(l)+" "
    if e: arguments += "-e "+str(e)+" "
    if c: arguments += "-c "+str(c)+" "
    if fr: arguments += "-fr "+fr+" "
    if fm: arguments += "-fm "+str(fm)+" "
    if of: arguments += "-of "+str(of)+" "
    if mf: arguments += "-mf "
    if rm: arguments += "-rm "+str(rm)+" "
    # arguments += "--bit-parallel "
    if par != "":
        par = par.replace('_'," ")
        par = "-" + par
        arguments += par+" "
    
    arguments += "-o "+ os.path.basename(output) + " -od " + os.path.dirname(output)

    arg_strings  = arguments.split(' ')
    arg_ptr      = (c_char_p * (len(arg_strings) + 1))()

    arg_ptr[0] = "triplexator"  # to simulate calling from cmd line
    for i, s in enumerate(arg_strings):
        arg_ptr[i + 1] = s
    
    triplex_lib.pyTriplexator(len(arg_strings) + 1, arg_ptr)
    os.remove(os.path.join(output + ".summary"))
    os.remove(os.path.join(output + ".log"))


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
                if not line: continue
                elif line[0][0] == "#": continue
                elif len(line) < 2: continue
                else:
                    v = line[1]
                    if v == "NA": 
                        access.append(0)
                    else: 
                        try: v = 2**(-float(v))
                        except: continue
                        if v >= cut_off:
                            access.append(1)
                        else:
                            access.append(0)
                    i += 1
    return access

def split_gene_name(gene_name, org):
    
    if gene_name == None:
        return ""
    if gene_name[0:2] == "chr":
        return gene_name

    if org=="hg19": ani = "human"
    elif org=="hg38": ani = "human"
    elif org=="mm9": ani = "mouse"
    else: ani = None

    if not ani:
        if org == "tair10":
            p1 = "".join(['<a href="https://www.arabidopsis.org/servlets/TairObject?name=', gene_name,
                          '&type=locus" target="_blank">', gene_name, '</a>' ])
            return p1
        else:
            return gene_name
    else:    
        p1 = '<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?org='+ani+\
             "&db="+org+"&singleSearch=knownCanonical&position="
        p2 = '" style="text-align:left" target="_blank">'
        p3 = '</a>'
        
        #print(gene_name)
        if ":" in gene_name:
            genes = gene_name.split(":")
            genes = list(set(genes))
            dlist = []
            glist = []
            for i, g in enumerate(genes):
                if "(" in g:
                    d = g.partition('(')[2].partition(')')[0]
                    dlist.append((int(d)))
                    g = g.partition('(')[0]
                    glist.append(g)

                else:
                    if i == 0:
                        result = p1+g+p2+g+p3
                    else:
                        result += ","+p1+g+p2+g+p3
            if dlist:
                i = dlist.index(min(dlist))
                if dlist[i] < 0: 
                    ds = str(dlist[i])
                    result = p1+glist[i]+p2+glist[i]+p3+"("+ds+")"
                    try:
                        j = dlist.index(max(dlist))
                        if dlist[j] < 0: rds = str(dlist[j])
                        else:
                            rds = "+"+ str(dlist[j])
                            result += ","+p1+glist[j]+p2+glist[j]+p3+"("+rds+")"
                    except:
                        pass
                else: 
                    ds = "+"+str(dlist[i])
                    result = p1+glist[i]+p2+glist[i]+p3+"("+str(dlist[i])+")"
                
                    
        elif gene_name == ".":
            result = "none"

        else:
            if "(" in gene_name:
                d = gene_name.partition('(')[2].partition(')')[0]
                g = gene_name.partition('(')[0]
                result = p1+g+p2+g+p3+"("+d+")"
            else:
                result = p1+gene_name+p2+gene_name+p3
        
        return result


def lineplot(txp, rnalen, rnaname, dirp, sig_region, cut_off, log, ylabel, linelabel, 
             filename, ac=None, showpa=False, exons=None):
    # Plotting
    f, ax = plt.subplots(1, 1, dpi=300, figsize=(6,4))
    
    # Extract data points
    x = range(rnalen)
    #print(rnalen)
    if log:
        all_y = [1] * rnalen
        p_y = [1] * rnalen
        a_y = [1] * rnalen
    else:
        all_y = [0] * rnalen
        p_y = [0] * rnalen
        a_y = [0] * rnalen

    txp.remove_duplicates_by_dbs()
    for rd in txp:
        #print(str(rd.rna.initial), str(rd.rna.final))
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
        max_y = max(all_y)+0.5
        min_y = 1
        ylabel += "(log10)"
    else:
        max_y = float(max(all_y) * 1.1)
        min_y = 0

    if ac:
        min_y = float(max_y*(-0.09))
    
    
    # Plotting
    for rbs in sig_region:
        rect = patches.Rectangle(xy=(rbs.initial,0), width=len(rbs), height=max_y, facecolor=sig_color, 
                                 edgecolor="none", alpha=0.5, lw=None, label="Significant DBD")
        ax.add_patch(rect)
    
    lw = 1.5
    if showpa:
        ax.plot(x, all_y, color=target_color, alpha=1, lw=lw, label="Parallel + Anti-parallel")
        ax.plot(x, p_y, color="purple", alpha=1, lw=lw, label="Parallel")
        ax.plot(x, a_y, color="dimgrey", alpha=.8, lw=lw, label="Anti-parallel")
    else:
        ax.plot(x, all_y, color="mediumblue", alpha=1, lw=lw, label=linelabel)

    # RNA accessbility
    if ac:
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
                pac = ax.add_patch(patches.Rectangle((last_i, min_y), i-last_i, -min_y,
                                   hatch='///', fill=False, snap=False, linewidth=0, label="RNA accessibility"))
                drawing = False
            else:
                continue

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    legend_h = []
    legend_l = []
    for uniqlabel in uniq(labels):
        legend_h.append(handles[labels.index(uniqlabel)])
        legend_l.append(uniqlabel)
    ax.legend(legend_h, legend_l, 
              bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
              prop={'size':9}, ncol=3)

    # XY axis
    ax.set_xlim(left=0, right=rnalen )
    ax.set_ylim( [min_y, max_y] ) 
    for tick in ax.xaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    for tick in ax.yaxis.get_major_ticks(): tick.label.set_fontsize(9) 
    ax.set_xlabel(rnaname+" sequence (bp)", fontsize=9)
    
    ax.set_ylabel(ylabel,fontsize=9, rotation=90)
    
    if None:
        if exons and len(exons) > 1:
            w = 0
            i = 0
            h = (max_y - min_y)*0.02

            for exon in exons:
                l = abs(exon[2] - exon[1])
                
                #print([i,l,w])
                #ax.axvline(x=w, color="gray", alpha=0.5, zorder=100)
                if i % 2 == 0:
                    rect = matplotlib.patches.Rectangle((w,max_y-h),l,h, color="moccasin")
                else:
                    rect = matplotlib.patches.Rectangle((w,max_y-h),l,h, color="gold")
                ax.add_patch(rect)
                i += 1
                w += l
            ax.text(rnalen*0.01, max_y-2*h, "exon boundaries", fontsize=5, color='black')

    f.tight_layout(pad=1.08, h_pad=None, w_pad=None)

    f.savefig(os.path.join(dirp, filename), facecolor='w', edgecolor='w',  
              bbox_extra_artists=(plt.gci()), bbox_inches='tight', dpi=300)
    # PDF
    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(12) 
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(12)
    ax.xaxis.label.set_size(14)
    ax.yaxis.label.set_size(14) 
    ax.legend(legend_h, legend_l, 
              bbox_to_anchor=(0., 1.02, 1., .102), loc=2, mode="expand", borderaxespad=0., 
              prop={'size':12}, ncol=3)
    pp = PdfPages(os.path.splitext(os.path.join(dirp,filename))[0] +'.pdf')
    pp.savefig(f,  bbox_inches='tight') # bbox_extra_artists=(plt.gci()),
    pp.close()

def load_dump(path, filename):
    file = open(os.path.join(path,filename),'r')
    object = pickle.load(file)
    file.close()
    print("\tLoading from file: "+filename)
    return object

def dump(object, path, filename):
    file = open(os.path.join(path,filename),'wb')
    pickle.dump(object,file)
    file.close()
    print("\tDump to file: "+filename)

def check_triplexator_path():
    try:
        cf = ConfigurationFile()
        with open(os.path.join(cf.data_dir,"triplexator_path.txt")) as f:
            l = f.readline()
            l = l.strip("\n")
            return l
    except:
        print("Please define the path to Triplexator by command: rgt-TDF triplexator -path <PATH>")
        sys.exit(1)
    
def rna_associated_gene(rna_regions, name, organism):
    if rna_regions:
        s = [ rna_regions[0][0], min([e[1] for e in rna_regions]), 
              max([e[2] for e in rna_regions]), rna_regions[0][3] ]
        g = GenomicRegionSet("RNA associated genes")
        g.add( GenomicRegion(chrom=s[0], initial=s[1], final=s[2], name=name, orientation=s[3]) )
        asso_genes = g.gene_association(organism=organism, promoterLength=1000, threshDist=1000000, show_dis=True)
        #print(name)
        #print( [ a.name for a in asso_genes ] )
        #print( [ a.proximity for a in asso_genes ])
        genes = asso_genes[0].name.split(":")
        #proxs = asso_genes[0].proximity.split(":")
        closest_genes = []
        for n in genes:
            if name not in n: closest_genes.append(n)
        closest_genes = set(closest_genes)
        #print(closest_genes)
        #for n in proxs:
        #    if name not in n: closest_genes.append(n)
        if len(closest_genes) == 0:
            return "."
        else:
            return ":".join(closest_genes)
    else:
        return "."

def rank_array(a):
    a = numpy.array(a)
    sa = numpy.searchsorted(numpy.sort(a), a)
    return sa

def dbd_regions(exons, sig_region, rna_name, output,out_file=False, temp=None, fasta=True):
    """Generate the BED file of significant DBD regions and FASTA file of the sequences"""
    if len(sig_region) == 0:
        return
    #print(self.rna_regions)
    if not exons:
        return
    else:
        dbd = GenomicRegionSet("DBD")
        dbdmap = {}
        if len(exons) == 1:
            print("## Warning: No information of exons in the given RNA sequence, the DBD position may be problematic. ")
        for rbs in sig_region:
            loop = True

            if exons[0][3] == "-":
                while loop:
                    cf = 0
                    for exon in exons:
                        #print(exon)

                        l = abs(exon[2] - exon[1])
                        tail = cf + l
                        #print("cf:   " + str(cf))
                        #print("tail: " + str(tail) )

                        if cf <= rbs.initial <=  tail:
                            dbdstart = exon[2] - rbs.initial + cf
                            
                            if rbs.final <= tail: 
                                #print("1")
                                dbdend = exon[2] - rbs.final + cf
                                if dbdstart > dbdend: dbdstart, dbdend = dbdend, dbdstart
                                dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                       initial=dbdstart, final=dbdend, 
                                                       orientation=exons[0][3], 
                                                       name=str(rbs.initial)+"-"+str(rbs.final) ) )
                                dbdmap[str(rbs)] = dbd[-1].toString() + " strand:-"
                                loop = False
                                break
                            elif rbs.final > tail:

                                subtract = l + cf - rbs.initial
                                #print("2")
                                #print("Subtract: "+str(subtract))
                                if dbdstart > exon[1]: dbdstart, exon[1] = exon[1], dbdstart
                                dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                       initial=dbdstart, final=exon[1], 
                                                       orientation=exons[0][3], 
                                                       name=str(rbs.initial)+"-"+str(rbs.initial+subtract)+"_split1" ) )
                        
                        elif rbs.initial < cf and rbs.final <= tail: 
                            #print("3")
                            dbdstart = exon[2]
                            dbdend = exon[2] - rbs.final + rbs.initial + subtract
                            if dbdstart > dbdend: dbdstart, dbdend = dbdend, dbdstart
                            dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                   initial=dbdstart, final=dbdend, 
                                                   orientation=exons[0][3], 
                                                   name=str(cf)+"-"+str(rbs.final)+"_split2" ) )
                            dbdmap[str(rbs)] = dbd[-2].toString() + " & " + dbd[-1].toString() + " strand:-"
                            loop = False
                            break

                        elif rbs.initial > tail:
                            pass

                        cf += l
                        
                    loop = False
            else:

                while loop:
                    cf = 0
                    for exon in exons:
                        #print(exon)
                        l = exon[2] - exon[1]
                        tail = cf + l
                        #print("cf:   " + str(cf))
                        #print("tail: " + str(tail) )
                        if cf <= rbs.initial <=  tail:
                            dbdstart = exon[1] + rbs.initial - cf
                            
                            if rbs.final <= tail: 
                                #print("1")
                                dbdend = exon[1] + rbs.final -cf
                                dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                       initial=dbdstart, final=dbdend, 
                                                       orientation=exons[0][3], 
                                                       name=str(rbs.initial)+"-"+str(rbs.final) ) )
                                dbdmap[str(rbs)] = dbd[-1].toString() + " strand:+"
                                loop = False
                                break
                            elif rbs.final > tail:

                                subtract = l + cf - rbs.initial
                                #print("2")
                                #print("Subtract: "+str(subtract))
                                dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                       initial=dbdstart, final=exon[2], 
                                                       orientation=exons[0][3], 
                                                       name=str(rbs.initial)+"-"+str(rbs.initial+subtract)+"_split1" ) )

                        elif rbs.initial < cf and rbs.final <= tail: 
                            #print("3")
                            dbdstart = exon[1]
                            dbdend = exon[1] + rbs.final - rbs.initial - subtract
                            dbd.add( GenomicRegion(chrom=exons[0][0], 
                                                   initial=dbdstart, final=dbdend, 
                                                   orientation=exons[0][3], 
                                                   name=str(cf)+"-"+str(rbs.final)+"_split2" ) )
                            dbdmap[str(rbs)] = dbd[-2].toString() + " & " + dbd[-1].toString() + " strand:+"
                            loop = False
                            break

                        elif rbs.initial > tail:
                            pass

                        cf += l
                        
                    loop = False
    if not out_file:        
        dbd.write_bed(filename=os.path.join(output, "DBD_"+rna_name+".bed"))
    else:
        # print(dbd)
        # print(dbd.sequences[0])
        dbd.write_bed(filename=output+".bed")
    # FASTA
    if fasta:
        #print(dbdmap)
        if not out_file:
            seq = pysam.Fastafile(os.path.join(output,"rna_temp.fa"))
            fasta_f = os.path.join(output, "DBD_"+rna_name+".fa")
        else:
            seq = pysam.Fastafile(os.path.join(temp,"rna_temp.fa"))
            fasta_f = output+".fa"

        with open(fasta_f, 'w') as fasta:
            
            for rbs in sig_region:
                try: info = dbdmap[str(rbs)]
                except: 
                    continue
                fasta.write(">"+ rna_name +":"+str(rbs.initial)+"-"+str(rbs.final)+ " "+ info +"\n")
                #print(seq.fetch(rbs.chrom, max(0, rbs.initial), rbs.final))
                if dbdmap[str(rbs)][-1] == "-":
                    fasta.write(seq.fetch(rbs.chrom, max(0, rbs.initial-1), rbs.final-1)+"\n" )
                else: 
                    fasta.write(seq.fetch(rbs.chrom, max(0, rbs.initial+1), rbs.final+1)+"\n" )

def connect_rna(rna, temp, rna_name):
    """Generate FASTA file merging all exons and return the number of exons and sequence length"""
    seq = ""
    exons = 0
    with open(rna) as f:
        for line in f:
            if line.startswith(">"):
                exons += 1
            else:
                line = line.strip()
                seq += line

    with open(os.path.join(temp,"rna_temp.fa"), "w") as r:
        print(">"+rna_name, file=r)
        for s in [seq[i:i + 80] for i in range(0, len(seq), 80)]:
            print(s, file=r)
    return [exons, len(seq)]

def get_dbss(input_BED,output_BED,rna_fasta,output_rbss,organism,l,e,c,fr,fm,of,mf,rm,temp,tp):
    regions = GenomicRegionSet("Target")
    regions.read_bed(input_BED)
    regions.gene_association(organism=organism)

    connect_rna(rna_fasta, temp=temp, rna_name="RNA")
    rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
    rnas.read_fasta(os.path.join(temp,"rna_temp.fa"))
    rna_regions = get_rna_region_str(os.path.join(temp,rna_fasta))
    # print(rna_regions)
    genome = GenomeData(organism)
    genome_path = genome.get_genome()
    txp = find_triplex(rna_fasta=rna_fasta, dna_region=regions, 
                       temp=temp, organism=organism, remove_temp=False, tp=tp,
                       l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, genome_path=genome_path,
                       prefix="targeted_region", dna_fine_posi=True)

    print("Total binding events:\t",str(len(txp)))
    txp.write_bed(output_BED)
    rbss = txp.get_rbs()
    dbd_regions(exons=rna_regions, sig_region=rbss, rna_name="rna", output=output_rbss, 
                out_file=True, temp=temp, fasta=False)
    # print(rbss.sequences)
    # print(len(rbss))
    # rbss.write_bed(output_rbss)

def get_rna_region_str(rna):
    """Getting the rna region from the information header with the pattern:
            REGION_chr3_51978050_51983935_-_
        or  chr3:51978050-51983935 -    """
    rna_regions = []
    with open(rna) as f:
        for line in f:
            if line[0] == ">":
                line = line.strip()
                if "REGION" in line:
                    line = line.split()
                    for i, e in enumerate(line):
                        if "REGION" in e:
                            e = e.split("_")
                            #print(e)
                            try:
                                rna_regions.append([e[1], int(e[2]), int(e[3]), e[4]])
                            except:
                                rna_regions.append([e[1], int(e[3]), int(e[4]), e[5]])
                
                elif "chr" in line:
                    line = line.partition("chr")[2]
                    chrom = "chr" + line.partition(":")[0]
                    start = int(line.partition(":")[2].partition("-")[0])
                    end = int(line.partition(":")[2].partition("-")[2].split()[0])
                    sign = line.partition(":")[2].partition("-")[2].split()[1]
                    if sign == "+" or sign == "-" or sign == ".":
                        rna_regions.append([chrom, start, end, sign])

                    else:
                        print(line)

                else:
                    rna_regions = None
                    break
    if rna_regions:
        rna_regions.sort(key=lambda x: x[1])
        if rna_regions[0][3] == "-":
            rna_regions = rna_regions[::-1]
    return rna_regions


def no_binding_response(args, rna_regions, rna_name, organism):
    print("*** Find no triple helices binding on the given RNA")

    pro_path = os.path.join(os.path.dirname(args.o), "profile.txt")
    exp = os.path.basename(args.o)
    if args.de:
        tar_reg = os.path.basename(args.de)
    else:
        tar_reg = os.path.basename(args.bed)
    r_genes = rna_associated_gene(rna_regions=rna_regions, name=rna_name, organism=organism)
    newlines = []
    if os.path.isfile(pro_path):
        with open(pro_path, 'r') as f:
            new_exp = True
            for line in f:
                line = line.strip()
                line = line.split("\t")
                if line[0] == exp:
                    newlines.append([exp, args.rn, args.o.split("_")[-1],
                                     args.organism, tar_reg, "0",
                                     "-", "1.0", r_genes, "No triplex found"])
                    new_exp = False
                else:
                    newlines.append(line)
            if new_exp:
                newlines.append([exp, args.rn, args.o.split("_")[-1],
                                 args.organism, tar_reg, "0",
                                 "-", "1.0", r_genes, "No triplex found"])
    else:
        newlines.append(["Experiment", "RNA_names", "Tag", "Organism", "Target_region", "No_sig_DBDs",
                         "Top_DBD", "p-value", "closest_genes"])
        newlines.append([exp, args.rn, args.o.split("_")[-1],
                         args.organism, tar_reg, "0",
                         "-", "1.0", r_genes, "No triplex found"])
    with open(pro_path, 'w') as f:
        for lines in newlines:
            print("\t".join(lines), file=f)

    # shutil.rmtree(args.o)
    # list_all_index(path=os.path.dirname(args.o), show_RNA_ass_gene=promoter.rna_regions)
    revise_index(root=os.path.dirname(os.path.dirname(args.o)), show_RNA_ass_gene=promoter.rna_regions)
    shutil.rmtree(args.o)
    sys.exit(1)

def write_stat(stat, filename):
    """Write the statistics into file"""
    order_stat = ["name", "genome", "exons", "seq_length",
                  "target_regions", "background_regions",
                  "DBD_all", "DBD_sig",
                  "DBSs_target_all", "DBSs_target_DBD_sig",
                  "DBSs_background_all", "DBSs_background_DBD_sig"]
    with open(filename, "w") as f:
        for k in order_stat:
            print("\t".join([k,stat[k]]), file=f)

def integrate_stat(path):
    """Integrate all statistics within a directory"""
    base = os.path.basename(path)
    order_stat = ["name", "genome", "exons", "seq_length",
                  "target_regions", "background_regions",
                  "DBD_all", "DBD_sig",
                  "DBSs_target_all", "DBSs_target_DBD_sig",
                  "DBSs_background_all", "DBSs_background_DBD_sig"]
    nested_dict = lambda: defaultdict(nested_dict)
    data = nested_dict()

    for item in os.listdir(path):
        pro = os.path.join(path, item, "stat.txt")
        if os.path.isfile(pro):
            with open(pro) as f:
                for line in f:
                    l = line.split()
                    data[item][l[0]] = l[1]
    with open(os.path.join(path,"statistics_"+base+".txt"), "w") as g:
        print("\t".join(order_stat), file=g)

        for item in data.keys():
            print("\t".join([data[item][o] for o in order_stat]), file=g)