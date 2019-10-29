# Python Libraries

import os
import re
import pysam
import shutil
import pickle
import errno
from ctypes import *
from collections import *
import natsort as natsort_ob

# Local Libraries
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
import matplotlib
from matplotlib import colors
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator


# Distal Libraries
from ..SequenceSet import SequenceSet
from ..viz.plotTools import output_array
from ..GenomicRegion import GenomicRegion
from .RNADNABindingSet import RNADNABindingSet
from ..GenomicRegionSet import GenomicRegionSet
from ..motifanalysis.Statistics import multiple_test_correction
from ..Util import SequenceType, Html, ConfigurationFile, GenomeData, LibraryPath

# Color code for all analysis
target_color = "mediumblue"
nontarget_color = "darkgrey"
sig_color = "powderblue"



order_stat = ["title", "name", "genome", "exons", "seq_length",
              "target_regions", "background_regions",
              "DBD_all", "DBD_sig",
              "DBSs_target_all", "DBSs_target_DBD_sig",
              "DBSs_background_all", "DBSs_background_DBD_sig", "p_value",
              "Norm_DBD", "Norm_DBS", "Norm_DBS_sig",
              "associated_gene", "expression", "loci", "autobinding",
              "MA_G","MA_T","MP_G","MP_T","RA_A","RA_G","YP_C","YP_T",
              "uniq_MA_G", "uniq_MA_T", "uniq_MP_G", "uniq_MP_T", "uniq_RA_A", "uniq_RA_G", "uniq_YP_C", "uniq_YP_T",
              "target_in_trans", "target_in_cis", "target_local",
              "background_in_trans", "background_in_cis", "background_local"]

              # "Mix_Antiparallel_A", "Mix_Antiparallel_G", "Mix_Antiparallel_T",
              # "Mix_Parallel_C", "Mix_Parallel_G", "Mix_Parallel_T",
              # "Purine_Antiparallel_A", "Purine_Antiparallel_G", "Purine_Antiparallel_T",
              # "Pyrimidine_Parallel_C", "Pyrimidine_Parallel_G", "Pyrimidine_Parallel_T"]
####################################################################################
####################################################################################

def print2(summary, string):
    """ Show the message on the console and also save in summary. """
    print(string)
    summary.append(string)


def output_summary(summary, directory, filename):
    """Save the summary log file into the defined directory"""
    pd = os.path.join(current_dir, directory)
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
    return list(map(try_int, re.findall(r'(\d+|\D+)', s)))


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


def list_all_index(path, link_d=None):
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
    type_list = 'sssssssssssssssssss'
    col_size_list = [20] * 20
    c = 0

    header_list = ["No.", "Experiments", "RNA", "Closest genes",
                   "Exon", "Length", "Expression*",
                   "Norm DBS*",
                   "Norm DBD*",  "No sig. DBD",
                   "Organism", "Target region",
                   "Rank*"]

    profile_f = open(os.path.join(path, "profile.txt"), 'r')
    profile = {}
    for line in profile_f:
        line = line.strip()
        line = line.split("\t")
        if line[0] == "Experiment": continue
        elif len(line) > 5: profile[line[0]] = line[1:]
    profile_f.close()

    # sig_list = []

    for i, exp in enumerate(profile):
        c += 1
        if profile[exp][10] == "-":
            new_line = [str(c), exp, profile[exp][0]]
        else:
            new_line = [str(c), '<a href="' + os.path.join(exp, "index.html") + '">' + exp + "</a>", profile[exp][0]]
        new_line += [ profile[exp][12],#3 close genes
                      profile[exp][1], #4 exon
                      profile[exp][2], #5 length
                      profile[exp][13] ]#6 exp

        if float(profile[exp][11]) < 0.05:
            new_line += [ profile[exp][6], #7 norm DBS
                          profile[exp][8], #8 norm DBD
                          profile[exp][9]] #9 sig DBD
                          # profile[exp][10], #10 Top DBD
                          # "<font color=\"red\">" + \
                          # profile[exp][11] + "</font>"]
            # sig_list.append(True)
        else:
            new_line += [str(0),  # 7 norm DBS
                         str(0),  # 8 norm DBD
                         profile[exp][9]]  # 9 sig DBD
                         # profile[exp][10],  # 10 Top DBD
                         # profile[exp][11]]
            # sig_list.append(False)

        new_line += [ profile[exp][4], profile[exp][5] ]

        data_table.append(new_line)

    rank_dbd = len(data_table) - rank_array([float(x[8]) for x in data_table])
    rank_dbs = len(data_table) - rank_array([float(x[7]) for x in data_table])

    rank_exp = len(data_table) - rank_array([0 if x[6] == "n.a." else float(x[6]) for x in data_table ])

    rank_sum = [x + y + z for x, y, z  in zip(rank_dbd, rank_dbs, rank_exp)]

    nd = [ d + [str(rank_sum[i])] for i, d in enumerate(data_table) ]

    nd = natsort_ob.natsorted(nd, key=lambda x: x[-1])
    html.add_zebra_table(header_list, col_size_list, type_list, nd,
                         align=10, cell_align="left", sortable=True)

    html.add_fixed_rank_sortable()
    html.write(os.path.join(path, "index.html"))


def revise_index(root):
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
    for d, p in plist.items():
        list_all_index(path=os.path.dirname(p),
                       link_d=dir_list)

def update_profile(dirpath, expression, name_list=None):
    header_list = ["Experiment", "RNA_names", "Tag", "Organism", "Target_region", "No_sig_DBDs",
                   "Top_DBD", "p-value", "closest_genes",
                   "Mix_Antiparallel_A", "Mix_Antiparallel_G", "Mix_Antiparallel_T",
                   "Mix_Parallel_C", "Mix_Parallel_G", "Mix_Parallel_T",
                   "Purine_Antiparallel_A", "Purine_Antiparallel_G", "Purine_Antiparallel_T",
                   "Pyrimidine_Parallel_C", "Pyrimidine_Parallel_G", "Pyrimidine_Parallel_T"]
    profiles = []
    pro = os.path.join(dirpath, "profile.txt")
    if not os.path.isfile(pro):
        print("There is no profile.txt in this directory.")
        return

    if expression:
        gene_exp = {}
        with open(expression) as f:
            for line in f:
                l = line.strip().split()
                gene_exp[l[0].partition(".")[0]] = l[1]
        profile_temp = []
        with open(pro) as f:
            for line in f:
                if not line.startswith("Experiment"):
                    l = line.strip().split("\t")
                    if len(l) == 12:
                        profile_temp.append(l + [gene_exp[l[0]]])
        with open(pro, "w") as f:
            h = ["Experiment", "RNA_names", "exon", "length",
                 "Tag", "Organism", "Target_region",
                 "Norm_DBS", "Norm_DBS_on_sig_DBD",
                 "Norm_DBD", "No_sig_DBDs", "Top_DBD",
                 "p-value", "closest_genes", "expression"]
            print("\t".join(h), file=f)
            for line in profile_temp:
                print("\t".join(line), file=f)
    if name_list:
        pass
    else:
        for item in os.listdir(dirpath):
            stat = os.path.join(dirpath, item, "stat.txt")
            summary = os.path.join(dirpath, item, "summary.txt")
            if os.path.isfile(stat) and os.path.isfile(summary):
                with open(stat) as f:
                    for line in f:
                        l = line.strip().split()
                        if l[0] == "name":
                            each_name = l[1]
                            each_tag = l[1].split("_")[-1]
                        if l[0] == "genome": each_organism = l[1]
                        if l[0] == "DBD_sig": each_DBD_sig = l[1]
                        if l[0] == "p_value": each_p_value = l[1]
                        if l[0] == "Mix_Antiparallel_A": MA_A = l[1]
                        if l[0] == "Mix_Antiparallel_G": MA_G = l[1]
                        if l[0] == "Mix_Antiparallel_T": MA_T = l[1]
                        if l[0] == "Mix_Parallel_C": MP_C = l[1]
                        if l[0] == "Mix_Parallel_G": MP_G = l[1]
                        if l[0] == "Mix_Parallel_T": MP_T = l[1]
                        if l[0] == "Purine_Antiparallel_A": RA_A = l[1]
                        if l[0] == "Purine_Antiparallel_G": RA_G = l[1]
                        if l[0] == "Purine_Antiparallel_T": RA_T = l[1]
                        if l[0] == "Pyrimidine_Parallel_C": YP_C = l[1]
                        if l[0] == "Pyrimidine_Parallel_G": YP_G = l[1]
                        if l[0] == "Pyrimidine_Parallel_T": YP_T = l[1]

                with open(summary) as g:
                    for line in g:
                        if "rgt-TDF" in line and " -de " in line:
                            l = line.strip().split()
                            each_target_region = os.path.basename(l[l.index("-de") + 1])
                profiles.append([item,each_name,each_tag,each_organism,each_target_region,
                                 each_DBD_sig,"n.a.",each_p_value,"-",
                                 MA_A, MA_G, MA_T, MP_C, MP_G, MP_T,
                                 RA_A, RA_G, RA_T, YP_C, YP_G, YP_T])

    with open(pro, "w") as f:
        print("\t".join(header_list), file=f)
        for line in profiles:
            print("\t".join(line), file=f)


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
    exps = natsorted(list(matrix.keys()))
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


    fig = plt.figure(figsize=(len(matrix) * 1.5, len(rnas) * 2.5))
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

    axmatrix.set_xticks(list(range(data.shape[1])))
    axmatrix.set_xticklabels(exps, minor=False, ha="left")
    axmatrix.xaxis.set_label_position('top')
    axmatrix.xaxis.tick_top()
    plt.xticks(rotation=40, fontsize=10)

    axmatrix.set_yticks(list(range(data.shape[0])))
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
    try:
        fig.savefig(os.path.join(path, 'condition_lncRNA_dendrogram.png'))
        # fig.savefig(os.path.join(path, 'condition_lncRNA_dendrogram.pdf'), format="pdf")
    except:
        pass


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
    if (isinstance(value,str)):
        try: value = float(value)
        except: return value
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

    7  8  9  10  11  12  13  14  15          16                 17
    l, e, c, fr, fm, of, mf, rm, filter_bed, self.genome_path,  par
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
                       par=input[17], genome_path=input[16],
                       dna_fine_posi=False)

    txp.merge_rbs(rbss=input[5], rm_duplicate=True)

    txpf = find_triplex(rna_fasta=input[1], dna_region=random, temp=input[3], 
                       organism=input[4], prefix=str(input[0]), remove_temp=True, 
                       l=int(input[7]), e=int(input[8]),  c=input[9], fr=input[10], 
                       fm=input[11], of=input[12], mf=input[13], rm=input[14], 
                       par=input[17], genome_path=input[16],
                       dna_fine_posi=True)

    txpf.merge_rbs(rbss=input[5], rm_duplicate=True)
    sys.stdout.flush()
    print("".join(["="]*int(input[6])), end="")

    return [ [len(tr) for tr in list(txp.merged_dict.values()) ], [len(dbss) for dbss in list(txpf.merged_dict.values())] ]


def save_sequence(dir_name, filename, regions, genome_path):
    """
    Fetch sequence into FASTA file according to the given BED file
    """
    genome = pysam.Fastafile(genome_path)
    # print(regions)
    with open(os.path.join(dir_name, filename), 'w') as output:
        for region in regions:
            if "_" not in region.chrom:
                print(">"+ region.toString(), file=output)
                print(genome.fetch(region.chrom, max(0, region.initial), region.final), file=output)


def find_triplex(rna_fasta, dna_region, temp, organism, l, e, dna_fine_posi, genome_path, prefix="", remove_temp=False, 
                 c=None, fr=None, fm=None, of=None, mf=None, rm=None, par="", autobinding=False, seq=False):
    """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
    
    # Generate FASTA 
    save_sequence(dir_name=temp, filename="dna_"+prefix+".fa", regions=dna_region, genome_path=genome_path)

    # Triplexator
    run_triplexator(ss=rna_fasta, ds=os.path.join(temp,"dna_"+prefix+".fa"), 
                    output=os.path.join(temp, "dna_"+prefix+".txp"), 
                    l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, rm=rm, par=par)
    # Autobinding
    if autobinding:
        run_triplexator(ss=rna_fasta, ds=os.path.join(temp, "dna_" + prefix + ".fa"),
                        output=os.path.join(temp, "autobinding_" + prefix + ".txp"),
                        l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, rm=rm, par=par + "_auto-binding-file")
    # Read txp
    txp = RNADNABindingSet("dna")
    txp.read_tpx(os.path.join(temp, "dna_"+prefix+".txp"), dna_fine_posi=dna_fine_posi, seq=True)
    txp.remove_duplicates()

    if remove_temp:
        silentremove(os.path.join(temp,"dna_"+prefix+".fa"))
        silentremove(os.path.join(temp,"dna_"+prefix+".txp"))

    return txp


def run_triplexator(ss, ds, output, l=None, e=None, c=None, fr=None, fm=None, of=None, mf=None, rm=None, par="", autobinding=None, summary_file=False):
    """Perform Triplexator"""
    #triplexator_path = check_triplexator_path()
    # triplexator -ss -ds -l 15 -e 20 -c 2 -fr off -fm 0 -of 1 -rm
    triclass = LibraryPath()
    triplex_lib_path = triclass.get_triplexator()
    triplex_lib  = cdll.LoadLibrary(triplex_lib_path)

    arguments = ""
    if not autobinding:
        if ss: arguments += "-ss "+ss+" "
        if ds: arguments += "-ds "+ds+" "
    else:
        arguments += "-as " + autobinding + " "

    if l: arguments += "-l "+str(l)+" "
    if e: arguments += "-e "+str(e)+" "
    if c: arguments += "-c "+str(c)+" "
    if fr: arguments += "-fr "+fr+" "
    if fm: arguments += "-fm "+str(fm)+" "
    if of: arguments += "-of "+str(of)+" "
    if mf: arguments += "-mf "
    if rm: arguments += "-rm "+str(rm)+" "
    arguments += "--bit-parallel -g 0 "
    if par != "":
        par = par.replace('_'," ")
        par = "-" + par
        arguments += par+" "
    
    arguments += "-o "+ os.path.basename(output) + " -od " + os.path.dirname(output)

    arg_strings  = arguments.split(' ')
    # print(arg_strings)
    arg_ptr      = (c_char_p * (len(arg_strings) + 1))()

    arg_ptr[0] = b"triplexator -bp"  # to simulate calling from cmd line
    for i, s in enumerate(arg_strings):
        arg_ptr[i + 1] = s.encode("utf-8")
    # print(arg_strings)
    triplex_lib.pyTriplexator(len(arg_strings) + 1, arg_ptr)
    if not summary_file:
        silentremove(os.path.join(output + ".summary"))
    silentremove(os.path.join(output + ".log"))

def run_triplexes(arguments):
    triclass = LibraryPath()
    triplex_lib_path = triclass.get_triplexator()
    triplex_lib = cdll.LoadLibrary(triplex_lib_path)
    arg_strings = arguments
    # print(arg_strings)
    arg_ptr = (c_char_p * (len(arg_strings) + 1))()

    arg_ptr[0] = b"triplexator"  # to simulate calling from cmd line
    for i, s in enumerate(arg_strings):
        arg_ptr[i + 1] = s.encode("utf-8")
    triplex_lib.pyTriplexator(len(arg_strings) + 1, arg_ptr)

def region_link_internet(organism, region):
    ani = None
    if organism == "hg19":
        ani = "human"
    elif organism == "hg38":
        ani = "human"
    elif organism == "mm9":
        ani = "mouse"
    if ani:
        region_link = "".join(['<a href="http://genome.ucsc.edu/cgi-bin/hgTracks?db=', organism,
                               "&position=", region.chrom, "%3A", str(region.initial), "-",
                               str(region.final), '" style="text-align:left" target="_blank">',
                               region.toString(space=True), '</a>'])
    else:
        if organism == "tair10":
            region_link = "".join(
                ['<a href="http://tairvm17.tacc.utexas.edu/cgi-bin/gb2/gbrowse/arabidopsis/?name=',
                 region.chrom, "%3A", str(region.initial), "..", str(region.final),
                 '" target="_blank">',
                 region.toString(space=True), '</a>'])
        else:
            region_link = region.toString(space=True)
    return region_link




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
        asso_genes = g.gene_association(organism=organism, promoter_length=1000, show_dis=True)

        # genes = asso_genes[0].name.split(":")
        closest_genes = []
        for g in asso_genes:
            gnames = g.name.split(":")
            for n in gnames:
                if name not in n: closest_genes.append(n)
        closest_genes = set(closest_genes)
        if len(closest_genes) == 0:
            return "."
        else:
            return ":".join(closest_genes)
    else:
        return "."

def rank_array(a):
    try:
        a = numpy.array(a)
    except:
        a = numpy.array([float(b) for b in a])
    sa = numpy.searchsorted(numpy.sort(a), a)
    return sa

def dbd_regions(exons, sig_region, rna_name, output,out_file=False, temp=None, fasta=True):
    """Generate the BED file of significant DBD regions and FASTA file of the sequences"""
    if len(sig_region) == 0:
        return
    #print(self.rna_regions)
    if not exons:
        pass
    else:
        dbd = GenomicRegionSet("DBD")
        dbdmap = {}
        if len(exons) == 1:
            print("## Warning: If more than 1 exon, the DBD position may be problematic. ")
        for rbs in sig_region:
            loop = True

            if exons[0][3] == "-":
                while loop:
                    cf = 0
                    for exon in exons:
                        #print(exon)

                        l = abs(exon[2] - exon[1])
                        tail = cf + l

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
            dbd.write(filename=os.path.join(output, rna_name+"_DBDs.bed"))
        else:
            # print(dbd)
            # print(dbd.sequences[0])
            dbd.write(filename=output)
    # FASTA
    if fasta:
        #print(dbdmap)
        if not out_file:
            seq = pysam.Fastafile(os.path.join(output,"rna_temp.fa"))
            fasta_f = os.path.join(output, rna_name+"_DBDs.fa")
        else:
            seq = pysam.Fastafile(os.path.join(temp,"rna_temp.fa"))
            fasta_f = output+".fa"

        with open(fasta_f, 'w') as fasta:
            for rbs in sig_region:
                print(">"+ rna_name +":"+str(rbs.initial)+"-"+str(rbs.final), file=fasta)
                s = seq.fetch(rna_name, max(0, rbs.initial), rbs.final)
                for ss in [s[i:i + 80] for i in range(0, len(s), 80)]:
                    print(ss, file=fasta)

        seq.close()

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

def get_dbss(input_BED,output_BED,rna_fasta,output_rbss,organism,l,e,c,fr,fm,of,mf,rm,temp):
    regions = GenomicRegionSet("Target")
    regions.read(input_BED)
    regions.gene_association(organism=organism, show_dis=True)

    connect_rna(rna_fasta, temp=temp, rna_name="RNA")
    rnas = SequenceSet(name="rna", seq_type=SequenceType.RNA)
    rnas.read_fasta(os.path.join(temp,"rna_temp.fa"))
    rna_regions = get_rna_region_str(os.path.join(temp,rna_fasta))
    # print(rna_regions)
    genome = GenomeData(organism)
    genome_path = genome.get_genome()
    txp = find_triplex(rna_fasta=rna_fasta, dna_region=regions, 
                       temp=temp, organism=organism, remove_temp=False,
                       l=l, e=e, c=c, fr=fr, fm=fm, of=of, mf=mf, genome_path=genome_path,
                       prefix="targeted_region", dna_fine_posi=True)

    print("Total binding events:\t",str(len(txp)))
    txp.write_bed(output_BED)
    txp.write_txp(filename=output_BED.replace(".bed",".txp"))
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

    def decode_loci(loci):
        l = loci.partition("chr")[2]
        chrom = "chr" + l.partition(":")[0]
        start = int(l.partition(":")[2].partition("-")[0])
        end = int(l.partition(":")[2].partition("-")[2].split()[0])
        return chrom, start, end

    rna_regions = []
    with open(rna) as f:
        for line in f:
            if line[0] == ">":
                line = line.strip()
                # Loci
                if "REGION" in line:
                    l = line.split()
                    for i, e in enumerate(l):
                        if "REGION" in e:
                            e = e.split("_")
                            #print(e)
                            try:
                                r = [e[1], int(e[2]), int(e[3]), e[4]]
                            except:
                                r = [e[1], int(e[3]), int(e[4]), e[5]]
                
                elif "chr" in line:
                    try:
                        line = line.split()
                        loci_ind = 0
                        sign = None
                        score = None
                        for i, l in enumerate(line):
                            if l.startswith("chr") and ":" in l and "-" in l:
                                c, s, e = decode_loci(l)
                                loci_ind = i
                            if "strand" in l:
                                sign = l.partition("strand=")[2][0]
                            if "score" in l:
                                score = float(l.partition("score=")[2])

                        if not sign:
                            sign = line[loci_ind+1]
                        if sign == "+" or sign == "-" or sign == ".":
                            r = [c, s, e, sign]
                        if score:
                            r.append(score)
                    except:
                        break

                else:
                    r = []
                # Score
                if r:
                    rna_regions.append(r)

    return rna_regions


def no_binding_response(args, stat):
    print("*** Find no DBD having DBS with cutoff = "+str(args.ccf))

    # pro_path = os.path.join(os.path.dirname(args.o), "profile.txt")
    # exp = os.path.basename(args.o)
    # try:
    #     if args.de:
    #         tar = args.de
    #     else:
    #         tar = args.bed
    # except:
    #     tar = args.bed

    # stat["DBD_all"] = 0
    # stat["DBD_sig"] = 0
    # stat["DBSs_target_all"] = 0
    # stat["DBSs_target_DBD_sig"] = 0
    # stat["DBSs_background_all"] = 0
    # stat["DBSs_background_DBD_sig"] = 0
    # stat["p_value"] = "-"

    # save_profile(rna_regions=rna_regions, rna_name=rna_name,
    #              organism=organism, output=args.o, bed=args.bed,
    #              geneset=tar, stat=stat, topDBD=["-", 1], sig_DBD=[],
    #              expression=expression)

    # revise_index(root=os.path.dirname(os.path.dirname(args.o)))
    # shutil.rmtree(args.o)
    for f in os.listdir(args.o):
        if os.path.isfile(os.path.join(args.o, f)):
            silentremove(os.path.join(args.o, f))

    write_stat(stat, os.path.join(args.o, "stat.txt"))
    # sys.exit(0)

def write_stat(stat, filename):
    """Write the statistics into file"""

    with open(filename, "w") as f:
        for k in order_stat:
            try:
                print("\t".join([k, str(stat[k])]), file=f)
            except:
                continue

def integrate_stat(path):
    """Integrate all statistics within a directory"""
    base = os.path.basename(path)

    data = {}

    for item in os.listdir(path):
        pro = os.path.join(path, item, "stat.txt")
        if os.path.isfile(pro):
            data[item] = {}
            with open(pro) as f:
                for line in f:
                    l = line.strip().split("\t")
                    data[item][l[0]] = l[1]
            # print(data[item])

            data[item]["Norm_DBD"] = value2str(float(data[item]["DBD_all"]) / int(data[item]["seq_length"]) * 1000)
            data[item]["Norm_DBS"] = value2str(float(data[item]["DBSs_target_all"]) / int(data[item]["seq_length"]) * 1000)
            data[item]["Norm_DBS_sig"] = value2str(float(data[item]["DBSs_target_DBD_sig"]) / int(data[item]["seq_length"]) * 1000)
            data[item]["title"] = item

    with open(os.path.join(path, "statistics_"+base+".txt"), "w") as g:
        print("\t".join(order_stat), file=g)
        for item in data:
            print("\t".join([data[item][o] for o in order_stat]), file=g)

def summerize_stat(target, link_d, score=True):
    base = os.path.basename(target)
    h = os.path.join(target, "index.html")
    stat = os.path.join(target, "statistics_" + base + ".txt")
    fp = "./style"
    html = Html(name=base, links_dict=link_d,
                fig_rpath=fp, homepage="../index.html",
                RGT_header=False, other_logo="TDF")
    html.add_heading(target)
    data_table = []
    type_list = 'sssssssssssss'
    col_size_list = [20] * 20
    c = 0
    if score:
        header_list = ["No.", "RNA", "Closest genes",
                       "Exon", "Length", "Score*",
                       "Norm DBS*", "Norm DBD*", "Number sig_DBD", "Autobinding",
                       "Organism", "Target region", "Rank*"]
    else:
        header_list = ["No.", "RNA", "Closest genes",
                       "Exon", "Length",
                       "Norm DBS*", "Norm DBD*", "Number sig_DBD", "Autobinding",
                       "Organism", "Target region", "Rank*"]

    with open(stat) as f_stat:
        for line in f_stat:
            if line.startswith("title"):
                continue
            elif not line.strip():
                continue
            else:
                c += 1
                l = line.strip().split()
                hh = "./" + l[0] + "/index.html"
                if l[14] == "0":
                    tag = l[0]
                else:
                    tag = '<a href="' + hh + '">' + l[0] + "</a>"
                if score:
                    data_table.append([str(c), tag, l[17],
                                       l[3], l[4], l[18],
                                       l[15], l[14], l[8], l[20],
                                       l[2], l[5]])
                else:
                    data_table.append([str(c), tag, l[17],
                                       l[3], l[4],
                                       l[15], l[14], l[8], l[20],
                                       l[2], l[5]])
    # print(data_table)
    rank_dbd = len(data_table) - rank_array([float(x[7]) for x in data_table])
    rank_dbs = len(data_table) - rank_array([float(x[6]) for x in data_table])
    if score:
        rank_exp = len(data_table) - rank_array([0 if x[5] == "n.a." else abs(float(x[5])) for x in data_table])
        rank_sum = [x + y + z for x, y, z in zip(rank_dbd, rank_dbs, rank_exp)]
    else:
        rank_sum = [x + y for x, y in zip(rank_dbd, rank_dbs)]

    for i, d in enumerate(data_table):
        d.append(str(rank_sum[i]))
    nd = natsort_ob.natsorted(data_table, key=lambda x: x[-1])
    html.add_zebra_table(header_list, col_size_list, type_list, nd,
                         sortable=True, clean=True)
    html.add_fixed_rank_sortable()
    html.write(h)


def merge_DBD_regions(path):
    """Merge all available DBD regions in BED format. """
    base = os.path.basename(path)
    dir_name = os.path.basename(os.path.dirname(path))
    dbd_pool = GenomicRegionSet(dir_name + "_" + base)
    for rna in os.listdir(path):
        if os.path.isdir(os.path.join(path, rna)):
            f = os.path.join(path, rna, rna+"_DBDs.bed")
            if os.path.exists(f):
                # print(f)
                dbd = GenomicRegionSet(rna)
                dbd.read(f)
                for r in dbd: r.name = rna+"_"+r.name
                dbd_pool.combine(dbd)
    # print(len(dbd_pool))
    dbd_pool.write(os.path.join(path, "DBD_"+dir_name + "_" + base +".bed"))

def merge_DBSs(path):
    """Merge all available DBD regions in BED format. """
    base = os.path.basename(path)
    dir_name = os.path.basename(os.path.dirname(path))
    dbss_pool = GenomicRegionSet(dir_name + "_" + base)
    for rna in os.listdir(path):
        if os.path.isdir(os.path.join(path, rna)):
            f = os.path.join(path, rna, rna+"_dbss.bed")
            if os.path.exists(f):
                # print(f)
                dbss = GenomicRegionSet(rna)
                dbss.read(f)
                for r in dbss: r.name = rna+"_"+r.name
                dbss_pool.combine(dbss)
    # print(len(dbd_pool))
    dbss_pool.write(os.path.join(path, "DBSs_"+dir_name + "_" + base +".bed"))

def merge_DNA_counts(path):
    """Merge all available DBD regions in BED format. """
    # base = os.path.basename(path)
    # dir_name = os.path.basename(os.path.dirname(path))
    m = []
    header = []
    for rna in os.listdir(path):
        if os.path.isdir(os.path.join(path, rna)):
            f = os.path.join(path, rna, "DNA_cov.txt")
            if os.path.exists(f):
                with open(f) as ff:
                    for i, line in enumerate(ff):
                        l = line.strip()
                        if i == 0:
                            header = ["RNA_name"] + l.split("\t")
                        else:
                            m.append([rna] + l.split("\t"))
    # print(len(m))
    with open(os.path.join(path, "RNA_DNA_matrix.txt"), "w") as fm:
        print("\t".join([str(x) for x in header]), file=fm)
        for l in m:
            print("\t".join([str(x) for x in l]), file=fm)

def save_profile(rna_regions, rna_name, organism, output, bed,\
                 stat, topDBD, sig_DBD, expression, geneset=None):
    """Save statistics for comparison with other results"""

    pro_path = os.path.join(os.path.dirname(output), "profile.txt")
    exp = os.path.basename(output)
    # tag = os.path.basename(os.path.dirname(rnafile))
    if geneset:
        tar_reg = os.path.basename(geneset)
    else:
        tar_reg = os.path.basename(bed)
    # RNA name with region
    # if rna_regions:
    #     exon = str(len(rna_regions))
    # else:
    #     exon = "-"
    # rna = self.rna_name
    # RNA associated genes
    r_genes = rna_associated_gene(rna_regions=rna_regions, name=rna_name, organism=organism)

    newlines = []

    this_rna = [exp, rna_name, stat["exons"], stat["seq_length"],
                output.split("_")[-1], organism, tar_reg,
                value2str(float(stat["DBSs_target_all"])/int(stat["seq_length"])*1000),
                value2str(float(stat["DBSs_target_DBD_sig"])/int(stat["seq_length"])*1000),
                value2str(float(stat["DBD_all"])/int(stat["seq_length"])*1000), str(len(sig_DBD)),
                topDBD[0], value2str(topDBD[1]), r_genes, value2str(expression)]
    # try:
    if os.path.isfile(pro_path):
        with open(pro_path, 'r') as f:
            new_exp = True
            for line in f:
                line = line.strip()
                line = line.split("\t")
                if line[0] == exp:
                    newlines.append(this_rna)
                    new_exp = False
                elif line[0] == "Experiment":
                    continue
                else:
                    newlines.append(line)
            if new_exp:
                newlines.append(this_rna)
    else:
        newlines.append(this_rna)

    try: newlines.sort(key=lambda x: float(x[12]))
    except: pass
    newlines = [["Experiment", "RNA_names", "exon", "length",
                 "Tag", "Organism", "Target_region",
                 "Norm_DBS", "Norm_DBS_on_sig_DBD",
                 "Norm_DBD", "No_sig_DBDs", "Top_DBD",
                 "p-value", "closest_genes", "expression"]] + newlines

    with open(pro_path, 'w') as f:
        for lines in newlines:
            print("\t".join(lines), file=f)

def silentremove(filename):
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

def integrate_html(target):
    # Project level index file
    condition_list = []  # name, link, no. tests, no. sig.

    for item in os.listdir(target):
        if item == "style": continue
        if os.path.isfile(os.path.join(target, item)):
            continue
        elif os.path.isdir(os.path.join(target, item)):
            h = os.path.join(item, "index.html")
            stat = os.path.join(target, item, "statistics_" + item + ".txt")

            if os.path.isfile(stat):
                nt = 0
                ns = 0
                with open(stat) as f:
                    for line in f:
                        line = line.strip().split("\t")
                        if line[0] == "title": continue
                        nt += 1
                        if line[13] == "-":
                            pass
                        elif float(line[13]) < 0.05:
                            ns += 1
                # print([item, h, str(nt), str(ns)])
                condition_list.append([item, h, str(nt), str(ns)])

    if len(condition_list) > 0:
        # print(condition_list)
        link_d = {}
        for item in os.listdir(os.path.dirname(target)):
            if os.path.isfile(os.path.dirname(target) + "/" + item + "/index.html"):
                link_d[item] = "../" + item + "/index.html"

        fp = condition_list[0][0] + "/style"
        html = Html(name="Directory: " + target, links_dict=link_d,
                    fig_rpath=fp,
                    RGT_header=False, other_logo="TDF")
        html.add_heading("All conditions in: " + target + "/")
        data_table = []
        type_list = 'sssssssssssss'
        col_size_list = [20] * 10
        c = 0
        header_list = ["No.", "Conditions", "No. tests", "No. sig. tests"]
        for i, exp in enumerate(condition_list):
            c += 1
            data_table.append([str(c),
                               '<a href="' + exp[1] + '">' + exp[0] + "</a>",
                               exp[2], exp[3]])
        html.add_zebra_table(header_list, col_size_list, type_list, data_table,
                             sortable=True, clean=True)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(target, "index.html"))

def shorten_dir(path):
    if path.count("/") < 3:
        return path
    else:
        n = path.count("/") - 3 + 1
        return path.split("/",n)[-1]

def purge(dir_name, pattern):
    for f in os.listdir(dir_name):
        print([f, pattern])
        if re.search(pattern, f):
            silentremove(os.path.join(dir_name, f))