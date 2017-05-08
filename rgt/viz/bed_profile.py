# Python Libraries
from __future__ import print_function
from __future__ import division
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy

# Local Libraries
# Distal Libraries
from rgt.Util import Html
from rgt.CoverageSet import *
from rgt.ExperimentalMatrix import *
from shared_function import output_array, group_refque, color_groupded_region, multiple_correction, value2str
# Local test
dir = os.getcwd()

###########################################################################################
#                    BED Profile
###########################################################################################

class BED_profile:
    def __init__(self, input_path, organism, args):

        if os.path.isdir(input_path):
            self.beds = []
            self.bednames = []
            for dirpath, dnames, fnames in os.walk(input_path):
                for f in fnames:
                    if f.endswith(".bed"):
                        name = os.path.basename(f).replace(".bed", "")
                        bed = GenomicRegionSet(name)
                        bed.read_bed(os.path.join(dirpath, f))
                        self.beds.append(bed)
                        self.bednames.append(name)

        elif os.path.isfile(input_path):
            if input_path.endswith(".bed"):
                name = os.path.basename(input_path).replace(".bed", "")
                bed = GenomicRegionSet(name)
                bed.read_bed(input_path)
                self.beds = [bed]
                self.bednames = [name]
            else:
                self.EM = ExperimentalMatrix()
                self.EM.read(input)
                self.beds = self.EM.get_regionsets()
                self.bednames = self.EM.get_regionsnames()

        self.organism = organism
        self.chromosomes = GenomicRegionSet(organism)
        self.chromosomes.get_genome_data(organism=organism, chrom_X=True)
        self.stats = OrderedDict()
        self.ind_col = {}
        # self.fig_size = (6, 6)
        size_panel = 6
        rows = len(self.beds)
        cols = 1
        if args.biotype:
            self.ind_col["Biotype"] = cols
            cols += 1
        if args.repeats:
            self.ind_col["Repeats"] = cols
            cols += 1
        self.fig_f, self.fig_axs = plt.subplots(rows, cols, dpi=300, figsize=(cols*size_panel, rows*size_panel))
        # print(len(self.fig_axs))
        # self.figures = OrderedDict()
        # for bed in self.bednames:
        #     self.figures[bed] = []


    def cal_statistics(self):
        for i, bed in enumerate(self.beds):
            self.stats[self.bednames[i]] = {}
            self.stats[self.bednames[i]]["Number of regions"] = len(bed)
            self.stats[self.bednames[i]]["Average length"] = bed.average_size()
            self.stats[self.bednames[i]]["s.d. of length"] = bed.size_variance()
            self.stats[self.bednames[i]]["max"] = bed.max_size()
            self.stats[self.bednames[i]]["min"] = bed.min_size()
            self.stats[self.bednames[i]]["Internal overlaps"] = len(bed) - len(bed.merge(w_return=True))
        # print(self.stats)

    def plot_distribution_length(self, outdir):
        for i, bed in enumerate(self.beds):
            dis = [len(r) for r in bed]
            try:
                ax = self.fig_axs[i, 0]
            except:
                try:
                    ax = self.fig_axs[i]
                except:
                    ax = self.fig_axs

            ax.hist(dis, bins=500)
            ax.set_title(self.bednames[i])
            ax.set_xlabel("Length (bp)")
            ax.set_ylabel("Frequency")

    # def plot_distribution_chromosome(self, outdir):
    #
    #     for i, bed in enumerate(self.beds):
    #         dis = [len(r) for r in bed]
    #         try:
    #             ax = self.fig_axs[i, 1]
    #         except:
    #             try:
    #                 ax = self.fig_axs[i]
    #             except:
    #                 ax = self.fig_axs
    #
    #         ax.hist(dis)
    #         ax.set_title(self.bednames[i])
    #         ax.set_xlabel("Length (bp)")
    #         ax.set_ylabel("Frequency")

    def plot_ref(self, ref_dir, tag):
        if os.path.isdir(ref_dir):
            refs = []
            refs_names = []
            for f in os.listdir(ref_dir):
                if f.endswith(".bed"):
                    name = os.path.basename(f).replace(".bed", "")
                    bed = GenomicRegionSet(name)
                    bed.read_bed(os.path.join(ref_dir, f))
                    refs.append(bed)
                    refs_names.append(name)
        else:
            print("*** Error: Not a valid directory: "+ref_dir)
            sys.exit(1)

        for i, bed in enumerate(self.beds):
            overlapping_counts = []
            for j, ref in enumerate(refs):
                overlapping_counts.append(bed.count_by_regionset(ref))
            # Plot
            try:
                ax = self.fig_axs[i, self.ind_col[tag]]
            except:
                try:
                    ax = self.fig_axs[i]
                except:
                    ax = self.fig_axs

            ax.bar(left=range(len(overlapping_counts)), height=overlapping_counts)

            ax.set_ylabel("Number")
            ax.set_xticks(range(len(overlapping_counts)))
            ax.set_xticklabels(refs_names, fontsize=7, rotation=20, ha="right")
        ax.set_xlabel(tag)




    def save_pdf(self, filename):

        pp = PdfPages(filename)
        pp.savefig(self.fig_f, bbox_inches='tight')
        pp.close()