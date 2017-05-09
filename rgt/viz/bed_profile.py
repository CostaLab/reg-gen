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
        size_panel = 8
        rows = len(self.beds)
        cols = 1
        if args.biotype:
            self.ind_col["Biotype"] = cols
            cols += 1
        if args.repeats:
            self.ind_col["Repeats"] = cols
            cols += 1
        if args.genposi:
            self.ind_col["Genetic position"] = cols
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

            ax.hist(dis, bins=200, color="cornflowerblue", linewidth=0)
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

    def plot_ref(self, ref_dir, tag, other=False):
        if os.path.isdir(ref_dir):
            refs = []
            refs_names = []
            for f in os.listdir(ref_dir):
                if f.endswith(".bed"):
                    name = os.path.basename(f).replace(".bed", "")
                    bed = GenomicRegionSet(name)
                    bed.read_bed(os.path.join(ref_dir, f))
                    # bed.merge()
                    refs.append(bed)
                    refs_names.append(name)
        else:
            print("*** Error: Not a valid directory: "+ref_dir)
            sys.exit(1)
        if other:
            refs_names.append("others")

        for i, bed in enumerate(self.beds):
            overlapping_counts = []
            for j, ref in enumerate(refs):
                overlapping_counts.append(bed.count_by_regionset(ref))
            if other:
                # b = bed
                # for r in refs:
                #     b = b.subtract(r, whole_region=True)
                # overlapping_counts.append(len(b))
                overlapping_counts.append(max(0, len(bed) - sum(overlapping_counts)))

            # Plot
            try:
                ax = self.fig_axs[i, self.ind_col[tag]]
            except:
                try:
                    ax = self.fig_axs[i]
                except:
                    ax = self.fig_axs

            ax.bar(left=range(len(overlapping_counts)), height=overlapping_counts, color="cornflowerblue", linewidth=0)
            ax.set_title(self.bednames[i])
            ax.set_ylabel("Number")
            ax.set_xticks([x + 0.4 for x in range(len(overlapping_counts))])
            ax.set_xticklabels(refs_names, fontsize=7, rotation=20, ha="right")

        ax.set_xlabel(tag)




    def save_fig(self, filename):
        self.fig_f.savefig(filename+".png", facecolor='w', edgecolor='w',
                  bbox_inches='tight', dpi=400)
        pp = PdfPages(filename+".pdf")
        pp.savefig(self.fig_f, bbox_inches='tight')
        pp.close()


    def gen_html(self, directory, title):
        dir_name = os.path.basename(directory)
        # check_dir(directory)
        html_header = dir_name + " / " + title
        link_d = OrderedDict()
        link_d["BED profile"] = "index.html"
        link_d["Parameters"] = "parameters.html"


        html = Html(name=html_header, links_dict=link_d, fig_dir=os.path.join(directory, "style"),
                    fig_rpath="../style", RGT_header=False, other_logo="viz", homepage="../index.html")

        type_list = 's'*10
        col_size_list = [10] * 10

        data_table = []
        header_list = ["No", "BED File", "Number", "Ave. Len.", "s.d. Len.", "Min", "Max", "Internal overlap"]
        c = 0
        for bed in self.bednames:
            c += 1
            data_table.append([ str(c), bed,
                                str(self.stats[bed]["Number of regions"]),
                                str(self.stats[bed]["Average length"]),
                                str(self.stats[bed]["s.d. of length"]),
                                str(self.stats[bed]["min"]),
                                str(self.stats[bed]["max"]),
                                str(self.stats[bed]["Internal overlaps"])])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table,
                             border_list=None, sortable=True, clean=True)
        html.add_figure("figure_"+title+".png", align=50)
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))