# Python Libraries
from __future__ import print_function
from __future__ import division
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy
import natsort

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

            index = natsort.index_natsorted(self.bednames)
            self.beds = natsort.order_by_index(self.beds, index)
            self.bednames = natsort.order_by_index(self.bednames, index)

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
        if args.genposi:
            self.ind_col["Genetic position"] = cols
            cols += 1
        if args.labels:
            for label in args.labels:
                self.ind_col[label] = cols
                cols += 1
        self.fig_f, self.fig_axs = plt.subplots(rows+1, cols, dpi=300, figsize=(cols*size_panel, rows*size_panel))
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

    def plot_distribution_length(self):
        dis = []
        for i, bed in enumerate(self.beds):
            dis.append([numpy.log10(len(r)) for r in bed])
        max_len = max([max(x) for x in dis])

        for i in range(len(self.beds)+1):
            try:
                ax = self.fig_axs[i, 0]
            except:
                try:
                    ax = self.fig_axs[i]
                except:
                    ax = self.fig_axs
            if i == 0:
                ax.set_title("Length")
                violin = ax.violinplot(dis, [x + 1 for x in range(len(dis))], widths=0.7, showmeans=False,
                                       showextrema=True, showmedians=True)

                for pc in violin['bodies']:
                    pc.set_color('cornflowerblue')
                    # pc.set_edgecolor('cornflowerblue')
                    pc.set_alpha(0.6)
                # violin['cmeans'].set_color('b')
                violin['cmins'].set_color('b')
                violin['cmaxes'].set_color('b')
                violin['cbars'].set_color('b')
                violin['cmedians'].set_color('b')
                # bp = ax.boxplot(dis, notch=False, sym='o', vert=True, whis=1.5, positions=None,
                #                 widths=None, patch_artist=True, bootstrap=None,
                #                 autorange=True)
                # z = 10
                # plt.setp(bp['whiskers'], color='black', linestyle='-', linewidth=0.8, zorder=z)
                # plt.setp(bp['fliers'], markerfacecolor='gray', color='white', alpha=0.3, markersize=1.8, zorder=z)
                # plt.setp(bp['caps'], color='gray', zorder=z)
                # plt.setp(bp['medians'], color='black', linewidth=1.5, zorder=z + 1)
                # plt.setp(bp['boxes'], linewidth=0, facecolor="cornflowerblue")

                ax.set_xticks([x + 1 for x in range(len(dis))])
                ax.set_xticklabels(self.bednames, fontsize=7, rotation=20, ha="right")
                ax.set_ylabel("Length (log10)")

            elif i > 0:
                ax.hist(dis[i-1], bins=100, color="cornflowerblue", linewidth=0)
                ax.set_title(self.bednames[i-1])
                ax.set_xlabel("Length (log10)")
                ax.set_ylabel("Frequency")
                ax.set_xlim([0, max_len])

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

    def plot_ref(self, ref_dir, tag, other=False, strand=False):
        print("Processing "+tag+" ....")
        refs = []
        refs_names = []
        if os.path.isdir(ref_dir):
            for f in os.listdir(ref_dir):
                if f.endswith(".bed"):
                    name = os.path.basename(f).replace(".bed", "")
                    bed = GenomicRegionSet(name)
                    bed.read_bed(os.path.join(ref_dir, f))
                    # bed.merge()
                    refs.append(bed)
                    refs_names.append(name)
        elif os.path.isfile(ref_dir) and ref_dir.endswith(".bed"):
            name = os.path.basename(ref_dir).replace(".bed", "")
            bed = GenomicRegionSet(name)
            bed.read_bed(ref_dir)
            # bed.merge()
            refs.append(bed)
            refs_names.append(name)
        else:
            print("*** Error: Not a valid directory: "+ref_dir)
            sys.exit(1)

        index = natsort.index_natsorted(refs_names)
        refs = natsort.order_by_index(refs, index)
        refs_names = natsort.order_by_index(refs_names, index)

        if other:
            refs_names.append("Else")
        if strand:
            ref_plus = []
            ref_minus = []
            for ref in refs:
                ref_plus.append(ref.filter_strand(strand="+"))
                ref_minus.append(ref.filter_strand(strand="-"))
        # Counting through all references
        overlapping_counts = []
        for i, bed in enumerate(self.beds):
            c = []
            if strand:
                bed_plus = bed.filter_strand(strand="+")
                bed_minus = bed.filter_strand(strand="-")
            for j, ref in enumerate(refs):
                if strand:
                    c.append(bed_plus.count_by_regionset(ref_plus[j]) + bed_minus.count_by_regionset(ref_minus[j]))
                else:
                    c.append(bed.count_by_regionset(ref))
            if other:

                c.append(max(0, len(bed) - sum(c)))
            overlapping_counts.append(c)

        if other:
            color_list = plt.cm.Set1(numpy.linspace(0, 1, len(refs_names))).tolist()
        else:
            color_list = plt.cm.Set1(numpy.linspace(0, 0.95, len(refs_names))).tolist()

        for i in range(len(self.beds) + 1):
            # Plot
            try:
                ax = self.fig_axs[i, self.ind_col[tag]]
            except:
                try:
                    ax = self.fig_axs[i]
                except:
                    ax = self.fig_axs
            if i == 0:
                # print(overlapping_counts)
                proportion = []
                for counts in overlapping_counts:
                    ss = sum(counts)
                    proportion.append([x/ss * 100 for x in counts])
                # print(proportion)
                ptable = []
                for j in range(len(refs_names)):
                    ptable.append([x[j] for x in proportion])
                # print(ptable)
                width = 0.6
                bottom = [0] * len(self.bednames)
                for j, y in enumerate(ptable):
                    ax.bar(range(len(self.bednames)), y, width=width, bottom=bottom, color=color_list[j],
                           edgecolor="none", align='center')
                    bottom = [x + y for x, y in zip(bottom, y)]
                ax.set_title(tag)
                ax.yaxis.tick_left()
                ax.set_xticks(range(len(self.bednames)))
                ax.set_xticklabels(self.bednames, fontsize=7, rotation=20, ha="right")
                ax.set_ylabel("Percentage %")
                # ax.tick_params(axis='x', which='both', top='off', bottom='off', labelbottom='on')
                ax.set_ylim([0, 100])
                ax.set_xlim([-0.5, len(self.bednames)-0.5])

            elif i > 0:
                x = [ x - 0.4 for x in range(len(overlapping_counts[i-1])) ]
                ax.bar(left=x, height=overlapping_counts[i-1],
                       color=color_list, linewidth=0)
                ax.set_title(self.bednames[i-1])
                ax.set_ylabel("Number")
                ax.set_xticks([x for x in range(len(overlapping_counts[i-1]))])
                ax.set_xticklabels(refs_names, fontsize=7, rotation=20, ha="right")
                ax.set_xlim([-0.5, len(overlapping_counts[i-1])-0.5])

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
        header_list = ["No", "BED File", "Number", "Length Ave.", "Length s.d.", "Min", "Max", "Internal overlap"]
        c = 0
        for bed in self.bednames:
            c += 1
            data_table.append([ str(c), bed,
                                str(self.stats[bed]["Number of regions"]),
                                "{0:.2f}".format(self.stats[bed]["Average length"]),
                                "{0:.2f}".format(self.stats[bed]["s.d. of length"]),
                                str(self.stats[bed]["min"]),
                                str(self.stats[bed]["max"]),
                                str(self.stats[bed]["Internal overlaps"])])

        html.add_zebra_table(header_list, col_size_list, type_list, data_table,
                             border_list=None, sortable=True, clean=True)
        html.add_figure("figure_"+title+".png", align=50, width=str(300*len(self.ind_col.keys())))
        html.add_fixed_rank_sortable()
        html.write(os.path.join(directory, title, "index.html"))