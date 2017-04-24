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
    def __init__(self, input, organism):
        self.EM = ExperimentalMatrix()
        self.EM.read(input)
        self.beds = self.EM.get_regionsets()
        self.bednames = self.EM.get_regionsnames()
        self.organism = organism
        self.chromosomes = GenomicRegionSet(organism)
        self.chromosomes.get_genome_data(organism=organism,chrom_X=True)
        self.stats = OrderedDict()
        self.fig_size = (6, 6)
        self.figures = OrderedDict()
        for bed in self.bednames:
            self.figures[bed] = []


    def cal_statistics(self):
        for i, bed in enumerate(self.beds):
            self.stats[self.bednames[i]] = {}
            self.stats[self.bednames[i]]["Number of regions"] = len(bed)
            self.stats[self.bednames[i]]["Average length"] = bed.average_size()
            self.stats[self.bednames[i]]["s.d. of length"] = bed.size_variance()
            self.stats[self.bednames[i]]["Internal overlaps"] = len(bed) - len(bed.merge())

    def plot_distribution_length(self, outdir):
        for i, bed in enumerate(self.beds):
            dis = [len(r) for r in bed]
            plot_f = "len_dis_" + self.bednames[i] + ".pdf"
            self.figures[self.bednames[i]].append(plot_f)
            # Plot the distribution

            with PdfPages(os.path.join(outdir, plot_f)) as pdf:
                plt.figure(figsize=self.fig_size)
                plt.hist(dis)
                plt.title(self.bednames[i])
                plt.xlabel("Length (bp)")
                plt.ylabel("Frequency")
                # Draw the mean

                pdf.savefig()  # saves the current figure into a pdf page
                plt.close()


    def plot_distribution_chromosome(self, outdir):

        # Convert genomic loci into coordinate

        for i, bed in enumerate(self.beds):
            # Plot the distribution
            # Save
            os.path.join(outdir, "dis_chrom_" + self.bednames[i] + ".pdf")

    def plot_biotype(self, outdir):


    def


