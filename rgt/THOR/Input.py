# Python
from __future__ import print_function
import os
import sys
from os.path import basename, join, isfile, isdir, exists
from optparse import OptionGroup, OptionParser
from datetime import datetime

# Internal

# from .input_parser import input_parser
from ..Util import which, npath
from .. import __version__



def _callback_list(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: int(x), value.split(',')))


def _callback_list_float(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: float(x), value.split(',')))



class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def InputParsers():
    parser = HelpfulOptionParser(usage=__doc__)

    parser.add_option("-n", "--name", default=None, dest="name", type="string",
                      help="Experiment's name and prefix for all files that are created.")
    parser.add_option("-m", "--merge", default=False, dest="merge", action="store_true",
                      help="Merge peaks which have a distance less than the estimated mean fragment size "
                           "(recommended for histone data). [default: do not merge]")
    parser.add_option("--no-merge-bin", default=True, dest="merge_bin", action="store_false",
                      help="Merge the overlapping bin before filtering by p-value."
                           "[default: Merging bins]")
    parser.add_option("--housekeeping-genes", default=None, dest="housekeeping_genes", type="str",
                      help="Define housekeeping genes (BED format) used for normalizing. [default: %default]")
    parser.add_option("--output-dir", dest="outputdir", default=None, type="string",
                      help="Store files in output directory. [default: %default]")
    parser.add_option("--report", dest="report", default=False, action="store_true",
                      help="Generate HTML report about experiment. [default: %default]")
    parser.add_option("--deadzones", dest="deadzones", default=None,
                      help="Define blacklisted genomic regions avoided for analysis (BED format). [default: %default]")
    parser.add_option("--no-correction", default=False, dest="no_correction", action="store_true",
                      help="Do not use multipe test correction for p-values (Benjamini/Hochberg). [default: %default]")
    parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.1, type="float",
                      help="P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. "
                           "[default: %default]")
    parser.add_option("--exts", default=None, dest="exts", type="str", action='callback', callback=_callback_list,
                      help="Read's extension size for BAM files (comma separated list for each BAM file in config "
                           "file). If option is not chosen, estimate extension sizes. [default: %default]")
    parser.add_option("--factors-inputs", default=None, dest="factors_inputs", type="str", action="callback",
                      callback=_callback_list_float,
                      help="Normalization factors for input-DNA (comma separated list for each BAM file in config "
                           "file). If option is not chosen, estimate factors. [default: %default]")
    parser.add_option("--scaling-factors", default=None, dest="scaling_factors_ip", type="str", action='callback',
                      callback=_callback_list_float,
                      help="Scaling factor for each BAM file (not control input-DNA) as comma separated list for "
                           "each BAM file in config file. If option is not chosen, follow normalization strategy "
                           "(TMM or HK approach) [default: %default]")
    parser.add_option("--save-input", dest="save_input", default=False, action="store_true",
                      help="Save input-DNA file if available. [default: %default]")
    parser.add_option("--version", dest="version", default=False, action="store_true",
                      help="Show script's version.")

    group = OptionGroup(parser, "Advanced options")
    group.add_option("--regions", dest="regions", default=None, type="string",
                     help="Define regions (BED format) to restrict the analysis, that is, where to train the HMM and "
                          "search for DPs. It is faster, but less precise.")
    group.add_option("-b", "--binsize", dest="binsize", default=100, type="int",
                     help="Size of underlying bins for creating the signal. [default: %default]")
    group.add_option("-s", "--step", dest="stepsize", default=50, type="int",
                     help="Stepsize with which the window consecutively slides across the genome to create the "
                          "signal. [default: %default]")
    group.add_option("--debug", default=False, dest="debug", action="store_true",
                     help="Output debug information. Warning: space consuming! [default: %default]")
    group.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true",
                     help="Do not normalize towards GC content. [default: %default]")
    group.add_option("--norm-regions", default=None, dest="norm_regions", type="str",
                     help="Restrict normalization to particular regions (BED format). [default: %default]")
    group.add_option("-f", "--foldchange", dest="foldchange", default=1.6, type="float",
                     help="Fold change parameter to define training set (t_1, see paper). [default: %default]")
    group.add_option("-t", "--threshold", dest="threshold", default=95, type="float",
                     help="Minimum signal support for differential peaks to define training set as percentage "
                          "(t_2, see paper). [default: %default]")
    group.add_option("--size", dest="size_ts", default=10000, type="int",
                     help="Number of bins the HMM's training set constists of. [default: %default]")
    group.add_option("--par", dest="par", default=1, type="int",
                     help="Percentile for p-value postprocessing filter. [default: %default]")
    group.add_option("--poisson", default=False, dest="poisson", action="store_true",
                     help="Use binomial distribution as emmission. [default: %default]")
    group.add_option("--single-strand", default=False, dest="singlestrand", action="store_true",
                     help="Allow single strand BAM file as input. [default: %default]")
    group.add_option("--m_threshold", default=80, dest="m_threshold", type="int",
                     help="Define the M threshold of percentile for training TMM. [default: %default]")
    group.add_option("--a_threshold", default=95, dest="a_threshold", type="int",
                     help="Define the A threshold of percentile for training TMM. [default: %default]")
    group.add_option("--rmdup", default=False, dest="rmdup", action="store_true",
                     help="Remove the duplicate reads [default: %default]")
    parser.add_option_group(group)

    (options, args) = parser.parse_args()
    options.save_wig = False
    options.exts_inputs = None
    options.verbose = False
    options.hmm_free_para = False

    if options.version:
        print("")
        print(__version__)
        sys.exit()

    if len(args) != 1:
        parser.error("Please give config file")

    config_path = npath(args[0])

    if not isfile(config_path):
        parser.error("Config file %s does not exist!" % config_path)

    bamfiles, genome, chrom_sizes, inputs, dims = input_files_parser(config_path)

    if not genome:
        options.no_gc_content = True

    if options.exts and len(options.exts) != len(bamfiles):
        parser.error("Number of Extension Sizes must equal number of bamfiles")

    if options.exts_inputs and len(options.exts_inputs) != len(inputs):
        parser.error("Number of Input Extension Sizes must equal number of input bamfiles")

    if options.scaling_factors_ip and len(options.scaling_factors_ip) != len(bamfiles):
        parser.error("Number of scaling factors for IP must equal number of bamfiles")

    for bamfile in bamfiles:
        if not isfile(bamfile):
            parser.error("BAM file %s does not exist!" % bamfile)

    if not inputs and options.factors_inputs:
        print("As no input-DNA, do not use input-DNA factors", file=sys.stderr)
        options.factors_inputs = None

    if options.factors_inputs and len(options.factors_inputs) != len(bamfiles):
        parser.error("factors for input-DNA must equal number of BAM files!")

    if inputs:
        for bamfile in inputs:
            if not isfile(bamfile):
                parser.error("BAM file %s does not exist!" % bamfile)

    if options.regions:
        if not isfile(options.regions):
            parser.error("Region file %s does not exist!" % options.regions)

    if genome and not isfile(genome):
        parser.error("Genome file %s does not exist!" % genome)

    if options.name is None:
        d = str(datetime.now()).replace("-", "_").replace(":", "_").replace(" ", "_").replace(".", "_").split("_")
        options.name = "THOR-exp" + "-" + "_".join(d[:len(d) - 1])

    if not which("wigToBigWig") or not which("bedGraphToBigWig") or not which("bigWigMerge"):
        print("Warning: wigToBigWig, bigWigMerge or bedGraphToBigWig not found! Signal will not be stored!",
              file=sys.stderr)

    if options.outputdir:
        options.outputdir = npath(options.outputdir)
        if isdir(options.outputdir) and sum(
                map(lambda x: x.startswith(options.name), os.listdir(options.outputdir))) > 0:
            parser.error("Output directory exists and contains files with names starting with your chosen experiment "
                         "name! Do nothing to prevent file overwriting!")
        if not exists(options.outputdir):
            os.mkdir(options.outputdir)
    else:
        options.outputdir = os.getcwd()

    options.name = join(options.outputdir, options.name)

    if options.report and isdir(join(options.outputdir, 'report_'+basename(options.name))):
        parser.error("Folder 'report_"+basename(options.name)+"' already exits in output directory!" 
                     "Do nothing to prevent file overwriting! "
                     "Please rename report folder or change working directory of THOR with the option --output-dir")

    if options.report:
        os.mkdir(join(options.outputdir, 'report_'+basename(options.name)+"/"))
        os.mkdir(join(options.outputdir, 'report_'+basename(options.name), 'pics/'))
        os.mkdir(join(options.outputdir, 'report_'+basename(options.name), 'pics/data/'))

    global FOLDER_REPORT
    global FOLDER_REPORT_PICS
    global FOLDER_REPORT_DATA
    global OUTPUTDIR
    global NAME

    FOLDER_REPORT = join(options.outputdir, 'report_'+basename(options.name)+"/")
    FOLDER_REPORT_PICS = join(options.outputdir, 'report_'+basename(options.name), 'pics/')
    FOLDER_REPORT_DATA = join(options.outputdir, 'report_'+basename(options.name), 'pics/data/')
    OUTPUTDIR = options.outputdir
    NAME = options.name

    if not inputs:
        print("Warning: Do not compute GC-content, as there is no input file", file=sys.stderr)

    if not genome:
        print("Warning: Do not compute GC-content, as there is no genome file", file=sys.stderr)

    if options.exts is None:
        options.exts = []

    if options.exts_inputs is None:
        options.exts_inputs = []

    return options, bamfiles, genome, chrom_sizes, dims, inputs


def get_data_block(filepath, feature):
    with open(filepath) as f:
        data = []
        read = False
        for line in f:
            line = line.strip()
            if line == '':
                continue
            if line.startswith("#") and line == "#" + str(feature):
                read = True

            if line.startswith("#") and line != "#" + str(feature):
                read = False

            if not line.startswith("#") and read:
                data.append(line)

    if len(data) == 1 and not (feature == "rep1" or feature == "rep2" or feature == "inputs1" or feature == "inputs2"):
        return data[0]
    else:
        return data


def input_files_parser(filepath):
    bamfiles_1 = get_data_block(filepath, "rep1")
    bamfiles_1 = map(npath, bamfiles_1)

    bamfiles_2 = get_data_block(filepath, "rep2")
    bamfiles_2 = map(npath, bamfiles_2)

    # genome is optional, so if we get an empty list
    # we set it to None, otherwise we normalise the path
    genome = get_data_block(filepath, "genome")
    genome = npath(genome) if genome else None

    # the chrom sizes are not optional, but right now it's undefined
    # what happens if the user doesn't specify them, or specifies more
    # than one. So we just relay whatever we got from the file.
    chrom_sizes = npath(get_data_block(filepath, "chrom_sizes"))
    chrom_sizes = npath(chrom_sizes) if chrom_sizes else chrom_sizes

    inputs1 = get_data_block(filepath, "inputs1")
    inputs1 = map(npath, inputs1)

    inputs2 = get_data_block(filepath, "inputs2")
    inputs2 = map(npath, inputs2)

    dims = [len(bamfiles_1), len(bamfiles_2)]

    if not inputs1 and not inputs2:
        inputs = None
    else:
        inputs = inputs1 + inputs2

    return bamfiles_1 + bamfiles_2, genome, chrom_sizes, inputs, dims

