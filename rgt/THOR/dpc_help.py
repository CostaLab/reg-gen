"""
%prog [options] CONFIG

THOR detects differential peaks in multiple ChIP-seq profiles associated
with two distinct biological conditions.

Copyright (C) 2014-2016  Manuel Allhoff (allhoff@aices.rwth-aachen.de)

This program comes with ABSOLUTELY NO WARRANTY. This is free 
software, and you are welcome to redistribute it under certain 
conditions. Please see LICENSE file for details.
"""

# Python
from __future__ import print_function
import os
import sys
import re
import time
import pysam
import numpy as np
from math import log, ceil
from operator import add
from os.path import splitext, basename, join, isfile, isdir, exists
from optparse import OptionParser, OptionGroup
from datetime import datetime
import shutil

# Internal
from rgt.THOR.postprocessing import merge_delete, filter_deadzones
from MultiCoverageSet import MultiCoverageSet
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.THOR.get_fast_gen_pvalue import get_log_pvalue_new
from input_parser import input_parser
from rgt.Util import which, npath
from rgt import __version__

import configuration

np.random.rand(42)

import matplotlib as mpl
#see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
mpl.use('Agg')
import matplotlib.pyplot as plt


def merge_output(signal_statics, options, no_bw_files, chrom_sizes):
    dim = signal_statics['dim']

    for i in range(dim[0]):
        for j in range(dim[1]):
            ## here are to output bed files for each signal file and output files
            temp_bed = npath(options.name + '-s%s-rep%s_temp.bed'% (i, j))

            files = [options.name + '-' + str(j) + '-s%s-rep%s.bw'%(i, j) for num in no_bw_files]
            if len(no_bw_files) > dim[0]*dim[1]:
                files = filter(lambda x: isfile(x), files)
                t = ['bigWigMerge'] + files + [temp_bed]
                c = " ".join(t)
                os.system(c)

                os.system("LC_COLLATE=C sort -k1,1 -k2,2n " + temp_bed + ' > ' + temp_bed +'.sort')

                t = ['bedGraphToBigWig', temp_bed + '.sort', chrom_sizes, options.name + '-s%s-rep%s.bw' % (i, j)]
                c = " ".join(t)
                os.system(c)

                for f in files:
                    os.remove(f)
                os.remove(temp_bed)
                os.remove(temp_bed + ".sort")
            else:
                ftarget = [options.name + '-s%s-rep%s.bw' %(i, j) for num in no_bw_files]
                for i in range(len(ftarget)):
                    c = ['mv', files[i], ftarget[i]]
                    c = " ".join(c)
                    os.system(c)


def dump_posteriors_and_viterbi(name, posteriors, DCS, states):
    print("Computing info...", file=sys.stderr)
    f = open(name + '-posts.bed', 'w')
    g = open(name + '-states-viterbi.bed', 'w')
    
    for i in range(len(DCS.indices_of_interest)):
        cov1, cov2 = DCS._get_sm_covs(i)
        p1, p2, p3 = posteriors[i][0], posteriors[i][1], posteriors[i][2]
        chrom, start, end = DCS._index2coordinates(DCS.indices_of_interest[i])
        
        print(chrom, start, end, states[i], cov1, cov2, sep='\t', file=g)
        print(chrom, start, end, max(p3, max(p1,p2)), p1, p2, p3, cov1, cov2, sep='\t', file=f)

    f.close()
    g.close()


def _compute_pvalue((x, y, side, distr)):
    a, b = int(np.mean(x)), int(np.mean(y))
    return -get_log_pvalue_new(a, b, side, distr)


def _get_log_ratio(l1, l2):
    l1, l2 = float(np.sum(np.array(l1))), float(np.sum(np.array(l2)))
    try:
        res = l1/l2
    except:
        return sys.maxint
    
    if res > 0:
        try:
            res = log(res)
            if np.isinf(res):
                return sys.maxint
            return res
        except:
            print('error to compute log ratio', l1, l2, file=sys.stderr)
            return sys.maxint
    else:
        return sys.maxint

def _merge_consecutive_bins(tmp_peaks, distr, merge=True):
    """Merge consecutive peaks and compute p-value. Return list
    <(chr, s, e, c1, c2, strand)> and <(pvalue)>"""
    peaks = []
    pvalues = []
    i, j, = 0, 0

    while i < len(tmp_peaks):
        j+=1
        c, s, e, c1, c2, strand, strand_pos, strand_neg = tmp_peaks[i]
        v1 = c1
        v2 = c2

        tmp_pos = [strand_pos]
        tmp_neg = [strand_neg]
        #merge bins
        while merge and i+1 < len(tmp_peaks) and e == tmp_peaks[i+1][1] and strand == tmp_peaks[i+1][5]:
            e = tmp_peaks[i+1][2]
            v1 = map(add, v1, tmp_peaks[i+1][3])
            v2 = map(add, v2, tmp_peaks[i+1][4])
            tmp_pos.append(tmp_peaks[i+1][6])
            tmp_neg.append(tmp_peaks[i+1][7])
            i += 1

        side = 'l' if strand == '+' else 'r'
        pvalues.append((v1, v2, side, distr))

        ratio = _get_log_ratio(tmp_pos, tmp_neg)
        peaks.append((c, s, e, v1, v2, strand, ratio))
        i += 1

    pvalues = map(_compute_pvalue, pvalues)
    assert len(pvalues) == len(peaks)

    return pvalues, peaks


def get_peaks(name, cov_set, states, exts, merge, distr, pcutoff, debug, no_correction, deadzones, merge_bin, p=70):
    """Merge Peaks, compute p-value and give out *.bed and *.narrowPeak"""
    exts = np.mean(exts)
    tmp_peaks = []
    tmp_data = []
    
    for i in range(len(cov_set.indices_of_interest)):
        if states[i] not in [1,2]:
            continue #ignore background states
        idx = cov_set.indices_of_interest[i]
        strand = '+' if states[i] == 1 else '-'
        covs_list, strand_covs_list = cov_set.get_sm_covs(idx, strand_cov=True)  # only return data not indices
        cov1 = covs_list[0]  #np.mean(covs_list[0])  # use avg to present covs
        cov2 = covs_list[1]  #np.mean(covs_list[1])
        
        # cov1_strand = np.sum(DCS.overall_coverage_strand[0][0][:,DCS.indices_of_interest[i]]) + np.sum(DCS.overall_coverage_strand[1][0][:,DCS.indices_of_interest[i]])
        # cov2_strand = np.sum(DCS.overall_coverage_strand[0][1][:,DCS.indices_of_interest[i]] + DCS.overall_coverage_strand[1][1][:,DCS.indices_of_interest[i]])
        ## here we need to set another function to get the cov_strand information
        cov1_strand = np.sum(strand_covs_list[0])
        cov2_strand = np.sum(strand_covs_list[1])
        ## get genome information from idx, to get chrom, end and start
        chrom, start, end = cov_set.sm_index2coordinates(idx)
        
        tmp_peaks.append((chrom, start, end, cov1, cov2, strand, cov1_strand, cov2_strand))
        side = 'l' if strand == '+' else 'r'
        tmp_data.append((np.sum(cov1), np.sum(cov2), side, distr))
    
    if not tmp_data:
        print('no data', file=sys.stderr)
        return [], [], []
    
    tmp_pvalues = map(_compute_pvalue, tmp_data)
    per = np.percentile(tmp_pvalues, p)
    
    tmp = []
    res = tmp_pvalues > per
    for j in range(len(res)):
        if res[j]:
            tmp.append(tmp_peaks[j])
    tmp_peaks = tmp

    pvalues, peaks, = _merge_consecutive_bins(tmp_peaks, distr, merge_bin) #merge consecutive peaks and compute p-value
    regions = merge_delete(exts, merge, peaks, pvalues) #postprocessing, returns GenomicRegionSet with merged regions
    if deadzones:
        regions = filter_deadzones(deadzones, regions)
    output = []
    pvalues = []
    ratios = []
    main_sep = ':' #sep <counts> main_sep <counts> main_sep <pvalue>
    int_sep = ';' #sep counts in <counts>
    
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        counts = ",".join(tmp[0:len(tmp)-1]).replace('], [', int_sep).replace('], ', int_sep).replace('([', '').replace(')', '').replace(', ', main_sep)
        pvalue = float(tmp[len(tmp)-2].replace(")", "").strip())
        ratio = float(tmp[len(tmp)-1].replace(")", "").strip())
        pvalues.append(pvalue)
        ratios.append(ratio)
        output.append((el.chrom, el.initial, el.final, el.orientation, counts))
    
    return ratios, pvalues, output


def _output_ext_data(ext_data_list, bamfiles):
    """Output textfile and png file of read size estimation"""
    names = [splitext(basename(bamfile))[0] for bamfile in bamfiles]

    for k, ext_data in enumerate(ext_data_list):
        f = open(configuration.FOLDER_REPORT_DATA + 'fragment_size_estimate_' + names[k] + '.data', 'w')
        for d in ext_data:
            print(d[0], d[1], sep='\t', file=f)
        f.close()
    
    for i, ext_data in enumerate(ext_data_list):
        d1 = map(lambda x: x[0], ext_data)
        d2 = map(lambda x: x[1], ext_data)
        ax = plt.subplot(111)
        plt.xlabel('shift')
        plt.ylabel('convolution')
        plt.title('Fragment Size Estimation')
        plt.plot(d2, d1, label=names[i])
    
    ax.legend()
    plt.savefig(configuration.FOLDER_REPORT_PICS + 'fragment_size_estimate.png')
    plt.close()


def get_all_chrom(bamfiles):
    chrom = set()
    for bamfile in bamfiles:
        bam = pysam.Samfile(bamfile, "rb" )
        for read in bam.fetch():
            c = bam.getrname(read.reference_id)
            if c not in chrom:
                chrom.add(c)
    return chrom


def initialize(options, strand_cov, genome_path, regionset, mask_file, signal_statics, inputs_statics,
               tracker, counter,test, end, output_bw=True):
    """Initialize the MultiCoverageSet
    Region_giver includes: regions to be analysed + regions to be masked + chrom_sizes file name + chrom_sizes_dict
    Use sampling methods to initialize certain part of data
    Get exp_data, we have one pattern, 0.1
    """
    if options.norm_regions:
        norm_regionset = GenomicRegionSet('norm_regions')
        norm_regionset.read(options.norm_regions)
    else:
        norm_regionset = None
    # options.binsize = 1000
    # options.stepsize = 500
    print("Begin reading", file=sys.stderr)
    start = time.time()
    cov_set = MultiCoverageSet(name=options.name, regionset=regionset,mask_file=mask_file, binsize=options.binsize, stepsize=options.stepsize, rmdup=options.rmdup, signal_statics=signal_statics, inputs_statics=inputs_statics,
                                     verbose=options.verbose, debug=options.debug, norm_regionset=norm_regionset, save_wig=options.save_wig, strand_cov=strand_cov,
                                     tracker=tracker, end=end, counter=counter, output_bw=output_bw,
                                     folder_report=configuration.FOLDER_REPORT, report=options.report, save_input=options.save_input, use_sm=True)

    elapsed_time = time.time() - start
    print("End reading using time %.3f s"%(elapsed_time), file=sys.stderr)

    options.no_gc_content = True
    if not options.no_gc_content and genome_path: # maybe we could use samples to get values not all data;; samples from indices, and around 1000 for it
        start = time.time()
        options.gc_hv = None
        options.gc_avg_T = None
        options.gc_hv, options.gc_avg_T = cov_set.normalization_by_gc_content(inputs_statics, genome_path, options.gc_hv, options.gc_avg_T, delta=0.01)
        elapsed_time = time.time() - start
        if configuration.VERBOSE:
            print("Compute GC-content using time %.3f s"%(elapsed_time), file=sys.stderr)

        # we need to save values for it for return ?? If we use another; [avg_T, hv] for each inputs files.

    cov_set.init_overall_coverage(strand_cov=strand_cov)

    if inputs_statics: # only inputs_statics exist we do it;
        start = time.time()
        options.factors_inputs = cov_set.normalization_by_input(signal_statics, inputs_statics, options.name, options.factors_inputs)
        elapsed_time = time.time() - start
        if configuration.VERBOSE:
            print("Normalize input-DNA using time %.3f s"%(elapsed_time), file=sys.stderr)

    if options.save_input: # !!! sth changes about parameters
        cov_set.output_input_bw(options.name, regionset, options.save_wig)
    # much complex, so we decay to change it
    start = time.time()
    options.scaling_factors_ip = cov_set.normalization_by_signal(options.name, options.scaling_factors_ip, signal_statics, options.housekeeping_genes, tracker, norm_regionset,
                                    options.report, options.m_threshold, options.a_threshold)
    elapsed_time = time.time() - start
    if configuration.VERBOSE:
        print('Normalize ChIP-seq profiles using time %.3f s' % (elapsed_time), file=sys.stderr)

    ## After this step, we have already normalized data, so we could output normalization data
    if output_bw:  # !!!!  sth changes about parameters
        cov_set._output_bw(options.name, regionset, options.save_wig, options.save_input)

    return cov_set


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def confirm(prompt=None, resp=False):
    """prompts for yes or no response from the user. Returns True for yes and
    False for no.

    'resp' should be set to the default value assumed by the caller when
    user simply types ENTER.

    >>> confirm(prompt='Create Directory?', resp=True)
    Create Directory? [y]|n:
    True
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y:
    False
    >>> confirm(prompt='Create Directory?', resp=False)
    Create Directory? [n]|y: y
    True

    """

    if prompt is None:
        prompt = 'Confirm'

    if resp:
        prompt = '%s [%s]|%s: ' % (prompt, 'y', 'n')
    else:
        prompt = '%s [%s]|%s: ' % (prompt, 'n', 'y')

    while True:
        ans = raw_input(prompt)
        if not ans:
            return resp
        if ans not in ['y', 'Y', 'n', 'N']:
            print('please enter y or n.')
            continue
        if ans == 'y' or ans == 'Y':
            return True
        if ans == 'n' or ans == 'N':
            return False


def _callback_list(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: int(x), value.split(',')))


def _callback_list_float(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: float(x), value.split(',')))


def handle_input():
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
    group.add_option("--size", dest="size_ts", default=1000, type="int",
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
    options.verbose = True
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

    bamfiles, genome, chrom_sizes, inputs, dims = input_parser(config_path)

    if not genome:
        options.no_gc_content = True

    if options.exts and len(options.exts) != len(bamfiles):
        parser.error("Number of Extension Sizes must equal number of bamfiles")

    if options.exts_inputs and len(options.exts_inputs) != len(inputs):
        parser.error("Number of Input Extension Sizes must equal number of input bamfiles")

    if options.scaling_factors_ip and len(options.scaling_factors_ip) != len(bamfiles):
        parser.error("Number of scaling factors for IP must equal number of bamfiles")

    for each_sample in bamfiles:
        for bamfile in each_sample:
            if not isfile(bamfile):
                parser.error("BAM file %s does not exist!" % bamfile)

    if not inputs and options.factors_inputs:
        print("As no input-DNA, do not use input-DNA factors", file=sys.stderr)
        options.factors_inputs = None

    if options.factors_inputs and len(options.factors_inputs) != len(bamfiles):
        parser.error("factors for input-DNA must equal number of BAM files!")

    if inputs:
        for each_sample in inputs:
            for bamfile in each_sample:
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
        # if exist then we judge if there exists one file with peak amd if it's then we save it;
        # else, we will delete the files
        if isdir(options.outputdir):
            if np.any(map(lambda x: re.search(r".*diffpeaks.bed$", x),os.listdir(options.outputdir))):
                if confirm(prompt="delete existing results?", resp=True):
                    shutil.rmtree(options.outputdir)
                else:
                    parser.error("Output directory exists and contains files with names starting with your chosen experiment name! "
                                 "Do nothing to prevent file overwriting!")
            else:
                shutil.rmtree(options.outputdir)

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

    configuration.FOLDER_REPORT = join(options.outputdir, 'report_'+basename(options.name)+"/")
    configuration.FOLDER_REPORT_PICS = join(options.outputdir, 'report_'+basename(options.name), 'pics/')
    configuration.FOLDER_REPORT_DATA = join(options.outputdir, 'report_'+basename(options.name), 'pics/data/')
    configuration.OUTPUTDIR = options.outputdir
    configuration.NAME = options.name

    if not inputs:
        print("Warning: Do not compute GC-content, as there is no input file", file=sys.stderr)

    if not genome:
        print("Warning: Do not compute GC-content, as there is no genome file", file=sys.stderr)

    if options.exts is None:
        options.exts = []

    if options.exts_inputs is None:
        options.exts_inputs = []

    return options, bamfiles, genome, chrom_sizes, dims, inputs
