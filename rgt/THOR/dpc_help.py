"""
%prog [CONFIG]

Find differential peaks in regions.

Caution: Pre-paper version

Author: Manuel Allhoff (allhoff@aices.rwth-aachen.de)

"""

from __future__ import print_function
from optparse import OptionParser, OptionGroup
from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.ODIN.get_extension_size import get_extension_size
import os.path
import sys
from MultiCoverageSet import MultiCoverageSet
#from get_gen_pvalue import get_log_pvalue
from rgt.ODIN.get_fast_gen_pvalue import get_log_pvalue_new
from math import log, log10
import multiprocessing
from input_parser import input_parser
#from rgt.ODIN import ODIN
import matplotlib as mpl #necessary to plot without x11 server (for cluster)
mpl.use('Agg')           #see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
import matplotlib.pyplot as plt
from random import sample
from scipy.optimize import curve_fit
import numpy as np
from numpy import linspace
from math import fabs
#from rgt.ODIN.postprocessing import merge_delete
from rgt.THOR.postprocessing import merge_delete
#from rgt.ODIN.dpc_help import _output_BED, _output_narrowPeak
from rgt.THOR.neg_bin import NegBin
from operator import add
from numpy import percentile
from norm_genelevel import norm_gene_level
from datetime import datetime
from rgt.ODIN.dpc_help import which
import pysam

def _func_quad_2p(x, a, c):
    """Return y-value of y=max(|a|*x^2 + x + |c|, 0),
    x may be an array or a single float"""
    res = []
    if type(x) is np.ndarray:
        for el in x:
            res.append(max(el, fabs(a) * el**2 + el + fabs(c)))
            
        return np.asarray(res)
    else:
        return max(x, fabs(a) * x**2 + x + fabs(c))

def _write_emp_func_data(data, name):
    """Write mean and variance data"""
    assert len(data[0]) == len(data[1])
    f = open(FOLDER_REPORT_DATA + name + '.data', 'w')
    for i in range(len(data[0])):
        print(data[0][i], data[1][i], sep='\t', file=f)
    f.close()
    

def _plot_func(plot_data, outputdir):
    """Plot estimated and empirical function"""

    maxs = [] #max for x (mean), max for y (var)
    for i in range(2): 
        tmp = np.concatenate((plot_data[0][i], plot_data[1][i])) #plot_data [(m, v, p)], 2 elements
        maxs.append(max(tmp[tmp < np.percentile(tmp, 90)]))

    for i in range(2):
        x = linspace(0, max(plot_data[i][0]), max(plot_data[i][0])+1)
        y = _func_quad_2p(x, plot_data[i][2][0], plot_data[i][2][1])
        
        for j in range(2):
            #use matplotlib to plot function and datapoints
            #and save datapoints to files
            ext = 'original'
            if j == 1:
                plt.xlim([0, maxs[0]])
                plt.ylim([0, maxs[1]])
                ext = 'norm'
            ax = plt.subplot(111)
            plt.plot(x, y, 'r', label = 'empirical datapoints') #plot polynom
            plt.scatter(plot_data[i][0], plot_data[i][1], label = 'fitted polynomial') #plot datapoints
            ax.legend()
            plt.xlabel('mean')
            plt.ylabel('variance')
            plt.title('Estimated Mean-Variance Function')
            name = "_".join(['mean', 'variance', 'func', 'cond', str(i), ext])
            _write_emp_func_data(plot_data[i], name)
            plt.savefig(FOLDER_REPORT_PICS + name + '.png')
            plt.close()

def _get_data_rep(overall_coverage, name, debug, sample_size):
    """Return list of (mean, var) points for samples 0 and 1"""
    data_rep = []
    for i in range(2):
        cov = np.asarray(overall_coverage[i]) #matrix: (#replicates X #bins)
        h = np.invert((cov==0).all(axis=0)) #assign True to columns != (0,..,0)
        cov = cov[:,h] #remove 0-columns
        r = np.random.randint(cov.shape[1], size=sample_size)
        r.sort()
        cov = cov[:,r]
             
        m = list(np.squeeze(np.asarray(np.mean(cov*1.0, axis=0))))
        n = list(np.squeeze(np.asarray(np.var(cov*1.0, axis=0))))
        
        assert len(m) == len(n)

        data_rep.append(zip(m, n))
        data_rep[i].append((0,0))
        data_rep[i] = np.asarray(data_rep[i])
        
    if debug:
        for i in range(2):
            np.save(str(name) + "-emp-data" + str(i) + ".npy", data_rep[i])
    
    for i in range(2):
        #print("percentile 99.75", np.percentile(data_rep[i][:,0], 99.75), file=sys.stderr)
        #print("percentile 99.75", np.percentile(data_rep[i][:,1], 99.75), file=sys.stderr)
        data_rep[i] = data_rep[i][data_rep[i][:,0] < np.percentile(data_rep[i][:,0], 99.75)]
        data_rep[i] = data_rep[i][data_rep[i][:,1] < np.percentile(data_rep[i][:,1], 99.75)]
    
    return data_rep
    
def _fit_mean_var_distr(overall_coverage, name, debug, verbose, outputdir, report, sample_size=10000):
    """Estimate empirical distribution (quadr.) based on empirical distribution"""
    done = False
    plot_data = [] #means, vars, paras
    
    while not done:
        data_rep = _get_data_rep(overall_coverage, name, debug, sample_size)
        res = []
        for i in range(2):
            try:
                m = np.asarray(map(lambda x: x[0], data_rep[i])) #means list
                v = np.asarray(map(lambda x: x[1], data_rep[i])) #vars list
                p, _ = curve_fit(_func_quad_2p, m, v) #fit quad. function to empirical data
                res.append(p)
                plot_data.append((m, v, p))
                if i == 1:
                    done = True
            except RuntimeError:
                print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                break #restart for loop
    
    if report:
        _plot_func(plot_data, outputdir)
                
    return lambda x: _func_quad_2p(x, p[0], p[1]), res

def dump_posteriors_and_viterbi(name, posteriors, DCS, states):
    print("Computing info...", file=sys.stderr)
    f = open(name + '-posts.bed', 'w')
    g = open(name + '-states-viterbi.bed', 'w')
    
    for i in range(len(DCS.indices_of_interest)):
        cov1, cov2 = _get_covs(DCS, i)
        p1, p2, p3 = posteriors[i][0], posteriors[i][1], posteriors[i][2]
        chrom, start, end = DCS._index2coordinates(DCS.indices_of_interest[i])
        
        print(chrom, start, end, states[i], cov1, cov2, sep = '\t', file=g)
        print(chrom, start, end, max(p3, max(p1,p2)), p1, p2, p3, cov1, cov2, sep = '\t', file=f)

    f.close()
    g.close()

lookup_pvalues = {}
def _compute_pvalue((x, y, side, distr)):
    #return 0.2
#     var =  np.var( x + y )
#     mu = np.mean( x + y )
#     alpha = max((var - mu) / np.square(mu), 0.00000000001)
#     m = NegBin(mu, alpha)
#     distr = {'distr_name': 'nb', 'distr': m}
#     
#     a, b = int(np.mean(x)), int(np.mean(y))
#     k = (a, b, mu, alpha)
#     if not lookup_pvalues.has_key(k):
#         lookup_pvalues[k] = -get_log_pvalue_new(a, b, side, distr)
#     
#     return lookup_pvalues[k]
    
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

def _merge_consecutive_bins(tmp_peaks, distr):
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
        while i+1 < len(tmp_peaks) and e == tmp_peaks[i+1][1] and strand == tmp_peaks[i+1][5]:
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
    

def _get_covs(DCS, i, as_list=False):
    """For a multivariant Coverageset, return mean coverage cov1 and cov2 at position i"""
    if not as_list:
        cov1 = int(np.mean(DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]))
    else:
        cov1 = DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]
        cov1 = map(lambda x: x[0], np.asarray((cov1)))
        cov2 = DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]
        cov2 = map(lambda x: x[0], np.asarray((cov2)))
    
    return cov1, cov2

def get_peaks(name, DCS, states, exts, merge, distr, pcutoff, debug, p=70):
    """Merge Peaks, compute p-value and give out *.bed and *.narrowPeak"""
    exts = np.mean(exts)
    tmp_peaks = []
    tmp_data = []
    
    for i in range(len(DCS.indices_of_interest)):
        if states[i] not in [1,2]:
            continue #ignore background states

        strand = '+' if states[i] == 1 else '-'
        cov1, cov2 = _get_covs(DCS, i, as_list=True)
        
        cov1_strand = np.sum(DCS.overall_coverage_strand[0][0][:,DCS.indices_of_interest[i]]) + np.sum(DCS.overall_coverage_strand[1][0][:,DCS.indices_of_interest[i]])
        cov2_strand = np.sum(DCS.overall_coverage_strand[0][1][:,DCS.indices_of_interest[i]] + DCS.overall_coverage_strand[1][1][:,DCS.indices_of_interest[i]])
        
        chrom, start, end = DCS._index2coordinates(DCS.indices_of_interest[i])
        
        tmp_peaks.append((chrom, start, end, cov1, cov2, strand, cov1_strand, cov2_strand))
        side = 'l' if strand == '+' else 'r'
        tmp_data.append((sum(cov1), sum(cov2), side, distr))
    
    tmp_pvalues = map(_compute_pvalue, tmp_data)
    per = np.percentile(tmp_pvalues, p)
    
    if debug:
        print('percentile for peak calling:', per, file=sys.stderr)
    
    tmp = []
    res = tmp_pvalues > per
    for j in range(len(res)):
        if res[j]:
            tmp.append(tmp_peaks[j])
    tmp_peaks = tmp

    #merge consecutive peaks and compute p-value
    pvalues, peaks, = _merge_consecutive_bins(tmp_peaks, distr)
    #postprocessing, returns GenomicRegionSet with merged regions
    regions = merge_delete(exts, merge, peaks, pvalues)
    #regions = merge_delete([0], False, peaks, pvalues) 
    output = []
    pvalues = []
    main_sep = ':' #sep <counts> main_sep <counts> main_sep <pvalue>
    int_sep = ';' #sep counts in <counts>
    
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        counts = ",".join(tmp[0:len(tmp)-1]).replace('], [', int_sep).replace('], ', int_sep).replace('([', '').replace(')', '').replace(', ', main_sep)
        pvalue = float(tmp[len(tmp)-2].replace(")", "").strip())
        ratio = float(tmp[len(tmp)-1].replace(")", "").strip())
        pvalues.append(pvalue)
        output.append((el.chrom, el.initial, el.final, el.orientation, counts, ratio))
    
    pcutoff = -log10(pcutoff)
    pv_pass = np.where(np.asarray(pvalues) >= pcutoff, True, False)

    output = np.array(output)
    output = output[pv_pass]
    pvalues = list(np.array(pvalues)[pv_pass])
    
    output = map(lambda x: tuple(x), list(output))
    
    assert(len(output) == len(pvalues))
    
    _output_BED(name, output, pvalues)
    #_output_narrowPeak(name, output, pvalues, strands_pos, strands_neg)

def _output_BED(name, output, pvalues):
    f = open(name + '-diffpeaks.bed', 'w')
     
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, strand, counts, ratio = output[i]
        print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, colors[strand], 0, counts, ratio, sep='\t', file=f)
    
    f.close()

def _output_narrowPeak(name, output, pvalues):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name + '-diffpeaks.narrowPeak', 'w')
    for i in range(len(pvalues)):
        c, s, e, strand, _ = output[i]
        print(c, s, e, 'Peak' + str(i), 0, strand, 0, pvalues[i], 0, -1, sep='\t', file=f)
    f.close()

#def _output_ext_data(ext_data, bamfile):
#    #write data
#    name = os.path.basename(os.path.splitext(a)[0])
#    f = open(FOLDER_REPORT_DATA +)
    
#    f = open(FOLDER_REPORT_DATA + + name + '.data', 'w')
#    for i in range(len(data[0])):
#        print(data[0][i], data[1][i], sep='\t', file=f)
#    f.close()
    

def _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, verbose):
    """Compute Extension sizes for bamfiles and input files"""
    start = 0
    end = 600
    ext_stepsize = 5

    #compute extension size
    if not exts:
        print("Computing read extension sizes for ChIP-DNA...", file=sys.stderr)
        for bamfile in bamfiles:
            e, ext_data = get_extension_size(bamfile, start=start, end=end, stepsize=ext_stepsize)
            #_output_ext_data(ext_data, bamfile)
            exts.append(e)
        #print(" ".join(exts), file=sys.stderr)

    if inputs and not exts_inputs:
        #print("Computing read extension sizes for input-DNA...", file=sys.stderr)
        #for inp in inputs:
        #    e, _ = get_extension_size(inp, start=start, end=end, stepsize=ext_stepsize)
        #    exts_inputs.append(e)
        exts_inputs = [5] * len(inputs)
        #print(" ".join(exts_inputs), file=sys.stderr)

    return exts, exts_inputs

def get_all_chrom(bamfiles):
    chrom = set()
    for bamfile in bamfiles:
        bam = pysam.Samfile(bamfile, "rb" )
        for read in bam.fetch():
            c = bam.getrname(read.reference_id)
            if c not in chrom:
                chrom.add(c)
    return chrom

def initialize(name, dims, genome_path, regions, stepsize, binsize, bamfiles, exts, \
               inputs, exts_inputs, factors_inputs, chrom_sizes, verbose, no_gc_content, \
               tracker, debug, norm_regions, scaling_factors_ip, save_wig, housekeeping_genes, test):
    """Initialize the MultiCoverageSet"""
    
    regionset = GenomicRegionSet(name)
    chrom_sizes_dict = {}
    #if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
#     if test:
#         contained_chrom = ['chr1', 'chr2']
#     else:
#         #contained_chrom = get_all_chrom(bamfiles)
#         contained_chrom = ['chr1', 'chr2']
    
    if regions is not None:
        print("Call DPs on specified regions.", file=sys.stderr)
        with open(regions) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                c, s, e = line[0], int(line[1]), int(line[2])
		#if c in contained_chrom:                
                regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
                chrom_sizes_dict[c] = e
    else:
        print("Call DPs on whole genome.", file=sys.stderr)
        with open(chrom_sizes) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                chrom, end = line[0], int(line[1])
		#if chrom in contained_chrom:
                regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
                chrom_sizes_dict[chrom] = end
    
    if not regionset.sequences:
        print('something wrong here', file=sys.stderr)
        sys.exit(2)
    
    if norm_regions:
        norm_regionset = GenomicRegionSet('norm_regions')
        norm_regionset.read_bed(norm_regions)
    else:
        norm_regionset = None
        
    
    regionset.sequences.sort()
    exts, exts_inputs = _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, verbose)
    tracker.write(text=str(exts).strip('[]'), header="Extension size (rep1, rep2, input1, input2)")
    
    multi_cov_set = MultiCoverageSet(name=name, regions=regionset, dims=dims, genome_path=genome_path, binsize=binsize, stepsize=stepsize,rmdup=True,\
                                  path_bamfiles = bamfiles, path_inputs = inputs, exts = exts, exts_inputs = exts_inputs, factors_inputs = factors_inputs, \
                                  chrom_sizes=chrom_sizes, verbose=verbose, no_gc_content=no_gc_content, chrom_sizes_dict=chrom_sizes_dict, debug=debug, \
                                  norm_regionset=norm_regionset, scaling_factors_ip=scaling_factors_ip, save_wig=save_wig, strand_cov=True,
                                  housekeeping_genes=housekeeping_genes, tracker=tracker)
    
    return multi_cov_set


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def _callback_list(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: int(x), value.split(',')))

def _callback_list_float(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: float(x), value.split(',')))

def input(laptop):
    parser = HelpfulOptionParser(usage=__doc__)
    if laptop:
        print("---------- TEST MODE ----------", file=sys.stderr)
        (options, args) = parser.parse_args()
        config_path = '/home/manuel/workspace/eclipse/office_share/blueprint/playground/input_test'
        args.append('')
        args[0] = config_path
        #config_path = '/home/manuel/workspace/eclipse/office_share/simulator/test.config'
        bamfiles, genome, chrom_sizes, inputs, dims = input_parser(config_path)
        options.regions = '/home/manuel/r.bed'
        options.exts = [200, 200, 200, 200, 200]
        options.exts_inputs = [200, 200, 200, 200, 200]
        options.pcutoff = 1
        options.name='test'
        options.merge=True
        options.stepsize=50
        options.binsize=100
        options.save_wig = False
        options.par = 70
        options.factors_inputs = None
        options.verbose = True
        options.no_gc_content = False
        options.debug = True
        options.norm_regions = None #'/home/manuel/data/testdata/norm_regions.bed'
        options.scaling_factors_ip = False
        options.housekeeping_genes = False
        options.distr='negbin'
        options.version = None
        options.outputdir = None #'/home/manuel/test/'
        options.report = False
    else:
        parser.add_option("-n", "--name", default=None, dest="name", type="string",\
                          help="Experiment's name and prefix for all files that are created.")
        parser.add_option("-m", "--merge", default=False, dest="merge", action="store_true", \
                          help="Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data). [default: %default]")
        parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.1, type="float",\
                          help="P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: %default]")
        parser.add_option("--exts", default=None, dest="exts", type="str", action='callback', callback=_callback_list,\
                          help="Read's extension size for BAM files. If option is not chosen, estimate extension sizes. [default: %default]")
        parser.add_option("--factors-inputs", default=None, dest="factors_inputs", type="str", action="callback", callback=_callback_list_float,\
                          help="Normalization factors for input-DNA. If option is not chosen, estimate factors. [default: %default]")
        parser.add_option("--par", dest="par", default=1, type="int",\
                          help="Percentile for p-value postprocessing filter. [default: %default]")
        parser.add_option("--scaling-factors", default=None, dest="scaling_factors_ip", type="str", action='callback', callback=_callback_list_float,\
                          help="Scaling factor for each IP-channel input BAM file [default: %default]")
        parser.add_option("--housekeeping-genes", default=None, dest="housekeeping_genes", type="str",\
                           help="Define housekeeping genes <BED> that are used for normalization [default: %default]")
        parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", \
                          help="Output among others initial state distribution, putative differential peaks, genomic signal and histograms (original and smoothed). [default: %default]")
        parser.add_option("--version", dest="version", default=False, action="store_true",\
                           help="Show script's version.")
        parser.add_option("--output-dir", dest="outputdir", default=None, type="string", \
                          help="All files are stored in output directory which is created if necessary.")
        parser.add_option("--report", dest="report", default=False, action="store_true", \
                          help="report.")
        parser.add_option("-f", "--foldchange", dest="foldchange", default=1.3, type="float",\
                          help="Foldchange for trainingsset [default: %default]")
        parser.add_option("-t", "--threshold", dest="threshold", default=20, type="float",\
                          help="Foldchange for trainingsset [default: %default]")
        parser.add_option("--size", dest="size_ts", default=10000, type="int",\
                          help="10000 of 2 free parameters for HMM, else 1000 [default: %default]")
        
        group = OptionGroup(parser, "Advanced options")
        group.add_option("--regions", dest="regions", default=None, type="string",\
                           help="Define regions (BED) where to call DPs.")
        group.add_option("-b", "--binsize", dest="binsize", default=100, type="int",\
                          help="Size of underlying bins for creating the signal.  [default: %default]")
        group.add_option("-s", "--step", dest="stepsize", default=50, type="int",\
                          help="Stepsize with which the window consecutively slides across the genome to create the signal.")
        group.add_option("--debug", default=False, dest="debug", action="store_true", \
                          help="Output debug information. Warning: space consuming! [default: %default]")
        group.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true", \
                          help="turn off GC content calculation")
        parser.add_option_group(group)
        parser.add_option("--norm-regions", default=None, dest="norm_regions", type="str", help="Define regions <BED> that are used for normalization [default: %default]")
        group.add_option("--three-parameter", default=False, dest="hmm_free_para", action="store_true", \
                          help="HMM with 3 free parameters[default: %default]")
        ##deprecated options
        #parser.add_option("--distr", dest="distr", default="negbin", type="str",\
        #                  help="HMM's emission distribution (negbin, binom). [default: %default]")
        #parser.add_option("--save-wig", dest="save_wig", default=False, action="store_true", help="save bw as well as wig files. Warning: space consuming! [default: %default]")
        #parser.add_option("--norm-regions", default=None, dest="norm_regions", type="str", help="Define regions <BED> that are used for normalization [default: %default]")
        #parser.add_option("--ext-inputs", default=None, dest="exts_inputs", type="str", action='callback', callback=_callback_list,\
        #                  help="Read's extension size for input files. If option is not chosen, estimate extension sizes. [default: %default]")
        
	(options, args) = parser.parse_args()

    options.distr = "negbin"
    options.save_wig = False
    #options.norm_regions = None
    options.exts_inputs = None
        
    if options.version:
        version = "version \"0.1alpha\""
        print("")
        print(version)
        sys.exit()
    
    if len(args) != 1:
        parser.error("Please give config file")
        
    config_path = args[0]

    if not os.path.isfile(config_path):
        parser.error("Config file %s does not exist!" %config_path)
        
    bamfiles, genome, chrom_sizes, inputs, dims = input_parser(config_path)
    
    if options.exts and len(options.exts) != len(bamfiles):
        parser.error("Number of Extension Sizes must equal number of bamfiles")
    
    if options.exts_inputs and len(options.exts_inputs) != len(inputs):
        parser.error("Number of Input Extension Sizes must equal number of input bamfiles")
        
    if options.scaling_factors_ip and len(options.scaling_factors_ip) != len(bamfiles):
        parser.error("Number of scaling factors for IP must equal number of bamfiles")
        
    for bamfile in bamfiles:
        if not os.path.isfile(bamfile):
            parser.error("BAM file %s does not exist!" %bamfile)
    
    if not inputs and options.factors_inputs:
        print("As no input-DNA, do not use input-DNA factors", file=sys.stderr)
        options.factors_inputs = None
    
    if options.factors_inputs and len(options.factors_inputs) != len(bamfiles):
        parser.error("factors for input-DNA must equal number of BAM files!")
    
    if inputs:
        for bamfile in inputs:
            if not os.path.isfile(bamfile):
                parser.error("BAM file %s does not exist!" %bamfile)
    
    if options.regions:
        if not os.path.isfile(options.regions):
            parser.error("Region file %s does not exist!" %options.regions)
    
    if not os.path.isfile(genome):
        parser.error("Genome file %s does not exist!" %genome)
    
    if options.name is None:
        d = str(datetime.now()).replace("-", "_").replace(":", "_").replace(" ", "_"). replace(".", "_").split("_")
        options.name = "THOR-exp" + "-" + "_".join(d[:len(d)-1])
    
    if not which("wigToBigWig"):
        print("Warning: wigToBigWig programm not found! Signal will not be stored!", file=sys.stderr)

    if options.outputdir:
        options.outputdir = os.path.expanduser(options.outputdir) #replace ~ with home path
        if os.path.isdir(options.outputdir) and sum(map(lambda x: x.startswith(options.name), os.listdir(options.outputdir))) > 0:
            parser.error("Output directory exists and contains files with names starting with your chosen experiment name! Do nothing to prevend file overwriting!")
        if not os.path.exists(options.outputdir):
            os.mkdir(options.outputdir)
    else:
        options.outputdir = os.getcwd() 
    
    options.name = os.path.join(options.outputdir, options.name)
    
    if options.report:
        os.mkdir(os.path.join(options.outputdir, 'report/'))
        os.mkdir(os.path.join(options.outputdir, 'report/pics/'))
        os.mkdir(os.path.join(options.outputdir, 'report/pics/data/'))
    
    global FOLDER_REPORT
    global FOLDER_REPORT_PICS
    global FOLDER_REPORT_DATA
    global OUTPUTDIR
    global NAME
    
    FOLDER_REPORT = os.path.join(options.outputdir, 'report/')
    FOLDER_REPORT_PICS = os.path.join(options.outputdir, 'report/pics/')
    FOLDER_REPORT_DATA = os.path.join(options.outputdir, 'report/pics/data/')
    OUTPUTDIR = options.outputdir
    NAME = options.name
    
    if options.exts is None:
        options.exts = []
    
    if options.exts_inputs is None:
        options.exts_inputs = []
    
    return options, bamfiles, genome, chrom_sizes, dims, inputs

