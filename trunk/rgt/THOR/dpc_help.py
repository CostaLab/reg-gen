from __future__ import print_function
from optparse import OptionParser
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
from matplotlib.pyplot import *
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

def _plot_func(m, v, p, name, xlimpara=None, ylimpara=None):
    """Plot estimated and empirical function"""
    x = linspace(0, max(m), max(m)+1)
    y = _func_quad_2p(x, p[0], p[1])
    
    if xlimpara is not None:
        xlim([0, xlimpara])
    if ylimpara is not None:
        ylim([0, ylimpara])
    
    plot(x, y)
    scatter(m, v)
    
    savefig(name + ".png")
    close()

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
    
def _fit_mean_var_distr(overall_coverage, name, debug, sample_size=10000):
    """Estimate empirical distribution (quadr.) based on empirical distribution"""
    done = False
    while not done:
        data_rep = _get_data_rep(overall_coverage, name, debug, sample_size)
        res = []
        for i in range(2):
            try:
                m = np.asarray(map(lambda x: x[0], data_rep[i])) #means list
                v = np.asarray(map(lambda x: x[1], data_rep[i])) #vars list
                p, _ = curve_fit(_func_quad_2p, m, v) #fit quad. function to empirical data
                res.append(p)
                _plot_func(m, v, p, str(name) + "-est-func" + str(i))
                _plot_func(m, v, p, str(name) + "-est-func-constraint" + str(i), xlimpara=200, ylimpara=3000)
                if i == 1:
                    done = True
            except RuntimeError:
                print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                break #restart for loop
                
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

def _merge_consecutive_bins(tmp_peaks, distr):
    """Merge consecutive peaks and compute p-value. Return list 
    <(chr, s, e, c1, c2, strand)> and <(pvalue)>"""
    peaks = []
    pvalues = []
    i, j, = 0, 0
    while i < len(tmp_peaks):
        j+=1
        c, s, e, c1, c2, strand = tmp_peaks[i]
        v1 = c1
        v2 = c2
        
        #merge bins
        while i+1 < len(tmp_peaks) and e == tmp_peaks[i+1][1] and strand == tmp_peaks[i+1][5]:
            e = tmp_peaks[i+1][2]
            v1 = map(add, v1, tmp_peaks[i+1][3])
            v2 = map(add, v2, tmp_peaks[i+1][4])
            i += 1
        s1 = v1
        s2 = v2
        side = 'l' if strand == '+' else 'r'
        pvalues.append((s1, s2, side, distr))
        peaks.append((c, s, e, s1, s2, strand))
        i += 1
    
    pvalues = map(_compute_pvalue, pvalues)
    assert len(pvalues) == len(peaks)
    
    return pvalues, peaks

def get_back(DCS, states):
    counts = []
#     print("H", file=sys.stderr)
    for i in range(len(DCS.indices_of_interest)):
        if states[i] == 0:
            cov1, cov2 = _get_covs(DCS, i)
            counts.append(cov1)
            counts.append(cov2)
#     print(len(counts), file=sys.stderr)
#     print(np.var(counts), file=sys.stderr)
#     print(np.mean(counts), file=sys.stderr)
    return np.var(counts), np.mean(counts)
        
    
def get_peaks(name, DCS, states, exts, merge, distr, pcutoff, p=70):
    """Merge Peaks, compute p-value and give out *.bed and *.narrowPeak"""
    exts = np.mean(exts)
    tmp_peaks = []
    tmp_data = []
    
    for i in range(len(DCS.indices_of_interest)):
        if states[i] not in [1,2]:
            continue #ignore background states
        
        strand = '+' if states[i] == 1 else '-'
        cov1, cov2 = _get_covs(DCS, i, as_list=True)
        c1, c2 = sum(cov1), sum(cov2)
        chrom, start, end = DCS._index2coordinates(DCS.indices_of_interest[i])
        tmp_peaks.append((chrom, start, end, cov1, cov2, strand))
        side = 'l' if strand == '+' else 'r'
        tmp_data.append((c1, c2, side, distr))
    
    tmp_pvalues = map(_compute_pvalue, tmp_data)
    per = np.percentile(tmp_pvalues, p)
    print('percentile', per, file=sys.stderr)
    
    tmp = []
    res = tmp_pvalues > per
    for j in range(len(res)):
        if res[j]:
            tmp.append(tmp_peaks[j])
    tmp_peaks = tmp

    #merge consecutive peaks and compute p-value
    pvalues, peaks = _merge_consecutive_bins(tmp_peaks, distr)
    #postprocessing, returns GenomicRegionSet with merged regions
    regions = merge_delete(exts, merge, peaks, pvalues) 
    #regions = merge_delete([0], False, peaks, pvalues) 
    output = []
    pvalues = []
    main_sep = ':' #sep <counts> main_sep <counts> main_sep <pvalue>
    int_sep = ';' #sep counts in <counts>
    
    for i, el in enumerate(regions):
        tmp = el.data.split(',')
        counts = ",".join(tmp).replace('], [', ';').replace('], ', int_sep).replace('([', '').replace(')', '').replace(', ', main_sep)
        pvalue = float(tmp[len(tmp)-1].replace(")", "").strip())
        pvalues.append(pvalue)
        output.append((el.chrom, el.initial, el.final, el.orientation, counts))
    
    pcutoff = -log10(pcutoff)
    pv_pass = np.where(np.asarray(pvalues) >= pcutoff, True, False)

    output = np.array(output)
    output = output[pv_pass]
    pvalues = list(np.array(pvalues)[pv_pass])
    output = map(lambda x: tuple(x), list(output))
    
    assert(len(output) == len(pvalues))
    
    _output_BED(name, output, pvalues)
    _output_narrowPeak(name, output, pvalues)

def _output_BED(name, output, pvalues):
    f = open(name + '-diffpeaks.bed', 'w')
     
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, strand, counts = output[i]
        print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, colors[strand], 0, counts, sep='\t', file=f)
    
    f.close()

def _output_narrowPeak(name, output, pvalues):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name + '-diffpeaks.narrowPeak', 'w')
    for i in range(len(pvalues)):
        c, s, e, strand, _ = output[i]
        print(c, s, e, 'Peak' + str(i), 0, strand, 0, pvalues[i], 0, -1, sep='\t', file=f)
    f.close()

def _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, verbose):
    """Compute Extension sizes for bamfiles and input files"""
    start = 0
    end = 600
    ext_stepsize = 5

    #compute extension size
    if not exts:
        print("Computing read extension sizes for ChIP-DNA...", file=sys.stderr)
        for bamfile in bamfiles:
            e, _ = get_extension_size(bamfile, start=start, end=end, stepsize=ext_stepsize)
            exts.append(e)
        print(exts, file=sys.stderr)

    if inputs and not exts_inputs:
        print("Computing read extension sizes for input-DNA...", file=sys.stderr)
        
        for inp in inputs:
            e, _ = get_extension_size(inp, start=start, end=end, stepsize=ext_stepsize)
            exts_inputs.append(e)
        print(exts_inputs, file=sys.stderr)

    return exts, exts_inputs

def initialize(name, dims, genome_path, regions, stepsize, binsize, bamfiles, exts, \
               inputs, exts_inputs, factors_inputs, chrom_sizes, verbose, no_gc_content, \
               tracker, debug, norm_regions, scaling_factors_ip, save_wig):
    """Initialize the MultiCoverageSet"""

    regionset = GenomicRegionSet(name)
    chrom_sizes_dict = {}
    #if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
    if regions is not None:
        print("Call DPs on specified regions.", file=sys.stderr)
        with open(regions) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                c, s, e = line[0], int(line[1]), int(line[2])
                regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
                chrom_sizes_dict[c] = e
    else:
        print("Call DPs on whole genome.", file=sys.stderr)
        with open(chrom_sizes) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                chrom, end = line[0], int(line[1])
                regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
                chrom_sizes_dict[chrom] = end
    
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
                                  norm_regionset=norm_regionset, scaling_factors_ip=scaling_factors_ip, save_wig=save_wig)
    
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
        #config_path = '/home/manuel/workspace/eclipse/office_share/simulator/test.config'
        bamfiles, regions, genome, chrom_sizes, inputs, dims = input_parser(config_path)
        options.exts = [200, 200, 200, 200, 200]
        options.exts_inputs = None #[200, 200, 200, 200, 200]
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
        options.norm_regions = '/home/manuel/data/testdata/norm_regions.bed'
        options.scaling_factors_ip = False
    else:
        parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.1, type="float",\
                          help="P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: %default]")
        parser.add_option("-m", "--merge", default=False, dest="merge", action="store_true", \
                          help="Merge peaks which have a distance less than the estimated fragment size (recommended for histone data). [default: %default]")
        parser.add_option("-b", "--binsize", dest="binsize", default=100, type="int",\
                          help="Size of underlying bins for creating the signal.  [default: %default]")
        parser.add_option("-s", "--step", dest="stepsize", default=50, type="int",\
                          help="Stepsize with which the window consecutively slides across the genome to create the signal.")
        parser.add_option("-n", "--name", default='run', dest="name", type="string",\
                          help="Experiment's name and prefix for all files that are created.")
        parser.add_option("--exts", default=None, dest="exts", type="str", action='callback', callback=_callback_list,\
                          help="Read's extension size for BAM files. If option is not chosen, estimate extension sizes. [default: %default]")
        parser.add_option("--ext-inputs", default=None, dest="exts_inputs", type="str", action='callback', callback=_callback_list,\
                          help="Read's extension size for input files. If option is not chosen, estimate extension sizes. [default: %default]")
        parser.add_option("--factors-inputs", default=None, dest="factors_inputs", type="float",\
                          help="Normalization factors for inputs. If option is not chosen, estimate factors. [default: %default]")
        parser.add_option("--par", dest="par", default=1, type="int",\
                          help="Percentile for p-value filter. [default: %default]")
        parser.add_option("--save-wig", dest="save_wig", default=False, action="store_true", help="save bw and wig. Warning: sapce consuming! [default: %default]")
        
        parser.add_option("--norm-regions", default=None, dest="norm_regions", type="str", help="Define regions <BED> that are used for normalization")
        parser.add_option("--scaling-factors", default=None, dest="scaling_factors_ip", type="str", action='callback', callback=_callback_list_float,\
                          help="Scaling factor for each BAM file [default: %default]")
        
        parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", \
                          help="Output among others initial state distribution, putative differential peaks, genomic signal and histograms (original and smoothed). [default: %default]")
        parser.add_option("--version", dest="version", default=False, action="store_true", help="Show script's version.")
        parser.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true", \
                          help="turn of GC content calculation")
        parser.add_option("--debug", default=False, dest="debug", action="store_true", \
                          help="Output debug information. Warning: space consuming! [default: %default]")
        (options, args) = parser.parse_args()

        if options.version:
            version = "version \"0.0.1alpha\""
            print("")
            print(version)
            sys.exit()
        
        if len(args) != 1:
            parser.error("Please give config file")
            
        config_path = args[0]
        bamfiles, regions, genome, chrom_sizes, inputs, dims = input_parser(config_path)
        
        if options.exts and len(options.exts) != len(bamfiles):
            parser.error("Number of Extension Sizes must equal number of bamfiles")
                
        
        if options.exts_inputs and len(options.exts_inputs) != len(inputs):
            parser.error("Number of Input Extension Sizes must equal number of input bamfiles")
            
        if options.scaling_factors_ip and len(options.scaling_factors_ip) != len(bamfiles):
            parser.error("Number of scaling factors must equal number of bamfiles")
#        bamfile_2 = args[1]
#        regions = args[2]
#        genome = args[3]
#        chrom_sizes = args[4]
#        
#        if not os.path.isfile(bamfile_1) or not os.path.isfile(bamfile_2) \
#            or not os.path.isfile(regions) or not os.path.isfile(genome):
#            parser.error("At least one input parameter is not a file")
#        
#        if options.name is None:
#            prefix = os.path.splitext(os.path.basename(bamfile_1))[0]
#            suffix = os.path.splitext(os.path.basename(bamfile_2))[0]
#            options.name = "-".join(['exp', prefix, suffix])
#        
#        if (options.input_1 is None and options.ext_input_1 is not None) \
#            or (options.input_2 is None and options.ext_input_2 is not None):
#            parser.error("Read extension size without input file (use -i)")
    
    
#    if options.input_1 is None and options.input_2 is None:
#        print("GC content is not calculated as there is no input file.", file=sys.stderr)
#        
#    if options.norm_strategy in [2, 4, 5] and (options.input_1 is None or options.input_2 is None):
#        parser.error("Please define input files for this normalization strategy!")
#        
#    if options.norm_strategy is not None and (options.input_factor_1 is not None or options.input_factor_1 is not None ):
#        parser.error("Input factors are not allowed for this normalization strategy!")
#    
#    if not options.no_gc_content and (options.input_1 is None or options.input_2 is None):
#        parser.error("GC content can only be computed with both input files.")
    
    if options.exts is None:
        options.exts = []
    if options.exts_inputs is None:
        options.exts_inputs = []
    
    return options, bamfiles, regions, genome, chrom_sizes, dims, inputs
