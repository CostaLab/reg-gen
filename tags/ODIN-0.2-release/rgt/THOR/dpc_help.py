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
from math import log
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

SIGNAL_CUTOFF = 10000

def _func_quad_2p(x, a, c):
    res = []
    if type(x) is np.ndarray:
        for el in x:
            res.append(max(el, fabs(a) * el**2 + el + fabs(c)))
            
        return np.asarray(res)
    else:
        return max(x, fabs(a) * x**2 + x + fabs(c))

def _fit_mean_var_distr(overall_coverage, name, verbose, cut=1.0, sample_size=10000):
        #list of (mean, var) points for samples 0 and 1
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
        
        print('cut', cut, file=sys.stderr)
#         for i in range(2): #shorten list
#             data_rep[i].sort()
#             data_rep[i] = data_rep[i][:int(len(data_rep[i]) * cut)]
        
        if verbose:
            for i in range(2):
                np.save(str(name) + "-data" + str(i) + ".npy", data_rep[i])
        
        res = []
        for i in range(2):
            m = np.asarray(map(lambda x: x[0], data_rep[i])) #means list
            v = np.asarray(map(lambda x: x[1], data_rep[i])) #vars list
            
            #p = np.polynomial.polynomial.polyfit(m, v, 2)
            p, _ = curve_fit(_func_quad_2p, m, v)
            print('popt ', p, file=sys.stderr)
            
            res.append(p)
            print('Length', len(m), file=sys.stderr)
            if verbose:
                print(p, file=sys.stderr)
                print(max(m), max(v), file=sys.stderr)
            
                x = linspace(0, max(m), max(m)+1)
                y = _func_quad_2p(x, p[0], p[1])
                plot(x, y)
                scatter(m, v)
                savefig(str(name) + "plot_original" + str(i) + ".png")
                close()
                
                plot(x, y)
                scatter(m, v)
                ylim([0, 3000])
                xlim([0, 200])
                savefig(str(name) + "plot" + str(i) + ".png")
                close()
        
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
    print("done...", file=sys.stderr)


def _compute_pvalue((x, y, side, distr)):
    if x == 'NA':
        return sys.maxint
    else:
        return -get_log_pvalue_new(x, y, side, distr)

def _get_covs(DCS, i):
    """For a multivariant Coverageset, return coverage cov1 and cov2 at position i"""
    cov1 = int(np.mean(DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]))
    cov2 = int(np.mean(DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]))
    
    return cov1, cov2

def get_peaks(name, DCS, states, distr):
    tmp_peaks = []
    for i in range(len(DCS.indices_of_interest)):
        if states[i] not in [1,2]:
            continue #ignore background states
        
        strand = '+' if states[i] == 1 else '-'
        
        cov1, cov2 = _get_covs(DCS, i)
        
        chrom, start, end = DCS._index2coordinates(DCS.indices_of_interest[i])
        
        tmp_peaks.append((chrom, start, end, cov1, cov2, strand))

    i, j, r = 0, 0, 0
    
    peaks = []
    pvalues = []
    
    while i < len(tmp_peaks):
        j+=1
        c, s, e, c1, c2, strand = tmp_peaks[i]
        v1 = [c1]
        v2 = [c2]
        
        #merge bins
        while i+1 < len(tmp_peaks) and e == tmp_peaks[i+1][1] and strand == tmp_peaks[i+1][5]:
            e = tmp_peaks[i+1][2]
            v1.append(tmp_peaks[i+1][3])
            v2.append(tmp_peaks[i+1][4])
            i += 1
        
        s1 = sum(v1)
        s2 = sum(v2)

        if s1 + s2 > SIGNAL_CUTOFF:
            pvalues.append(('NA', 'NA', 'NA', 'NA'))
        else:
            if strand == '+':
                pvalues.append((s1, s2, 'l', distr))
            else:
                pvalues.append((s1, s2, 'r', distr))
        
        peaks.append((c, s, e, s1, s2, strand))
        i += 1
    
    print('Number of Peaks where p-value is not calculated: ', pvalues.count(('NA', 'NA', 'NA', 'NA')), file=sys.stderr)
    
    #pool = multiprocessing.Pool(processes=2)#multiprocessing.cpu_count() * 3/2)
    pvalues = map(_compute_pvalue, pvalues)
    
    assert len(pvalues) == len(peaks)
    
    f = open(name + '-diffpeaks.bed', 'w')
    
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, c1, c2, strand = peaks[i]
        color = colors[strand]

        print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, \
              color, 0, str(c1) + ',' + str(c2) + ',' + str(pvalues[i]), sep='\t', file=f)

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
        #print('inputs', inputs, file=sys.stderr)
        #print('exts_inputs', exts_inputs, file=sys.stderr)
        
        for inp in inputs:
            e, b = get_extension_size(inp, start=start, end=end, stepsize=ext_stepsize)
            #print(b, file=sys.stderr)
            #print(inp, e, file=sys.stderr)
            exts_inputs.append(e)
            #print(exts_inputs, file=sys.stderr)
        print(exts_inputs, file=sys.stderr)
    
    return exts, exts_inputs
#    if verbose:
#        if 'values_1' in locals() and values_1 is not None:
#            with open(name + '-read-ext-1', 'w') as f:
#                for v, i in values_1:
#                    print(i, v, sep='\t', file=f)


def initialize(name, dims, genome_path, regions, stepsize, binsize, bamfiles, exts, \
               inputs, exts_inputs, factors_inputs, chrom_sizes, verbose, no_gc_content, tracker):
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
                
    regionset.sequences.sort()
    exts, exts_inputs = _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, verbose)
    tracker.write(text=str(exts).strip('[]'), header="Extension size (rep1, rep2, input1, input2)")
    
    multi_cov_set = MultiCoverageSet(name=name, regions=regionset, dims=dims, genome_path=genome_path, binsize=binsize, stepsize=stepsize,rmdup=True,\
                                  path_bamfiles = bamfiles, path_inputs = inputs, exts = exts, exts_inputs = exts_inputs, factors_inputs = factors_inputs, \
                                  chrom_sizes=chrom_sizes, verbose=verbose, no_gc_content=no_gc_content, chrom_sizes_dict=chrom_sizes_dict)
    
    return multi_cov_set


def _get_chrom_list(file):
    l = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            if line[0] not in l:
                l.append(line[0])
    return l

def _check_order(deadzones, regions, parser):
    chrom_dz = _get_chrom_list(deadzones)
    chrom_regions= _get_chrom_list(regions)
    #chrom_regions may be subset of chrom_dz, but in same ordering
    pos_old = -1
    tmp_dz = []
    #the list should have the same element
    for r_chrom in chrom_dz:
        if r_chrom in chrom_regions:
            tmp_dz.append(r_chrom)
    #they should be in the same order
    if not tmp_dz == chrom_regions:
        parser.error("Please make sure the deadzone file has the same order as the region file.")
    

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def _callback_list(option, opt, value, parser):
    setattr(parser.values, option.dest, map(lambda x: int(x), value.split(',')))

def input(laptop):
    parser = HelpfulOptionParser(usage=__doc__)
    if laptop:
        print("---------- TEST MODE ----------", file=sys.stderr)
        (options, args) = parser.parse_args()
        config_path = '/home/manuel/workspace/eclipse/office_share/blueprint/playground/input_test'
        bamfiles, regions, genome, chrom_sizes, inputs, dims = input_parser(config_path)
        options.exts = [200, 200, 200, 200, 200]
        options.exts_inputs = [200, 200, 200, 200, 200]
        options.pcutoff = 1
        options.name='test'
        options.stepsize=50
        options.binsize=100
        options.factors_inputs = None
        options.verbose = True
        options.no_gc_content = True
        options.cut_obs = 1.0
    else:
#        parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.01, type="float",\
#                          help="P-value cutoff for peak detection. Call only peaks with p-value lower than cutoff. [default: %default]")
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
        parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", \
                          help="Output among others initial state distribution, putative differential peaks, genomic signal and histograms (original and smoothed). [default: %default]")
        parser.add_option("--version", dest="version", default=False, action="store_true", help="Show script's version.")
        parser.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true", \
                          help="turn of GC content calculation")
        parser.add_option("-c", "--cut", dest="cut_obs", default=1.0, type="float",\
                          help="Cut for observation.")
        
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
#    
#    if options.deadzones is not None:
#        #check the ordering of deadzones and region
#        _check_order(options.deadzones, regions, parser)
    
    if options.exts is None:
        options.exts = []
    if options.exts_inputs is None:
        options.exts_inputs = []
    
    return options, bamfiles, regions, genome, chrom_sizes, dims, inputs
