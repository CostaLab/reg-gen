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
import pysam
import numpy as np
from math import fabs, log, ceil
from operator import add
from os.path import splitext, basename, join, isfile, isdir, exists


# Internal
from ..THOR.postprocessing import merge_delete, filter_deadzones
from .MultiCoverageSet import MultiCoverageSet
from ..GenomicRegionSet import GenomicRegionSet
from ..THOR.get_extension_size import get_extension_size
from ..THOR.get_fast_gen_pvalue import get_log_pvalue_new
from ..Util import which, npath

# External
from numpy import linspace
from scipy.optimize import curve_fit
import matplotlib as mpl
#see http://stackoverflow.com/questions/4931376/generating-matplotlib-graphs-without-a-running-x-server
mpl.use('Agg')
import matplotlib.pyplot as plt

FOLDER_REPORT = None

np.random.seed(42)


def merge_output(bamfiles, dims, options, no_bw_files, chrom_sizes):
    for i in range(len(bamfiles)):
        rep = i if i < dims[0] else i - dims[0]
        sig = 1 if i < dims[0] else 2

        temp_bed = npath(options.name + '-s%s-rep%s_temp.bed' % (sig, rep))

        files = [options.name + '-' + str(j) + '-s%s-rep%s.bw' %(sig, rep) for j in no_bw_files]
        if len(no_bw_files) > len(bamfiles):
            files = filter(lambda x: isfile(x), files)
            t = ['bigWigMerge'] + files + [temp_bed]
            c = " ".join(t)
            os.system(c)

            os.system("LC_COLLATE=C sort -k1,1 -k2,2n " + temp_bed + ' > ' + temp_bed +'.sort')

            t = ['bedGraphToBigWig', temp_bed + '.sort', chrom_sizes, options.name + '-s%s-rep%s.bw' % (sig, rep)]
            c = " ".join(t)
            os.system(c)

            for f in files:
                os.remove(f)
            os.remove(temp_bed)
            os.remove(temp_bed + ".sort")
        else:
            ftarget = [options.name + '-s%s-rep%s.bw' %(sig, rep) for j in no_bw_files]
            for i in range(len(ftarget)):
                c = ['mv', files[i], ftarget[i]]
                c = " ".join(c)
                os.system(c)


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
        x = linspace(0, max(plot_data[i][0]), int(ceil(max(plot_data[i][0]))))
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
            plt.plot(x, y, 'r', label = 'fitted polynomial') #plot polynom
            plt.scatter(plot_data[i][0], plot_data[i][1], label = 'empirical datapoints') #plot datapoints
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
        data_rep[i] = data_rep[i][data_rep[i][:,0] < np.percentile(data_rep[i][:,0], 99.75)]
        data_rep[i] = data_rep[i][data_rep[i][:,1] < np.percentile(data_rep[i][:,1], 99.75)]
    
    return data_rep


def _fit_mean_var_distr(overall_coverage, name, debug, verbose, outputdir, report, poisson, sample_size=5000):
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
                
                if len(m) > 0 and len(v) > 0: 
                    try:
                        p, _ = curve_fit(_func_quad_2p, m, v) #fit quad. function to empirical data
                    except:
                        print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                        break #restart for loop
                else:
                    p = np.array([0, 1])
                
                res.append(p)
                plot_data.append((m, v, p))
                if i == 1:
                    done = True
            except RuntimeError:
                print("Optimal parameters for mu-var-function not found, get new datapoints", file=sys.stderr)
                break #restart for loop
    
    if report:
        _plot_func(plot_data, outputdir)
    
    if poisson:
        print("Use Poisson distribution as emission", file=sys.stderr)
        p[0] = 0
        p[1] = 0
        res = [np.array([0, 0]), np.array([0, 0])]
    
    return lambda x: _func_quad_2p(x, p[0], p[1]), res

    
def dump_posteriors_and_viterbi(name, posteriors, DCS, states):
    print("Computing info...", file=sys.stderr)
    f = open(name + '-posts.bed', 'w')
    g = open(name + '-states-viterbi.bed', 'w')
    
    for i in range(len(DCS.indices_of_interest)):
        cov1, cov2 = _get_covs(DCS, i)
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
    

def _get_covs(DCS, i, as_list=False):
    """For a multivariant Coverageset, return mean coverage cov1 and cov2 at position i"""
    if not as_list:
        cov1 = int(np.mean(DCS.overall_coverage[0][:, DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:, DCS.indices_of_interest[i]]))
    else:
        cov1 = DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]
        cov1 = map(lambda x: x[0], np.asarray((cov1)))
        cov2 = DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]
        cov2 = map(lambda x: x[0], np.asarray((cov2)))
    
    return cov1, cov2


def get_peaks(name, DCS, states, exts, merge, distr, pcutoff, debug, no_correction, deadzones, merge_bin, p=70):
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
        f = open(FOLDER_REPORT_DATA + 'fragment_size_estimate_' + names[k] + '.data', 'w')
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
    plt.savefig(FOLDER_REPORT_PICS + 'fragment_size_estimate.png')
    plt.close()


def _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, report):
    """Compute Extension sizes for bamfiles and input files"""
    start = 0
    end = 600
    ext_stepsize = 5
    
    ext_data_list = []
    #compute extension size
    if not exts:
        print("Computing read extension sizes for ChIP-seq profiles", file=sys.stderr)
        for bamfile in bamfiles:
            e, ext_data = get_extension_size(bamfile, start=start, end=end, stepsize=ext_stepsize)
            exts.append(e)
            ext_data_list.append(ext_data)
    
    if report and ext_data_list:
        _output_ext_data(ext_data_list, bamfiles)
    
    if inputs and not exts_inputs:
        exts_inputs = [5] * len(inputs)

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
               tracker, debug, norm_regions, scaling_factors_ip, save_wig, housekeeping_genes, \
               test, report, chrom_sizes_dict, counter, end, gc_content_cov=None, avg_gc_content=None, \
               gc_hist=None, output_bw=True, save_input=False, m_threshold=80, a_threshold=95, rmdup=False):
    """Initialize the MultiCoverageSet"""
    regionset = regions
    regionset.sequences.sort()
    
    if norm_regions:
        norm_regionset = GenomicRegionSet('norm_regions')
        norm_regionset.read(norm_regions)
    else:
        norm_regionset = None
        
    exts, exts_inputs = _compute_extension_sizes(bamfiles, exts, inputs, exts_inputs, report)
    
    multi_cov_set = MultiCoverageSet(name=name, regions=regionset, dims=dims, genome_path=genome_path,
                                     binsize=binsize, stepsize=stepsize, rmdup=rmdup, path_bamfiles=bamfiles,
                                     path_inputs=inputs, exts=exts, exts_inputs=exts_inputs,
                                     factors_inputs=factors_inputs, chrom_sizes=chrom_sizes, verbose=verbose,
                                     no_gc_content=no_gc_content, chrom_sizes_dict=chrom_sizes_dict, debug=debug,
                                     norm_regionset=norm_regionset, scaling_factors_ip=scaling_factors_ip,
                                     save_wig=save_wig, strand_cov=True, housekeeping_genes=housekeeping_genes,
                                     tracker=tracker, gc_content_cov=gc_content_cov, avg_gc_content=avg_gc_content,
                                     gc_hist=gc_hist, end=end, counter=counter, output_bw=output_bw,
                                     folder_report=FOLDER_REPORT, report=report, save_input=save_input,
                                     m_threshold=m_threshold, a_threshold=a_threshold)
    return multi_cov_set


