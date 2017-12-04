#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Copyright (C) 2014-2016 Manuel Allhoff (allhoff@aices.rwth-aachen.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

%prog <BED,BAM> <BED,BAM> <NAME> <STEP_WIDTH> <Zero-counts 1/0> <hg19/mm9>

Computes p-q list, k and a from Diaz et al., 2012.

Features must have same length in both input files.

@author: Manuel Allhoff
"""

## This file could be used to write all normalization methods
## 1. normalization by inputs
## 2. normalization by signal files
### In MultiCoverageSet, only accept the factors and process w.r.t the factors

from __future__ import print_function
from optparse import OptionParser
from math import fabs
import sys, operator
import numpy as np
from random import sample
from scipy import sparse
import pysam


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""

    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))


def _accumulate(iterable):
    """Return running totals"""
    it = iter(iterable)
    total = next(it)
    yield [total[0], total[1]]
    for e1, e2 in it:
        total[0] = total[0] + e1
        total[1] = total[1] + e2
        yield [total[0], total[1]]


def write_pq_list(pq_list, max_index, max_value, factor1, factor2, filename):
    """Write p,q-list to file"""
    if pq_list:
        with open(filename, 'w') as f:
            print('#max index', 'max value', 'factor1', 'factor2', sep='\t', file=f)
            print('#' + str(max_index), str(max_value), str(factor1), str(factor2), sep='\t', file=f)
            for p, q in pq_list:
                print(p, q, file=f)


def get_bin_from_covs(cov, step_times=10):
    """
    We need to check the requirmetns for Diaz method; It's for non-overlapping bins
    Arguments:
        cov: [[chr1] , [chr2], [chr3]...]  [chr1]: [ 2,3,7,10,....]; If we store it as matrix, then we could get dimension information
    # here how could we get the bins for coverage if we use sm_coverage representation??
    Still we could use the todense function to tansform it, or to use sparse matrix generate new bins;
    Can we do it ???
    """
    # we use list to store them then we transfer them into array...
    cov_counts = []
    # so first we need to dig into cov and then get the that sizes
    # here if we use sparse matrix, could we also do this ??
    tmp_cov = np.ravel(cov.sm_overall_cov.todense()) # into one-dimension view
    step_size = (cov.binsize/cov.stepsize) * step_times
    # one problem is range of it, at end, we also want to get the last bin sizes then,
    for bin_idx in range(0, len(tmp_cov), step_size):
        bin_sum = sum(tmp_cov[bin_idx:(bin_idx+step_size):(cov.binsize/cov.stepsize)])
        cov_counts.append(bin_sum)
    return cov_counts


def _get_lists_from_cov(count_list, zero_counts, two_sample=False):
    if two_sample:
        count_list.sort(key=lambda x: x[0] + x[1])
        if not zero_counts:
            count_list = filter(lambda x: x[0] != 0.0 or x[1] != 0.0, count_list)
        pre_pq_list = list(_accumulate(count_list))
        pq_list = pre_pq_list
        max_index = int(len(pq_list) * 0.5)
        max_value = fabs(pq_list[max_index][1] - pq_list[max_index][0])
    else:
        count_list.sort(key=lambda a: a[0])
        if not zero_counts:
            count_list = filter(lambda x: x[0] != 0.0, count_list)
        pre_pq_list = list(_accumulate(count_list))
        # get k, a from Diaz et al., 2012, compute p, q from Diaz et al., 2012
        pq_list = map(lambda x: [x[0] / float(pre_pq_list[-1][0]), x[1] / float(pre_pq_list[-1][1])], pre_pq_list)
        max_index, max_value = max(enumerate(map(lambda x: fabs(x[0] - x[1]), pq_list)), key=operator.itemgetter(1))

    return pq_list, max_index, max_value, \
           float(pq_list[max_index][0]) / pq_list[max_index][1], \
           float(pq_list[max_index][1]) / pq_list[max_index][0]


def get_normalization_factor_by_cov(cov, inputs_cov,zero_counts,filename,debug,step_times=10,two_samples=False):
    """ cov and inputs_cov are coverage for signal and input, step_times is the times of original step width of coverage;
    zero_counts are 0, but the purpose???
    New method: to generate new bin statics from cov and input coves, cause they are corresponding before;
     pq_list, max_index,max_value and factor1 for each inputs file
        covs, inputs_cov should include binsize and step_size such information and then combine with step_times to get more information
    normalization to inputs and get factors
        """

    # first get new bin from cov
    cov_counts = get_bin_from_covs(cov, step_times=step_times)
    inputs_cov_counts = get_bin_from_covs(inputs_cov, step_times=step_times)
    # combine them together by zip
    count_list = zip(cov_counts, inputs_cov_counts)
    count_list = map(list, count_list)
    # after using zip, it creats a list but we can't change it, it's sad; So we need to change it
    # Maybe we don't need to use zip...
    # then count pq_list after it we get other features
    pq_list, max_index, max_value, factor1, factor2 = _get_lists_from_cov(count_list, zero_counts, two_samples)
    if debug:
        write_pq_list(pq_list, max_index, max_value, factor1, factor2, filename + '-pqlist')

    if not pq_list:
        return None, None

    if two_samples:
        l = 0.5
        s1 = sum([i for (i, j) in pq_list[:int(len(pq_list) * l)]])
        s2 = sum([j for (i, j) in pq_list[:int(len(pq_list) * l)]])

        # print("norm sums half: ", s1, s2, file=sys.stderr)

        if s1 > s2:
            return 2, factor1
        else:
            return 1, factor1
    else:
        return -1, factor1


## get signal normalization factors
def trim4TMM( m_values, a_values, m_threshold=80, a_threshold=95):
    """q=20 or q=5"""
    assert len(m_values) == len(a_values)

    perc_m_l = np.percentile(m_values, 100 - m_threshold)
    perc_m_h = np.percentile(m_values, m_threshold)
    perc_a_l = np.percentile(a_values, 100 - a_threshold)
    perc_a_h = np.percentile(a_values, a_threshold)

    try:
        res = filter(lambda x: not (x[0] > perc_m_h or x[0] < perc_m_l), \
                     filter(lambda x: not (x[1] > perc_a_h or x[1] < perc_a_l),
                            zip(list(m_values.squeeze()), list(a_values.squeeze()))))
    except:
        print('something wrong %s %s' % (len(m_values), len(a_values)), file=sys.stderr)
        return np.asarray(m_values), np.asarray(a_values)

    if res:
        return np.asarray(map(lambda x: x[0], res)), np.asarray(map(lambda x: x[1], res))
    else:
        print('TMM normalization: nothing trimmed...', file=sys.stderr)
        return np.asarray(m_values), np.asarray(a_values)


## overall_coverage should be the all coverage for signal files, but signal file are now in format sparse matrix ...
#  cuase we need to use sum and duplicate then we need to use sparse matrix, better ?? Then we could change it here
def get_sm_norm_TMM_factor(overall_coverage, m_threshold, a_threshold):
    """Normalize with TMM approach, based on PePr
      ref, we use the mean of two samples and compare to it
      data_rep present the data..
      we accept overall_coverage as one parameter, and deal with it;
      coverall_coverage [ signal_0[ file_0, file_1, file_2], signal_1[ file_0, file_1, file_2] ]
      But each file is sparse_matrix, and now we need to get non-zeros columns

    """
    factors_signal = []
    dim = overall_coverage['dim']
    # at least one file has one read in it...
    # or we could get mask_ref with all data are over 0
    # non_zeros columns for all coverage from each file
    # but we should change the value of overall_coverage; we just need to get factors
    valid_indices = None
    for i in range(dim[0]):
        for j in range(dim[1]):
            if valid_indices is None:
                valid_indices = set(overall_coverage['data'][i][j].indices)
            else:
                valid_indices &= set(overall_coverage['data'][i][j].indices)
    #mask_ref = reduce(lambda x, y: set(x) & set(y),
    #                         [overall_coverage[i][j].indices for j in range(dim[1]) for i in range(dim[0])])
    mask_ref = list(valid_indices)
    sm_ref = reduce(lambda x,y: x+y, [overall_coverage['data'][i][j][:, mask_ref] for i in range(dim[0]) for j in range(dim[1])])
    ref = sm_ref.data/float(dim[0]*dim[1])

    for i in range(dim[0]):
        factors_signal.append([])
        for j in range(dim[1]):  # normalize all replicates
            # get the data for each sample under each condition
            data_rep = overall_coverage['data'][i][j][:, mask_ref].data
            # here we sample data but in Manuel method, he uses the biggest ones...
            # we sort data by order which data_rep + ref descends and then choose the highest ones
            tmp_idx = np.argsort(data_rep + ref)[-min(len(data_rep), 50000):]
            # tmp_idx = sample(range(len(data_rep)), min(len(data_rep), 100000))
            tmp_ref = ref[tmp_idx]  # use index to make ref and data correspond
            data_rep = data_rep[tmp_idx]
            # calculate m_values and a_values
            m_values = np.log(tmp_ref / data_rep)
            a_values = 0.5 * np.log(data_rep * tmp_ref)
            try:  # assume they have a relations and then plot them to get scale factor.
                m_values, a_values = trim4TMM(m_values, a_values, m_threshold, a_threshold)
                f = 2 ** (np.sum(m_values * a_values) / np.sum(a_values))
                factors_signal[i].append(f)
            except:
                print('TMM normalization not successfully performed, do not normalize data', file=sys.stderr)
                factors_signal[i].append(1)

    return factors_signal


def get_gc_factor(cov, delta, genome_fpath, sample_size=10000):
    """ compute gc_content_factor from one input coverage(control coverage) just to get hv and T
    Return:
        for each bins ?? No, we can compute and then change it , no need to store it !!
            # compute_gc_proportion we have one bin, and we want to get gc_content_factor:
            # coverage region, chrom, genome it uses but we could get before we do this step;;
            # it means if we make sure it is from chrom-sizes files and then references to it
        h(v) for range from delta  [0, delta, 2 delta, ... , 1- delta]
        T  T for sum of h(v)

    Arguments:
        coverage: input coverage in sm_format, only coverage without chroms information??? We should have it !!
        delta: to control gc-proportion
        genome: we need to get the genome from chroms-sizes
    """
    # initialize gc-proportion array including 1.0;; and then we need to
    hv, lens = {}, {}
    for prop in np.arange(0, 1.0, delta, dtype=float):
        hv[prop] = []
        lens[prop] = 0

    # reading genome from file
    genome_fasta = pysam.Fastafile(genome_fpath)
    # fetch genome w.r.t. chromsome
    chroms= cov.genomicRegions.get_chrom()
    # chroms=['chr1','chr2','chr3']

    # coverage separated by chroms regions; so we need to find out first what it belongs to
    for i in range(len(chroms)):
        # here we need to count all, too slow;; we need sample
        for bin_idx in sample(cov.sm_coverage[i].indices, k=min(len(cov.sm_coverage[i].indices), sample_size)): # for only bins with non zeros reads
        # for bin_idx in range(cov.sm_coverage[i].shape[-1]): # for all bins
            s = bin_idx * cov.stepsize
            e = s + cov.binsize
            genome_seq = genome_fasta.fetch(reference=chroms[i], start=s, end=e)
            prop = get_gc_content_proportion(genome_seq)
            # how to make prop proper to access hv?? We could change the length and other staff;
            if prop is not None:
                # if prop is zeros, it's fine
                prop_key = int(prop / delta) * delta
                lens[prop_key] += 1
                if cov.sm_coverage[i][:,bin_idx].getnnz() > 0 :
                    hv[prop_key].append(cov.sm_coverage[i][:,bin_idx].data) # put one bin num into hv w.r.t props

    avg_T = 0.0
    for prop, nums in hv.items():
        # problem if hv.item is empty; and also, if
        if nums:
            hv[prop] = float(sum(nums))/lens[prop]
            avg_T += hv[prop]
    avg_T *= delta

    return hv, avg_T


def get_gc_content_proportion(seq):
    """according to one genome we get the gc-content for one bin """
    seq = seq.upper()
    count = seq.count("C") + seq.count("G")
    if len(seq) > 0:
        return float(count) / len(seq)
    else:
        return None


if __name__ == "__main__":
    ## test correctness of this function compared to old method
    # old method

    # new method

    class Cov(object):
        pass
    cov = Cov()
    cov.binsize = 100
    cov.stepsize = 50
    s = [0, 0, 1, 3, 5, 7, 0, 3, 2, 0, 0, 0]
    # cov.overall_cov = sparse.csr_matrix(s)

    inputs_cov = Cov()
    inputs_cov.binsize = 100
    inputs_cov.stepsize = 50
    t = [0, 0, 1, 3, 0, 1, 0, 3, 2, 0, 6, 0]
    # inputs_cov.overall_cov = sparse.csr_matrix(t)

    ## test count gc-content
    cov.coverage = [s[:],t[:],s[:5]]
    inputs_cov.coverage = [t[:], s[:], t[:5]]

    cov.sm_coverage = []
    for i in range(len(cov.coverage)):
        cov.sm_coverage.append(sparse.csr_matrix(cov.coverage[i]))

    inputs_cov.sm_coverage = []
    for i in range(len(inputs_cov.coverage)):
        inputs_cov.sm_coverage.append(sparse.csr_matrix(inputs_cov.coverage[i]))

    genome_fpath= "/home/kefang/programs/THOR_example_data/refactor/genome_hg19.fa"

    hv, avg_T = get_gc_factor(inputs_cov, 0.01, genome_fpath)
    print(avg_T)

    #cov = {'binsize':100, 'stepsize':50,'overall_cov':[0,0,1,3,5,7,0,3,2,0,0,0]}
    #inputs_cov = {'binsize':100, 'stepsize':50,'overall_cov':[0,0,1,3,0,1,0,3,2,0,6,0]}
    # get_normalization_factor_by_cov(cov, inputs_cov, 0, 'test', True, step_times=2, two_samples=False)
    """
    # test normalization factor by signal
    orig_cov = np.asarray([[[0,0,4,3,5,7,0,3,2,0,0,0], [0,0,1,2,5,5,0,3,0,0,0,0]],[[0,0,1,3,0,1,0,3,2,0,6,0], [2,0,1,3,0,1,0,1,2,0,0,3]]])
    #z = get_norm_TMM_factor(orig_cov, m_threshold=80, a_threshold=95)
    #print(z)
    print('new method')
    overall_cov = []
    for i in range(orig_cov.shape[0]):
        overall_cov.append([])
        for j in range(orig_cov.shape[1]):
            tmp = sparse.csr_matrix(orig_cov[i][j])
            overall_cov[i].append(tmp)

    #s = get_sm_norm_TMM_factor(overall_cov, m_threshold=80,a_threshold=95)
    #print(s)
    cov = {'dim':[2,2], 'data':overall_cov}
    z1 = dpc_help.fit_sm_mean_var_distr(cov, 'fun-var-test', True, False, 'test', False, False, sample_size=5000)
    #z2 = dpc_help.fit_mean_var_distr(orig_cov, 'fun-var-test', True, False, 'test', False, False, sample_size=5000)
    print(z1)
    # print(z2)
    """