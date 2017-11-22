"""
THOR detects differential peaks in multiple ChIP-seq profiles between
two distinct biological conditions.

Copyright (C) 2014-2016 Manuel Allhof (allhoff@aices.rwth-aachen.de)

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

@author: Manuel Allhoff
"""

from __future__ import print_function
import sys
import gc
from random import sample
import numpy as np
from normalize import get_normalization_factor_by_cov, get_sm_norm_TMM_factor,get_gc_factor
from norm_genelevel import norm_gene_level
from rgt.CoverageSet import CoverageSet
import configuration


class MultiCoverageSet(): # DualCoverageSet
    """Not inherit from DualCoveragSet, instead we can use it to represent DualCoverageSet"""

    def _init_coverage(self, statics, region_giver,rmdup, binsize, stepsize, strand_cov,use_sm=False):
        """Return covs in list as CoverageSet for one statics
        self._help_init(signal_statics, inputs_statics,region_giver, norm_regionset, rmdup, binsize, stepsize, strand_cov = strand_cov)
        But before we need to do statistics about the file, get information, how much read for this regions are in this fields..
        Better we need to do is get all info of all fields, and extract data from it to combine the training fields.
        Later we use the paras to estimate peaks.. Maybe according to regions
        mask_file is used to filter data before we do other operations.
        """
        if statics:
            covs = []
            for i in range(statics['dim'][0]):
                covs.append([])
                for j in range(statics['dim'][1]):
                    covs[i].append(CoverageSet('file_' + str(i)+'_'+str(j), region_giver.valid_regionset))
                    covs[i][j].coverage_from_bam(bam_file=statics['data'][i][j]['fname'], extension_size=statics['data'][i][j]['extension_size'], rmdup=rmdup, binsize=binsize, \
                                    stepsize=stepsize, mask_file=region_giver.mask_file, get_strand_info=strand_cov, use_sm=use_sm)
            return covs
        else:
            return None

    def init_overall_coverage(self, strand_cov=True):
        """Convert coverage data (and optionally strand data) to matrix list"""
        # overall_coverage format are [[coverage_0_0 , coverage_0_1 ], [ coverage_1_0, coverag_1_1 ]]
        self.overall_coverage = {'dim':self.dim, 'data':[]} # here we could add bins size into dim
        for i in range(self.dim[0]):
            self.overall_coverage['data'].append([])
            for j in range(self.dim[1]):
                # make it use sparse matrix; but here we need to consider how to achieve it??
                # also, what we can do is normalization_by_inputs;; here what we can do is ??
                self.overall_coverage['data'][i].append(self.covs[i][j].sm_overall_cov)
        if strand_cov:
            self.overall_coverage_strand = {'dim':self.dim, 'data':[]}
            for i in range(self.dim[0]):
                self.overall_coverage_strand['data'].append([])
                for j in range(self.dim[1]):
                    self.overall_coverage_strand['data'][i].append(self.covs[i][j].cov_strand_all)

    def _is_cov_valid(self, statics, covs):
        """test if data coverage valid
            # statics is not None but covs are None or not big enough, return False
            # How to judge that data are big enough?? We could get chroms_length from statics
            # count positive signal and compare it to the length of chroms
            # statics is not None and covs are big enough, return True
            # static is None and covs are None, return True
        """
        if statics:
            if not covs:
                return False
            for i in range(self.dim[0]):
                for j in range(self.dim[1]):
                    if covs[i][j].sm_overall_cov.sum() < covs[i][j].sm_overall_cov.shape[-1] * 0.00001:
                        return False
        return True

    def __init__(self, name, region_giver, binsize, stepsize, norm_regionset, \
                 verbose, debug, rmdup, signal_statics, inputs_statics,  save_wig, strand_cov,
                 tracker, end, counter, hv=None, avg_T=None, gc_hist=None, output_bw=True,\
                 folder_report=None, report=None, save_input=False, ignored_regions=None, use_sm=False):
        """Compute CoverageSets, GC-content and normalize input-DNA and IP-channel"""
        # one improvement is to make the left one key_word parameter and we parse it, not like this, all in a list
        self.region_giver = region_giver
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.dim = signal_statics['dim']
        self.hv = hv
        self.avg_T = avg_T
        self.gc_hist = gc_hist
        self.end = end
        self.counter = counter # use of counter ???
        self.no_data = False
        self.FOLDER_REPORT = folder_report

        configuration.DEBUG = debug
        configuration.VERBOSE = verbose

        # here we need to judge if the coverage fine, or not; Actually before we have read statitics data,so it's fine
        # could give info about doing what; reading signal files, reading inputs files
        if signal_statics:
            print("Reading CHIP bamfiles", file=sys.stderr)
            self.covs = self._init_coverage(signal_statics, region_giver, rmdup, binsize, stepsize, strand_cov=strand_cov, use_sm=use_sm)

        if inputs_statics:
            print("Reading inputs bamfiles", file=sys.stderr)
            self.inputs_covs = self._init_coverage(inputs_statics, region_giver, rmdup, binsize, stepsize, strand_cov=strand_cov, use_sm=use_sm)

        if not self._is_cov_valid(signal_statics, self.covs) or not self.covs:
            self.data_valid = False
            return None
        if not self._is_cov_valid(inputs_statics, self.inputs_covs):
            self.data_valid = False
            return None
        # init overall_coverage for signal fiels
        # self._init_overall_coverage(strand_cov=strand_cov)
        self.indices_of_interest = []

    def _get_sm_covs(self, idx):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i in integer?? or in float?? """
        cov1 = int(sum([self.overall_coverage['data'][0][j][:,idx] for j in range(self.dim[1])]).toarray()/self.dim[1])
        cov2 = int(sum([self.overall_coverage['data'][1][j][:,idx] for j in range(self.dim[1])]).toarray()/self.dim[1])
        return cov1, cov2

    def output_input_bw(self, name, chrom_sizes, save_wig):
        """print inputs bw"""
        for sig in range(self.dim[0]):
            for rep in range(self.dim[1]):
                if self.inputs_covs:
                    self.inputs_covs[sig][rep].write_bigwig(name + '-' + str(self.counter) + '-input-s%s-rep%s.bw' %(sig, rep), chrom_sizes, save_wig=save_wig, end=self.end)

    def _output_bw(self, name, chrom_sizes, save_wig, save_input):
        """Output bigwig files"""
        for sig in range(self.dim[0]):
            for rep in range(self.dim[1]):
            
                self.covs[sig][rep].write_bigwig(name + '-' + str(self.counter) + '-s%s-rep%s.bw' %(sig, rep), chrom_sizes, save_wig=save_wig, end=self.end)
        
        #ra = [self.covs_avg, self.input_avg] if self.inputs else [self.covs_avg]
        #for k, d in enumerate(ra):
        #    g = self.covs if k == 0 else self.inputs
        #    for j in range(2):
        #        d[j] = deepcopy(g[0]) if j == 0 else deepcopy(g[self.dim_1])
        #        r = range(1, self.dim_1) if j == 0 else range(self.dim_1 + 1, self.dim_1 + self.dim_2)
        #        f = 1./self.dim_1 if j == 0 else 1./self.dim_2
        #        for i in r:
        #            d[j].add(g[i])
        #        d[j].scale(f)
        #        n = name + '-s%s.bw' %(j+1) if k == 0 else name + '-s%s-input.bw' %(j+1)
        #        d[j].write_bigwig(n, chrom_sizes, save_wig=save_wig, end=self.end)
        
        self.covs_avg = None
        self.input_avg = None
        if self.inputs_covs:
            for i in range(len(self.covs)):
                self.inputs_covs[i] = None #last time that we need this information, delete it
        gc.collect()

    def get_max_colsum(self):
        """Sum over all columns and add maximum"""
        # how to change it to self.cov ???
        return self.overall_coverage[0].sum(axis=0).max() + self.overall_coverage[1].sum(axis=0).max()
    
    def output_overall_coverage(self, path):
        for j in range(2):
            f = open(path + str(j), 'w')
            for i in range(self.overall_coverage[j].shape[1]):
                print(self.overall_coverage[j][:,i].T, file=f)

    def normalization_by_gc_content(self, inputs_statics, genome_fpath, gc_hv, gc_avg_T, delta):
        """normalization by gc-content, applied on both inputs and output data"""
        if inputs_statics and gc_hv is None:
            print("Compute gc factors including gv_hv and gc_avg_T", file=sys.stderr)
            gc_hv, gc_avg_T = [], []
            for i in range(self.dim[0]):
                gc_hv.append([])
                gc_avg_T.append([])
                for j in range(self.dim[1]):
                    hv, avg_T = get_gc_factor(self.inputs_covs[i][j], delta, genome_fpath)
                    gc_hv[i].append(hv)
                    gc_avg_T[i].append(avg_T)
        for i in range(self.dim[0]):
            for j in range(self.dim[1]):
                self.covs[i][j].normalization_by_gc_content(gc_hv[i][j], gc_avg_T[i][j], genome_fpath, delta)
                self.inputs_covs[i][j].normalization_by_gc_content(gc_hv[i][j], gc_avg_T[i][j], genome_fpath, delta)
        self.hv = gc_hv
        self.avg_T = gc_avg_T
        return gc_hv, gc_avg_T

    def normalization_by_input(self, signal_statics, inputs_statics, name, factors_inputs):
        """Normalize input-DNA. Use predefined factors or follow Diaz et al, 2012"""
        if factors_inputs:
            if configuration.VERBOSE:
                print("Use with predefined factors", file=sys.stderr)
            for i in range(signal_statics['dim'][0]):
                for j in range(signal_statics['dim'][1]):
                    self.inputs_covs[i][j].sm_scale(factors_inputs[i][j])
                    self.covs[i][j].sm_subtract(self.inputs_covs[i][j])
        elif inputs_statics:
            factors_inputs = []
            print("Compute factors", file=sys.stderr)

            for i in range(signal_statics['dim'][0]):
                factors_inputs.append([])
                for j in range(signal_statics['dim'][1]):
                # get_normalization_factor is 1-to-1 so, we use only covs, step_times is the times of original coverage bin size
                    _, factor = get_normalization_factor_by_cov(self.covs[i][j], self.inputs_covs[i][j], step_times=3, zero_counts=0, \
                                                    filename=name + '-norm' + str(i), debug=configuration.DEBUG,two_samples=False)
                    if factor:
                        print("Normalize input of Signal %s, Rep %s with factor %s" \
                              % (i, j, round(factor, configuration.ROUND_PRECISION)), file=sys.stderr)
                        self.inputs_covs[i][j].sm_scale(factor)
                        ## this is where we should look into the codes.... If after doing inputs, all data turn into zeros...
                        self.covs[i][j].sm_subtract(self.inputs_covs[i][j])
                        factors_inputs[i].append(factor)

        self.factors_inputs = factors_inputs
        return factors_inputs

    def normalization_by_signal(self, name, factors_ip, signal_statics, housekeeping_genes, tracker, norm_regionset, report,
                                 m_threshold, a_threshold):
        """Normalize signal. comparing data to data"""
        if not factors_ip and housekeeping_genes:
            print('Use housekeeping gene approach', file=sys.stderr)
            factors_ip, _ = norm_gene_level(signal_statics, housekeeping_genes, name, verbose=True, folder = self.FOLDER_REPORT, report=report)
        elif not factors_ip:
            if norm_regionset:
                print('Use TMM approach based on peaks', file=sys.stderr)
                # what to do with norm_regionset??
                # norm_regionset_coverage = self.init_overall_coverage(cov_strand=True) #TMM approach based on peaks
                # factors_ip = get_sm_norm_TMM_factor(norm_regionset_coverage,m_threshold, a_threshold)
            else:
                print('Use global TMM approach ', file=sys.stderr)
                factors_ip = get_sm_norm_TMM_factor(self.overall_coverage, m_threshold, a_threshold) #TMM approach
        
        if factors_ip:
            for i in range(self.dim[0]):
                for j in range(self.dim[1]):
                    self.covs[i][j].sm_scale(factors_ip[i][j])
                    self.overall_coverage['data'][i][j].data *= factors_ip[i][j]
                    if configuration.DEBUG:
                        print('Use scaling factor %s' %round(factors_ip[i][j], configuration.ROUND_PRECISION), file=sys.stderr)
        self.factors_ip = factors_ip
        return factors_ip

    def _index2coordinates(self, index):
        """Translate index within coverage array to genomic coordinates."""
        iter = self.genomicRegions.__iter__()
        r = iter.next()
        sum = r.final
        last = 0
        i = 0
        while sum <= index * self.stepsize:
            last += len(self.covs[0].coverage[i])
            try:
                r = iter.next()
            except StopIteration:
                sum += r.final
                i += 1
                break
            sum += r.final
            i += 1
        
        return r.chrom, (index-last) * self.stepsize, \
            min((index-last) * self.stepsize + self.stepsize, r.final)
                              
    def __len__(self):
        """Return number of observations."""
        return len(self.indices_of_interest)
    
    def get_observation(self, mask=np.array([])):
        """Return indices of observations. Do not consider indices contained in <mask> array"""
        mask = np.asarray(mask)
        if not mask.size:
            mask = np.array([True]*self._get_bin_number())
        return np.asarray(np.concatenate((self.overall_coverage[0][:,mask].T, self.overall_coverage[1][:,mask].T), axis=1))

    def _compute_sm_score(self):
        """Compute score for each observation (based on Xu et al.), but firstly only to get non-zeros numbers, which we have already done it"""
        # we use sparse-matrix to get scores
        # return a sum list of rates for each bin compared with average bins..
        # when it is is for whole overall_coverage..
        # Firstly to get average of two sparse matrix for one signal
        # we see overall_coverage is alreday in saprse matrix
        self.scores = None
        for i in range(self.dim[0]):
            # get signal_avg for all data
            signal_avg = sum([self.overall_coverage['data'][i][j].mean() for j in range(self.dim[1])])
            # use signal_avg to compare to each element in overall_coverage columns
            signal_rate = sum(self.overall_coverage['data'][i])/(self.dim[1] * (signal_avg))
            if self.scores is None:
                self.scores = signal_rate
            else:
                self.scores += signal_rate

    def _get_bin_number(self):
        """Return number of bins"""
        return self.overall_coverage[0].shape[1]

    def compute_sm_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows:
        - score must be > 0, i.e. everthing
        - overall coverage in library 1 and 2 must be > 3"""
        try:
            self._compute_sm_score()
            # threshold = 2.0 / (self.scores.shape[0])  # before it's considered if it's zero, now add some thresholds.
            threshold = 0.0
            # if we use the method before , we could see, all scores are from data greater than 0
            # second condition, we need average reads in one column bigger than l
            signal_sum = 0
            for i in range(self.dim[0]):
                signal_sum += sum(self.overall_coverage['data'][i])

            # signal_avg_sum /= (self.dim[1])
            self.indices_of_interest = np.intersect1d((self.scores > threshold).indices, (signal_sum > l*self.dim[1]).indices)  # 2/(m*n) thres = 2 /(self.scores.shape[0])
        except:
            self.indices_of_interest = None

    def write_test_samples(self, name, l):
        f = open(name, 'w')
        
        for el1, el2 in l:
            print(el1, el2, sep='\t', file=f)
        f.close()
    
    def output_training_set(self, name, training_set, s0_v, s1_v, s2_v):
        """Output debug info for training_set computation."""
        f=open(name + '-trainingset.bed', 'w')
        for l in training_set:
            chrom, s, e = self._index2coordinates(l)
            print(chrom, s, e, sep ='\t', file=f)
        f.close()
        
        self.write_test_samples(name + '-s0', s0_v)
        self.write_test_samples(name + '-s1', s1_v)
        self.write_test_samples(name + '-s2', s2_v)
    
    def get_training_set(self, test, name, threshold, min_t, y=1000, ex=2):
        """Return HMM's training set (max <y> positions). Enlarge each contained bin by <ex>.
           If first sample can't represent data, we need to resample it from population, self.indices_of_interest..
           Other way, we could build the samples for training, but then it will cause other troubles, maybe...
           we need to filter data, make diff_cov firstly non zeros and then get half part from it..s

           If we have exp_data, it's actually the same here, so we don't need more information
        """
        ## whatever we get, we need to transform it into integer value??
        #  then we estimate HMM distribution; if it is float value; could we get it ?
        # firstly to get the differences of cov_avg values and then use it to count diff...
        # so here, we maybe need the cov_avg to provide convenience; One is for score, indices_of_interest
        signal_diff = None
        for i in range(self.dim[0]):
            if signal_diff is None:
                signal_diff = sum(self.overall_coverage['data'][i])
            else:
                signal_diff = signal_diff - sum(self.overall_coverage['data'][i])
                signal_diff.data = np.abs(signal_diff.data)

        # order of list doesn't matter
        diff_cov = int(np.percentile(filter(lambda x: x > 0, signal_diff.data), min_t)/self.dim[1]) # get the mean of each

        if test:
            diff_cov, threshold = 2, 1.5
        
        if configuration.DEBUG:
            print('Training set parameters: threshold: %s, diff_cov: %s' %(threshold, diff_cov), file=sys.stderr)

        # here we append s0, s1 and s2;; So if once failure, but it append it into that; which is not right!!
        # compute training set parameters, re-compute training set if criteria do not hold

        done =False
        while not done:
            s0, s1, s2, tmp = [], [], [], []
            if diff_cov == 1 and threshold == 1.1:
                print("No differential peaks detected", file=sys.stderr)
                sys.exit()
            steps = 0
            for idx in sample(self.indices_of_interest, min(y, len(self.indices_of_interest))):
                cov1, cov2 = self._get_sm_covs(idx)
                steps += 1
                #apply criteria for initial peak calling
                if ((cov1 +1 ) / (float(cov2) + 1) > threshold and cov1+cov2 > diff_cov/2) or cov1-cov2 > diff_cov:
                    s1.append((idx, cov1, cov2))
                elif ((cov1 + 1) / (float(cov2) + 1) < 1/threshold and cov1+cov2 > diff_cov/2) or cov2-cov1 > diff_cov:
                    s2.append((idx, cov1, cov2))
                else:
                    s0.append((idx, cov1, cov2))
            
                if steps % 500 == 0 and len(s0) > y and len(s1) > y and len(s2) > y:
                    tmp = []
                    for el in [s0, s1, s2]:
                        el = np.asarray(el)
                        if not test:
                            el = el[
                                np.logical_and(
                                    el[:, 1] < np.percentile(el[:, 1], 97.5) * (el[:, 1] > np.percentile(el[:, 1], 2.5)),
                                    el[:, 2] < np.percentile(el[:, 2], 97.5) * (el[:, 2] > np.percentile(el[:, 2], 2.5))),:]

                        tmp.append(el)

                    l = np.min([len(tmp[0]), len(tmp[1]), len(tmp[2])])
                    if l >= y:
                        break
            
            if len(s1) < 100/2 and len(s2) > 2*100:
                s1 = map(lambda x: (x[0], x[2], x[1]), s2)
            if len(s2) < 100/2 and len(s1) > 2*100:
                s2 = map(lambda x: (x[0], x[2], x[1]), s1)
            
            if len(s1) < 100 or len(s2) < 100:
                diff_cov -= 15
                threshold -= 0.1
                diff_cov = max(diff_cov, 1)
                threshold = max(threshold, 1.1)
            else:
                done = True
        
        if configuration.DEBUG:
            print('Final training set parameters: threshold: %s, diff_cov: %s' %(threshold, diff_cov), file=sys.stderr)
        
        #optimize training set, extend each bin
        # here we need to combine data to sample all data, if they are not meet our requirements, we need to sample from all population data
        # and the population data are directly from interest_of_points..
        if tmp == [] :
            for el in [s0, s1, s2]:
                el = np.asarray(el)
                if not test:
                    el = el[np.where(
                        np.logical_and(el[:, 1] < np.percentile(el[:, 1], 97.5) * (el[:, 1] > np.percentile(el[:, 1], 2.5)),
                                       el[:, 2] < np.percentile(el[:, 2], 97.5) * (el[:, 2] > np.percentile(el[:, 2], 2.5))))]

                tmp.append(el)

            l = np.min([len(tmp[0]), len(tmp[1]), len(tmp[2]), y])

        print('the smallest length l is %d between %d, %d, %d, %d'%(l,y, len(tmp[0]),len(tmp[1]),len(tmp[2])))
        s0 = sample(tmp[0], l)
        s1 = sample(tmp[1], l)
        s2 = sample(tmp[2], l)

        tmp2 = []
        for i, ss in enumerate([s0, s1, s2]):
            while np.any(np.sum(ss, axis=0) < len(ss)):
                print('resample because data is not spatial')
                ss = sample(tmp[i], l)
            tmp2.append(ss)

        s0_v = map(lambda x: (x[1], x[2]), tmp2[0])
        s1_v = map(lambda x: (x[1], x[2]), tmp2[1])
        s2_v = map(lambda x: (x[1], x[2]), tmp2[2])

        
        extension_set = set()
        for i, _, _ in s0 + s1 + s2:
            for j in range(max(0, i - ex), i + ex + 1): #extend bins
                extension_set.add(j)
         
        tmp = s0 + s1 + s2
        training_set = map(lambda x: x[0], tmp) + list(extension_set)
         
        training_set = list(training_set)
        training_set.sort()
        
        if configuration.DEBUG:
            self.output_training_set(name, training_set, s0_v, s1_v, s2_v)
        
        return training_set, s0_v, s1_v, s2_v