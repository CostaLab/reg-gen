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
from normalize import get_normalization_factor_by_cov, get_norm_TMM_factor
from DualCoverageSet import DualCoverageSet
from norm_genelevel import norm_gene_level
from rgt.CoverageSet import CoverageSet, get_gc_context
import configuration


class MultiCoverageSet(DualCoverageSet):

    def _help_init(self, signal_statics, inputs_statics, region_giver, norm_regionset, rmdup, binsize, stepsize, strand_cov):
        """Return self.covs and self.inputs as CoverageSet
        self._help_init(signal_statics, inputs_statics,region_giver, norm_regionset, rmdup, binsize, stepsize, strand_cov = strand_cov)
        But before we need to do statistics about the file, get information, how much read for this regions are in this fields..
        Better we need to do is get all info of all fields, and extract data from it to combine the training fields.
        Later we use the paras to estimate peaks.. Maybe according to regions
        mask_file is used to filter data before we do other operations.
        """
        # here we make covs into 2-D list
        self.covs = []
        self.covs_avg = []
        for i in range(signal_statics['dim'][0]):
            self.covs.append([])
            self.covs_avg.append(CoverageSet('cov_avg' + str(i), region_giver.valid_regionset))
            for j in range(signal_statics['dim'][1]):
                self.covs[i].append(CoverageSet('file_' + str(i)+'_'+str(j), region_giver.valid_regionset))
                self.covs[i][j].coverage_from_bam(bam_file=signal_statics['data'][i][j]['fname'], extension_size=signal_statics['data'][i][j]['extension_size'], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize, mask_file=region_giver.mask_file, get_strand_info=strand_cov)

        if inputs_statics:
            self.inputs_covs = []
            self.inputs_covs_avg = []
            for i in range(inputs_statics['dim'][0]):
                self.inputs_covs.append([])
                self.inputs_covs_avg.append(CoverageSet('inputs_cov_avg' + str(i), region_giver.valid_regionset))
                for j in range(inputs_statics['dim'][1]):
                    self.inputs_covs[i].append(CoverageSet('file_' + str(i) + '_' + str(j), region_giver.valid_regionset))
                    self.inputs_covs[i][j].coverage_from_bam(bam_file=inputs_statics['data'][i][j]['fname'],extension_size=inputs_statics['data'][i][j]['extension_size'],
                                                      rmdup=rmdup, binsize=binsize, stepsize=stepsize, mask_file=region_giver.mask_file, get_strand_info=strand_cov)
        else:
            self.inputs_covs = []
        """
        if norm_regionset:
            self.norm_regions = [CoverageSet('norm_region' + str(i), norm_regionset) for i in range(dim)]
            for i, c in enumerate(self.norm_regions):
                c.coverage_from_bam(bam_file=path_bamfiles[i], extension_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                    stepsize=stepsize, mask_file=region_giver.mask_file, get_strand_info = strand_cov)
            self.input_avg = [CoverageSet('input_avg'  + str(i), region_giver.regionset) for i in range(2)]
        else:
            self.norm_regions = None
        """

    def _get_covs(self, DCS, i):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i"""
        cov1 = int(np.mean(DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]))

        return cov1, cov2
    
    def _compute_gc_content(self, no_gc_content, inputs_statics, stepsize, binsize, genome_path, name, region_giver):
        """Compute GC-content, please pay attension to dimension changes!!!"""
        if not no_gc_content and inputs_statics and self.gc_content_cov is None:
            print("Compute GC-content", file=sys.stderr)
            for i in range(self.dim[0]):
                for j in range(self.dim[1]):
                    inputs_cov = self.inputs_covs[i][j] #1 to 1 mapping between input and cov
                    self.gc_content_cov, self.avg_gc_content, self.gc_hist = get_gc_context(stepsize, binsize, genome_path, inputs_cov.coverage, region_giver.valid_chrom_sizes)
                    self._norm_gc_content(inputs_cov[i][j].coverage, self.gc_content_cov, self.avg_gc_content)
                    self._norm_gc_content(inputs_cov.coverage, self.gc_content_cov, self.avg_gc_content)

                #if VERBOSE:
                #    self.print_gc_hist(name + '-s%s-rep%s-' %(sig, rep), gc_hist)
                #    cov.write_bigwig(name + '-s%s-rep%s-gc.bw' %(sig, rep), chrom_sizes)

    def _output_input_bw(self, name, chrom_sizes, save_wig):
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

    def _help_get_data(self, i, type):
        if type != 'normregion':
            for j in range(len(self.covs[i].genomicRegions)):
                if type == 'cov':
                    yield self.covs[i].coverage[j]
                elif type == 'strand':
                    yield self.covs[i].cov_strand_all[j]
        elif type == 'normregion':
            for j in range(len(self.norm_regions[i].genomicRegions)):
                yield self.norm_regions[i].coverage[j]
    
    def _help_init_overall_coverage(self, cov_strand=True):
        """Convert coverage data (and optionally strand data) to matrix list"""
        # overall_coverage format are [[coverage_0_0 , coverage_0_1 ], [ coverage_1_0, coverag_1_1 ]]
        overall_coverage = []
        for i in range(self.dim[0]):
            overall_coverage.append([])
            for j in range(self.dim[1]):
                # make it use sparse matrix
                overall_coverage[i].append(self.covs[i][j].sm_overall_cov)
        if cov_strand:
            overall_coverage_strand = []
            for i in range(self.dim[0]):
                overall_coverage_strand.append([])
                for j in range(self.dim[1]):
                    overall_coverage_strand[i].append(self.covs[i][j].cov_strand_all)
            return np.asarray(overall_coverage), np.asarray(overall_coverage_strand)
        else:
            return np.asarray(overall_coverage)

    def count_positive_signal(self):
        """get num of all read in coverage; to judge the data validation, if it's too small then we dispose it for training,
        maybe we don't need it; if wa can be sure that training samples are over it, from statistics data!!!"""
        return np.sum([self.covs[i][j].coverage for j in range(self.dim[1]) for i in range(self.dim[0])])
    
    def __init__(self, name, region_giver, genome_path, binsize, stepsize, norm_regionset, \
                 verbose, debug, no_gc_content, rmdup, signal_statics, inputs_statics, \
                 factors_inputs, scaling_factors_ip, save_wig, strand_cov, housekeeping_genes,\
                 tracker, end, counter, gc_content_cov=None, avg_gc_content=None, gc_hist=None, output_bw=True,\
                 folder_report=None, report=None, save_input=False, m_threshold=80, a_threshold=95, ignored_regions=None):
        """Compute CoverageSets, GC-content and normalize input-DNA and IP-channel"""
        # one improvement is to make the left one key_word parameter and we parse it, not like this, all in a list
        """    
        regionset = region_giver.regionset
        chrom_sizes = region_giver.chrom_sizes_file
        chrom_sizes_dict = region_giver.get_chrom_dict()
        """
        self.region_giver = region_giver
        binsize = 1000
        stepsize = 500
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.dim = signal_statics['dim']
        self.gc_content_cov = gc_content_cov
        self.avg_gc_content = avg_gc_content
        self.gc_hist = gc_hist
        self.scaling_factors_ip = scaling_factors_ip
        self.factors_inputs = factors_inputs
        self.end = end
        self.counter = counter # use of counter ???
        self.no_data = False
        self.FOLDER_REPORT = folder_report

        configuration.DEBUG = debug
        configuration.VERBOSE = verbose
        
        #make data nice
        # we can put help_init codes here and for compute_gc_content or normalization we see it as a process and do it outside the class
        # This is only data; MCV model
        self._help_init(signal_statics, inputs_statics,region_giver, norm_regionset, rmdup, binsize, stepsize, strand_cov=strand_cov)
        #if self.count_positive_signal() < 1:
        #    self.no_data = True
        #    return None
        self._compute_gc_content(no_gc_content, inputs_statics, stepsize, binsize, genome_path, name, region_giver)
        self._normalization_by_input(signal_statics, inputs_statics, name, factors_inputs, save_input)
        if save_input:
            self._output_input_bw(name, region_giver.chrom, save_wig)
            
        self.overall_coverage, self.overall_coverage_strand = self._help_init_overall_coverage(cov_strand=True)

        # much complex, so we decay to change it
        self._normalization_by_signal(name, scaling_factors_ip, signal_statics, housekeeping_genes, tracker, norm_regionset, report,
                                      m_threshold, a_threshold)

        ## After this step, we have already normalized data, so we could output normalization data
        if output_bw:
            self._output_bw(name, region_giver.chrom_sizes_file, save_wig, save_input)

        # this step we could use sparse matrix, it shows automatically the indices of non-zeros bins
        # but one thing is indices_of_interest needs more processes;
        # self.scores = np.zeros(len(self.covs[0][0].overall_cov))

        self.indices_of_interest = []
    
    def get_max_colsum(self):
        """Sum over all columns and add maximum"""
        # how to change it to self.cov ???
        return self.overall_coverage[0].sum(axis=0).max() + self.overall_coverage[1].sum(axis=0).max()
    
    def output_overall_coverage(self, path):
        for j in range(2):
            f = open(path + str(j), 'w')
            for i in range(self.overall_coverage[j].shape[1]):
                print(self.overall_coverage[j][:,i].T, file=f)
    
    def _normalization_by_input(self, signal_statics, inputs_statics, name, factors_inputs, save_input):
        """Normalize input-DNA. Use predefined factors or follow Diaz et al, 2012"""
        
        if configuration.VERBOSE:
            print("Normalize input-DNA", file=sys.stderr)
        
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

    def _normalization_by_signal(self, name, factors_ip, signal_statics, housekeeping_genes, tracker, norm_regionset, report,
                                 m_threshold, a_threshold):
        """Normalize signal. comparing data to data"""
        
        if configuration.VERBOSE:
            print('Normalize ChIP-seq profiles', file=sys.stderr)
        
        if not factors_ip and housekeeping_genes:
            print('Use housekeeping gene approach', file=sys.stderr)
            factors_ip, _ = norm_gene_level(signal_statics, housekeeping_genes, name, verbose=True, folder = self.FOLDER_REPORT, report=report)
        elif not factors_ip:
            if norm_regionset:
                print('Use TMM approach based on peaks', file=sys.stderr)
                norm_regionset_coverage = self._help_init_overall_coverage(cov_strand=False) #TMM approach based on peaks
                factors_ip = get_norm_TMM_factor(norm_regionset_coverage,m_threshold, a_threshold)
            else:
                print('Use global TMM approach ', file=sys.stderr)
                factors_ip = get_norm_TMM_factor(self.overall_coverage, m_threshold, a_threshold) #TMM approach
        
        if factors_ip:
            for i in range(self.dim[0]):
                for j in range(self.dim[1]):
                    self.covs[i][j].sm_scale(factors_ip[i][j])
                    # this is actually one redundant step, if overall_coverage[i][j] from covs[i][j]
                    # so if we change covs[i][j], it should change.. cause overall_coverage is not Coverage Set, so we can't use method sm_scale
                    self.overall_coverage[i][j].data *= factors_ip[i][j]
                    if configuration.DEBUG:
                        print('Use scaling factor %s' %round(factors_ip[i][j], configuration.ROUND_PRECISION), file=sys.stderr)
        
        self.factors_ip = factors_ip

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
    
    def _compute_score(self):
        """Compute score for each observation (based on Xu et al.)"""
        # after np.squeeze, we remove single-dimensional entry.. What does it make ??? seems nothing about process
        # interest_region is column scores are greater than one values...
        # self.scores = np.sum(np.asarray([np.squeeze(np.asarray(self.overall_coverage[i][j])) /float(np.sum(self.overall_coverage[i][j])) for j in xrange(self.dim_2) for i in range(self.dim_1)]), axis=0)/self.dim_2
        # old methods to count interest regions
        self.scores = sum([np.squeeze(np.asarray(np.mean(self.overall_coverage[i], axis=0))) / float(np.mean(self.overall_coverage[i])) for i in range(2)])

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
            signal_avg = sum([self.overall_coverage[i][j].mean() for j in range(self.dim[1])])
            # use signal_avg to compare to each element in overall_coverage columns
            signal_rate = sum(self.overall_coverage[i])/(self.dim[1] * (signal_avg))
            if self.scores is None:
                self.scores = signal_rate
            else:
                self.scores += signal_rate


    def _get_bin_number(self):
        """Return number of bins"""
        return self.overall_coverage[0].shape[1]
    
    def compute_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows: 
        - score must be > 0, i.e. everthing
        - overall coverage in library 1 and 2 must be > 3"""
        
        try:
            self._compute_score()
            # threshold = 2.0 / (self.scores.shape[0])  # before it's considered if it's zero, now add some thresholds.
            threshold = 0.0
            self.indices_of_interest = np.where(self.scores > threshold)[0] #2/(m*n) thres = 2 /(self.scores.shape[0])
            tmp = np.where(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) + np.squeeze(np.asarray(np.mean(self.overall_coverage[1], axis=0))) > 10)[0]
            tmp2 = np.intersect1d(self.indices_of_interest, tmp)
            self.indices_of_interest = tmp2
        except:
            self.indices_of_interest = None

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
            signal_avg = sum([self.overall_coverage[i]/self.dim[1] for i in range(self.dim[0])])
            self.indices_of_interest = np.intersect1d((self.scores > threshold).indices, (signal_avg > l).indices)  # 2/(m*n) thres = 2 /(self.scores.shape[0])
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
    
    def get_training_set(self, test, exp_data, name, foldchange, min_t, y=1000, ex=2):
        """Return HMM's training set (max <y> positions). Enlarge each contained bin by <ex>.
           If first sample can't represent data, we need to resample it from population, self.indices_of_interest..
           Other way, we could build the samples for training, but then it will cause other troubles, maybe...
           we need to filter data, make diff_cov firstly non zeros and then get half part from it..s
        """
        threshold = foldchange
        diff_cov = int(np.percentile(filter(lambda x: x>0, np.abs(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) - \
                                            np.squeeze(np.asarray(np.mean(self.overall_coverage[1], axis=0))))), min_t))

        if test:
            diff_cov, threshold = 2, 1.5
        
        if configuration.DEBUG:
            print('Training set parameters: threshold: %s, diff_cov: %s' %(threshold, diff_cov), file=sys.stderr)
        
        s0, s1, s2 , tmp = [], [], [], []
        
        # compute training set parameters, re-compute training set if criteria do not hold
        print('The length of indices_of_interest')
        print(len(self.indices_of_interest))
        rep=True
        while rep:

            if diff_cov == 1 and threshold == 1.1:
                print("No differential peaks detected", file=sys.stderr)
                sys.exit()
            steps = 0
            for i in sample(range(len(self.indices_of_interest)), len(self.indices_of_interest)):
                cov1, cov2 = self._get_covs(exp_data, i)
                steps += 1
                #apply criteria for initial peak calling
                if ((cov1 +1 ) / (float(cov2) + 1) > threshold and cov1+cov2 > diff_cov/2) or cov1-cov2 > diff_cov:
                    s1.append((self.indices_of_interest[i], cov1, cov2))
                elif ((cov1 + 1) / (float(cov2) + 1) < 1/threshold and cov1+cov2 > diff_cov/2) or cov2-cov1 > diff_cov:
                    s2.append((self.indices_of_interest[i], cov1, cov2))
                else:
                    s0.append((self.indices_of_interest[i], cov1, cov2))
            
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
                rep = False
        
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
