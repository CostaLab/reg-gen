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
import os
from normalize import *
from norm_genelevel import norm_gene_level
from rgt.CoverageSet import CoverageSet
import configuration


class MultiCoverageSet(): # DualCoverageSet
    """Not inherit from DualCoveragSet, instead we can use it to represent DualCoverageSet"""

    def _init_coverage(self, statics, valid_regionset, mask_file, rmdup, binsize, stepsize, strand_cov,use_sm=False):
        """Return covs in list as CoverageSet for one statics
        self._help_init(signal_statics, inputs_statics,regionset, norm_regionset, rmdup, binsize, stepsize, strand_cov = strand_cov)
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
                    covs[i].append(CoverageSet('file_' + str(i)+'_'+str(j), valid_regionset))
                    covs[i][j].coverage_from_bam(bam_file=statics['data'][i][j]['fname'], extension_size=statics['data'][i][j]['extension_size'], rmdup=rmdup, binsize=binsize, \
                                    stepsize=stepsize, mask_file=mask_file, get_strand_info=strand_cov)
                    if use_sm: # maybe one thiing is how to deal with different operation on sparse matrix
                        tmp = [cov_to_smatrix(chrom_cov) for chrom_cov in covs[i][j].coverage]
                        covs[i][j].sm_cov_strand_all = [cov_to_smatrix(chrom_cov_strand) for chrom_cov_strand in
                                                        covs[i][j].cov_strand_all]

                        del covs[i][j].coverage, covs[i][j].coverageorig, covs[i][j].overall_cov, covs[i][j].cov_strand_all
                        # want to save memory spaces but should we also delete coveragerig??
                        covs[i][j].coverage = tmp
                        covs[i][j].sm_overall_cov = sparse.hstack(covs[i][j].coverage, format='csr')

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
                self.overall_coverage['data'][i].append(self.covs[i][j].sm_overall_cov)
        if strand_cov:
            self.overall_strand_coverage = {'dim':self.dim, 'data':[]}
            for i in range(self.dim[0]):
                self.overall_strand_coverage['data'].append([])
                for j in range(self.dim[1]):
                    self.overall_strand_coverage['data'][i].append(self.covs[i][j].sm_cov_strand_all)

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

    def __init__(self, name, regionset, mask_file, binsize, stepsize, rmdup, signal_statics, inputs_statics, strand_cov, hv=None, avg_T=None, use_sm=False):
        """Compute CoverageSets, GC-content and normalize input-DNA and IP-channel"""
        # one improvement is to make the left one key_word parameter and we parse it, not like this, all in a list
        self.name = name
        self.regionset = regionset
        self.maskfile = mask_file
        self.binsize = binsize
        self.stepsize = stepsize
        self.rmdup = rmdup # bool if remove duplicate

        self.dim = signal_statics['dim']
        self.hv = hv
        self.avg_T = avg_T
        self.data_valid = True

        # here we need to judge if the coverage fine, or not; Actually before we have read statitics data,so it's fine
        # could give info about doing what; reading signal files, reading inputs files
        if signal_statics:
            print("Begin reading CHIP bamfiles", file=sys.stderr)
            self.covs = self._init_coverage(signal_statics, regionset, mask_file, rmdup, binsize, stepsize, strand_cov=strand_cov, use_sm=use_sm)
            print("End reading CHIP bamfiles", file=sys.stderr)

        if inputs_statics:
            print("Reading inputs bamfiles", file=sys.stderr)
            self.inputs_covs = self._init_coverage(inputs_statics, regionset, mask_file, rmdup, binsize, stepsize, strand_cov=strand_cov, use_sm=use_sm)
            print("End reading CHIP bamfiles", file=sys.stderr)
        else:
            self.inputs_covs = None

        if not self._is_cov_valid(signal_statics, self.covs) or not self.covs:
            self.data_valid = False
            return None
        if not self._is_cov_valid(inputs_statics, self.inputs_covs):
            self.data_valid = False
            return None

    def normalization_by_gc_content(self, inputs_statics, genome_fpath, gc_hv, gc_avg_T, delta):
        """normalization by gc-content, applied on both inputs and output data
        actually we need to do it in class CoverageSet , but now we need to build sth on it, so we need to change methods
        """
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
                cov_normalization_by_gc_content(self.covs[i][j], gc_hv[i][j], gc_avg_T[i][j], genome_fpath, delta)
                cov_normalization_by_gc_content(self.inputs_covs[i][j], gc_hv[i][j], gc_avg_T[i][j], genome_fpath, delta)
        self.hv = gc_hv
        self.avg_T = gc_avg_T
        return gc_hv, gc_avg_T

    def normalization_by_input(self, signal_statics, inputs_statics, name, factors_inputs):
        """Normalize input-DNA. Use predefined factors or follow Diaz et al, 2012"""

        if inputs_statics and not factors_inputs:
            factors_inputs = []
            print("Compute factors", file=sys.stderr)

            for i in range(signal_statics['dim'][0]):
                factors_inputs.append([])
                for j in range(signal_statics['dim'][1]):
                # get_normalization_factor is 1-to-1 so, we use only covs, step_times is the times of original coverage bin size
                    _, factor = get_normalization_factor_by_cov(self.covs[i][j], self.inputs_covs[i][j], step_times=3, zero_counts=0, \
                                                    filename=name + '-norm' + str(i), debug=configuration.DEBUG,two_samples=False)
                    factors_inputs[i].append(factor)

        if factors_inputs:
            for i in range(signal_statics['dim'][0]):
                for j in range(signal_statics['dim'][1]):
                    if configuration.VERBOSE:
                        print("Normalize input of Signal %s, Rep %s with factor %s" \
                              % (i, j, round(factors_inputs[i][j], configuration.ROUND_PRECISION)), file=sys.stderr)
                    sm_scale(self.inputs_covs[i][j], factors_inputs[i][j])
                    sm_subtract(self.covs[i][j],self.inputs_covs[i][j])

        return factors_inputs

    def normalization_by_signal(self, name, factors_ip, signal_statics, housekeeping_genes, report,
                                 m_threshold, a_threshold):
        """Normalize signal. comparing data to data"""
        if not factors_ip and housekeeping_genes:
            print('Use housekeeping gene approach', file=sys.stderr)
            factors_ip, _ = norm_gene_level(signal_statics, housekeeping_genes, name, verbose=True, folder=configuration.FOLDER_REPORT, report=report)

        elif not factors_ip:
            print('Use global TMM approach ', file=sys.stderr)
            factors_ip = get_sm_norm_TMM_factor(self.overall_coverage, m_threshold, a_threshold) #TMM approach
        
        if factors_ip:
            for i in range(self.dim[0]):
                for j in range(self.dim[1]):
                    if configuration.VERBOSE:
                        print('Use scaling factor %s' %round(factors_ip[i][j], configuration.ROUND_PRECISION), file=sys.stderr)

                    sm_scale(self.covs[i][j], factors_ip[i][j])
                    self.overall_coverage['data'][i][j].data *= factors_ip[i][j]
        self.factors_ip = factors_ip
        return factors_ip

    def sm_index2coordinates(self, index):
        """how to get genome information only depends one index of overall_coverage_information
        what we can depend on is self.regionset to build connection to it;
        how to test it ??
        """
        # judge order and len of chroms in regionset or we can use information of covs to get it
        # judge idx in which covs field, result is order; less than or greater than covs.shape[-1]
        order = 0 # len(self.covs[0][0].coverage)
        index -= self.covs[0][0].coverage[order].shape[-1]
        while index >= 0:
            order += 1
            index -= self.covs[0][0].coverage[order].shape[-1]
        # when it stops index < 0
        index += self.covs[0][0].coverage[order].shape[-1]
        # get regionset[order] and then use idx to subtract it len, we get initial idx for it
        chrm_name = self.regionset[order].chrom
        # use binsize and step.size to get i
        start = index * self.stepsize
        end = min(start + self.binsize, self.regionset[order].final)
        return chrm_name, start, end

    def get_sm_covs(self, indices, strand_cov=False):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i
        strand_cov : True also return strand_cov information, False, only coverage information
        _get_sm_covs: return all covs according to one idx, [0_sample_0, 0_sample_1, 1_sample_0,1_sample_1]
        """
        cov1 = [np.squeeze(self.overall_coverage['data'][0][j][:, indices].toarray()) for j in range(self.dim[1])]
        cov2 = [np.squeeze(self.overall_coverage['data'][1][j][:, indices].toarray()) for j in range(self.dim[1])]
        if strand_cov:
            strand_cov1 = [np.squeeze(self.overall_strand_coverage['data'][0][j][indices,:].toarray()) for j in range(self.dim[1])]
            strand_cov2 = [np.squeeze(self.overall_strand_coverage['data'][1][j][indices,:].toarray()) for j in range(self.dim[1])]
            return [cov1, cov2], [strand_cov1, strand_cov2]
        else:
            return [cov1, cov2]

    def _compute_sm_score(self):
        """Compute score for each observation (based on Xu et al.), but firstly only to get non-zeros numbers, which we have already done it
            # return a sum list of rates for each bin compared with average bins..
            # when it is is for whole overall_coverage..
            # Firstly to get average of two sparse matrix for one signal
        """
        self.scores = None
        for i in range(self.dim[0]):
            # get signal_avg for all data in one sample # here maybe sth wrong
            signal_sum = sum([self.overall_coverage['data'][i][j].sum() for j in range(self.dim[1])])
            # use signal_avg to compare to each element in overall_coverage columns
            signal_rate = sum(self.overall_coverage['data'][i])/signal_sum
            if self.scores is None:
                self.scores = signal_rate
            else:
                self.scores += signal_rate

    def compute_sm_putative_region_index(self, l=3, eta=1):
        """Compute putative differential peak regions as follows:
        - score must be > 0, i.e. everthing
        - overall coverage in library 1 and 2 must be > 3"""
        try:
            self._compute_sm_score()
            threshold = 2.0 / (eta*self.overall_coverage['data'][0][0].shape[-1])  # before it's considered if it's zero, now add some thresholds.
            # threshold = 0.0
            # if we use the method before , we could see, all scores are from data greater than 0
            # second condition, we need average reads in one column bigger than l
            signal_sum = 0
            for i in range(self.dim[0]):
                signal_sum += sum(self.overall_coverage['data'][i])
            self.indices_of_interest = np.intersect1d((self.scores > threshold).indices, (signal_sum > l*self.dim[1]).indices)  # 2/(m*n) thres = 2 /(self.scores.shape[0])
        except:
            self.indices_of_interest = None

    def output_input_bw(self, filename, chrom_sizes_file, save_wig):
        """output inputs as bigwig file"""
        for sig in range(self.dim[0]):
            for rep in range(self.dim[1]):

                tmp_path = filename + '-s%s-rep%s.bw' %(sig, rep) + '.wig'
                write_wig(self.inputs_covs[sig][rep], tmp_path)
                t = ['wigToBigWig', "-clip", tmp_path, chrom_sizes_file, filename]
                c = " ".join(t)
                os.system(c)
                if not save_wig:
                    os.remove(tmp_path)

    def output_signal_bw(self, filename, chrom_sizes_file, save_wig):
        """Output signal files as bigwig files"""
        for sig in range(self.dim[0]):
            for rep in range(self.dim[1]):
                # self.covs[sig][rep].write_bigwig(name + '-s%s-rep%s.bw' %(sig, rep), chrom_sizes_file, save_wig=save_wig)
                tmp_path = filename + '-s%s-rep%s.bw' %(sig, rep) + '.wig'
                write_wig(self.covs[sig][rep], tmp_path)
                t = ['wigToBigWig', "-clip", tmp_path, chrom_sizes_file, filename]
                c = " ".join(t)
                os.system(c)
                if not save_wig:
                    os.remove(tmp_path)


def cov_to_smatrix(cov):
    """change list of coverage into sparse matrix by using scipy"""
    cov_mtx = sparse.csr_matrix(cov, dtype=float)
    return cov_mtx


def sm_scale(cov, factor):
    """we scale sparse matrix representation with factor"""
    cov.sm_overall_cov *= factor


def sm_add(cov, sm_cs):
    """subtract coverage set in sparse matrix format"""
    # firstly test if they have same region; else not do it!!
    assert set(cov.genomicRegions.get_chrom()) == set(sm_cs.genomicRegions.get_chrom())
    # we use overall_cov, so not bother for one chrom and then another, we assume they are in the same order.
    # self.sm_overall_cov is 1-D array
    cov.sm_overall_cov += sm_cs.sm_overall_cov  # we could assume all values are non-negative


def sm_subtract(cov, sm_cs):
    """subtract coverage set in sparse matrix format"""
    # firstly test if they have same region; else not do it!!
    assert set(cov.genomicRegions.get_chrom()) == set(sm_cs.genomicRegions.get_chrom())
    cov.sm_overall_cov -= sm_cs.sm_overall_cov
    cov.sm_overall_cov.data = np.clip(cov.sm_overall_cov.data, 0,
                                      cov.sm_overall_cov.max())  # make negative into zeros
    cov.sm_overall_cov.eliminate_zeros()  # eliminate zeros elements in matrix


def cov_normalization_by_gc_content(cov, hv, avg_T, genome_fpath,delta=0.2):
    """After we get gc-factor from one input-coverage, and then we apply it on coverage from one signal file
    """
    genome_fasta = pysam.Fastafile(genome_fpath)
    # fetch genome w.r.t. chromsome
    chroms = cov.genomicRegions.get_chrom()

    # coverage separated by chroms regions; so we need to find out first what it belongs to
    for i in range(len(chroms)):
        for bin_idx in cov.coverage[i].indices:
            s = bin_idx * cov.stepsize
            e = s + cov.binsize
            genome_seq = genome_fasta.fetch(reference=chroms[i], start=s, end=e)
            prop = get_gc_content_proportion(genome_seq)
            # but we need to make prop in range of delta range
            prop_key = int(prop / delta) * delta
            if hv[prop_key]: # sometimes zeros occur, do not consider
                cov.coverage[i][:,bin_idx] *= (avg_T/hv[prop_key])


def isvalid_training_data(state_data, threshold):
    """test if data are valid or not valid for training set
    Arguments : ss is one state data, including cov1 and cov2
     Threshold: if the column sum of state_data is less than threshold * num(state_data)
    Results: return True or False w.r.t threshold
    """
    if np.any(np.sum(state_data, axis=0) < len(state_data) * threshold):
        return False
    else:
        return True


def transform_data_for_HMM(training_cov):
    """change training cov into trainig data format for HMM;
    training_cov is in format [[ array(), array()], [array() , array()]]
    training_data is in format [[0_sample_0, 0_sample_1, 1_sample_0,1_sample_1],[], [] ....]
    return training data
    """
    tmp = [zip(*training_cov[i]) for i in range(len(training_cov))]
    return np.concatenate(tmp,axis=1)


def get_training_set(exp_data, test, name, threshold, min_t, y=1000, ex=0):
    """Return HMM's training set (max <y> positions). Enlarge each contained bin by <ex>.
       If first sample can't represent data, we need to resample it from population, self.indices_of_interest..
       Other way, we could build the samples for training, but then it will cause other troubles, maybe...
       we need to filter data, make diff_cov firstly non zeros and then get half part from it..s
    ## return :
    training data sepcially for s0, s1 and s2
    training data directly for next step
    """
    # firstly to get the differences of cov_avg values and then use it to count diff...
    signal_diff = None
    for i in range(exp_data.dim[0]):
        if signal_diff is None:
            signal_diff = sum(exp_data.overall_coverage['data'][i])
        else:
            signal_diff = signal_diff - sum(exp_data.overall_coverage['data'][i])
            signal_diff.data = np.abs(signal_diff.data)

    # order of list doesn't matter
    diff_cov = np.percentile(filter(lambda x: x > 0, signal_diff.data), min_t) / exp_data.dim[1]  # get the mean of each

    if test:
        diff_cov, threshold = 2, 1.5
    # here we append s0, s1 and s2;; So if once failure, but it append it into that; which is not right!!
    # compute training set parameters, re-compute training set if criteria do not hold

    done = False
    while not done:
        s0, s1, s2, tmp = [], [], [], []
        if diff_cov == 1 and threshold == 1.1:
            print("No differential peaks detected", file=sys.stderr)
            sys.exit()

        for idx in sample(exp_data.indices_of_interest, min(y, len(exp_data.indices_of_interest))):
            covs_list = exp_data.get_sm_covs(idx)  # only return data not indices
            cov1 = np.mean(covs_list[0])  # use avg to present covs
            cov2 = np.mean(covs_list[1])
            # apply criteria for initial peak calling
            if ((cov1 + 1) / (float(cov2) + 1) > threshold and cov1 + cov2 > diff_cov / 2) or cov1 - cov2 > diff_cov:
                s1.append((idx, cov1, cov2))
            elif ((cov1 + 1) / (float(cov2) + 1) < 1 / threshold and cov1 + cov2 > diff_cov / 2) or cov2 - cov1 > diff_cov:
                s2.append((idx, cov1, cov2))
            else:
                s0.append((idx, cov1, cov2))

        # on side data is smaller , then we use reverse data to build training data
        if len(s1) < 100 / 2 and len(s2) > 2 * 100:
            # should append data until same length
            # s1 = map(lambda x: (x[0], x[2], x[1]), s2)
            append_data = sample(s2, k=200-len(s1))
            append_data = map(lambda x: (x[0], x[2], x[1]), append_data)
            s1.extend(append_data)
        if len(s2) < 100 / 2 and len(s1) > 2 * 100:
            # s2 = map(lambda x: (x[0], x[2], x[1]), s1)
            append_data = sample(s1, k=200-len(s1))
            append_data = map(lambda x: (x[0], x[2], x[1]), append_data)
            s2.extend(append_data)

        # not enough data for state training data , then we change parameter
        if not test and (len(s1) < 100 or len(s2) < 100):
            diff_cov = max(diff_cov - 15, 1)
            threshold = max(threshold - 0.1, 1.1)
        else:
            done = True

    if configuration.DEBUG:
        print('Final training set parameters: threshold: %s, diff_cov: %s' % (threshold, diff_cov), file=sys.stderr)

    all_data = []
    for el in [s0, s1, s2]:
        el = np.asarray(el)
        if not test:
            # here the condition we can't judge if it's good
            """
            el = el[np.where(
                np.logical_and(el[:, 1] < np.percentile(el[:, 1], 98.5) * (el[:, 1] > np.percentile(el[:, 1], 1.5)),
                               el[:, 2] < np.percentile(el[:, 2], 98.5) * (el[:, 2] > np.percentile(el[:, 2], 1.5))))]
            """
            el = el[np.where(np.logical_and(el[:, 1]<np.percentile(el[:, 1], 98), el[:, 2]<np.percentile(el[:, 2], 98)))]

        all_data.append(el)

    training_size = min(min([len(all_data[i]) for i in range(3)]),y)
    # get good enough training data sets
    for i in range(len(all_data)):
        if isvalid_training_data(all_data[i], 0.5):
            ss = sample(all_data[i], training_size)
            limits = 10
            while not isvalid_training_data(ss, 0.3) and limits:
                ss = sample(all_data[i], training_size)
                limits -= 1
            all_data[i]=ss
        else:
            ## need to go back to sample_methods again
            ## if we could also judge it before we do it, it would be great
            pass

    # first column is idx information
    training_indices = set() #set([all_data[i][:, 0] for i in range(len(all_data))])
    for i in range(len(all_data)):
        for data in all_data[i]:
            idx = range(max(0, int(data[0]) - ex), int(data[0]) + ex + 1)
            training_indices |= set(idx)
    training_indices = list(training_indices)

    training_cov = exp_data.get_sm_covs(training_indices)

    training_data = transform_data_for_HMM(training_cov)
    s0_v = map(lambda x: (x[1], x[2]), all_data[0])
    s1_v = map(lambda x: (x[1], x[2]), all_data[1])
    s2_v = map(lambda x: (x[1], x[2]), all_data[2])

    if configuration.DEBUG:
        output_training_set(name, s0_v, s1_v, s2_v)

    return training_data, s0_v, s1_v, s2_v


def write_test_samples( name, l):
    f = open(name, 'w')

    for el1, el2 in l:
        print(el1, el2, sep='\t', file=f)
    f.close()


def output_training_set(name, s0_v, s1_v, s2_v):
    """Output debug info for training_set computation."""
    write_test_samples(name + '-s0', s0_v)
    write_test_samples(name + '-s1', s1_v)
    write_test_samples(name + '-s2', s2_v)


def write_wig(sm_cov, filename):
    """Output coverage in filename of wig format."""
    f = open(filename, 'w')
    i = 0
    for region in sm_cov.genomicRegions:
        print('variableStep chrom=' + str(region.chrom) + ' span=' + str(sm_cov.stepsize), file=f)
        c = sm_cov.coverage[i]
        i += 1
        for j in c.indices:
            print(j * sm_cov.stepsize + ((sm_cov.binsize-sm_cov.stepsize)/2), c[:,j], file=f)
    f.close()