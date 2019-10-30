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

import sys
import gc
from random import sample
import numpy as np
from .normalize import get_normalization_factor
from .DualCoverageSet import DualCoverageSet
from .norm_genelevel import norm_gene_level
from ..CoverageSet import CoverageSet, get_gc_context
from functools import reduce

EPSILON = 1 ** -320
ROUND_PRECISION = 3
DEBUG = None
VERBOSE = None


class MultiCoverageSet(DualCoverageSet):
    def _help_init(self, path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, dim, regions,
                   norm_regionset, strand_cov):
        """Return self.covs and self.inputs as CoverageSet"""
        self.exts = exts
        self.covs = [CoverageSet('file' + str(i), regions) for i in range(dim)]
        for i, c in enumerate(self.covs):
            c.coverage_from_bam(bam_file=path_bamfiles[i], extension_size=exts[i], rmdup=rmdup, binsize=binsize,
                                stepsize=stepsize, get_strand_info=strand_cov)
        self.covs_avg = [CoverageSet('cov_avg' + str(i), regions) for i in range(2)]
        if path_inputs:
            self.inputs = [CoverageSet('input' + str(i), regions) for i in range(len(path_inputs))]
            for i, c in enumerate(self.inputs):
                c.coverage_from_bam(bam_file=path_inputs[i], extension_size=exts_inputs[i], rmdup=rmdup,
                                    binsize=binsize, stepsize=stepsize, get_strand_info=strand_cov)
            self.input_avg = [CoverageSet('input_avg' + str(i), regions) for i in range(2)]
        else:
            self.inputs = []

        if norm_regionset:
            self.norm_regions = [CoverageSet('norm_region' + str(i), norm_regionset) for i in range(dim)]
            for i, c in enumerate(self.norm_regions):
                c.coverage_from_bam(bam_file=path_bamfiles[i], extension_size=exts[i], rmdup=rmdup, binsize=binsize,
                                    stepsize=stepsize, get_strand_info=strand_cov)
            self.input_avg = [CoverageSet('input_avg' + str(i), regions) for i in range(2)]
        else:
            self.norm_regions = None

    def _get_covs(self, DCS, i):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i"""
        cov1 = int(np.mean(DCS.overall_coverage[0][:, DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:, DCS.indices_of_interest[i]]))

        return cov1, cov2

    def _compute_gc_content(self, no_gc_content, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes,
                            chrom_sizes_dict):
        """Compute GC-content"""
        if not no_gc_content and path_inputs and self.gc_content_cov is None:
            print("Compute GC-content", file=sys.stderr)
            for i, cov in enumerate(self.covs):
                inputfile = self.inputs[i]  # 1 to 1 mapping between input and cov
                rep = i if i < self.dim_1 else i - self.dim_1
                sig = 1 if i < self.dim_1 else 2
                self.gc_content_cov, self.avg_gc_content, self.gc_hist = get_gc_context(stepsize, binsize, genome_path,
                                                                                        inputfile.coverage,
                                                                                        chrom_sizes_dict)
                self._norm_gc_content(cov.coverage, self.gc_content_cov, self.avg_gc_content)
                self._norm_gc_content(inputfile.coverage, self.gc_content_cov, self.avg_gc_content)

                # if VERBOSE:
                #    self.print_gc_hist(name + '-s%s-rep%s-' %(sig, rep), gc_hist)
                #    cov.write_bigwig(name + '-s%s-rep%s-gc.bw' %(sig, rep), chrom_sizes)

    def _output_input_bw(self, name, chrom_sizes, save_wig):
        for i in range(len(self.covs)):
            rep = i if i < self.dim_1 else i - self.dim_1
            sig = 1 if i < self.dim_1 else 2
            if self.inputs:
                self.inputs[i].write_bigwig(name + '-' + str(self.counter) + '-input-s%s-rep%s.bw' % (sig, rep),
                                            chrom_sizes, save_wig=save_wig, end=self.end)

    def _output_bw(self, name, chrom_sizes, save_wig, save_input):
        """Output bigwig files"""
        for i in range(len(self.covs)):
            rep = i if i < self.dim_1 else i - self.dim_1
            sig = 1 if i < self.dim_1 else 2

            self.covs[i].write_bigwig(name + '-' + str(self.counter) + '-s%s-rep%s.bw' % (sig, rep), chrom_sizes,
                                      save_wig=save_wig, end=self.end)

        # ra = [self.covs_avg, self.input_avg] if self.inputs else [self.covs_avg]
        # for k, d in enumerate(ra):
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
        if self.inputs:
            for i in range(len(self.covs)):
                self.inputs[i] = None  # last time that we need this information, delete it
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

        tmp = [[], []]
        tmp2 = [[[], []], [[], []]]

        for k in range(2):
            it = list(range(self.dim_1)) if k == 0 else list(range(self.dim_1, self.dim_1 + self.dim_2))
            for i in it:
                if cov_strand:
                    tmp_el = reduce(lambda x, y: np.concatenate((x, y)), self._help_get_data(i, 'cov'))
                    tmp[k].append(tmp_el)

                    tmp_el = [(x[0], x[1]) for x in
                              reduce(lambda x, y: np.concatenate((x, y)), self._help_get_data(i, 'strand'))]
                    tmp2[k][0].append([x[0] for x in tmp_el])
                    tmp2[k][1].append([x[1] for x in tmp_el])
                else:
                    tmp_el = reduce(lambda x, y: np.concatenate((x, y)), self._help_get_data(i, 'normregion'))
                    tmp[k].append(tmp_el)

        if cov_strand:
            # 1. or 2. signal -> pos/neg strand -> matrix with rep x bins
            overall_coverage_strand = [[np.matrix(tmp2[0][0]), np.matrix(tmp2[0][1])],
                                       [np.matrix(tmp2[1][0]), np.matrix(tmp2[0][1])]]
            # list of matrices: #replicates (row) x #bins (columns)
            overall_coverage = [np.matrix(tmp[0]), np.matrix(tmp[1])]

            return overall_coverage, overall_coverage_strand
        else:
            return [np.matrix(tmp[0]), np.matrix(tmp[1])]

    def count_positive_signal(self):
        return np.sum([self.covs[i].coverage for i in range(self.dim_1 + self.dim_2)])

    def __init__(self, name, dims, regions, genome_path, binsize, stepsize, chrom_sizes, norm_regionset,
                 verbose, debug, no_gc_content, rmdup, path_bamfiles, exts, path_inputs, exts_inputs,
                 factors_inputs, chrom_sizes_dict, scaling_factors_ip, save_wig, strand_cov, housekeeping_genes,
                 tracker, end, counter, gc_content_cov=None, avg_gc_content=None, gc_hist=None, output_bw=True,
                 folder_report=None, report=None, save_input=False, m_threshold=80, a_threshold=95):
        """Compute CoverageSets, GC-content and normalize input-DNA and IP-channel"""
        self.genomicRegions = regions
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.chrom_sizes_dict = chrom_sizes_dict
        self.dim_1, self.dim_2 = dims
        self.exts = exts
        self.exts_inputs = exts_inputs
        self.gc_content_cov = gc_content_cov
        self.avg_gc_content = avg_gc_content
        self.gc_hist = gc_hist
        self.scaling_factors_ip = scaling_factors_ip
        self.factors_inputs = factors_inputs
        self.end = end
        self.counter = counter
        self.no_data = False
        self.FOLDER_REPORT = folder_report
        global DEBUG, VERBOSE
        DEBUG = debug
        VERBOSE = verbose

        # make data nice
        self._help_init(path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, sum(dims), regions,
                        norm_regionset, strand_cov=strand_cov)
        if self.count_positive_signal() < 1:
            self.no_data = True
            return
        self._compute_gc_content(no_gc_content, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes,
                                 chrom_sizes_dict)
        self._normalization_by_input(path_bamfiles, path_inputs, name, factors_inputs, save_input)
        if save_input:
            self._output_input_bw(name, chrom_sizes, save_wig)

        self.overall_coverage, self.overall_coverage_strand = self._help_init_overall_coverage(cov_strand=True)

        self._normalization_by_signal(name, scaling_factors_ip, path_bamfiles, housekeeping_genes, tracker,
                                      norm_regionset, report,
                                      m_threshold, a_threshold)

        if output_bw:
            self._output_bw(name, chrom_sizes, save_wig, save_input)

        self.scores = np.zeros(len(self.overall_coverage[0]))
        self.indices_of_interest = []

    def get_max_colsum(self):
        """Sum over all columns and add maximum"""
        return self.overall_coverage[0].sum(axis=0).max() + self.overall_coverage[1].sum(axis=0).max()

    def output_overall_coverage(self, path):
        for j in range(2):
            f = open(path + str(j), 'w')
            for i in range(self.overall_coverage[j].shape[1]):
                print(self.overall_coverage[j][:, i].T, file=f)

    def _normalization_by_input(self, path_bamfiles, path_inputs, name, factors_inputs, save_input):
        """Normalize input-DNA. Use predefined factors or follow Diaz et al, 2012"""

        if VERBOSE:
            print("Normalize input-DNA", file=sys.stderr)

        if factors_inputs:
            if VERBOSE:
                print("Use with predefined factors", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                self.inputs[i].scale(factors_inputs[i])
                self.covs[i].subtract(self.inputs[i])
        elif path_inputs:
            factors_inputs = []
            print("Compute factors", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                rep = i if i < self.dim_1 else i - self.dim_1
                sig = 0 if i < self.dim_1 else 1
                j = 0 if i < self.dim_1 else 1
                _, n = get_normalization_factor(path_bamfiles[i], path_inputs[i], step_width=1000, zero_counts=0, \
                                                filename=name + '-norm' + str(i), debug=DEBUG,
                                                chrom_sizes_dict=self.chrom_sizes_dict, two_sample=False, stop=True)
                if n is not None:
                    print("Normalize input of Signal %s, Rep %s with factor %s" \
                          % (sig, rep, round(n, ROUND_PRECISION)), file=sys.stderr)
                    self.inputs[i].scale(n)

                    self.covs[i].subtract(self.inputs[i])
                    factors_inputs.append(n)

        self.factors_inputs = factors_inputs

    def _trim4TMM(self, m_values, a_values, m_threshold=80, a_threshold=95):
        """q=20 or q=5"""
        assert len(m_values) == len(a_values)

        mask = np.asarray([not x for x in np.isinf(m_values) + np.isinf(a_values)])

        m_values = m_values[mask]
        a_values = a_values[mask]

        perc_m_l = np.percentile(m_values, 100 - m_threshold)
        perc_m_h = np.percentile(m_values, m_threshold)
        perc_a_l = np.percentile(a_values, 100 - a_threshold)
        perc_a_h = np.percentile(a_values, a_threshold)

        try:
            res = [x for x in [x for x in zip(list(m_values.squeeze()), list(a_values.squeeze())) if
                               not (x[1] > perc_a_h or x[1] < perc_a_l)] if not (x[0] > perc_m_h or x[0] < perc_m_l)]
        except:
            print('something wrong %s %s' % (len(m_values), len(a_values)), file=sys.stderr)
            return np.asarray(m_values), np.asarray(a_values)

        if res:
            return np.asarray([x[0] for x in res]), np.asarray([x[1] for x in res])
        else:
            print('TMM normalization: nothing trimmed...', file=sys.stderr)
            return np.asarray(m_values), np.asarray(a_values)

    def _norm_TMM(self, overall_coverage, m_threshold, a_threshold):
        """Normalize with TMM approach, based on PePr"""
        scaling_factors_ip = []
        for j, cond_max in enumerate([self.dim_1, self.dim_2]):
            for i in range(cond_max):  # normalize all replicates
                ref = np.asarray(np.sum(overall_coverage[0], axis=0) + np.sum(overall_coverage[1], axis=0),
                                 dtype='float') / (self.dim_1 + self.dim_2)

                mask_ref = ref > 0
                ref = ref[mask_ref]
                data_rep = np.asarray(overall_coverage[j][i, :])[mask_ref]
                tmp = list(zip(data_rep, ref, data_rep + ref))
                tmp.sort(key=lambda x: x[2], reverse=True)
                tmp = tmp[:min(len(tmp), 10000)]

                data_rep = np.asarray([x[0] for x in tmp])
                ref = np.asarray([x[1] for x in tmp])
                assert len(data_rep) == len(ref)
                m = data_rep > 0
                data_rep = data_rep[m]
                ref = ref[m]

                m_values = np.log(ref / data_rep)
                a_values = 0.5 * np.log(data_rep * ref)
                try:
                    m_values, a_values = self._trim4TMM(m_values, a_values, m_threshold, a_threshold)
                    f = 2 ** (np.sum(m_values * a_values) / np.sum(a_values))
                    scaling_factors_ip.append(f)
                except:
                    print('TMM normalization not successfully performed, do not normalize data', file=sys.stderr)
                    scaling_factors_ip.append(1)

        return scaling_factors_ip

    def _normalization_by_signal(self, name, scaling_factors_ip, bamfiles, housekeeping_genes, tracker, norm_regionset,
                                 report,
                                 m_threshold, a_threshold):
        """Normalize signal"""

        if VERBOSE:
            print('Normalize ChIP-seq profiles', file=sys.stderr)

        if not scaling_factors_ip and housekeeping_genes:
            print('Use housekeeping gene approach', file=sys.stderr)
            scaling_factors_ip, _ = norm_gene_level(bamfiles, housekeeping_genes, name, verbose=True,
                                                    folder=self.FOLDER_REPORT, report=report)
        elif not scaling_factors_ip:
            if norm_regionset:
                print('Use TMM approach based on peaks', file=sys.stderr)
                norm_regionset_coverage = self._help_init_overall_coverage(
                    cov_strand=False)  # TMM approach based on peaks
                scaling_factors_ip = self._norm_TMM(norm_regionset_coverage, m_threshold, a_threshold)
            else:
                print('Use global TMM approach ', file=sys.stderr)
                scaling_factors_ip = self._norm_TMM(self.overall_coverage, m_threshold, a_threshold)  # TMM approach

        for i in range(len(scaling_factors_ip)):
            self.covs[i].scale(scaling_factors_ip[i])

        if scaling_factors_ip:
            for j, cond in enumerate([self.dim_1, self.dim_2]):
                for i in range(cond):  # normalize all replicates
                    k = i if j == 0 else i + self.dim_1
                    self.overall_coverage[j][i, :] *= scaling_factors_ip[k]
                    if DEBUG:
                        print('Use scaling factor %s' % round(scaling_factors_ip[k], ROUND_PRECISION), file=sys.stderr)

        self.scaling_factors_ip = scaling_factors_ip

    def _index2coordinates(self, index):
        """Translate index within coverage array to genomic coordinates."""
        iter = self.genomicRegions.__iter__()
        r = next(iter)
        sum = r.final
        last = 0
        i = 0
        while sum <= index * self.stepsize:
            last += len(self.covs[0].coverage[i])
            try:
                r = next(iter)
            except StopIteration:
                sum += r.final
                i += 1
                break
            sum += r.final
            i += 1

        return r.chrom, (index - last) * self.stepsize, \
               min((index - last) * self.stepsize + self.stepsize, r.final)

    def __len__(self):
        """Return number of observations."""
        return len(self.indices_of_interest)

    def get_observation(self, mask=np.array([])):
        """Return indices of observations. Do not consider indices contained in <mask> array"""
        mask = np.asarray(mask)
        if not mask.size:
            mask = np.array([True] * self._get_bin_number())
        return np.asarray(
            np.concatenate((self.overall_coverage[0][:, mask].T, self.overall_coverage[1][:, mask].T), axis=1))

    def _compute_score(self):
        """Compute score for each observation (based on Xu et al.)"""
        self.scores = sum([np.squeeze(np.asarray(np.mean(self.overall_coverage[i], axis=0))) / float(
            np.mean(self.overall_coverage[i])) for i in range(2)])

    def _get_bin_number(self):
        """Return number of bins"""
        return self.overall_coverage[0].shape[1]

    def compute_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows: 
        - score must be > 0, i.e. everthing
        - overall coverage in library 1 and 2 must be > 3"""

        try:
            self._compute_score()
            self.indices_of_interest = np.where(self.scores > 0)[0]  # 2/(m*n)
            tmp = np.where(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) + np.squeeze(
                np.asarray(np.mean(self.overall_coverage[1], axis=0))) > 10)[0]
            tmp2 = np.intersect1d(self.indices_of_interest, tmp)
            self.indices_of_interest = tmp2
        except:
            self.indices_of_interest = None
        # print(len(self.indices_of_interest), file=sys.stderr)
        # tmp = set()
        # for i in self.indices_of_interest:
        #    for j in range(max(0, i-l), i+l+1):
        #        tmp.add(j)
        # tmp = list(tmp)
        # tmp.sort()
        # self.indices_of_interest = np.array(tmp)

    def write_test_samples(self, name, l):
        f = open(name, 'w')

        for el1, el2 in l:
            print(el1, el2, sep='\t', file=f)
        f.close()

    def output_training_set(self, name, training_set, s0_v, s1_v, s2_v):
        """Output debug info for training_set computation."""
        f = open(name + '-trainingset.bed', 'w')
        for l in training_set:
            chrom, s, e = self._index2coordinates(l)
            print(chrom, s, e, sep='\t', file=f)
        f.close()

        self.write_test_samples(name + '-s0', s0_v)
        self.write_test_samples(name + '-s1', s1_v)
        self.write_test_samples(name + '-s2', s2_v)

    def get_training_set(self, test, exp_data, name, foldchange, min_t, y=5000, ex=2):
        """Return HMM's training set (max <y> positions). Enlarge each contained bin by <ex>."""
        threshold = foldchange
        diff_cov = int(np.percentile(np.abs(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) - \
                                            np.squeeze(np.asarray(np.mean(self.overall_coverage[1], axis=0)))), min_t))

        if test:
            diff_cov, threshold = 2, 1.5

        if DEBUG:
            print('Training set parameters: threshold: %s, diff_cov: %s' % (threshold, diff_cov), file=sys.stderr)

        s0, s1, s2 = [], [], []

        # compute training set parameters, re-compute training set if criteria do not hold
        rep = True
        while rep:
            for i in sample(list(range(len(self.indices_of_interest))), min(y, len(self.indices_of_interest))):
                cov1, cov2 = self._get_covs(exp_data, i)

                # apply criteria for initial peak calling
                if (cov1 / max(float(cov2), 1) > threshold and cov1 + cov2 > diff_cov / 2) or cov1 - cov2 > diff_cov:
                    s1.append((self.indices_of_interest[i], cov1, cov2))
                elif (cov1 / max(float(cov2),
                                 1) < 1 / threshold and cov1 + cov2 > diff_cov / 2) or cov2 - cov1 > diff_cov:
                    s2.append((self.indices_of_interest[i], cov1, cov2))
                else:
                    s0.append((self.indices_of_interest[i], cov1, cov2))

                if len(s0) > y and len(s1) > y and len(s2) > y:
                    break

            if diff_cov == 1 and threshold == 1.1:
                print("No differential peaks detected", file=sys.stderr)
                sys.exit()

            if len(s1) < 100 / 2 and len(s2) > 2 * 100:
                s1 = [(x[0], x[2], x[1]) for x in s2]
            if len(s2) < 100 / 2 and len(s1) > 2 * 100:
                s2 = [(x[0], x[2], x[1]) for x in s1]

            if len(s1) < 100 or len(s2) < 100:
                diff_cov -= 15
                threshold -= 0.1
                diff_cov = max(diff_cov, 1)
                threshold = max(threshold, 1.1)
            else:
                rep = False

        if DEBUG:
            print('Final training set parameters: threshold: %s, diff_cov: %s' % (threshold, diff_cov), file=sys.stderr)

        # optimize training set, extend each bin
        tmp = []
        for i, el in enumerate([s0, s1, s2]):
            el = np.asarray(el)
            if not test:
                el = el[el[:, 1] < np.percentile(el[:, 1], 90)]
                el = el[el[:, 2] < np.percentile(el[:, 2], 90)]
            tmp.append(el)

        s0 = tmp[0]
        s1 = tmp[1]
        s2 = tmp[2]

        l = np.min([len(s1), len(s2), len(s0), y])

        s0 = sample(list(s0), l)
        s1 = sample(list(s1), l)
        s2 = sample(list(s2), l)

        s0_v = [(x[1], x[2]) for x in s0]
        s1_v = [(x[1], x[2]) for x in s1]
        s2_v = [(x[1], x[2]) for x in s2]

        extension_set = set()
        for i, _, _ in s0 + s1 + s2:
            for j in range(max(0, i - ex), i + ex + 1):  # extend bins
                extension_set.add(j)

        tmp = s0 + s1 + s2
        training_set = [x[0] for x in tmp] + list(extension_set)

        training_set = list(training_set)
        training_set.sort()

        if DEBUG:
            self.output_training_set(name, training_set, s0_v, s1_v, s2_v)

        return training_set, s0_v, s1_v, s2_v
