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

@author: Manuel Allhoff
"""


import sys
import numpy as np
from os import path
from random import sample
from rgt.CoverageSet import CoverageSet
from rgt.CoverageSet import get_gc_context
from .normalize import get_normalization_factor
from functools import reduce

EPSILON=1e-320

class DualCoverageSet():
    def _get_BAM_names(self, input, ip):
        # get names of BAM files for bw files
        name_bam = path.splitext(path.basename(ip))[0]
        if input is not None:
            name_input = path.splitext(path.basename(input))[0]
        else:
            name_input = None

        return name_bam, name_input

    def __init__(self, name, region, genome_path, binsize, stepsize, rmdup, file_1, ext_1, file_2, ext_2, \
                 input_1, ext_input_1, input_factor_1, input_2, ext_input_2, input_factor_2, chrom_sizes, verbose,
                 norm_strategy, no_gc_content, deadzones, \
                 factor_input_1, factor_input_2, chrom_sizes_dict, debug, tracker):
        self.genomicRegions = region
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.cov1 = CoverageSet('first file', region)
        self.cov2 = CoverageSet('second file', region)

        print("Loading reads...", file=sys.stderr)
        self.cov1.coverage_from_bam(bam_file=file_1, extension_size=ext_1, rmdup=rmdup, binsize=binsize,
                                    stepsize=stepsize, mask_file=deadzones)
        self.cov2.coverage_from_bam(bam_file=file_2, extension_size=ext_2, rmdup=rmdup, binsize=binsize,
                                    stepsize=stepsize, mask_file=deadzones)

        map_input = {1: {'input': input_1, 'input_factor': input_factor_1, 'ext': ext_input_1, 'cov-ip': self.cov1,
                         'ip': file_1},
                     2: {'input': input_2, 'input_factor': input_factor_2, 'ext': ext_input_2, 'cov-ip': self.cov2,
                         'ip': file_2}}

        if not no_gc_content and input_1 is not None and input_2 is not None:
            print("Computing GC content", file=sys.stderr)
        else:
            print("Do not compute GC content", file=sys.stderr)

        for i in [1, 2]:
            input = map_input[i]
            name_bam, name_input = self._get_BAM_names(input['input'], input['ip'])

            if debug:  # 0: output raw IP
                input['cov-ip'].write_bigwig(name + '-debug-0-' + name_bam + '.bw', chrom_sizes)

            if input['input'] is not None:
                input['cov-input'] = CoverageSet('%s file' % input['input'], region)
                input['cov-input'].coverage_from_bam(bam_file=input['input'], extension_size=input['ext'], rmdup=rmdup,
                                                     binsize=binsize, stepsize=stepsize)
                map_input[i]['cov-input'] = input['cov-input']

            if not no_gc_content and input['input'] is not None:
                gc_content_cov, avg_gc_content, gc_hist = get_gc_context(stepsize, binsize, genome_path,
                                                                         input['cov-input'].coverage, chrom_sizes_dict)

                self._norm_gc_content(input['cov-ip'].coverage, gc_content_cov, avg_gc_content)
                self._norm_gc_content(input['cov-input'].coverage, gc_content_cov, avg_gc_content)

                if debug:  # 1: output after GC
                    self.print_gc_hist(name + '-' + name_input, gc_hist)  # print hist data
                    input['cov-input'].write_bigwig(name + '-debug-1-' + name_input + '.bw', chrom_sizes)
                    input['cov-ip'].write_bigwig(name + '-debug-1-' + name_bam + '.bw', chrom_sizes)

        norm_done = False
        print("Normalizing signals", file=sys.stderr)
        for i in [1, 2]:
            input = map_input[i]
            name_bam, name_input = self._get_BAM_names(input['input'], input['ip'])

            # TODO: uncomment here!
            norm_done = self.normalization(map_input, i, norm_strategy, norm_done, name, debug, factor_input_1,
                                           factor_input_2, chrom_sizes_dict, tracker)

            if input['input'] is not None:
                input['cov-input'].write_bigwig(name + '-' + name_input + '-normalized.bw', chrom_sizes)
            input['cov-ip'].write_bigwig(name + '-' + name_bam + '-normalized.bw', chrom_sizes)

        # make one array for the coverage
        self.first_overall_coverage = reduce(lambda x, y: np.concatenate((x, y)),
                                             [self.cov1.coverage[i] for i in range(len(self.cov1.genomicRegions))])
        self.second_overall_coverage = reduce(lambda x, y: np.concatenate((x, y)),
                                              [self.cov2.coverage[i] for i in range(len(self.cov2.genomicRegions))])
        assert (len(self.first_overall_coverage) == len(self.second_overall_coverage))

        self.scores = np.zeros(len(self.first_overall_coverage))
        self.indices_of_interest = []

    def normalization(self, map_input, i, norm_strategy, norm_done, name, debug, factor_input_1, factor_input_2,
                      chrom_sizes_dict, tracker):
        input = map_input[i]

        # compute normalization factor
        # pre-defined values
        if input['input_factor'] is not None and i != 1:
            print("Normalize by Diaz and pre-defined values...", input['input_factor'], file=sys.stderr)
            print("Normalize file 1 with input normalization factor %s" % (map_input[1]['input_factor']),
                  file=sys.stderr)
            print("Normalize file 2 with input normalization factor %s" % (map_input[2]['input_factor']),
                  file=sys.stderr)
            tracker.write(text=str(map_input[1]['input_factor']) + ',' + str(map_input[2]['input_factor']),
                          header="Predefined Normalization factor of Input")

            map_input[1]['cov-input'].scale(map_input[1]['input_factor'])
            map_input[2]['cov-input'].scale(map_input[2]['input_factor'])
            map_input[1]['cov-ip'].subtract(map_input[1]['cov-input'])
            map_input[2]['cov-ip'].subtract(map_input[2]['cov-input'])

        # naive norm.
        if not norm_done and norm_strategy == 1:
            if factor_input_1 is None or factor_input_2 is None:
                s1 = sum([sum(map_input[1]['cov-ip'].coverage[i]) for i in
                          range(len(map_input[1]['cov-ip'].genomicRegions))])
                s2 = sum([sum(map_input[2]['cov-ip'].coverage[i]) for i in
                          range(len(map_input[2]['cov-ip'].genomicRegions))])
                if s1 > s2:
                    map_input[2]['cov-ip'].scale(s1 / float(s2))
                    print("Normalize file 2 by signal with estimated factor %s " % (round(s1 / float(s2), 3)),
                          file=sys.stderr)
                    tracker.write(text=str(round(s1 / float(s2), 3)), header="Normalization factor of signal 2")
                elif s2 >= s1:
                    print("Normalize file 1 by signal with estimated factor %s " % (round(s2 / float(s1), 3)),
                          file=sys.stderr)
                    tracker.write(text=str(round(s2 / float(s1), 3)), header="Normalization factor of signal 1")
                    map_input[1]['cov-ip'].scale(s2 / float(s1))

                norm_done = True
            else:
                map_input[1]['cov-ip'].scale(factor_input_1)
                print("Normalize file 1 by signal with given factor %s " % round(factor_input_1, 3), file=sys.stderr)
                tracker.write(text=str(round(factor_input_1, 3)), header="Predefined Normalization factor of signal 1")

                map_input[2]['cov-ip'].scale(factor_input_2)
                print("Normalize file 2 by signal with given factor %s " % round(factor_input_2, 3), file=sys.stderr)
                tracker.write(text=str(round(factor_input_2, 3)), header="Predefined Normalization factor of signal 2")
                norm_done = True

        # diaz and naive
        if i != 1 and norm_strategy == 5:
            # apply diaz
            _, map_input[1]['input_factor'] = get_normalization_factor(map_input[1]['ip'], map_input[1]['input'],
                                                                       step_width=1000, zero_counts=0, \
                                                                       filename=name + '-norm' + str(i), debug=debug,
                                                                       chrom_sizes_dict=chrom_sizes_dict,
                                                                       two_sample=False)
            _, map_input[2]['input_factor'] = get_normalization_factor(map_input[2]['ip'], map_input[2]['input'],
                                                                       step_width=1000, zero_counts=0, \
                                                                       filename=name + '-norm' + str(i), debug=debug,
                                                                       chrom_sizes_dict=chrom_sizes_dict,
                                                                       two_sample=False)

            print("Normalize input with factor %s and %s" % (
            round(map_input[1]['input_factor'], 3), round(map_input[2]['input_factor'], 3)), file=sys.stderr)
            tracker.write(
                text=str(round(map_input[1]['input_factor'], 3)) + ',' + str(round(map_input[2]['input_factor'], 3)),
                header="Input Normalization factors")

            map_input[1]['cov-input'].scale(map_input[1]['input_factor'])
            map_input[2]['cov-input'].scale(map_input[2]['input_factor'])

            map_input[1]['cov-ip'].subtract(map_input[1]['cov-input'])
            map_input[2]['cov-ip'].subtract(map_input[2]['cov-input'])

            if factor_input_1 is None or factor_input_2 is None:
                # apply naive method
                s1 = sum([sum(map_input[1]['cov-ip'].coverage[i]) for i in
                          range(len(map_input[1]['cov-ip'].genomicRegions))])
                s2 = sum([sum(map_input[2]['cov-ip'].coverage[i]) for i in
                          range(len(map_input[2]['cov-ip'].genomicRegions))])

                if s1 > s2:
                    map_input[2]['cov-ip'].scale(s1 / float(s2))
                    print("Normalize file 2 by signal with estimated factor %s " % (round(s1 / float(s2), 3)),
                          file=sys.stderr)
                    tracker.write(text=str(round(s1 / float(s2), 3)), header="Normalization factor of signal 2")
                elif s2 >= s1:
                    print("Normalize file 1 by signal with estimated factor %s " % (round(s2 / float(s1), 3)),
                          file=sys.stderr)
                    map_input[1]['cov-ip'].scale(s2 / float(s1))
                    tracker.write(text=str(round(s2 / float(s1), 3)), header="Normalization factor of signal 1")
            else:
                map_input[1]['cov-ip'].scale(factor_input_1)
                print("Normalize file 1 by signal with given factor %s " % round(factor_input_1, 3), file=sys.stderr)
                tracker.write(text=str(round(factor_input_1, 3)), header="Normalization factor of signal 1")
                map_input[2]['cov-ip'].scale(factor_input_2)
                print("Normalize file 2 by signal with given factor %s " % round(factor_input_2, 3), file=sys.stderr)
                tracker.write(text=str(round(factor_input_2, 3)), header="Normalization factor of signal 2")
        return norm_done

    def print_gc_hist(self, name, gc_hist):
        f = open(name + '-gc-content.data', 'w')
        for i in range(len(gc_hist)):
            print(i, gc_hist[i], file=f)
        f.close()

    def _norm_gc_content(self, cov, gc_cov, gc_avg):
        for i in range(len(cov)):
            assert len(cov[i]) == len(gc_cov[i])
            #            cov[i] = gc_cov[i]
            cov[i] = np.array(cov[i])
            gc_cov[i] = np.array(gc_cov[i])
            gc_cov[i][gc_cov[i] < EPSILON] = gc_avg  # sometimes zeros occur, do not consider
            cov[i] = cov[i] * gc_avg / gc_cov[i]
            cov[i] = cov[i].clip(0, max(max(cov[i]), 0))  # neg. values to 0
            cov[i] = cov[i].astype(int)

    def _index2coordinates(self, index):
        """Translate index within coverage array to genomic coordinates."""
        iter = self.genomicRegions.__iter__()
        r = next(iter)
        sum = r.final
        last = 0
        i = 0
        while sum <= index * self.stepsize:
            last += len(self.cov1.coverage[i])
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
        if not mask.size:
            mask = np.array([True] * len(self.first_overall_coverage))
        return np.array([self.first_overall_coverage[mask], self.second_overall_coverage[mask]]).T

    def _compute_score(self):
        """Compute score for each observation (based on Xu et al.)"""
        self.scores = self.first_overall_coverage / float(sum(self.first_overall_coverage)) + \
                      self.second_overall_coverage / float(sum(self.second_overall_coverage))

    def compute_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows: 
        - score must be > 2/(m*n) (m=#obs, n=0.9 (default) )
        - overall coverage in library 1 and 2 must be > 3
        - extend resulting sites by l steps in both directions. """
        m = len(self.first_overall_coverage)
        n = 0.9
        self._compute_score()
        #        print('before filter step:', len(self.scores), file=sys.stderr)
        # print(self.first_overall_coverage, self.second_overall_coverage, file=sys.stderr)
        self.indices_of_interest = np.where(self.scores > 2 / (m * n))[0]
        #        print('after first filter step: ', len(self.indices_of_interest), file=sys.stderr)
        tmp = np.where(self.first_overall_coverage + self.second_overall_coverage > 3)[0]
        tmp2 = np.intersect1d(self.indices_of_interest, tmp)
        #        print('length of intersection set: ', len(tmp), file=sys.stderr)
        self.indices_of_interest = tmp2
        #        print('after second filter step: ', len(self.indices_of_interest), file=sys.stderr)
        # extend regions by l steps
        tmp = set()
        for i in self.indices_of_interest:
            for j in range(max(0, i - l), i + l + 1):
                tmp.add(j)
        tmp = list(tmp)
        tmp.sort()
        self.indices_of_interest = np.array(tmp)

    def get_initial_dist(self, filename):
        """Write BED file with initial state distribution"""
        states = []
        threshold = 2.0
        for i in self.indices_of_interest:
            c1 = self.first_overall_coverage[i]
            c2 = self.second_overall_coverage[i]

            if c1 + c2 <= 3:
                state = 0
            elif c1 / max(float(c2), 1) > threshold or c1 - c2 > 10:
                state = 1
            elif c1 / max(float(c2), 1) < 1 / threshold or c2 - c1 > 10:
                state = 2
            else:
                state = 0

            states.append(state)

        f = open(filename, 'w')
        for j in range(len(states)):
            i = self.indices_of_interest[j]
            chrom, start, end = self._index2coordinates(i)
            s = states[j]
            print(chrom, start, end, s, self.first_overall_coverage[i], self.second_overall_coverage[i], sep='\t',
                  file=f)

        f.close()

    def write_putative_regions(self, path):
        """Write putative regions (defined by criteria mentioned in method) as BED file."""
        with open(path, 'w') as f:
            for i in self.indices_of_interest:
                chrom, start, end = self._index2coordinates(i)
                print(chrom, start, end, file=f)

    def get_training_set(self, exp_data, x, verbose, name, debug, constraint_chrom, min_fc=2.0):
        """Return linked genomic positions (at least <x> positions) to train HMM.
        Grep randomly a position within a putative region, and take then the entire region."""
        training_set = set()
        ts1 = set()
        ts2 = set()
        threshold = min_fc
        diff_cov = 10

        if constraint_chrom is not None:
            print("HMM training set based on %s" % constraint_chrom, file=sys.stderr)

        for i in range(len(self.indices_of_interest)):
            chrom, start, end = self._index2coordinates(i)
            if constraint_chrom is not None and chrom != constraint_chrom:
                continue
            cov1 = exp_data.first_overall_coverage[self.indices_of_interest[i]]
            cov2 = exp_data.second_overall_coverage[self.indices_of_interest[i]]

            if cov1 / max(float(cov2), 1) > threshold or cov1 - cov2 > diff_cov:
                ts1.add(i)
            if cov1 / max(float(cov2), 1) < 1 / threshold or cov2 - cov1 > diff_cov:
                ts2.add(i)

        l = min(min(len(ts1), len(ts2)), x)
        tmp = set(sample(ts1, l)) | set(sample(ts2, l))

        for i in tmp:
            training_set.add(self.indices_of_interest[i])
            # search up
            while i + 1 < len(self.indices_of_interest) and self.indices_of_interest[i + 1] == self.indices_of_interest[
                i] + 1:
                training_set.add(self.indices_of_interest[i + 1])
                i += 1
            # search down
            while i - 1 > 0 and self.indices_of_interest[i - 1] == self.indices_of_interest[i] - 1:
                training_set.add(self.indices_of_interest[i - 1])
                i -= 1


                #        while len(training_set) < x:
                #            i = randrange(1, len(self.indices_of_interest)-1)
                #
                #            if i in used:
                #                continue #todo: super ugly
                #            used.add(i)
                #            training_set.add(self.indices_of_interest[i])
                #            #search up
                #            while i+1 < len(self.indices_of_interest) and self.indices_of_interest[i+1] == self.indices_of_interest[i]+1:
                #                training_set.add(self.indices_of_interest[i+1])
                #                i += 1
                #            #search down
                #            while i-1 > 0 and self.indices_of_interest[i-1] == self.indices_of_interest[i]-1:
                #                training_set.add(self.indices_of_interest[i-1])
                #                i -= 1

        training_set = list(training_set)
        training_set.sort()
        if debug:
            f = open(name + '-trainingset.bed', 'w')
            for l in training_set:
                chrom, s, e = self._index2coordinates(l)
                print(chrom, s, e, sep='\t', file=f)
            f.close()

        return np.array(training_set)
