from __future__ import print_function
from rgt.CoverageSet import CoverageSet
import numpy as np
from random import sample, randrange
from time import time
from rgt.ODIN.gc_content import get_gc_context
import sys
from rgt.ODIN.normalize import get_normalization_factor


EPSILON = 1**-320

class MultiCoverageSet():
    def _help_init(self, path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, dim, regions):
        """Return self.covs and self.inputs as CoverageSet"""
        
        self.covs = [CoverageSet('file'+str(i), regions) for i in range(dim)]
        for i, c in enumerate(self.covs):
            c.coverage_from_bam(bam_file=path_bamfiles[i], read_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize)
        
        if path_inputs:
            self.inputs = [CoverageSet('input' + str(i), regions) for i in range(len(path_inputs))]
            for i, c in enumerate(self.inputs):
                c.coverage_from_bam(bam_file=path_inputs[i], read_size=exts_inputs[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize)
        else:
            self.inputs = []
    
    def _get_covs(self, DCS, i):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i"""
        cov1 = int(np.mean(DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]))
    
        return cov1, cov2
    
    def _compute_gc_content(self, no_gc_content, verbose, path_inputs, stepsize, binsize, genome_path, input, name, chrom_sizes):
        """Compute GC-content"""
        if not no_gc_content and path_inputs:
            print("Compute GC-content", file=sys.stderr)
            for i, cov in enumerate(self.covs):
                input = self.inputs[i] #1 to 1 mapping between input and cov
                
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 1 if i < self.dim_1 else 2
                  
                gc_content_cov, avg_gc_content, gc_hist = get_gc_context(stepsize, binsize, genome_path, input.coverage)
                self._norm_gc_content(cov.coverage, gc_content_cov, avg_gc_content)
                self._norm_gc_content(input.coverage, gc_content_cov, avg_gc_content)
            
                if verbose:
                    self.print_gc_hist(name + '-s%s-rep%s-' %(sig, rep), gc_hist)
                    input.write_bigwig(name + '-s%s-rep%s-input-gc.bw' %(sig, rep), chrom_sizes)
                    cov.write_bigwig(name + '-s%s-rep%s-gc.bw' %(sig, rep), chrom_sizes)
        else:
            print("Do not compute GC-content", file=sys.stderr)
    
    
    def __init__(self, name, dims, regions, genome_path, binsize, stepsize, chrom_sizes, \
                 verbose, no_gc_content, rmdup, path_bamfiles, exts, path_inputs, exts_inputs, factors_inputs):
        """Compute CoverageSets, GC-content and normalization"""
        self.genomicRegions = regions
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        
        self.dim_1, self.dim_2 = dims
        
        #make data nice
        self._help_init(path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, sum(dims), regions)
        self._compute_gc_content(no_gc_content, verbose, path_inputs, stepsize, binsize, genome_path, input, name, chrom_sizes)
        self._normalization_by_input(path_bamfiles, path_inputs, name, verbose)
        self._normalization_by_signal(name, verbose)
        
        for i in range(len(self.covs)):
            rep = i if i < self.dim_1 else i-self.dim_1
            sig = 1 if i < self.dim_1 else 2
            
            self.covs[i].write_bigwig(name + '-s%s-rep%s.bw' %(sig, rep), chrom_sizes)
            
        #make data in nice list of two matrices
        tmp = [[], []]
        for k in range(2):
            it = range(self.dim_1) if k == 0 else range(self.dim_1, self.dim_1 + self.dim_2)
            for i in it:
                tmp_el = reduce(lambda x,y: np.concatenate((x,y)), [self.covs[i].coverage[j] for j in range(len(self.covs[i].genomicRegions))])
                tmp[k].append(tmp_el)
            
        self.overall_coverage = [np.matrix(tmp[0]), np.matrix(tmp[1])] #list of matrices: #replicates (row) x #bins (columns)
        
        self.scores = np.zeros(len(self.overall_coverage[0]))
        self.indices_of_interest = []
    
    
    def _normalization_by_input(self, path_bamfiles, path_inputs, name, verbose):
        """Normalize with regard to input file"""
        if path_inputs:
            print("Normalize", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 1 if i < self.dim_1 else 2
                j = 0 if i < self.dim_1 else 1
                
                _, n = get_normalization_factor(path_bamfiles[i], path_inputs[i], step_width=1000, zero_counts=0, \
                                                                  genome='mm9', filename=name + '-norm' + str(i), verbose=False, two_sample=False)
                
                print("Factor: normalize input with input factor %s (Signal %s, Rep %s)"\
                       %(round(n, 3), sig, rep) , file=sys.stderr)
                self.inputs[i].scale(n)
                self.covs[i].subtract(self.inputs[i])
    
#     def _get_signal_sums(self):
#         s1 = sum([sum([sum(self.covs[k].coverage[i]) for i in range(len(self.covs[k].genomicRegions))]) for k in range(self.dim_1)])
#         s2 = sum([sum([sum(self.covs[k].coverage[i]) for i in range(len(self.covs[k].genomicRegions))]) for k in range(self.dim_1, self.dim_1+self.dim_2)])
#         
#         return s1, s2
    
    def _normalization_by_signal(self, name, verbose):
        """Normalize signal"""
        signals = [sum([sum(self.covs[k].coverage[i]) for i in range(len(self.covs[k].genomicRegions))]) for k in range(self.dim_1 + self.dim_2)]
        print("Normalize by signal", file=sys.stderr)
        
        means_signal = [np.mean(signals[:self.dim_1]), np.mean(signals[self.dim_1:])]
        max_index = means_signal.index(max(means_signal))
        
        if max_index == 1:
            r = range(self.dim_1)
            f = means_signal[1] / means_signal[0]
        if max_index == 0:
            r = range(self.dim_1, self.dim_1 + self.dim_2)
            f = means_signal[0] / means_signal[1]
        
        print(r, f, file=sys.stderr)
        
        for i in r:
            self.covs[i].scale(f)
        
        
    def print_gc_hist(self, name, gc_hist):
        f = open(name + 'gc-content.data', 'w')
        for i in range(len(gc_hist)):
            print(i, gc_hist[i], file=f)
        f.close()
    
    def _norm_gc_content(self, cov, gc_cov, gc_avg):
        for i in range(len(cov)):
            assert len(cov[i]) == len(gc_cov[i])
#            cov[i] = gc_cov[i]
            cov[i] = np.array(cov[i])
            gc_cov[i] = np.array(gc_cov[i])
            gc_cov[i][gc_cov[i] < EPSILON] = gc_avg #sometimes zeros occur, do not consider
            cov[i] = cov[i] * gc_avg / gc_cov[i]
            cov[i] = cov[i].clip(0, max(max(cov[i]), 0)) #neg. values to 0
            cov[i] = cov[i].astype(int)
        
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
        
        return r.chrom, (index-last) * self.stepsize + ((self.binsize-self.stepsize)/2), \
            min((index-last) * self.stepsize + ((self.binsize+self.stepsize)/2), r.final)
                              
    def __len__(self):
        """Return number of observations."""
        return len(self.indices_of_interest)
    
    def get_observation(self, mask=np.array([])):
        """Return indices of observations. Do not consider indices contained in <mask> array"""
        if not mask.size:
            mask = np.array([True]*self._get_bin_number())
            
        return np.asarray(np.concatenate((self.overall_coverage[0][:,mask].T, self.overall_coverage[1][:,mask].T), axis=1))
    
    def _compute_score(self):
        """Compute score for each observation (based on Xu et al.)"""
        self.scores = sum([np.squeeze(np.asarray(np.mean(self.overall_coverage[i], axis=0))) / float(np.mean(self.overall_coverage[i])) for i in range(2)])
    
    def _get_bin_number(self):
        return self.overall_coverage[0].shape[1]
    
    def compute_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows: 
        - score must be > 2/(m*n) (m=#obs, n=0.9 (default) )
        - overall coverage in library 1 and 2 must be > 3
        - extend resulting sites by l steps in both directions. """
        m = self._get_bin_number()
        n = 0.9
        self._compute_score()
        print('before filter step:', len(self.scores), file=sys.stderr)
        self.indices_of_interest = np.where(self.scores > 2/(m*n))[0]
        print('after first filter step: ', len(self.indices_of_interest), file=sys.stderr)
        tmp = np.where(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) + np.squeeze(np.asarray(np.mean(self.overall_coverage[1], axis=0))) > 3)[0]
        tmp2 = np.intersect1d(self.indices_of_interest, tmp)
        print('length of intersection set: ', len(tmp), file=sys.stderr)
        self.indices_of_interest = tmp2
        print('after second filter step: ', len(self.indices_of_interest), file=sys.stderr)

        tmp = set()
        for i in self.indices_of_interest:
            for j in range(max(0, i-l), i+l+1):
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
            elif c1 / max(float(c2), 1) > threshold or c1-c2>10:
                state = 1
            elif c1 / max(float(c2), 1) < 1/threshold or c2-c1>10:
                state = 2
            else:
                state = 0
            
            states.append(state)
        
        f = open(filename, 'w')
        for j in range(len(states)):
            i = self.indices_of_interest[j]
            chrom, start, end = self._index2coordinates(i)
            s = states[j]
            print(chrom, start, end, s, self.first_overall_coverage[i], self.second_overall_coverage[i], sep='\t', file=f)
            
        f.close()
   
    
        
    def write_putative_regions(self, path):
        """Write putative regions (defined by criteria mentioned in method) as BED file."""
        with open(path, 'w') as f:
            for i in self.indices_of_interest:
                chrom, start, end = self._index2coordinates(i)
                print(chrom, start, end, file=f)
            
    def write_test_samples(self, name, l):
        f = open(name, 'w')
        
        for el in l:
            print(el, file=f)
        f.close()
    
    def get_training_set(self, exp_data, x, verbose, name, y):
        """Return linked genomic positions (at least <x> positions) to train HMM.
        Grep randomly a position within a putative region, and take then the entire region."""
        training_set = set()
        ts1, ts2 = set(), set()
        threshold = 2.0
        diff_cov = 10
        s0, s1, s2 = set(), set(), set()

        for i in range(len(self.indices_of_interest)):
            cov1, cov2 = self._get_covs(exp_data, i)
            
            #for parameter fitting for function
            if (cov1 / max(float(cov2), 1) > threshold and cov1+cov2 > diff_cov/2) or cov1-cov2 > diff_cov:
                s1.add((cov1, cov2))
            elif (cov1 / max(float(cov2), 1) < 1/threshold and cov1+cov2 > diff_cov/2) or cov2-cov1 > diff_cov:
                s2.add((cov1, cov2))
            else:
                s0.add((cov1, cov2))
            
            #for training set
            if cov1 / max(float(cov2), 1) > threshold or cov1-cov2 > diff_cov:
                ts1.add(i)
            if cov1 / max(float(cov2), 1) < 1/threshold or cov2-cov1 > diff_cov:
                ts2.add(i)
        
        if verbose:
            self.write_test_samples(name + '-s1', s1)
            self.write_test_samples(name + '-s2', s2)
        
        s0 = sample(s0, min(y, len(s0)))
        s1 = sample(s1, min(y, len(s1)))
        s2 = sample(s2, min(y, len(s2)))
        
        l = min(min(len(ts1), len(ts2)), x)
        tmp = set(sample(ts1, l)) | set(sample(ts2, l))
        
        for i in tmp:
            training_set.add(self.indices_of_interest[i])
            #search up
            while i+1 < len(self.indices_of_interest) and self.indices_of_interest[i+1] == self.indices_of_interest[i]+1:
                training_set.add(self.indices_of_interest[i+1])
                i += 1
            #search down
            while i-1 > 0 and self.indices_of_interest[i-1] == self.indices_of_interest[i]-1:
                training_set.add(self.indices_of_interest[i-1])
                i -= 1
        
        training_set = list(training_set)
        training_set.sort()
        if verbose:
            f=open(name + '-trainingset.bed', 'w')
            for l in training_set:
                chrom, s, e = self._index2coordinates(l)
                print(chrom, s, e, sep ='\t', file=f)
            f.close()
            
        return np.array(training_set), s0, s1, s2