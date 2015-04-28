from __future__ import print_function
from rgt.CoverageSet import CoverageSet
import numpy as np
from random import sample
from rgt.ODIN.gc_content import get_gc_context
import sys
from rgt.ODIN.normalize import get_normalization_factor
from math import fabs
from rgt.ODIN.DualCoverageSet import DualCoverageSet
from rgt.GenomicRegionSet import GenomicRegionSet
from copy import deepcopy

EPSILON = 1**-320

class MultiCoverageSet(DualCoverageSet):
    def _help_init(self, path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, dim, regions, norm_regionset):
        """Return self.covs and self.inputs as CoverageSet"""
        self.exts = exts
        self.covs = [CoverageSet('file' + str(i), regions) for i in range(dim)]
        for i, c in enumerate(self.covs):
            c.coverage_from_bam(bam_file=path_bamfiles[i], read_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize)
        self.covs_avg = [CoverageSet('cov_avg'  + str(i) , regions) for i in range(2)]
        if path_inputs:
            self.inputs = [CoverageSet('input' + str(i), regions) for i in range(len(path_inputs))]
            for i, c in enumerate(self.inputs):
                c.coverage_from_bam(bam_file=path_inputs[i], read_size=exts_inputs[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize)
            self.input_avg = [CoverageSet('input_avg'  + str(i), regions) for i in range(2)]
        else:
            self.inputs = []
            
        if norm_regionset:
            self.norm_regions = [CoverageSet('norm_region' + str(i), norm_regionset) for i in range(dim)]
            for i, c in enumerate(self.norm_regions):
                c.coverage_from_bam(bam_file=path_bamfiles[i], read_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                    stepsize=stepsize)
            self.input_avg = [CoverageSet('input_avg'  + str(i), regions) for i in range(2)]
        else:
            self.norm_regions = None
    
    def _get_covs(self, DCS, i):
        """For a multivariant Coverageset, return coverage cov1 and cov2 at position i"""
        cov1 = int(np.mean(DCS.overall_coverage[0][:,DCS.indices_of_interest[i]]))
        cov2 = int(np.mean(DCS.overall_coverage[1][:,DCS.indices_of_interest[i]]))
    
        return cov1, cov2
    
    def _compute_gc_content(self, no_gc_content, verbose, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes, chrom_sizes_dict):
        """Compute GC-content"""
        if not no_gc_content and path_inputs:
            print("Compute GC-content", file=sys.stderr)
            for i, cov in enumerate(self.covs):
                inputfile = self.inputs[i] #1 to 1 mapping between input and cov
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 1 if i < self.dim_1 else 2
                gc_content_cov, avg_gc_content, gc_hist = get_gc_context(stepsize, binsize, genome_path, inputfile.coverage, chrom_sizes_dict)
                self._norm_gc_content(cov.coverage, gc_content_cov, avg_gc_content)
                self._norm_gc_content(inputfile.coverage, gc_content_cov, avg_gc_content)
            
                if verbose:
                    self.print_gc_hist(name + '-s%s-rep%s-' %(sig, rep), gc_hist)
                    cov.write_bigwig(name + '-s%s-rep%s-gc.bw' %(sig, rep), chrom_sizes)
        else:
            print("Do not compute GC-content, as there is no input.", file=sys.stderr)
    
    
    def _output_bw(self, name, chrom_sizes, save_wig):
        """Output bigwig files"""
        for i in range(len(self.covs)):
            rep = i if i < self.dim_1 else i-self.dim_1
            sig = 1 if i < self.dim_1 else 2
            if self.inputs:
                self.inputs[i].write_bigwig(name + '-s%s-rep%s-input.bw' %(sig, rep), chrom_sizes, save_wig)
            self.covs[i].write_bigwig(name + '-s%s-rep%s.bw' %(sig, rep), chrom_sizes, save_wig)
        
        ra = [self.covs_avg, self.input_avg] if self.inputs else [self.covs_avg]
        for k, d in enumerate(ra):
            g = self.covs if k == 0 else self.inputs
            for j in range(2):
                d[j] = deepcopy(g[0]) if j == 0 else deepcopy(g[self.dim_1])
                r = range(1, self.dim_1) if j == 0 else range(self.dim_1 + 1, self.dim_1 + self.dim_2)
                f = 1./self.dim_1 if j == 0 else 1./self.dim_2
                for i in r:
                    d[j].add(g[i])
                d[j].scale(f)
                n = name + '-s%s.bw' %(j+1) if k == 0 else name + '-s%s-input.bw' %(j+1)
                d[j].write_bigwig(n, chrom_sizes, save_wig)

    def __init__(self, name, dims, regions, genome_path, binsize, stepsize, chrom_sizes, norm_regionset, \
                 verbose, debug, no_gc_content, rmdup, path_bamfiles, exts, path_inputs, exts_inputs, \
                 factors_inputs, chrom_sizes_dict, scaling_factors_ip, save_wig):
        """Compute CoverageSets, GC-content and normalization"""
        self.genomicRegions = regions
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.chrom_sizes_dict = chrom_sizes_dict
        self.dim_1, self.dim_2 = dims
        
        #make data nice
        self._help_init(path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, sum(dims), regions, norm_regionset)
        self._compute_gc_content(no_gc_content, verbose, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes, chrom_sizes_dict)
        self._normalization_by_input(path_bamfiles, path_inputs, name, debug)
        self._normalization_by_signal(name, scaling_factors_ip)
        
        self._output_bw(name, chrom_sizes, save_wig) 
        
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
    
    def get_max_colsum(self):
        """Sum over all columns and add maximum"""
        return self.overall_coverage[0].sum(axis=0).max() + self.overall_coverage[1].sum(axis=0).max()
    
    def output_overall_coverage(self, path):
        for j in range(2):
            f = open(path + str(j), 'w')
            for i in range(self.overall_coverage[j].shape[1]):
                print(self.overall_coverage[j][:,i].T, file=f)
    
    def _normalization_by_input(self, path_bamfiles, path_inputs, name, debug):
        """Normalize with regard to input file"""
        if path_inputs:
            print("Normalize by input-DNA", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 0 if i < self.dim_1 else 1
                j = 0 if i < self.dim_1 else 1
                _, n = get_normalization_factor(path_bamfiles[i], path_inputs[i], step_width=1000, zero_counts=0, \
                                                filename=name + '-norm' + str(i), debug=debug, chrom_sizes_dict=self.chrom_sizes_dict, two_sample=False, stop=True)
                
                if n is not None:
                    print("Normalize input of Signal %s, Rep %s with factor %s"\
                           %(sig, rep, round(n, 3)) , file=sys.stderr)
                    self.inputs[i].scale(n)
                    self.covs[i].subtract(self.inputs[i])
    
    def _normalization_by_signal(self, name, scaling_factors_ip):
        """Normalize signal"""
        if scaling_factors_ip:
            print("Normalize signal by predefined scaling factors", file=sys.stderr)
            assert len(scaling_factors_ip) == len(self.covs)
            for i in range(len(scaling_factors_ip)):
                self.covs[i].scale(scaling_factors_ip[i])
        else:
        
            if self.norm_regions:
                print("Normalize by signal (on specified regions)", file=sys.stderr)
                signals = [sum([sum(self.norm_regions[k].coverage[i]) for i in range(len(self.norm_regions[k].genomicRegions))]) for k in range(self.dim_1 + self.dim_2)]
            else:
                print("Normalize by signal (genomewide)", file=sys.stderr)
                signals = [sum([sum(self.covs[k].coverage[i]) for i in range(len(self.covs[k].genomicRegions))]) for k in range(self.dim_1 + self.dim_2)]
            
            means_signal = [np.mean(signals[:self.dim_1]), np.mean(signals[self.dim_1:])]
            max_index = means_signal.index(max(means_signal))
            
            if max_index == 1:
                r = range(self.dim_1)
                f = means_signal[1] / means_signal[0]
                print("Normalize first set of replicates with factor %s" %(round(f, 2)), file=sys.stderr)
            if max_index == 0:
                r = range(self.dim_1, self.dim_1 + self.dim_2)
                f = means_signal[0] / means_signal[1]
                print("Normalize second set of replicates with factor %s" %(round(f, 2)), file=sys.stderr)
            
            for i in r:
                self.covs[i].scale(f)
           
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
        if not mask.size:
            mask = np.array([True]*self._get_bin_number())
            
        return np.asarray(np.concatenate((self.overall_coverage[0][:,mask].T, self.overall_coverage[1][:,mask].T), axis=1))
    
    def _compute_score(self):
        """Compute score for each observation (based on Xu et al.)"""
        self.scores = sum([np.squeeze(np.asarray(np.mean(self.overall_coverage[i], axis=0))) / float(np.mean(self.overall_coverage[i])) for i in range(2)])
    
    def _get_bin_number(self):
        """Return number of bins"""
        return self.overall_coverage[0].shape[1]
    
    def compute_putative_region_index(self, l=5):
        """Compute putative differential peak regions as follows: 
        - score must be > 0, i.e. everthing
        - overall coverage in library 1 and 2 must be > 3"""
        
        self._compute_score()
        self.indices_of_interest = np.where(self.scores > 0)[0] #2/(m*n)
        tmp = np.where(np.squeeze(np.asarray(np.mean(self.overall_coverage[0], axis=0))) + np.squeeze(np.asarray(np.mean(self.overall_coverage[1], axis=0))) > 3)[0]
        tmp2 = np.intersect1d(self.indices_of_interest, tmp)
        self.indices_of_interest = tmp2

        #tmp = set()
        #for i in self.indices_of_interest:
        #    for j in range(max(0, i-l), i+l+1):
        #        tmp.add(j)
        #tmp = list(tmp)
        #tmp.sort()
        #self.indices_of_interest = np.array(tmp)
         
    def write_test_samples(self, name, l):
        f = open(name, 'w')
        
        for el1, el2 in l:
            print(el1, el2, sep='\t', file=f)
        f.close()
    
    def debug_output_get_training_set(self, name, training_set, s0_v, s1_v, s2_v):
        """Output debug info for training_set computation."""
        f=open(name + '-trainingset.bed', 'w')
        for l in training_set:
            chrom, s, e = self._index2coordinates(l)
            print(chrom, s, e, sep ='\t', file=f)
        f.close()
        
        self.write_test_samples(name + '-s0', s0_v)
        self.write_test_samples(name + '-s1', s1_v)
        self.write_test_samples(name + '-s2', s2_v)
    
    def get_training_set(self, test, exp_data, debug, name, y=5000, ex=2):
        """Return genomic positions (max <y> positions) and enlarge them by <ex> bins to train HMM."""
        threshold = 1.3
        #diff_cov = 20
        
        diff_cov = max(20, np.percentile(np.append(np.asarray(self.overall_coverage[0].flatten())[0], np.asarray(self.overall_coverage[1].flatten())[0]), 95))
        t = np.percentile(np.append(np.asarray(self.overall_coverage[0].flatten())[0], np.asarray(self.overall_coverage[1].flatten())[0]), 95)
        print('training diff_cov: %s (%s)' %(diff_cov, t), file=sys.stderr)
        
        if test:
            diff_cov = 2
            threshold = 1.5
        
        s0, s1, s2 = [], [], []
        
        #if debug:
        #    self.output_overall_coverage('signal')
        
        rep=True
        while rep:
            for i in range(len(self.indices_of_interest)):
                cov1, cov2 = self._get_covs(exp_data, i)
    
                #apply criteria for initial peak calling
                if (cov1 / max(float(cov2), 1) > threshold and cov1+cov2 > diff_cov/2) or cov1-cov2 > diff_cov:
                    s1.append((i, cov1, cov2))
                elif (cov1 / max(float(cov2), 1) < 1/threshold and cov1+cov2 > diff_cov/2) or cov2-cov1 > diff_cov:
                    s2.append((i, cov1, cov2))
                elif fabs(cov1 - cov2) < diff_cov/2 and cov1 + cov2 > diff_cov/4:
                    s0.append((i, cov1, cov2))
            
            print(diff_cov, threshold, len(s0), len(s1), len(s2), file=sys.stderr)
            
            if diff_cov == 1 and threshold == 1.1:
                print("No differential peaks detected", file=sys.stderr)
                sys.exit()
            
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
            
        
        
        
        tmp = []
        for i, el in enumerate([s0, s1, s2]):
            el = np.asarray(el)
            if not test:
                el = el[el[:,1] < np.percentile(el[:,1], 90)]
                el = el[el[:,2] < np.percentile(el[:,2], 90)]
            tmp.append(el)
        
        s0 = tmp[0]
        s1 = tmp[1]
        s2 = tmp[2]
        
        l = np.min([len(s1), len(s2), len(s0), y])
        
        s0 = sample(s0, l)
        s1 = sample(s1, l)
        s2 = sample(s2, l)
        
        s0_v = map(lambda x: (x[1], x[2]), s0)
        s1_v = map(lambda x: (x[1], x[2]), s1)
        s2_v = map(lambda x: (x[1], x[2]), s2)
        
        #enlarge training set(assumption everything is in indices_of_interest)
        extension_set = set()
        for i, _, _ in s0 + s1 + s2:
            for j in range(max(0, i - ex), i + ex + 1):
                extension_set.add(j)
        
        tmp = s0 + s1 + s2
        training_set = map(lambda x: x[0], tmp) + list(extension_set)
        
        training_set = list(training_set)
        training_set.sort()
        
        if debug:
            self.debug_output_get_training_set(name, training_set, s0_v, s1_v, s2_v)
            
        return np.array(training_set), s0_v, s1_v, s2_v
