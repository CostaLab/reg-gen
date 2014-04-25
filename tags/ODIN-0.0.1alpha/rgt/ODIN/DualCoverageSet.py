from __future__ import print_function
from .. CoverageSet import CoverageSet
import numpy as np
from random import sample, randrange
from time import time
from help_hmm import EPSILON
from gc_content import get_gc_context
import sys
from normalize import get_normalization_factor

class DualCoverageSet():
    def __init__(self, name, region, genome_path, binsize, stepsize, rmdup, file_1, ext_1, file_2, ext_2,\
                 input_1, ext_input_1, input_factor_1, input_2, ext_input_2, input_factor_2, chrom_sizes, verbose, norm_strategy, no_gc_content, deadzones):
        self.genomicRegions = region
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        
        self.cov1 = CoverageSet('first file', region)
        self.cov2 = CoverageSet('second file', region)
        
        print("Loading reads...", file=sys.stderr)
        self.cov1.coverage_from_bam(bam_file=file_1, read_size=ext_1, rmdup=rmdup, binsize=binsize, stepsize=stepsize, mask_file=deadzones)
        self.cov2.coverage_from_bam(bam_file=file_2, read_size=ext_2, rmdup=rmdup, binsize=binsize, stepsize=stepsize, mask_file=deadzones)
        
        map_input = {1: {'input': input_1, 'input_factor': input_factor_1, 'ext': ext_input_1, 'cov-ip': self.cov1, 'ip': file_1}, 
                     2: {'input': input_2, 'input_factor': input_factor_2, 'ext': ext_input_2, 'cov-ip': self.cov2, 'ip': file_2}}
        
        if no_gc_content:
            print("Do not consider GC content model", file=sys.stderr)
        else:
            print("Computing GC content", file=sys.stderr)
            
        for i in [1, 2]:
            #get GC-content
            input = map_input[i]
            
            if verbose:
                input['cov-ip'].write_bigwig(name + '-gc-s%s-1.bw' %i, chrom_sizes)
                #input['cov-ip'].write_bed(name + '-gc-s%s-1.bed' %i)
            
            if input['input'] is not None:
                input['cov-input'] = CoverageSet('%s file' %input['input'], region)
                input['cov-input'].coverage_from_bam(bam_file=input['input'], read_size=input['ext'], rmdup=rmdup, binsize=binsize, stepsize=stepsize)
                map_input[i]['cov-input'] = input['cov-input']
            
            if not no_gc_content:
                if input['input'] is not None:
                    gc_content_cov, avg_gc_content, gc_hist = get_gc_context(stepsize, binsize, genome_path, input['cov-input'].coverage)
                    #print("Gamma: %s" %avg_gc_content, file=sys.stderr)
                    self._norm_gc_content(input['cov-ip'].coverage, gc_content_cov, avg_gc_content)
                    self._norm_gc_content(input['cov-input'].coverage, gc_content_cov, avg_gc_content)
                    
                if verbose:
                    if not no_gc_content:
                        self.print_gc_hist(name + '-s%s-' %i, gc_hist)
                    if input['input'] is not None:
                        input['cov-input'].write_bigwig(name + '-gc-s%s-input-2.bw' %i, chrom_sizes)
                        #input['cov-input'].write_bed(name + '-gc-s%s-input-2.bed' %i)
                    input['cov-ip'].write_bigwig(name + '-gc-s%s-2.bw' %i, chrom_sizes)
                    #input['cov-ip'].write_bed(name + '-gc-s%s-2.bed' %i)
        
        norm_done = False
        for i in [1, 2]:
            input = map_input[i]
            #print(sum([sum(self.cov1.coverage[j]) for j in range(len(self.cov1.genomicRegions))]), file=sys.stderr)
            #print(sum([sum(self.cov2.coverage[j]) for j in range(len(self.cov2.genomicRegions))]), file=sys.stderr)
            norm_done = self.normalization(map_input, i, norm_strategy, norm_done, name, verbose)
            #print(sum([sum(self.cov1.coverage[j]) for j in range(len(self.cov1.genomicRegions))]), file=sys.stderr)
            #print(sum([sum(self.cov2.coverage[j]) for j in range(len(self.cov2.genomicRegions))]), file=sys.stderr)
            
            if verbose:
                if input['input'] is not None:
                    input['cov-input'].write_bigwig(name + '-gc-s%s-input-3'%i, chrom_sizes)
                    #input['cov-input'].write_bed(name + '-gc-s%s-input-3.bed'%i)
                input['cov-ip'].write_bigwig(name + '-gc-s%s-3.bw'%i, chrom_sizes)
                #input['cov-ip'].write_bed(name + '-gc-s%s-3.bed'%i)
                
        #make one array for the coverage
        self.first_overall_coverage = reduce(lambda x,y: np.concatenate((x,y)), [self.cov1.coverage[i] for i in range(len(self.cov1.genomicRegions))])
        self.second_overall_coverage = reduce(lambda x,y: np.concatenate((x,y)), [self.cov2.coverage[i] for i in range(len(self.cov2.genomicRegions))])
        assert(len(self.first_overall_coverage) == len(self.second_overall_coverage))
        
        self.scores = np.zeros(len(self.first_overall_coverage))
        self.indices_of_interest = []
    
    
    def normalization(self, map_input, i, norm_strategy, norm_done, name, verbose):
        input = map_input[i]
        
        #compute normalization factor
        #pre-defined values
        if input['input_factor'] is not None and i != 1:
            print("Normalize by Diaz and pre-defined values...", input['input_factor'], file=sys.stderr)
            print("Normalize file 1 with input normalization factor %s" %(map_input[1]['input_factor']), file=sys.stderr)
            print("Normalize file 2 with input normalization factor %s" %(map_input[2]['input_factor']), file=sys.stderr)
            
            map_input[1]['cov-input'].scale(map_input[1]['input_factor'])
            map_input[2]['cov-input'].scale(map_input[2]['input_factor'])
            map_input[1]['cov-ip'].subtract(map_input[1]['cov-input'])
            map_input[2]['cov-ip'].subtract(map_input[2]['cov-input'])

        #naive norm.
        if not norm_done and norm_strategy == 1: 
            
            s1 = sum([sum(map_input[1]['cov-ip'].coverage[i]) for i in range(len(map_input[1]['cov-ip'].genomicRegions))])
            s2 = sum([sum(map_input[2]['cov-ip'].coverage[i]) for i in range(len(map_input[2]['cov-ip'].genomicRegions))])
            
            if s1 > s2:
                map_input[2]['cov-ip'].scale(s1/float(s2))
                print("Factor: normalize file 2 by signal with factor %s: " %(s1/float(s2)), file=sys.stderr)
            elif s2 >= s1:
                print("Factor: normalize file 1 by signal with factor %s: " %(s2/float(s1)), file=sys.stderr)
                map_input[1]['cov-ip'].scale(s2/float(s1))

            norm_done = True
        
        #Diaz
        if norm_strategy == 2:
            print("Normalize by input (Diaz)", file=sys.stderr)
            _, input['input_factor'] = get_normalization_factor(input['ip'], input['input'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm' + str(i), verbose=verbose, two_sample=False) #todo: solve the genome parameter
            print("Factor: normalize file %s with input normalization factor %s" %(i, input['input_factor']), file=sys.stderr)
            input['cov-input'].scale(input['input_factor']) #scale input for normalization
            input['cov-ip'].subtract(input['cov-input']) #TODO: case: multiple input files
        
        #own
        if not norm_done and norm_strategy == 3:
            print("Normalize by own method", file=sys.stderr)
            smaller_sample, factor = get_normalization_factor(map_input[1]['ip'], map_input[2]['ip'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm'+ str(i), verbose=verbose, two_sample=True)
            factor = max(1/factor, factor) #increase current
            print("Factor: normalize file %s with factor %s" %(smaller_sample, factor), file=sys.stderr)
            map_input[smaller_sample]['cov-ip'].scale(factor)
            norm_done = True
        
        #diaz and own
        if i != 1 and norm_strategy == 4:
            print("Normalize by input (Diaz) and our own method", file=sys.stderr)
            #apply diaz
            _, map_input[1]['input_factor'] = get_normalization_factor(map_input[1]['ip'], map_input[1]['input'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm' + str(i), verbose=verbose, two_sample=False)
            _, map_input[2]['input_factor'] = get_normalization_factor(map_input[2]['ip'], map_input[2]['input'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm' + str(i), verbose=verbose, two_sample=False)
            
            print("Factor: normalize input with factor %s and %s" %(map_input[1]['input_factor'], map_input[2]['input_factor']), file=sys.stderr)
            map_input[1]['cov-input'].scale(map_input[1]['input_factor'])
            map_input[2]['cov-input'].scale(map_input[2]['input_factor'])
            
            map_input[1]['cov-ip'].subtract(map_input[1]['cov-input'])
            map_input[2]['cov-ip'].subtract(map_input[2]['cov-input'])
            
            #apply our own
            smaller_sample, factor = get_normalization_factor(map_input[1]['ip'], map_input[2]['ip'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm-sample' + str(i), verbose=verbose, two_sample=True)
            
            factor = max(1/factor, factor) #increase current
            print("Factor: normalize file %s with factor %s" %(smaller_sample, factor), file=sys.stderr)
            map_input[smaller_sample]['cov-ip'].scale(factor)
        
        
        #diaz and naive
        if i != 1 and norm_strategy == 5:
            print("Normalizing...", file=sys.stderr)
            #apply diaz
            _, map_input[1]['input_factor'] = get_normalization_factor(map_input[1]['ip'], map_input[1]['input'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm' + str(i), verbose=verbose, two_sample=False)
            _, map_input[2]['input_factor'] = get_normalization_factor(map_input[2]['ip'], map_input[2]['input'], step_width=1000, zero_counts=0, \
                                                              genome='mm9', filename=name + '-norm' + str(i), verbose=verbose, two_sample=False)
            
            print("Normalize input with factor %s and %s" %(map_input[1]['input_factor'], map_input[2]['input_factor']), file=sys.stderr)
            map_input[1]['cov-input'].scale(map_input[1]['input_factor'])
            map_input[2]['cov-input'].scale(map_input[2]['input_factor'])
            
            map_input[1]['cov-ip'].subtract(map_input[1]['cov-input'])
            map_input[2]['cov-ip'].subtract(map_input[2]['cov-input'])
            
            #apply naive method
            s1 = sum([sum(map_input[1]['cov-ip'].coverage[i]) for i in range(len(map_input[1]['cov-ip'].genomicRegions))])
            s2 = sum([sum(map_input[2]['cov-ip'].coverage[i]) for i in range(len(map_input[2]['cov-ip'].genomicRegions))])
            
            if s1 > s2:
                map_input[2]['cov-ip'].scale(s1/float(s2))
                print("Normalize file 2 by signal with factor %s: " %(s1/float(s2)), file=sys.stderr)
            elif s2 >= s1:
                print("Normalize file 1 by signal with factor %s: " %(s2/float(s1)), file=sys.stderr)
                map_input[1]['cov-ip'].scale(s2/float(s1))

        
        return norm_done

    
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
            last += len(self.cov1.coverage[i])
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
            mask = np.array([True]*len(self.first_overall_coverage))
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
        #print('before filter step:', len(self.scores), file=sys.stderr)
        self.indices_of_interest = np.where(self.scores > 2/(m*n))[0]
        #print('after first filter step: ', len(self.indices_of_interest), file=sys.stderr)
        tmp = np.where(self.first_overall_coverage + self.second_overall_coverage > 3)[0]
        tmp2 = np.intersect1d(self.indices_of_interest, tmp)
        #print('length of intersection set: ', len(tmp), file=sys.stderr)
        self.indices_of_interest = tmp2
        #print('after second filter step: ', len(self.indices_of_interest), file=sys.stderr)
        #extend regions by l steps
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
            
    def get_training_set(self, exp_data, x, verbose, name):
        """Return linked genomic positions (at least <x> positions) to train HMM.
        Grep randomly a position within a putative region, and take then the entire region."""
        training_set = set()
        ts1 = set()
        ts2 = set()
        threshold = 2.0
        diff_cov = 10

        for i in range(len(self.indices_of_interest)):
            cov1 = exp_data.first_overall_coverage[self.indices_of_interest[i]]
            cov2 = exp_data.second_overall_coverage[self.indices_of_interest[i]]
            
            if cov1 / max(float(cov2), 1) > threshold or cov1-cov2 > diff_cov:
                ts1.add(i)
            if cov1 / max(float(cov2), 1) < 1/threshold or cov2-cov1 > diff_cov:
                ts2.add(i)
        
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
        if verbose:
            f=open(name + '-trainingset.bed', 'w')
            for l in training_set:
                chrom, s, e = self._index2coordinates(l)
                print(chrom, s, e, sep ='\t', file=f)
            f.close()
            
        return np.array(training_set)
