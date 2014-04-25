#!/usr/bin/env python
# -*- coding: utf-8 -*-
from __future__ import print_function
import numpy as np
import sys
from math import log, fabs
from dpc_help import initialize
from dpc_help import dump_posteriors_and_viterbi
from dpc_help import get_peaks
from dpc_help import input
from hmm_binom_2d3s import BinomialHMM2d3s
from random import sample
from .. Util import GenomeData

def write(name, l):
    f = open(name, 'w')
    
    for el in l:
        print(el, file=f)
    f.close()

def get_init_parameters(name, indices_of_interest, first_overall_coverage, second_overall_coverage, \
                        verbose, x = 10000, threshold = 2.0, diff_cov = 10):

    tmp = sum( [ first_overall_coverage[i] + second_overall_coverage[i] for i in indices_of_interest]) / 2
    n_ = np.array([tmp, tmp])
    #print('n_: ', n_, file=sys.stderr)
    
    s0 = []
    s1 = []
    s2 = []
    so = []
    i = 1
    while i < x:
        state = None
        j = sample(indices_of_interest, 1)[0]
        
        cov1 = first_overall_coverage[j]
        cov2 = second_overall_coverage[j]
        
        if cov1 + cov2 <= 3:
            s0.append((cov1,cov2))
            so.append((cov1,cov2))
        elif (cov1 / max(float(cov2), 1) > threshold and cov1+cov2 > diff_cov/2) or cov1-cov2 > diff_cov:
            s1.append((cov1,cov2))
            so.append((cov1,cov2))
        elif (cov1 / max(float(cov2), 1) < 1/threshold and cov1+cov2 > diff_cov/2) or cov2-cov1 > diff_cov:
            s2.append((cov1,cov2))
            so.append((cov1,cov2))
        i += 1
    
    if verbose:
        write(name + 's0', s0)
        write(name + 's1', s1)
        write(name + 's2', s2)
        write(name + 'soverall', so)
    
    
    #get observation that occurs most often:
    m_ =[float(np.argmax(np.bincount(map(lambda x: x[0], s1)))), float(np.argmax(np.bincount(map(lambda x: x[1], s2)))) ]
    #print('m_', m_, file=sys.stderr)
    
    p_ = [[-1,-1,-1],[-1,-1,-1]] #first: 1. or 2. emission, second: state
    
    p_[0][0] = 1. / n_[0]
    p_[1][0] = 1. / n_[1]
       
    p_[0][1] = m_[0] / n_[0]
    p_[1][1] = p_[1][0]
    
    p_[0][2] = p_[0][0]
    p_[1][2] = m_[1] / n_[1]
    
    #print('p_', p_, file=sys.stderr)
    
    return n_, p_
 
def main():
    test = False
    options, bamfile_1, bamfile_2, genome = input(test)
    
    ######### WORK! ##########
    exp_data = initialize(name=options.name, genome_path=genome, regions=options.regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bam_file_1 = bamfile_1, ext_1=options.ext_1, \
                          bam_file_2 = bamfile_2, ext_2=options.ext_2, \
                          input_1=options.input_1, input_factor_1=options.input_factor_1, ext_input_1=options.ext_input_1, \
                          input_2=options.input_2, input_factor_2=options.input_factor_2, ext_input_2=options.ext_input_2, \
                          chrom_sizes = options.chrom_sizes,
                          verbose = options.verbose, norm_strategy=options.norm_strategy, no_gc_content=options.no_gc_content, deadzones=None)

    #print('Number of regions to be considered by the HMM:', len(exp_data), file=sys.stderr)
    exp_data.compute_putative_region_index()
    print('Number of regions with putative differential peaks:', len(exp_data.indices_of_interest), file=sys.stderr)
    
    if options.verbose:
        exp_data.write_putative_regions(options.name + '-putative-peaks.bed')
    
    print('Computing training set...',file=sys.stderr)
    training_set = exp_data.get_training_set(exp_data, min(len(exp_data.indices_of_interest) / 3, 600000), options.verbose, options.name)
    training_set_obs = exp_data.get_observation(training_set)
         
    n_, p_ = get_init_parameters(options.name, exp_data.indices_of_interest, exp_data.first_overall_coverage, exp_data.second_overall_coverage, options.verbose)
     
    print('Training HMM...', file=sys.stderr)
    m = BinomialHMM2d3s(n_components=3, n=n_, p=p_)
      
    if options.verbose:
        exp_data.get_initial_dist(options.name + '-initial-states.bed')
    
    m.fit([training_set_obs])
     
    if options.verbose:
        #print('p', m.p, file=sys.stderr)
        print("Final HMM's transition matrix: ", file=sys.stderr)
        print(m._get_transmat(), file=sys.stderr)
         
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
     
    if options.verbose: 
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
 
    get_peaks(name=options.name, states=states, DCS=exp_data, distr={'distr_name': "binomial", 'n': m.n[0], 'p': m.p[0][1]})
    
if __name__ == '__main__':
    main()
    