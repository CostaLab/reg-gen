#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog

Find differential peaks in regions.

Author: Manuel Allhoff (allhoff@aices.rwth-aachen.de)

"""

from __future__ import print_function
import numpy as np
import sys
from math import log, fabs
from dpc_help import initialize
from dpc_help import dump_posteriors_and_viterbi
from dpc_help import get_peaks
from dpc_help import input
from dpc_help import _fit_mean_var_distr
from neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters
from random import sample
import multiprocessing

def write(name, l):
    f = open(name, 'w')
    
    for el in l:
        print(el, file=f)
    f.close()

def _get_training_sets(indices_of_interest, overall_coverage, name, verbose, x = 10000, threshold = 2.0, diff_cov = 10):
    """Return s0,s1,s2 (list of tuples) with initial peaks"""
    s0 = []
    s1 = []
    s2 = []
    so = []
    i = 1
    while i < x:
        j = sample(indices_of_interest, 1)[0]
        
        cov1 = exp_data.overall_coverage[0][:,indices_of_interest[j]].sum()
        cov2 = exp_data.overall_coverage[1][:,indices_of_interest[j]].sum()
        
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
        write(name + '-s0', s0)
        write(name + '-s1', s1)
        write(name + '-s2', s2)
        write(name + '-soverall', so)
        
    return s0, s1, s2

if __name__ == '__main__':
    test = True
    options, bamfiles, regions, genome, chrom_sizes, dims, inputs = input(test)
    
    ######### WORK! ##########
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, \
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes)
    
    func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, verbose = True, cut=options.cut_obs, sample_size=200)
    print('Number of regions to be considered by the HMM:', len(exp_data), file=sys.stderr)
    exp_data.compute_putative_region_index()
    print('Number of regions with putative differential peaks:', len(exp_data.indices_of_interest), file=sys.stderr)

    if options.verbose:
        exp_data.write_putative_regions(options.name + '-putative-peaks.bed')
    print('Compute training set...',file=sys.stderr)
    
    training_set = exp_data.get_training_set(exp_data, min(len(exp_data.indices_of_interest) / 3, 600000), options.verbose, options.name)
    training_set_obs = exp_data.get_observation(training_set)
    
    _, s1, s2 = _get_training_sets(exp_data.indices_of_interest, exp_data.overall_coverage, options.name, options.verbose)
    init_alpha, init_mu, init_para_func = get_init_parameters(s1, s2)
      
    print('Training HMM...', file=sys.stderr)
    m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], para_func = init_para_func)

#     if options.verbose:
#         exp_data.get_initial_dist(options.name + '-initial-states.bed')
       
    m.fit([training_set_obs])
      
    print('...done', file=sys.stderr)
    
#     if options.verbose:
#         print('p', m.p, file=sys.stderr)
#         print("Final HMM's transition matrix: ", file=sys.stderr)
#         print(m._get_transmat(), file=sys.stderr)
#          
#     print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
#     posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
#     states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
#     print("...done", file=sys.stderr)
#      
#     if options.verbose: 
#         dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
#  
#     print('t', m.n, m.p, file=sys.stderr)
#     get_peaks(name=options.name, states=states, DCS=exp_data, distr={'distr_name': "binomial", 'n': m.n[0], 'p': m.p[0][1]})
     
     
    
