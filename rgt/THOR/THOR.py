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
from tracker import Tracker

if __name__ == '__main__':
    test = False
    options, bamfiles, regions, genome, chrom_sizes, dims, inputs = input(test)
    
    ######### WORK! ##########
    tracker = Tracker(options.name + '-setup.info')
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, \
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes, tracker=tracker)
    
    func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, verbose = True, cut=options.cut_obs, sample_size=200)
    print("func para", func_para, file=sys.stderr)
    #func_para=[[1, 2, 0.3], [1, 2, 0.3]]
    print('Number of regions to be considered by the HMM:', len(exp_data), file=sys.stderr)
    exp_data.compute_putative_region_index()
    print('Number of regions with putative differential peaks:', len(exp_data.indices_of_interest), file=sys.stderr)

    if options.verbose:
        exp_data.write_putative_regions(options.name + '-putative-peaks.bed')
    print('Compute training set...',file=sys.stderr)
    
    training_set, s1, s2 = exp_data.get_training_set(exp_data, min(len(exp_data.indices_of_interest) / 3, 600000), options.verbose, options.name, 10000)
    training_set_obs = exp_data.get_observation(training_set)
    
    init_alpha, init_mu = get_init_parameters(s1, s2)
      
    print('Training HMM...', file=sys.stderr)
    m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], para_func = func_para)

#     if options.verbose:
#         exp_data.get_initial_dist(options.name + '-initial-states.bed')
       
    m.fit([training_set_obs])
    
    print('...done', file=sys.stderr)
    tracker.write(text=str(m.mu) + "\n" + str(m.alpha), header="HMM's Neg. Bin. Emission (mu,alpha)")
    
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    print("...done", file=sys.stderr)
    
    if options.verbose:
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)

    n = np.mean([exp_data.overall_coverage[i][:,exp_data.indices_of_interest].sum() for i in range(2)])
    p = m.mu[0,0] / n
    
    tracker.write(text=str(n) + " " + str(p), header="Binomial distr for p-value estimates (n,p)")
    
    get_peaks(name=options.name, states=states, DCS=exp_data, distr={'distr_name': "binomial", 'n': n, 'p': p})
     
     
    
