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
from rgt.THOR.neg_bin import NegBin
import cProfile

def _get_pvalue_distr(exp_data, mu, alpha, tracker):
    """Derive NB1 parameters for p-value calculation"""
    mu = mu[0,0]
    alpha = alpha[0,0]
    tracker.write(text=str(mu) + " " + str(alpha), header="Neg. Bin. distribution for p-value estimates (mu, alpha)")
    nb = NegBin(mu, alpha)
    print(mu, alpha, file=sys.stderr)
    n = np.mean([np.sum(exp_data.overall_coverage[i]) for i in range(2)])
    print(n, file=sys.stderr)
    p = mu / n
    print(p, file=sys.stderr)
    return {'distr_name': 'binomial', 'n': n, 'p': p}
    #return {'distr_name': 'nb', 'distr': nb}

def main():
    test = False
    options, bamfiles, regions, genome, chrom_sizes, dims, inputs = input(test)
    print(options.debug, file=sys.stderr)
    ######### WORK! ##########
    tracker = Tracker(options.name + '-setup.info')
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, debug=options.debug,\
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes, tracker=tracker)
    
    func, func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, options.debug, sample_size=20000)
    tracker.write(text=str(func_para), header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) ")
    exp_data.compute_putative_region_index()

    print('Compute training set...',file=sys.stderr)

    training_set, s0, s1, s2 = exp_data.get_training_set(test, exp_data, options.debug, options.name, 10000, 3)
    training_set_obs = exp_data.get_observation(training_set)
    
    init_alpha, init_mu = get_init_parameters(s0, s1, s2)
    tracker.write(text=str(init_mu)  + str('\n') + str(init_alpha), header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu,alpha)")
    
    print('Training HMM...', file=sys.stderr)
    m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
    m.fit([training_set_obs])

    tracker.write(text=str(m.mu)  + str('\n') + str(m.alpha), header="Final HMM's Neg. Bin. Emission distribution (mu,alpha)")
    tracker.write(text=str(m._get_transmat()), header="Transmission matrix")
    
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    
    if options.debug:
        posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
    
    distr = _get_pvalue_distr(exp_data, m.mu, m.alpha, tracker)
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr, merge=options.merge, exts=exp_data.exts, pcutoff=options.pcutoff)
    
if __name__ == '__main__':
    main() 
    
