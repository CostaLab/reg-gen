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
from dpc_help import get_back
from neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters
from random import sample
import multiprocessing
from tracker import Tracker
from rgt.THOR.neg_bin import NegBin
import cProfile

def _get_pvalue_distr(exp_data, mu, alpha, tracker):
    """Derive NB1 parameters for p-value calculation"""
    mu = mu[0,0]
    #alpha = max(3, alpha[0,0])
    alpha = alpha[0,0] / 10000.
    tracker.write(text=str(mu) + " " + str(alpha), header="Neg. Bin. distribution for p-value estimates (mu, alpha)")
    nb = NegBin(mu, alpha)
    return {'distr_name': 'nb', 'distr': nb}

    #print(mu, alpha, file=sys.stderr)
    #n = np.mean([np.sum(exp_data.overall_coverage[i]) for i in range(2)])
    #print(n, file=sys.stderr)
    #p = mu / n
    #print(p, file=sys.stderr)
    #tracker.write(text=str(n) + " " + str(p), header="Bin. distribution for p-value estimates (n, p)")
    #return {'distr_name': 'binomial', 'n': n, 'p': p}

def main():
    test = False
    options, bamfiles, regions, genome, chrom_sizes, dims, inputs = input(test)

    ######### WORK! ##########
    tracker = Tracker(options.name + '-setup.info')
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, debug=options.debug,\
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes, \
                          tracker=tracker, norm_regions=options.norm_regions, scaling_factors_ip = options.scaling_factors_ip, save_wig=options.save_wig)
    
    func, func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, options.debug, sample_size=20000)
    tracker.write(text=func_para[0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) ")
    tracker.write(text=func_para[1])
    
    exp_data.compute_putative_region_index()

    print('Compute training set...',file=sys.stderr)

    training_set, s0, s1, s2 = exp_data.get_training_set(test, exp_data, options.debug, options.name, 10000, 3)
    training_set_obs = exp_data.get_observation(training_set)
    
    init_alpha, init_mu = get_init_parameters(s0, s1, s2)
    tracker.write(text=init_mu, header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu,alpha)")
    tracker.write(text=init_alpha)
    
    print('Training HMM...', file=sys.stderr)
    m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
    m.fit([training_set_obs])

    tracker.write(text=m.mu, header="Final HMM's Neg. Bin. Emission distribution (mu,alpha)")
    tracker.write(text=m.alpha)
    tracker.write(text=m._get_transmat(), header="Transmission matrix")
    
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    back_var, back_mean = get_back(exp_data, states)
    tracker.write(text=back_var, header="background variance")
    tracker.write(text=back_mean, header="background mean")
    a = (back_var - back_mean) / (back_mean ** 2)
    tracker.write(text=a, header="new alpha")
    
    if options.debug:
        posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
    
    distr = _get_pvalue_distr(exp_data, m.mu, m.alpha, tracker)
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr, merge=options.merge, exts=exp_data.exts, pcutoff=options.pcutoff, p=options.par)
    
if __name__ == '__main__':
    main() 
    
