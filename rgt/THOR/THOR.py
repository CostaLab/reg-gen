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
from dpc_help import _fit_mean_var_distr
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
    
    options, bamfiles, genome, chrom_sizes, dims, inputs = input(test)

    ######### WORK! ##########
    tracker = Tracker(options.name + '-setup.info')
    
    #tracker.write('test')
    #sys.exit()
    
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=options.regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, debug=options.debug,\
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes, \
                          tracker=tracker, norm_regions=options.norm_regions, scaling_factors_ip = options.scaling_factors_ip, save_wig=options.save_wig, \
                          housekeeping_genes=options.housekeeping_genes, test=test)
    
    func, func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, options.debug, sample_size=5000, verbose=options.debug, outputdir = options.outputdir, report=options.report)
    tracker.write(text=func_para[0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) ")
    tracker.write(text=func_para[1])
    sys.stderr.flush()
    exp_data.compute_putative_region_index()

    print('Compute training set...',file=sys.stderr)
    sys.stderr.flush()
    l, s0, s1, s2 = exp_data.get_training_set(test, exp_data, options.debug, options.name, 10000, 1)
    
    sys.stderr.flush()
    #print(training_set_obs[training_set_obs>0])
    if options.distr == "negbin":
        from rgt.THOR.neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters
        print('Training HMM...', file=sys.stderr)
        sys.stderr.flush()
        init_alpha, init_mu = get_init_parameters(s0, s1, s2)
        m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
        training_set_obs = exp_data.get_observation(sample(range(exp_data.overall_coverage[0].shape[1]), l))
        m.fit([training_set_obs])
        for i in range(5):
            print(i, file=sys.stderr)
            init_alpha, init_mu = get_init_parameters(s0, s1, s2)
            m_new = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
             
            training_set_obs = exp_data.get_observation(sample(range(exp_data.overall_coverage[0].shape[1]), l))
            m_new.fit([training_set_obs])
            if m_new.em_prob > m.em_prob:
                print('take HMM no. %s, new: %s, old: %s ' %(i, m_new.em_prob, m.em_prob), file=sys.stderr)
                m = m_new
                 
        tracker.write(text=init_mu, header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu,alpha)")
        tracker.write(text=init_alpha)
        tracker.write(text=m.mu, header="Final HMM's Neg. Bin. Emission distribution (mu,alpha)")
        tracker.write(text=m.alpha)
         
    tracker.write(text=m._get_transmat(), header="Transmission matrix")
    
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
     
    if options.distr == 'negbin':
        distr = _get_pvalue_distr(exp_data, m.mu, m.alpha, tracker)
    else:
        distr = {'distr_name': 'binomial', 'n': m.n[0], 'p': m.p[0][1]}
     
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr, merge=options.merge, exts=exp_data.exts, pcutoff=options.pcutoff, debug=options.debug, p=options.par)
    
if __name__ == '__main__':
    main() 
