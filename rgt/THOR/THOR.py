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
from dpc_help import get_back
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
                          housekeeping_genes=options.housekeeping_genes)
    
    func, func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, options.debug, sample_size=20000, verbose=options.debug)
    tracker.write(text=func_para[0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) ")
    tracker.write(text=func_para[1])
    
    exp_data.compute_putative_region_index()

    print('Compute training set...',file=sys.stderr)
    
    if options.distr == "binom":
        cov0, cov1 = [], []
        for i in range(exp_data.overall_coverage[0].shape[1]):
            if i % 100000 == 0:
                print(i, exp_data.overall_coverage[0].shape[1])
            #if i > 80000:
            #    break
            cov0.append(np.sum(exp_data.overall_coverage[0][:,i])) #np.sum(b, axis=1) ???
            cov1.append(np.sum(exp_data.overall_coverage[1][:,i]))

        exp_data.dim_1, exp_data.dim_2 = 1, 1
        exp_data.overall_coverage[0] = np.matrix(cov0)
        exp_data.overall_coverage[1] = np.matrix(cov1)

    training_set, s0, s1, s2 = exp_data.get_training_set(test, exp_data, options.debug, options.name, 10000, 3)
    training_set_obs = exp_data.get_observation(training_set)
    
    if options.distr == "negbin":
        from rgt.THOR.neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters
        init_alpha, init_mu = get_init_parameters(s0, s1, s2)
        tracker.write(text=init_mu, header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu,alpha)")
        tracker.write(text=init_alpha)
        m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
        print('Training HMM...', file=sys.stderr)
        m.fit([training_set_obs])
        tracker.write(text=m.mu, header="Final HMM's Neg. Bin. Emission distribution (mu,alpha)")
        tracker.write(text=m.alpha)
    elif options.distr == "binom":
        print('use binom distr', file=sys.stderr)
        #from rgt.ODIN.binom_hmm_2d_3s import BinomialHMM2d3s, get_init_parameters
        from rgt.ODIN.hmm_binom_2d3s import BinomialHMM2d3s, get_init_parameters
        tmp = 0
        for i in range(len(exp_data.indices_of_interest)):
            c1, c2 = exp_data._get_covs(exp_data, i)
            tmp += sum([c1, c2])
        n_, p_ = get_init_parameters(s1, s2, count=tmp)
        print(n_, p_, file=sys.stderr)
        #tracker.write(text=n_, header="Inital parameter estimate for HMM's Bin. Emission distribution (n, p)")
        #tracker.write(text=np.asarray(p_))
        #m = BinomialHMM(n_components=3, p = p_, startprob=[1,0,0], n = n_, dim_cond_1=dims[0], dim_cond_2=dims[1])
        #print('Training HMM...', file=sys.stderr)
        #m.fit([training_set_obs])
        #tracker.write(text=m.n, header="Final HMM's Neg. Bin. Emission distribution (mu,alpha)")
        #tracker.write(text=m.p)
        #tmp = sum( [ exp_data.first_overall_coverage[i] + exp_data.second_overall_coverage[i] for i in exp_data.indices_of_interest]) / 2
        #n_, p_ = get_init_parameters(s1, s2, count=tmp)
        m = BinomialHMM2d3s(n_components=3, n=n_, p=p_)
        m.fit([training_set_obs])
        #m.save_setup(tracker)
        distr_pvalue={'distr_name': "binomial", 'n': m.n[0], 'p': m.p[0][1]}
        
    tracker.write(text=m._get_transmat(), header="Transmission matrix")
    
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    
    if options.debug:
        posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
    
    if options.distr == 'negbin':
        distr = _get_pvalue_distr(exp_data, m.mu, m.alpha, tracker)
    else:
        distr = {'distr_name': 'binomial', 'n': m.n[0], 'p': m.p[0][1]}
    
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr, merge=options.merge, exts=exp_data.exts, pcutoff=options.pcutoff, debug=options.debug, p=options.par)
    
if __name__ == '__main__':
    main() 
    
