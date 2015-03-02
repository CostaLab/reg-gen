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
import numpy as np
from rgt.THOR.tracker import Tracker
from dpc_help import get_bibtex_entry

from random import sample

def write(name, l):
    f = open(name, 'w')
    
    for el in l:
        print(el, file=f)
    f.close()

def _get_training_sets(indices_of_interest, first_overall_coverage, second_overall_coverage, name, debug, x = 10000, threshold = 2.0, diff_cov = 10):
    """Return s0,s1,s2 (list of tuples) with initial peaks"""
    s0 = []
    s1 = []
    s2 = []
    so = []
    i = 1
    while i < x:
        state = None
        #print(indices_of_interest, file=sys.stderr)
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
    
    if len(s1) == 0:
        s1 = map(lambda x: (x[1], x[0]), s2)
    if len(s2) == 0:
        s2 = map(lambda x: (x[1], x[0]), s1)
    
    if debug:
        write(name + '-s0', s0)
        write(name + '-s1', s1)
        write(name + '-s2', s2)
        write(name + '-soverall', so)
    return s0, s1, s2


def main():
    test = False
    if test:
        print("---------- TEST MODE ----------", file=sys.stderr)
    
    options, bamfile_1, bamfile_2, genome, chrom_sizes = input(test)
    
    ######### WORK! ##########
    tracker = Tracker(options.name + '-setup.info')
    
    if options.no_correction:
        tracker.write(text=options.pcutoff, header="p-value cutoff, no p-value correction")
    else:
        tracker.write(text=options.pcutoff, header="p-value cutoff, with p-value correction (Benjamini/Hochberg)")
    
    exp_data, ext_sizes = initialize(name=options.name, genome_path=genome, regions=options.regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bam_file_1 = bamfile_1, ext_1=options.ext_1, bam_file_2 = bamfile_2, ext_2=options.ext_2, \
                          input_1=options.input_1, input_factor_1=options.input_factor_1, ext_input_1=options.ext_input_1, \
                          input_2=options.input_2, input_factor_2=options.input_factor_2, ext_input_2=options.ext_input_2, \
                          chrom_sizes = chrom_sizes,
                          verbose = options.verbose, norm_strategy=options.norm_strategy, no_gc_content=options.no_gc_content, deadzones=options.deadzones, \
                          factor_input_1 = options.factor_input_1, factor_input_2 = options.factor_input_2, debug=options.debug, tracker=tracker)
    exp_data.compute_putative_region_index()
    print('Number of regions with putative differential peaks:', len(exp_data.indices_of_interest), file=sys.stderr)
    
    if options.debug:
        exp_data.write_putative_regions(options.name + '-putative-peaks.bed')
    print('Compute training set...',file=sys.stderr)
    training_set = exp_data.get_training_set(exp_data, min(len(exp_data.indices_of_interest) / 3, 600000), options.verbose, options.name, options.debug, options.constchrom)
    training_set_obs = exp_data.get_observation(training_set)
    
    _, s1, s2 = _get_training_sets(exp_data.indices_of_interest, exp_data.first_overall_coverage, exp_data.second_overall_coverage, options.name, options.debug)
    
    #choose proper HMM
    print('Training HMM...', file=sys.stderr)
    if options.distr == 'binom':
        #for binomial distribution
        print("Use binomial HMM", file=sys.stderr)
        from hmm_binom_2d3s import BinomialHMM2d3s, get_init_parameters
        tmp = sum( [ exp_data.first_overall_coverage[i] + exp_data.second_overall_coverage[i] for i in exp_data.indices_of_interest]) / 2
        n_, p_ = get_init_parameters(s1, s2, count=tmp)
        m = BinomialHMM2d3s(n_components=3, n=n_, p=p_)
        m.fit([training_set_obs])
        m.save_setup(tracker)
        distr_pvalue={'distr_name': "binomial", 'n': m.n[0], 'p': m.p[0][1]}
    elif options.distr == 'poisson':
        print("Use poisson mixture HMM", file=sys.stderr)
        from hmm_mixture_poisson_2d3s import PoissonHMM2d3s, get_init_parameters
        distr_magnitude = options.mag
        n_components = 3
        n_features = 2
        initial_c, initial_p = get_init_parameters(s1, s2, distr_magnitude=distr_magnitude, n_components=n_components, n_features=n_features)
        
        f = map(lambda x: x+1, range(distr_magnitude))
        #g = map(lambda x: x+1, range(distr_magnitude))
        #f = map(lambda x: x/(float(distr_magnitude)), g)
        m = PoissonHMM2d3s(c=initial_c, distr_magnitude=distr_magnitude, factors=f, p=initial_p)
        
        m.fit([training_set_obs])
        n = sum( [ exp_data.first_overall_coverage[i] + exp_data.second_overall_coverage[i] for i in exp_data.indices_of_interest]) / 2
        mean = np.mean([m.get_mean(0,0), m.get_mean(0,1)])
        p = mean / float(n)
        m.save_setup(tracker, n, p)
        distr_pvalue={'distr_name': "binomial", 'n': n, 'p': p}
    elif options.distr == 'poisson-c':
        print("Use poisson constrained mixture HMM", file=sys.stderr)
        from hmm_mixture_constpoisson_2d3s import PoissonHMM2d3s, get_init_parameters
        distr_magnitude = options.mag
        n_components = 3
        n_features = 2
        initial_c, initial_p = get_init_parameters(s1, s2, distr_magnitude=distr_magnitude, n_components=n_components, n_features=n_features)
        
        f = map(lambda x: x+1, range(distr_magnitude))
        m = PoissonHMM2d3s(c=initial_c, distr_magnitude=distr_magnitude, factors=f, p=initial_p)
        
        m.fit([training_set_obs])
        n = sum( [ exp_data.first_overall_coverage[i] + exp_data.second_overall_coverage[i] for i in exp_data.indices_of_interest]) / 2
        mean = np.mean([m.get_mean(0,0), m.get_mean(0,1)])
        p = mean / float(n)
        m.save_setup(tracker, n, p)
        distr_pvalue={'distr_name': "binomial", 'n': n, 'p': p}
        
    if options.debug:
        exp_data.get_initial_dist(options.name + '-initial-states.bed')
      
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
     
    if options.debug: 
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
 
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr_pvalue, ext_size=np.mean(ext_sizes), merge=options.merge, pcutoff=options.pcutoff, no_correction=options.no_correction)


if __name__ == '__main__':
    main()