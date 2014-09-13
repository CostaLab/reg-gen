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

from random import sample

def write(name, l):
    f = open(name, 'w')
    
    for el in l:
        print(el, file=f)
    f.close()

def _get_training_sets(indices_of_interest, first_overall_coverage, second_overall_coverage, name, verbose, x = 10000, threshold = 2.0, diff_cov = 10):
    """Return s0,s1,s2 (list of tuples) with initial peaks"""
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
        write(name + '-s0', s0)
        write(name + '-s1', s1)
        write(name + '-s2', s2)
        write(name + '-soverall', so)
    return s0, s1, s2


def main():
    test = False
    options, bamfile_1, bamfile_2, genome, chrom_sizes = input(test)
    #print(options.verbose, file=sys.stderr)
    ######### WORK! ##########
    exp_data, ext_sizes = initialize(name=options.name, genome_path=genome, regions=options.regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bam_file_1 = bamfile_1, ext_1=options.ext_1, \
                          bam_file_2 = bamfile_2, ext_2=options.ext_2, \
                          input_1=options.input_1, input_factor_1=options.input_factor_1, ext_input_1=options.ext_input_1, \
                          input_2=options.input_2, input_factor_2=options.input_factor_2, ext_input_2=options.ext_input_2, \
                          chrom_sizes = chrom_sizes,
                          verbose = options.verbose, norm_strategy=options.norm_strategy, no_gc_content=options.no_gc_content, deadzones=options.deadzones, \
                          factor_input_1 = options.factor_input_1, factor_input_2 = options.factor_input_2)
    print('done', file=sys.stderr)
    print('Number of regions to be considered by the HMM:', len(exp_data), file=sys.stderr)
    exp_data.compute_putative_region_index()
    print('Number of regions with putative differential peaks:', len(exp_data.indices_of_interest), file=sys.stderr)
    
    if options.verbose:
        exp_data.write_putative_regions(options.name + '-putative-peaks.bed')
    print('Compute training set...',file=sys.stderr)
    training_set = exp_data.get_training_set(exp_data, min(len(exp_data.indices_of_interest) / 3, 600000), options.verbose, options.name)
    training_set_obs = exp_data.get_observation(training_set)
    
    _, s1, s2 = _get_training_sets(exp_data.indices_of_interest, exp_data.first_overall_coverage, exp_data.second_overall_coverage, options.name, options.verbose)
    
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
        #if options.verbose:
        f = open(options.name + '-setup.info', 'w')
        f.write('Binomial n, p\n')
        f.write("%s %s\n" %(m.n, m.p))
        f.write('p-value settings\n')
        f.write("%s %s\n" %(m.n[0], m.p[0][1]))
        f.close()
        distr_pvalue={'distr_name': "binomial", 'n':m.n[0], 'p':m.p[0][1]}
    elif options.distr == 'poisson':
        print("Use poisson mixture HMM", file=sys.stderr)
        from hmm_mixture_constpoisson_2d3s import PoissonHMM2d3s, get_init_parameters
        distr_magnitude = options.mag
        n_components = 3
        n_features = 2
        initial_c, initial_p = get_init_parameters(s1, s2, distr_magnitude=distr_magnitude, n_components=n_components, n_features=n_features)
        m = PoissonHMM2d3s(c=initial_c, distr_magnitude=distr_magnitude, factors=map(lambda x: x+1, range(distr_magnitude)), p=initial_p)
        
        m.fit([training_set_obs])
        n = sum( [ exp_data.first_overall_coverage[i] + exp_data.second_overall_coverage[i] for i in exp_data.indices_of_interest]) / 2
        mean = np.mean([m.get_mean(0,0), m.get_mean(0,1)])
        p = mean / float(n)
        #if options.verbose:
        f = open(options.name + '-setup.info', 'w')
        f.write('Poisson P\n')
        f.write("%s \n" %m.p)
        f.write('Poisson C\n')
        f.write("%s \n" %m.c)
        f.write('Poisson p-value settings\n')
        f.write("%s %s\n" %(n, p))
        f.close()
        distr_pvalue={'distr_name': "binomial", 'n': n, 'p': p}
        
    if options.verbose:
        exp_data.get_initial_dist(options.name + '-initial-states.bed')
      
    print('...done', file=sys.stderr)
    if options.verbose:
        print('p', m.p, file=sys.stderr)
        print("Final HMM's transition matrix: ", file=sys.stderr)
        print(m._get_transmat(), file=sys.stderr)
         
    print("Computing HMM's posterior probabilities and Viterbi path", file=sys.stderr)
    posteriors = m.predict_proba(exp_data.get_observation(exp_data.indices_of_interest))
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    print("...done", file=sys.stderr)
     
    if options.verbose: 
        dump_posteriors_and_viterbi(name=options.name, posteriors=posteriors, states=states, DCS=exp_data)
 
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr_pvalue, ext_size=np.mean(ext_sizes), merge=options.merge)


if __name__ == '__main__':
    main()