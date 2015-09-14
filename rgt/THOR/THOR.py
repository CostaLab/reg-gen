#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
from dpc_help import get_peaks, input, _fit_mean_var_distr, initialize
from tracker import Tracker
from rgt.THOR.neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters, _get_pvalue_distr
TEST = False #enable to test THOR locally

def _write_info(tracker, report, **data):
    """Write information to tracker"""
    tracker.write(text=data['func_para'][0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (a)")
    tracker.write(text=data['func_para'][1], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (c)")
    tracker.write(text=data['init_mu'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu)")
    tracker.write(text=data['init_alpha'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (alpha)")
    tracker.write(text=data['m'].mu, header="Final HMM's Neg. Bin. Emission distribution (mu)")
    tracker.write(text=data['m'].alpha, header="Final HMM's Neg. Bin. Emission distribution (alpha)")
    tracker.write(text=data['m']._get_transmat(), header="Transmission matrix")
    
    if report:
        tracker.make_html()
    
def main():
    options, bamfiles, genome, chrom_sizes, dims, inputs = input(TEST)
    tracker = Tracker(options.name + '-setup.info')
    #print(FOLDER_REPORT, file=sys.stderr)
    exp_data = initialize(name=options.name, dims=dims, genome_path=genome, regions=options.regions, stepsize=options.stepsize, binsize=options.binsize, \
                          bamfiles = bamfiles, exts=options.exts, inputs=inputs, exts_inputs=options.exts_inputs, debug=options.debug,\
                          verbose = options.verbose, no_gc_content=options.no_gc_content, factors_inputs=options.factors_inputs, chrom_sizes=chrom_sizes, \
                          tracker=tracker, norm_regions=options.norm_regions, scaling_factors_ip = options.scaling_factors_ip, save_wig=options.save_wig, \
                          housekeeping_genes=options.housekeeping_genes, test=TEST, report=options.report)
    func, func_para = _fit_mean_var_distr(exp_data.overall_coverage, options.name, options.debug, verbose=options.verbose, \
                                          outputdir = options.outputdir, report=options.report)
    exp_data.compute_putative_region_index()

    print('Compute HMM\'s training set', file=sys.stderr)
    training_set, s0, s1, s2 = exp_data.get_training_set(TEST, exp_data, options.name, options.foldchange, options.threshold, options.size_ts, 3)
    init_alpha, init_mu = get_init_parameters(s0, s1, s2)
    m = NegBinRepHMM(alpha = init_alpha, mu = init_mu, dim_cond_1 = dims[0], dim_cond_2 = dims[1], func = func)
    training_set_obs = exp_data.get_observation(training_set)
    
    print('Train HMM', file=sys.stderr)
    m.fit([training_set_obs], options.hmm_free_para)
    
    print("Compute HMM's posterior probabilities and Viterbi path to call differential peaks", file=sys.stderr)
    states = m.predict(exp_data.get_observation(exp_data.indices_of_interest))
    
    distr = _get_pvalue_distr(exp_data, m.mu, m.alpha, tracker)
    get_peaks(name=options.name, states=states, DCS=exp_data, distr=distr, merge=options.merge, \
              exts=exp_data.exts, pcutoff=options.pcutoff, debug=options.debug, p=options.par,\
              no_correction=options.no_correction, deadzones=options.deadzones)
    
    _write_info(tracker, options.report, func_para=func_para, init_mu=init_mu, init_alpha=init_alpha, m=m)
    