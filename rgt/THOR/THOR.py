#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
THOR detects differential peaks in multiple ChIP-seq profiles associated
with two distinct biological conditions.

Copyright (C) 2014-2016 Manuel Allhoff (allhoff@aices.rwth-aachen.de)

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

@author: Manuel Allhoff
"""

# Python
from __future__ import print_function
import sys

# Internal
from dpc_help import get_peaks, initialize, merge_output, handle_input
from help_hmm import fit_sm_mean_var_distr
from tracker import Tracker
from postprocessing import _output_BED, _output_narrowPeak
from rgt.THOR.neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters, _get_pvalue_distr
from rgt.THOR.RegionGiver import RegionGiver
from rgt.THOR.postprocessing import filter_by_pvalue_strand_lag
from rgt import __version__

import configuration
from get_statistics import get_file_statistics, compute_extension_sizes, update_statics_extension_sizes, is_stats_valid
from MultiCoverageSet import get_training_set, transform_data_for_HMM
# External


TEST = False #enable to test THOR locally


def _write_info(tracker, report, **data):
    """Write information to tracker"""
    tracker.write(text=data['func_para'][0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (a)")
    tracker.write(text=data['func_para'][1], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (c)")
    #tracker.write(text=data['init_mu'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu)")
    #tracker.write(text=data['init_alpha'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (alpha)")
    #tracker.write(text=data['m'].mu, header="Final HMM's Neg. Bin. Emission distribution (mu)")
    #tracker.write(text=data['m'].alpha, header="Final HMM's Neg. Bin. Emission distribution (alpha)")
    #tracker.write(text=data['m']._get_transmat(), header="Transmission matrix")
    # print(configuration.FOLDER_REPORT)
    if report:
        tracker.make_html(configuration.FOLDER_REPORT)


def train_HMM(region_giver, options, signal_statics, inputs_statics, genome, tracker):
    """Train HMM
    Firstly, we get training set: exp_data, and then estimate parameter
    get training data: we need to consider all file data;
        -- get estimate extension size for all files w.r.t to whole data: bt on the other hand, we see it's an attribute for bamfile
        But actually we could subtract it from train_HMM, because we need it also other parts
        -- choose samples w.r.t to sample distribution to get exp_data, but according to data, we could choose different len of samples
    estimate parameters: all data included

    Return: distribution parameters
    """
    # stats_total, stats_data, isspatial = get_read_statistics(signal_files[i][j], chrom_fname)

    exp_data = initialize(options=options,strand_cov=True, genome_path=genome, regionset=region_giver.valid_regionset, mask_file=region_giver.mask_file,
                          signal_statics=signal_statics, inputs_statics=inputs_statics)
    # test if exp_data satisfies conditions to be good enough
    #if exp_data.count_positive_signal() > len(train_regions.sequences[0]) * 0.00001:
    #    tracker.write(text=" ".join(map(lambda x: str(x), exp_data.exts)), header="Extension size (rep1, rep2, input1, input2)")
    #    tracker.write(text=map(lambda x: str(x), exp_data.scaling_factors_ip), header="Scaling factors")
        # break
    ## here we need to use tracker to record
    ## per bins mean, var
    func, func_para = fit_sm_mean_var_distr(exp_data.overall_coverage, options.name, options.debug,
                                          verbose=options.verbose, outputdir=options.outputdir,
                                          report=options.report, poisson=options.poisson)
    exp_data.compute_sm_putative_region_index()
     
    print('Compute HMM\'s training set', file=sys.stderr)

    training_data, s0, s1, s2 = get_training_set(exp_data, True, options.name, options.foldchange,
                                                       options.threshold, options.size_ts, 0)
    init_alpha, init_mu = get_init_parameters(s0, s1, s2)
    m = NegBinRepHMM(alpha=init_alpha, mu=init_mu, dim_cond_1=signal_statics['dim'][0], dim_cond_2=signal_statics['dim'][1], func=func)

    print('Train HMM', file=sys.stderr)
    m.fit([training_data], options.hmm_free_para)
    distr = _get_pvalue_distr(m.mu, m.alpha, tracker)
         
    return m, exp_data, func_para, init_mu, init_alpha, distr


def run_HMM(region_giver, options, signal_statics, inputs_statics, genome, tracker, m, distr):
    """Run trained HMM chromosome-wise on genomic signal and call differential peaks"""
    output, pvalues, ratios, no_bw_files = [], [], [], []
    print("Compute HMM's posterior probabilities and Viterbi path to call differential peaks", file=sys.stderr)
    
    for i, regionset in enumerate(region_giver):
        chrom = regionset.sequences[0].chrom
        ## before we test if data is big enough
        if not is_stats_valid(signal_statics,chrom):
            print("- not taking into account %s, since too less data" % chrom, file=sys.stderr)
            continue

        print("- taking into account %s" %chrom , file=sys.stderr)
        exp_data = initialize(options=options, strand_cov=True, genome_path=genome, regionset=regionset, mask_file=region_giver.mask_file,
                              signal_statics=signal_statics, inputs_statics=inputs_statics)
        if not exp_data.data_valid:
            print('less data and then consider next chromasome',file=sys.stderr)
            continue

        ## After this step, we have already normalized data, so we could output normalization data
        if options.save_input:
            exp_data.output_input_bw(options.name + '-' + str(i), region_giver.chrom_sizes_file, options.save_wig)
        exp_data.output_signal_bw(options.name + '-' + str(i), region_giver.chrom_sizes_file, options.save_wig)

        no_bw_files.append(i)
        exp_data.compute_sm_putative_region_index()

        if exp_data.indices_of_interest is None:
            continue
        cov_data = exp_data.get_sm_covs(exp_data.indices_of_interest)
        states = m.predict(transform_data_for_HMM(cov_data))
        
        inst_ratios, inst_pvalues, inst_output = get_peaks(name=options.name, states=states, cov_set=exp_data,
                                                           distr=distr, merge=options.merge, exts=options.exts,
                                                           pcutoff=options.pcutoff, debug=options.debug, p=options.par,
                                                           no_correction=options.no_correction,
                                                           merge_bin=options.merge_bin, deadzones=options.deadzones)

        output += inst_output
        pvalues += inst_pvalues
        ratios += inst_ratios
    
    res_output, res_pvalues, res_filter_pass = filter_by_pvalue_strand_lag(ratios, options.pcutoff, pvalues, output,
                                                                           options.no_correction, options.name,
                                                                           options.singlestrand)
    
    _output_BED(options.name, res_output, res_pvalues, res_filter_pass)
    _output_narrowPeak(options.name, res_output, res_pvalues, res_filter_pass)
    
    merge_output(signal_statics,  options, no_bw_files, region_giver.chrom_sizes_file)


def main():
    options, bamfiles, genome, chrom_sizes_file, dims, inputs_files = handle_input()

    tracker = Tracker(options.name + '-setup.info', bamfiles, genome, chrom_sizes_file, dims, inputs_files, options, __version__)
    # get statistical data for all files and compare them to create one big files..
    region_giver = RegionGiver(chrom_sizes_file, options.regions)
    # get statistic information for each file
    # stats_total, stats_data, read_size, and other parts..

    signal_statics = get_file_statistics(bamfiles, region_giver)
    region_giver.update_regions(signal_statics)

    inputs_statics = None
    # but how about input files, if we want extension size, then they are connected..But we could extract them outside
    # inputs_statics = get_file_statistics(inputs_files, region_giver)
    # region_giver.update_regions(inputs_statics)
    # compute extension size if option.ext are not given
    # for testing..
    options.exts = [72,96,114,204]
    # options.inputs_exts = [228,228,231,231]
    if options.exts:
        update_statics_extension_sizes(signal_statics, options.exts)
    else:
        options.exts, _ = compute_extension_sizes(signal_statics)

    if options.exts_inputs:
        # less one step how to get inputs size and then give it
        update_statics_extension_sizes(inputs_statics,options.exts_inputs)
    else:
        options.exts_inputs, _ = compute_extension_sizes(inputs_statics)
    # one function to transform these parameters, after we read and do it into callback function??
    # options.factors_inputs = [[0.692, 0.719], [0.726,0.708]]
    # options.scaling_factors_ip = [[1.0537533317434074, 0.9986756833314534], [0.99775712024261831, 1.1606581998669872]]
    # pass stats_total, stats_data, extension sizes to train_HMM
    m, exp_data, func_para, init_mu, init_alpha, distr = train_HMM(region_giver, options, signal_statics, inputs_statics, genome,  tracker)

    # we need to change region_giver and update the valid_regions
    region_giver.reset_regions()
    run_HMM(region_giver, options, signal_statics, inputs_statics, genome,  tracker, m, distr)
    
    _write_info(tracker, options.report, func_para=func_para, init_mu=init_mu, init_alpha=init_alpha, m=m)


if __name__ == "__main__":
    main()