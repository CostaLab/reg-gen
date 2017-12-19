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

# Internal
import configuration
from input_parser import handle_input
from peak_calling import get_peaks, initialize
from help_hmm import fit_sm_mean_var_distr
from tracker import Tracker
from postprocessing import output_BED, output_narrowPeak, merge_bw_output
from rgt.THOR.neg_bin_rep_hmm import NegBinRepHMM, get_init_parameters, _get_pvalue_distr
from rgt.THOR.RegionGiver import RegionGiver
from rgt.THOR.postprocessing import filter_by_pvalue_strand_lag
from rgt import __version__

from get_statistics import *
from MultiCoverageSet import get_training_set, transform_data_for_HMM
# External


TEST = False #enable to test THOR locally


def _write_info(tracker, report, **data):
    """Write information to tracker"""
    tracker.write(text=data['func_para'][0], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (a)")
    tracker.write(text=data['func_para'][1], header="Parameters for both estimated quadr. function y=max(|a|*x^2 + x + |c|, 0) (c)")
    tracker.write(text=data['init_mu'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (mu)")
    tracker.write(text=data['init_alpha'], header="Inital parameter estimate for HMM's Neg. Bin. Emission distribution (alpha)")
    tracker.write(text=data['m'].mu, header="Final HMM's Neg. Bin. Emission distribution (mu)")
    tracker.write(text=data['m'].alpha, header="Final HMM's Neg. Bin. Emission distribution (alpha)")
    # tracker.write(text=data['m']._get_transmat(), header="Transmission matrix")
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
    exp_data = initialize(options=options,strand_cov=True, genome_path=genome, regionset=region_giver.valid_regionset, mask_file=region_giver.mask_file,
                          signal_statics=signal_statics, inputs_statics=inputs_statics, verbose=True)

    tracker.write(text=map(lambda x: str(x), options.scaling_factors_ip), header="Scaling factors")

    data_valid = exp_data.compute_sm_putative_region_index()
    if not data_valid: # less than 50 peaks we consider
        print('putative region is not enough, less than 50, so no peak calling', file=sys.stderr)
        sys.exit()

    func, func_para = fit_sm_mean_var_distr(exp_data.overall_coverage, options.name, options.debug,
                                            verbose=options.verbose, outputdir=options.outputdir,
                                            report=True, poisson=options.poisson)

    print('Compute HMM\'s training set', file=sys.stderr)
    training_data, s0, s1, s2 = get_training_set(exp_data, options.name, options.foldchange,
                                                       options.threshold, options.size_ts, 3, test=False)
    init_alpha, init_mu = get_init_parameters(s0, s1, s2, report=True)

    m = NegBinRepHMM(alpha=init_alpha, mu=init_mu, dim_cond_1=signal_statics['dim'][0], dim_cond_2=signal_statics['dim'][1], func=func)

    print('Train HMM', file=sys.stderr)
    m.fit([training_data], options.hmm_free_para)
    distr = _get_pvalue_distr(m.mu, m.alpha, tracker)
    return m, func_para, init_mu, init_alpha, distr


def run_HMM(region_giver, options, signal_statics, inputs_statics, genome, tracker, m, distr):
    """Run trained HMM chromosome-wise on genomic signal and call differential peaks"""
    no_bw_files = []
    output, pvalues, ratios = [], [], []
    # pcutoff_output, pcutoff_pvalues, pcutoff_ratios = {}, {}, {}
    print("Compute HMM's posterior probabilities and Viterbi path to call differential peaks", file=sys.stderr)
    
    for i, regionset in enumerate(region_giver):
        chrom = regionset.sequences[0].chrom
        ## before we test if data is big enough
        if not is_stats_valid(signal_statics,chrom):
            print("- not taking into account %s, since too less data" % chrom, file=sys.stderr)
            continue
        else:
            print("- taking into account %s" %chrom , file=sys.stderr)
        exp_data = initialize(options=options, strand_cov=True, genome_path=genome, regionset=regionset, mask_file=region_giver.mask_file,
                              signal_statics=signal_statics, inputs_statics=inputs_statics)

        options.save_bw = False
        ## After this step, we have already normalized data, so we could output normalization data
        if options.save_input:
            print("Begin: output nomalized inputs read data into file", file=sys.stderr)
            exp_data.output_input_bw(options.name + '-' + str(i), region_giver.chrom_sizes_file)
            print("End: output nomalized inputs read data into file", file=sys.stderr)
        if options.save_bw:
            print("Begin : output nomalized signal read data into file", file=sys.stderr)
            exp_data.output_signal_bw(options.name + '-' + str(i), region_giver.chrom_sizes_file)
            print("End: output nomalized signal read data into file", file=sys.stderr)
        no_bw_files.append(i)
        # if we accept command to stop here, then we don't call diff-peaks, but only output normalized files
        if not options.call_peaks:
            continue
        # after training we don't need to verify the number of putative_region
        exp_data.compute_sm_putative_region_index()

        cov_data = exp_data.get_sm_covs(exp_data.indices_of_interest)
        states = m.predict(transform_data_for_HMM(cov_data))
        inst_ratios, inst_pvalues, inst_output = get_peaks(states=states, cov_set=exp_data,
                                                           distr=distr, merge=options.merge, exts=options.exts,
                                                           pcutoff=options.pcutoff, debug=options.debug, p=options.par,
                                                           no_correction=options.no_correction,
                                                           merge_bin=options.merge_bin)
        """
        for pi_value in inst_output.keys():
            if pi_value not in pcutoff_output.keys():
                pcutoff_output[pi_value] = []
                pcutoff_pvalues[pi_value] = []
                pcutoff_ratios[pi_value] = []

            pcutoff_output[pi_value] += inst_output[pi_value]
            pcutoff_pvalues[pi_value] += inst_pvalues[pi_value]
            pcutoff_ratios[pi_value] += inst_ratios[pi_value]
        """
        output += inst_output
        pvalues += inst_pvalues
        ratios += inst_ratios

    if options.save_bw:
        merge_bw_output(signal_statics,  options, no_bw_files, region_giver.chrom_sizes_file)

    if options.call_peaks: # if we don't have call_peaks, we only output nomalized data
        res_output, res_pvalues, res_filter_pass = filter_by_pvalue_strand_lag(ratios, options.pcutoff, pvalues, output,
                                                                               options.no_correction, options.name,
                                                                               options.singlestrand)

        output_BED(options.name, res_output, res_pvalues, res_filter_pass)
        output_narrowPeak(options.name, res_output, res_pvalues, res_filter_pass)


def main():
    options, signal_files, genome, chrom_sizes_file, dims, inputs_files = handle_input()

    tracker = Tracker(options.name + '-setup.info', signal_files, genome, chrom_sizes_file, dims, inputs_files, options, __version__)
    region_giver = RegionGiver(chrom_sizes_file, options.regions)

    # get statistic information for each file
    # stats_total, stats_data, read_size, and other parts..
    signal_statics = get_file_statistics(signal_files, region_giver)
    region_giver.update_regions(signal_statics)

    if options.exts:
        update_statics_extension_sizes(signal_statics, options.exts)
    else:
        options.exts = compute_extension_sizes(signal_statics, options.report)
    tracker.write(text=" ".join(map(lambda x: str(x), options.exts)),
                  header="Extension size for signal files (rep1, rep2, input1, input2)")

    # compute extension size for signal files if inputs_files exist
    if inputs_files:
        inputs_statics = get_file_statistics(inputs_files, region_giver)  # None
        region_giver.update_regions(inputs_statics)
        if options.exts_inputs:
            update_statics_extension_sizes(inputs_statics,options.exts_inputs)
        else:
            options.exts_inputs= compute_extension_sizes(inputs_statics, options.report)
            # Foe given data we don't need to adjust data;; Only what we have we need to do
            options.exts = adjust_extension_sizes(inputs_statics, signal_statics)
            options.exts_inputs = adjust_extension_sizes(signal_statics, inputs_statics)
        tracker.write(text=" ".join(map(lambda x: str(x), options.exts_inputs)),
                  header="Extension size for inputs files (rep1, rep2, input1, input2)")
    else:
        inputs_statics = None

    m, func_para, init_mu, init_alpha, distr = train_HMM(region_giver, options, signal_statics, inputs_statics, genome, tracker)

    # we need to change region_giver and update the valid_regions
    region_giver.reset_regions()
    run_HMM(region_giver, options, signal_statics, inputs_statics, genome,  tracker, m, distr)
    
    _write_info(tracker, options.report, func_para=func_para, init_mu=init_mu, init_alpha=init_alpha, m=m)


if __name__ == "__main__":
    main()