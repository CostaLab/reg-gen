###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import numpy as np
from collections import Counter

from sklearn.externals import joblib
from hmmlearn import hmm

# Internal
from ..Util import GenomeData
from signalProcessing import GenomicSignal
from rgt.GenomicRegionSet import GenomicRegionSet
from hmm import HMM
from biasTable import BiasTable

"""
Train a hidden Markov model (HMM) based on the annotation data

Authors: Eduardo G. Gusmao, Zhijian Li
"""


class TrainHMM:
    """
    Contains methods used to train a hidden Markov model
    """

    def __init__(self, reads_file, annotate_file, print_bed_file, output_location,
                 output_prefix, bias_table, organism, k_nb, chrom, start, end,
                 downstream_ext, upstream_ext, forward_shift, reverse_shift):
        self.reads_file = reads_file
        self.annotate_file = annotate_file
        self.print_bed_file = print_bed_file
        self.output_location = output_location
        self.output_prefix = output_prefix
        self.bias_table = bias_table
        self.organism = organism
        self.k_nb = k_nb
        self.chrom = chrom
        self.start = start
        self.end = end
        self.downstream_ext = downstream_ext
        self.upstream_ext = upstream_ext
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift

    def read_states_signals(self):
        # Read states from the annotation file
        states = ""
        with open(self.annotate_file) as f:
            for line in f:
                if len(line) < 2 or "#" in line or "=" in line:
                    continue
                ll = line.strip().split(" ")
                for state in ll[1:-1]:
                    states += state

        # If need to estimate bias table
        genome_data = GenomeData(self.organism)
        table = None

        # If the bias table is provided
        if self.bias_table:
            bias_table = BiasTable()
            bias_table_list = self.bias_table.split(",")
            table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                          table_file_name_R=bias_table_list[1])

        # Get the normalization and slope signal from the raw bam file
        raw_signal = GenomicSignal(self.reads_file)
        raw_signal.load_sg_coefs(slope_window_size=9)
        norm_signal, slope_signal = \
            raw_signal.get_signal2(self.chrom, self.start, self.end,
                                   self.downstream_ext, self.upstream_ext,
                                   self.forward_shift, self.reverse_shift,
                                   bias_table=table, genome_file_name=genome_data.get_genome())
        if self.print_bed_file:
            self.output_bed_file(states)

        return states, norm_signal, slope_signal

    def train(self):
        # Estimate the HMM parameters using the maximum likelihood method.
        states, norm_signal, slope_signal = self.read_states_signals()
        hmm_model = HMM()

        hmm_model.dim = 2
        # States number
        state_list = [int(state) for state in list(states)]
        hmm_model.states = len(np.unique(np.array(state_list)))

        # Initial state probabilities vector
        init_state = state_list[0]
        hmm_model.pi = [0.0] * hmm_model.states
        hmm_model.pi[init_state] = 1.0

        # Transition
        trans_matrix = np.zeros((hmm_model.states, hmm_model.states))
        for (x, y), c in Counter(zip(state_list, state_list[1:])).iteritems():
            trans_matrix[x, y] = c

        for i in range(hmm_model.states):
            trans_list = list()
            for j in range(hmm_model.states):
                trans_list.append(trans_matrix[i][j])
            trans_sum = sum(trans_list) * 1.0
            prob_list = [e / trans_sum for e in trans_list]

            # make sure that the sum of this line converge to 1
            total_prob = sum(prob_list)
            if total_prob != 1.0:
                prob_list[0] += (1.0 - total_prob)

            hmm_model.A.append(prob_list)

        # Emission
        for i in range(hmm_model.states):
            norm = list()
            slope = list()
            for j in range(len(state_list)):
                if state_list[j] == i:
                    norm.append(norm_signal[j])
                    slope.append(slope_signal[j])
            # Compute the mean of norm and slope signal
            means_list = list()
            means_list.append(np.mean(norm))
            means_list.append(np.mean(slope))
            hmm_model.means.append(means_list)

            # Compute covariance matrix of norm and slope signal
            covs_list = list()
            covs_matrix = np.cov(norm, slope)
            for j in range(hmm_model.dim):
                for k in range(hmm_model.dim):
                    covs_list.append(covs_matrix[j][k] + 0.000001) # covariance must be symmetric, positive-definite
            hmm_model.covs.append(covs_list)

        output_fname = os.path.join(self.output_location, self.output_prefix)
        hmm_model.save_hmm(output_fname)

    def train1(self):
        # Estimate the HMM parameters using the maximum likelihood method.
        states, norm_signal, slope_signal = self.read_states_signals()
        hmm_model = HMM()

        hmm_model.dim = 3
        # States number
        state_list = [int(state) for state in list(states)]
        hmm_model.states = len(np.unique(np.array(state_list)))

        # Initial state probabilities vector
        init_state = state_list[0]
        hmm_model.pi = [0.0] * hmm_model.states
        hmm_model.pi[init_state] = 1.0

        # Transition
        trans_matrix = np.zeros((hmm_model.states, hmm_model.states))
        for (x, y), c in Counter(zip(state_list, state_list[1:])).iteritems():
            trans_matrix[x, y] = c

        for i in range(hmm_model.states):
            trans_list = list()
            for j in range(hmm_model.states):
                trans_list.append(trans_matrix[i][j])
            trans_sum = sum(trans_list) * 1.0
            prob_list = [e / trans_sum for e in trans_list]

            # make sure that the sum of this line converge to 1
            total_prob = sum(prob_list)
            if total_prob != 1.0:
                prob_list[0] += (1.0 - total_prob)

            hmm_model.A.append(prob_list)

        # Emission
        for i in range(hmm_model.states):
            state_norm = list()
            state_slope = list()
            state_diff = list()
            for j in range(len(state_list)):
                if state_list[j] == i:
                    state_norm.append(norm_signal[j])
                    state_slope.append(slope_signal[j])
                    state_diff.append(diff_signal[j])
            # Compute the mean of norm and slope signal
            means_list = list()
            means_list.append(np.mean(state_norm))
            means_list.append(np.mean(state_slope))
            means_list.append(np.mean(state_diff))
            hmm_model.means.append(means_list)

            # Compute covariance matrix of norm, slope and subtrac signal
            signal = list()
            signal.append(state_norm)
            signal.append(state_slope)
            signal.append(state_diff)

            covs_list = list()
            covs_matrix = np.cov(np.array(signal))
            for j in range(hmm_model.dim):
                for k in range(hmm_model.dim):
                    covs_list.append(covs_matrix[j][k] + 0.000001) # covariance must be symmetric, positive-definite
            hmm_model.covs.append(covs_list)

        if self.estimate_bias_correction:
            model_fname = os.path.join(self.output_locaiton, "{}_{}".format(self.k_nb, self.atac_forward_shift))
        else:
            model_fname = os.path.join(self.output_locaiton, self.output_fname)
        hmm_model.save_hmm(model_fname)

    def output_bed_file(self, states):
        bed_fname = os.path.join(self.output_locaiton, self.output_fname)
        bed_fname += ".bed"

        state_dict = dict([(0, "BACK"), (1, "UPD"), (2, "TOPD"), (3, "DOWND"), (4, "FP")])
        color_dict = dict([(0, "50,50,50"), (1, "10,80,0"), (2, "20,40,150"), (3, "150,20,40"), (4, "198,150,0")])

        state_list = [int(state) for state in list(states)]
        current_state = state_list[0]
        start_postion = self.start
        is_print = False
        with open(bed_fname, "w") as bed_file:
            for i in range(len(state_list)):
                if state_list[i] != current_state:
                    end_position = self.start + i
                    is_print = True
                elif i == len(state_list) - 1:
                    end_position = self.end
                    is_print = True

                if is_print:
                    bed_file.write(self.chrom + " " + str(start_postion) + " " + str(end_position) + " "
                                   + state_dict[current_state] + " " + str(1000) + " . "
                                   + str(start_postion) + " " + str(end_position) + " "
                                   + color_dict[current_state] + "\n")
                    start_postion = end_position
                    current_state = state_list[i]
                    is_print = False

    def train2(self, observe_sequence):
        model = hmm.GaussianHMM(n_components=9, covariance_type="full", n_iter=100)
        model.fit(observe_sequence)
        joblib.dump(model, "filename.pkl")


