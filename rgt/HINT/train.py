###################################################################################################
# Libraries
###################################################################################################

# Python
import numpy as np
from collections import Counter

# Internal
from signalProcessing import GenomicSignal
from hmm import HMM

"""
Train a hidden Markov model (HMM) based on the annotation data

Authors: Eduardo G. Gusmao, Zhijian Li
"""


class TrainHMM:
    """
    Contains methods used to train a hidden Markov model
    """

    def __init__(self, bam_file, annotate_file, print_bed_file, output_fname):
        self.bam_file = bam_file
        self.annotate_fname = annotate_file
        self.print_bed_file = print_bed_file
        self.output_fname = output_fname
        self.chrom = "chr1"
        self.start = 211428000
        self.end = 211438000

    def read_states_signals(self):
        # Read states from the annotation file
        states = ""
        with open(self.annotate_fname) as annotate_file:
            for line in annotate_file:
                if (len(line) < 2 or "#" in line or "=" in line):
                    continue
                ll = line.strip().split(" ")
                for state in ll[1:-1]:
                    states += state

        # Get the normalization and slope signal from the raw bam file
        raw_signal = GenomicSignal(self.bam_file)
        raw_signal.load_sg_coefs(slope_window_size=9)
        norm_signal, slope_signal = raw_signal.get_signal(ref=self.chrom, start=self.start, end=self.end,
                                                          downstream_ext=1, upstream_ext=0,
                                                          forward_shift=0, reverse_shift=0)

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
            trans_sum = sum(trans_list)*1.0
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
            print(covs_matrix)
            for j in range(hmm_model.dim):
                for k in range(hmm_model.dim):
                    covs_list.append(covs_matrix[j][k])
            hmm_model.covs.append(covs_list)

        hmm_model.save_hmm(self.output_fname)