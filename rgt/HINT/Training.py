import os
import numpy as np
from collections import Counter
import joblib
from scipy import linalg
from argparse import SUPPRESS

# Internal
from rgt.Util import GenomeData
from rgt.HINT.signalProcessing import GenomicSignal
from rgt.HINT.hmm import HMM, SemiSupervisedGaussianHMM
from rgt.HINT.biasTable import BiasTable

"""
Train a hidden Markov model (HMM) based on the annotation data

Authors: Eduardo G. Gusmao, Zhijian Li
"""


def training_args(parser):
    # Parameters Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19")
    parser.add_argument("--reads-file", type=str, metavar="FILE", default=None,
                        help="The BAM file containing the DNase-seq or ATAC-seq reads")
    parser.add_argument("--chrom", type=str, metavar="STRING", default="chr1",
                        help="The name of chromosome used to train HMM")
    parser.add_argument("--start", type=int, metavar="INT", default=211428000)
    parser.add_argument("--end", type=int, metavar="INT", default=211438000)
    parser.add_argument("--annotate-file", type=str, metavar="STRING",
                        default=None, help="A annotate file containing all the states.")
    parser.add_argument("--bias-table", type=str, metavar="FILE1_F,FILE1_R", default=None)

    # Hidden Options
    parser.add_argument("--downstream-ext", type=int, metavar="INT", default=1, help=SUPPRESS)
    parser.add_argument("--upstream-ext", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--forward-shift", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=0, help=SUPPRESS)
    parser.add_argument("--k-nb", type=int, metavar="INT", default=6, help=SUPPRESS)

    parser.add_argument("--semi-supervised", action="store_true", default=False,
                        help="If used, HMM model will be trained using semi-supervised learning.")
    parser.add_argument("--signal-file", type=str, metavar="FILE", default=None,
                        help="The txt file containing the DNase-seq or ATAC-seq signals used to train HMM model.")
    parser.add_argument("--num-states", type=int, metavar="INT", default=7,
                        help="The states number of HMM model.")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written.")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default=None,
                        help="The prefix for results files.")


def training_run(args):
    if args.semi_supervised:
        semi_supervised(args)
    else:
        supervised(args)


def semi_supervised(args):
    signal_list = list()
    with open(args.signal_file) as f:
        for line in f.readlines():
            signal_list.append(list(map(float, line.strip().split("\t"))))

    states_prior = np.zeros(len(signal_list[0]))
    for i in range(495, 505):
        states_prior[i] = 1

    hmm_model = SemiSupervisedGaussianHMM(n_components=args.num_states, random_state=42, n_iter=10,
                                          covariance_type="full", states_prior=states_prior,
                                          fp_state=args.num_states - 1)

    sequene = np.array(signal_list).T
    hmm_model.fit(sequene)

    # make sure covariance is symmetric and positive-definite
    for i in range(hmm_model.n_components):
        while np.any(np.array(linalg.eigvalsh(hmm_model.covars_[i])) <= 0):
            hmm_model.covars_[i] += 0.000001 * np.eye(hmm_model.covars_[i].shape[0])

    output_fname = os.path.join(args.output_location, "{}.pkl".format(args.output_prefix))
    joblib.dump(hmm_model, output_fname)


def read_states_signals(args):
    # Read states from the annotation file
    states = ""
    with open(args.annotate_file) as f:
        for line in f:
            if len(line) < 2 or "#" in line or "=" in line:
                continue
            ll = line.strip().split(" ")
            for state in ll[1:-1]:
                states += state

    # If need to estimate bias table
    genome_data = GenomeData(args.organism)
    table = None

    # If the bias table is provided
    if args.bias_table:
        bias_table = BiasTable()
        bias_table_list = args.bias_table.split(",")
        table = bias_table.load_table(table_file_name_F=bias_table_list[0],
                                      table_file_name_R=bias_table_list[1])

    # Get the normalization and slope signal from the raw bam file
    raw_signal = GenomicSignal(args.reads_file)
    raw_signal.load_sg_coefs(slope_window_size=9)
    norm_signal, slope_signal = \
        raw_signal.get_signal(args.chrom, args.start, args.end,
                              args.downstream_ext, args.upstream_ext,
                              args.forward_shift, args.reverse_shift,
                              bias_table=table, genome_file_name=genome_data.get_genome())
    if args.print_bed_file:
        args.output_bed_file(states)

    return states, norm_signal, slope_signal


def supervised(args):
    # Estimate the HMM parameters using the maximum likelihood method.
    states, norm_signal, slope_signal = read_states_signals(args)
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
    for (x, y), c in list(Counter(list(zip(state_list, state_list[1:]))).items()):
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
                covs_list.append(covs_matrix[j][k] + 0.000001)  # covariance must be symmetric, positive-definite
        hmm_model.covs.append(covs_list)

    output_fname = os.path.join(args.output_location, args.output_prefix)
    hmm_model.save_hmm(output_fname)


def output_bed_file(args, states):
    bed_fname = os.path.join(args.output_location, "{}.bed".format(args.output_prefix))

    state_dict = dict([(0, "BACK"), (1, "UPD"), (2, "TOPD"), (3, "DOWND"), (4, "FP")])
    color_dict = dict([(0, "50,50,50"), (1, "10,80,0"), (2, "20,40,150"), (3, "150,20,40"), (4, "198,150,0")])

    state_list = [int(state) for state in list(states)]
    current_state = state_list[0]
    start_postion = args.start
    is_print = False
    with open(bed_fname, "w") as bed_file:
        for i in range(len(state_list)):
            if state_list[i] != current_state:
                end_position = args.start + i
                is_print = True
            elif i == len(state_list) - 1:
                end_position = args.end
                is_print = True

            if is_print:
                bed_file.write(args.chrom + " " + str(start_postion) + " " + str(end_position) + " "
                               + state_dict[current_state] + " " + str(1000) + " . "
                               + str(start_postion) + " " + str(end_position) + " "
                               + color_dict[current_state] + "\n")
                start_postion = end_position
                current_state = state_list[i]
                is_print = False
