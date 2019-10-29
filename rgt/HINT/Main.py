
import sys
import time
from random import seed
from argparse import ArgumentParser

import warnings

if not sys.warnoptions:
    warnings.simplefilter("ignore")

# Internal
from rgt import __version__
from rgt.HINT.Training import training_args, training_run
from rgt.HINT.Plotting import plotting_args, plotting_run
from rgt.HINT.DifferentialAnalysis import diff_analysis_args, diff_analysis_run
from rgt.HINT.Footprinting import footprinting_args, footprinting_run
from rgt.HINT.Estimation import estimation_args, estimation_run
from rgt.HINT.Evaluation import evaluation_args, evaluation_run
from rgt.HINT.Evidence import evidence_args, evidence_run
from rgt.HINT.Tracks import tracks_args, tracks_run

"""
HINT - HMM-based Identification of TF Footprints.
Finds transcription factor footprints given open chromatin data.

Basic Input:
- Regions (bed) in which to find footprints (i.e. enriched regions or hypersensitivity regions).
- Reads (bam) containing the open chromatin signal for DNase/ATAC and 0 <= N <= 3 histone modifications.

Authors: Eduardo G. Gusmao, Manuel Allhoff, Joseph Kuo and Ivan G. Costa.
"""


def main():
    """
    Main function that performs footprint analysis.
    Keyword arguments: None
    Return: None
    """
    start = time.time()

    seed(42)

    version_message = "HINT - Regulatory Analysis Toolbox (RGT) - v" + str(__version__)

    parser = ArgumentParser(prog='rgt-hint')
    parser.add_argument('--version', action='version', version=version_message)

    subparsers = parser.add_subparsers(help='Commands:')

    footprinting_parser = subparsers.add_parser('footprinting',
                                                help='detect footprints based on reads.bam and regions.bed')
    footprinting_args(footprinting_parser)
    footprinting_parser.set_defaults(func=footprinting_run)

    diff_analysis_parser = subparsers.add_parser('differential',
                                                 help='perform differential analysis based on footprints of two '
                                                      'conditions')
    diff_analysis_args(diff_analysis_parser)
    diff_analysis_parser.set_defaults(func=diff_analysis_run)

    plotting_parser = subparsers.add_parser('plotting', help='generate plots based on input')
    plotting_args(plotting_parser)
    plotting_parser.set_defaults(func=plotting_run)

    training_parser = subparsers.add_parser('training', help='train Hidden Markov models')
    training_args(training_parser)
    training_parser.set_defaults(func=training_run)

    estimation_parser = subparsers.add_parser('estimation', help='estimate sequence specific bias')
    estimation_args(estimation_parser)
    estimation_parser.set_defaults(func=estimation_run)

    evaluation_parser = subparsers.add_parser('evaluation',
                                              help='evaluate the predicted footprints based on ChIP-seq data')
    evaluation_args(evaluation_parser)
    evaluation_parser.set_defaults(func=evaluation_run)

    evidence_parser = subparsers.add_parser('evidence',
                                            help='create evidence file based on motif matching results and '
                                                 'ChIP-seq data')
    evidence_args(evidence_parser)
    evidence_parser.set_defaults(func=evidence_run)

    tracks_parser = subparsers.add_parser('tracks', help='create wig track file for visualization')
    tracks_args(tracks_parser)
    tracks_parser.set_defaults(func=tracks_run)

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    args.func(args)

    secs = time.time() - start
    m, s = divmod(secs, 60)
    h, m = divmod(m, 60)

    print()
    print("[total time: ", "%dh %dm %ds" % (h, m, s), "]", sep="")
