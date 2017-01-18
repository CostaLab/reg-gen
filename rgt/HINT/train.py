###################################################################################################
# Libraries
###################################################################################################

#Python
import os
import sys


"""
Train a hidden Markov model (HMM) based on the annotation data

Authors: Eduardo G. Gusmao, Zhijian Li
"""

class TrainHMM:
    """
    Contains methods used to train a hidden Markov model
    """

    def __init__(self, annotate_file, print_bed_file, output_location):
        self.annotate_fname = annotate_file
        self.print_bed_file = print_bed_file
        self.output_location = output_location

    def read_states(self):
        states = ""
        with open(self.annotate_fname) as annotate_file:
            for line in annotate_file:
                if (len(line) < 2 or "#" in line or "=" in line): continue
                ll = line.strip().split(" ")
                states += ll[1:-1]
        return states
