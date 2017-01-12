
###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import sys
import math
import operator
import numpy as np
from pickle import load
from sklearn.metrics import auc
from scipy.integrate import simps, trapz
from optparse import OptionParser,BadOptionError,AmbiguousOptionError

"""
Evaluate the footprints prediction using TF ChIP-seq or expression data.

Authors: Eduardo G. Gusmao.
"""