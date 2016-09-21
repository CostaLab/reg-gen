""" 
Binding for calling Triplexator as Python function.  

Author: Barna Zajzon
Revisions: 0.1 01.07.2016
"""

from ctypes import *
import os
from rgt.Util import Triplexator

#some defines

class runTriplexator():
    Triplexator = Triplexator()
    triplex_lib_path = Triplexator.get_path()

    def run(args):
        """
        Provides an interface for calling the Triplexator C++ shared library. TRIPLEXATOR_LIBRARY environment variable
        must be set to the location of the Triplexator shared library.
        :param args: string with command line options that would be used to directly call triplexator
        :return: nothing
        """
        # global TRIPLEXATOR_LIBRARY_ENV
        # triplex_lib_path = os.environ.get(TRIPLEXATOR_LIBRARY_ENV)

        if os.environ.get(TRIPLEXATOR_LIBRARY_ENV) is None:
            print "Please set the environment variable for the Triplexator library (" + TRIPLEXATOR_LIBRARY_ENV + ")."
        else:
            triplex_lib  = cdll.LoadLibrary(triplex_lib_path)
            arg_strings  = args.split(' ')
            arg_ptr      = (c_char_p * (len(arg_strings) + 1))()

            arg_ptr[0] = "triplexator"  # to simulate calling from cmd line
            for i, s in enumerate(arg_strings):
                arg_ptr[i + 1] = s

            triplex_lib.pyTriplexator(len(arg_strings) + 1, arg_ptr)

