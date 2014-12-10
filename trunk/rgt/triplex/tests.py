# Python Libraries
from __future__ import print_function
from __future__ import division
import sys
import os.path
lib_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
#sys.path.append(lib_path)

# Local Libraries
# Distal Libraries
from .. GenomicRegionSet import *
from .. ExperimentalMatrix import *
from rgt.Util import GenomeData, OverlapType, Html

dir = os.getcwd()
"""
Statistical tests and plotting tools for triplex binding site analysis

Author: Joseph Kuo

"""

#################################################################################################
##### FUNCTIONS #################################################################################
#################################################################################################
    

######### Universal functions



def main():
    #################################################################################################
    ##### PARAMETERS ################################################################################
    #################################################################################################
    
    
    
    parser = argparse.ArgumentParser(description='Provides various Statistical tests for triplex binding site analysis \
    \nAuthor: Joseph Kuo, Ivan Gesteira Costa Filho', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')
    
    ################### Projection test ##########################################
    parser_projection = subparsers.add_parser('randomtest',
                                              help='Test the association of')
    parser_projection.add_argument('output', help=helpoutput) 
    parser_projection.add_argument('-r', '--reference',help=helpreference)
    parser_projection.add_argument('-q', '--query', help=helpquery)
    parser_projection.add_argument('-t', '--title', default='projection_test', help=helptitle)
    parser_projection.add_argument('-g', default=None, help=helpgroupbb +" (Default:None)")
    parser_projection.add_argument('-c', default="regions", help=helpcolorbb +' (Default: regions)')
    parser_projection.add_argument('-bg', help="Define a BED file as background. If not defined, the background is whole genome according to the given organism.")
    parser_projection.add_argument('-union', action="store_true", help='Take the union of references as background for binominal test.')
    parser_projection.add_argument('-organism',default='hg19', help='Define the organism. (Default: hg19)')
    parser_projection.add_argument('-log', action="store_true", help='Set y axis of the plot in log scale.')
    parser_projection.add_argument('-color', action="store_true", help=helpDefinedColot)
    #parser_projection.add_argument('-pdf', action="store_true", help='Save the plot in pdf format.')
    #parser_projection.add_argument('-html', action="store_true", help='Save the figure in html format.')
    parser_projection.add_argument('-show', action="store_true", help='Show the figure in the screen.')
    parser_projection.add_argument('-table', action="store_true", help='Store the tables of the figure in text format.')
