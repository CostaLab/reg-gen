###################################################################################################
# Setup Logo Data Script
#   This script will create all motif logos. Be aware that this is a time consuming procedure
#   depending on the total number of PWMs.
# Dependencies:
#   Biopython -> http://biopython.org/wiki/Main_Page
###################################################################################################

# Import
from __future__ import print_function

from optparse import OptionParser
from os import path, walk, mkdir
from sys import platform

from Bio import motifs

# Optional Input
usage_message = "python setupLogoData.py [options]"

# Initializing Option Parser
parser = OptionParser(usage=usage_message)

# Parameter: RGT Data Location
parser.add_option("--all", dest="all", action="store_true", default=False,
                  help="Fetch all data sets.")
parser.add_option("--hocomoco", dest="hocomoco", action="store_true", default=False,
                  help="Creates logos only for HOCOMOCO repository.")
parser.add_option("--jaspar-vertebrates", dest="jaspar_vertebrates", action="store_true", default=False,
                  help="Creates logos only for Jaspar Vertebrates repository.")
parser.add_option("--uniprobe-primary", dest="uniprobe_primary", action="store_true", default=False,
                  help="Creates logos only for Uniprobe (primary motifs) repository.")
parser.add_option("--uniprobe-secondary", dest="uniprobe_secondary", action="store_true", default=False,
                  help="Creates logos only for Uniprobe (secondary motifs) repository.")
options, arguments = parser.parse_args()
if options.all:
    options.hocomoco = True
    options.jaspar_vertebrates = True
    options.uniprobe_primary = True
    options.uniprobe_secondary = True

###################################################################################################
# Parameters
###################################################################################################

# Current rgt data path
curr_dir = path.dirname(path.realpath(__file__))

# Platform
supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Logo Graphs
###################################################################################################

# Creating logos
output_logos_dir = path.join(curr_dir, "logos")
if not path.exists(output_logos_dir):
    mkdir(output_logos_dir)
for dir_name, subdir_list, file_list in walk(path.join(curr_dir, "motifs")):
    base_name = path.basename(dir_name)
    if ((options.hocomoco and base_name == "hocomoco") or
            (options.jaspar_vertebrates and base_name == "jaspar_vertebrates") or
            (options.uniprobe_primary and base_name == "uniprobe_primary") or
            (options.uniprobe_secondary and base_name == "uniprobe_secondary")):
        output_dir = path.join(curr_dir, "logos", base_name)
        if not path.exists(output_dir):
            mkdir(output_dir)
        else:
            continue
        print("Creating logos for " + base_name)
        for pwm_file_name in file_list:
            pwm_full_file_name = path.join(dir_name, pwm_file_name)
            if pwm_file_name.split(".")[-1] != "pwm":
                continue
            pwm_file = open(pwm_full_file_name, "r")
            logo_file_name = path.join(output_dir, ".".join(pwm_file_name.split(".")[:-1]) + ".png")
            pwm = motifs.read(pwm_file, "pfm")
            pwm.weblogo(logo_file_name, format="png_print", stack_width="medium", color_scheme="color_classic")
            pwm_file.close()
        print("OK")
