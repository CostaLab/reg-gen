###################################################################################################
# Setup Logo Data Script
#   This script will create all motif logos. Be aware that this is a time consuming procedure
#   depending on the total number of PWMs.
# Dependencies:
#   Biopython -> http://biopython.org/wiki/Main_Page
###################################################################################################

# Import
import argparse
from os import path, mkdir, listdir, walk
from sys import platform, exit

from Bio import motifs

# Initializing Option Parser
parser = argparse.ArgumentParser()

g = parser.add_mutually_exclusive_group()

g.add_argument("-a", "--all", dest="all", action="store_true", default=False, help="Fetch all data sets.")
g.add_argument('folders', metavar="repository_folder", nargs='*', default=[])

args = parser.parse_args()

if not args.all and not args.folders:
    parser.print_help()
    exit(1)

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

repositories = [r for r in listdir("motifs") if path.isdir(path.join(curr_dir, "motifs", r))]

if not repositories:
    print("ERROR: the motifs directory is empty")
    exit(1)

if not args.all:
    repositories = set(repositories)
    query = set(args.folders)

    if not repositories.issuperset(query):
        print("ERROR: query repositories %s do not exist" % str(list(query.difference(repositories))))
        exit(1)

    repositories = args.folders

print(">>> CREATING logos for", repositories)

for repo in repositories:
    dir_name = path.join(curr_dir, "motifs", repo)
    for _, _, file_list in walk(dir_name):
        output_dir = path.join(curr_dir, "logos", repo)

        if not path.exists(output_dir):
            mkdir(output_dir)

        print(">>", repo)

        for pwm_file_name in file_list:
            pwm_full_file_name = path.join(dir_name, pwm_file_name)
            if pwm_file_name.split(".")[-1] != "pwm":
                continue
            pwm_file = open(pwm_full_file_name, "r")
            logo_file_name = path.join(output_dir, ".".join(pwm_file_name.split(".")[:-1]) + ".png")
            pwm = motifs.read(pwm_file, "pfm")
            pwm.weblogo(logo_file_name, format="png_print", stack_width="medium", color_scheme="color_classic")
            pwm_file.close()
