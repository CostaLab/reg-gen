import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from os import path, mkdir, listdir, walk
from sys import platform, exit
from tqdm import tqdm
from Bio import motifs
import logomaker
import seaborn

# Initializing Option Parser
parser = argparse.ArgumentParser()

g = parser.add_mutually_exclusive_group()

g.add_argument(
    "-a",
    "--all",
    dest="all",
    action="store_true",
    default=False,
    help="Fetch all data sets.",
)
g.add_argument("folders", metavar="repository_folder", nargs="*", default=[])

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
    print(
        "ERROR: This package currently supports only unix-based systems (Linux and MAC OS X)."
    )
    exit(0)

###################################################################################################
# Logo Graphs
###################################################################################################
# Creating logos
output_logos_dir = path.join(curr_dir, "logos")

if not path.exists(output_logos_dir):
    mkdir(output_logos_dir)

repositories = [
    r for r in listdir("motifs") if path.isdir(path.join(curr_dir, "motifs", r))
]

if not repositories:
    print("ERROR: the motifs directory is empty")
    exit(1)

if not args.all:
    repositories = set(repositories)
    query = set(args.folders)

    if not repositories.issuperset(query):
        print(
            "ERROR: query repositories %s do not exist"
            % str(list(query.difference(repositories)))
        )
        exit(1)

    repositories = args.folders

print(">>> CREATING logos for", repositories)


def compute_icm(counts):
    pwm = {k: counts[k] for k in ("A", "C", "G", "T")}
    pwm = pd.DataFrame(data=pwm)
    pwm = pwm.add(1)
    pwm_prob = (pwm.T / pwm.T.sum()).T
    pwm_prob_log = np.log2(pwm_prob)
    pwm_prob_log = pwm_prob_log.mul(pwm_prob)
    info_content = pwm_prob_log.T.sum() + 2
    icm = pwm_prob.mul(info_content, axis=0)

    return icm


def plot_logo(icm, logo_file_name):
    plt.close()

    fig = plt.figure()
    fig.set_size_inches(len(icm), 3)
    ax = fig.add_subplot(111)
    logomaker.Logo(icm, ax=ax, show_spines=False, baseline_width=0)
    ax.set_xticks(range(len(icm)))
    ax.set_yticks(range(0, 3))

    ax.set_xticklabels(range(1, len(icm) + 1))
    ax.set_yticklabels(np.arange(0, 3, 1))

    seaborn.despine(ax=ax, offset=10, trim=True)
    
    fig.tight_layout()
    plt.savefig(logo_file_name, dpi=300)


for repo in repositories:
    dir_name = path.join(curr_dir, "motifs", repo)
    for _, _, file_list in walk(dir_name):
        output_dir = path.join(curr_dir, "logos", repo)

        if not path.exists(output_dir):
            mkdir(output_dir)

        print(">>", repo)

        for pwm_file_name in tqdm(file_list):
            pwm_full_file_name = path.join(dir_name, pwm_file_name)

            if pwm_file_name.split(".")[-1] != "pwm":
                continue

            with open(pwm_full_file_name) as handle:
                counts = motifs.read(handle, "pfm").counts
                icm = compute_icm(counts)
                logo_file_name = path.join(
                    output_dir, ".".join(pwm_file_name.split(".")[:-1]) + ".png"
                )
                plot_logo(icm=icm, logo_file_name=logo_file_name)
