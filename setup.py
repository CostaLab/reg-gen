import os
import sys
import io
import re
from shutil import copy
from pwd import getpwnam
from sys import platform, exit
from distutils import dir_util
from setuptools import setup, find_packages
from os import walk, chown, chmod, path, getenv, makedirs, remove
from optparse import OptionParser, BadOptionError, AmbiguousOptionError

if not sys.version_info[0] == 2:
    sys.exit("Sorry, Python 3 is not supported (yet)")

"""
Installs the RGT tool with standard setuptools options and additional
options specific for RGT.

Authors: Manuel Allhoff, Ivan G. Costa, Eduardo G. Gusmao, Joseph C.C. Kuo, Fabio Ticconi.

Installs the RGT tool with standard setuptools options and additional
options specific for RGT.
"""


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8")
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


current_version = find_version("rgt", "__version__.py")

###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Parameters
###################################################################################################

"""
Tools Dictionary:
  * Insert in the dictionary bellow a key+tuple X: (Y,Z,W,K) representing:
    - X: A string representing the name the user should provide for this script in order to install the tool.
    - Y: A string representing the name of the program which the user should type in the terminal in order to execute the tool.
    - Z: A string representing the path to the main function that executes the program.
    - W: A list with the package requirements for that program.
    - K: A list with external binary files that will be copied by the setup installation function to the user's bin folder.
Tools Dictionary Standard:
  * All programs should start with "rgt-" followed by the program name.
  * The main function called within the script must be termed "main".
"""

# TODO: sometime in the future, we should use environment markers for platform-specific dependencies.
# We currently can't, because the change is recent and many people have older python versions.
# We don't like crashes at installation time..

common_deps = ["cython",
               "numpy>=1.4.0",
               "scipy>=0.7.0",
               "pysam==0.11.1",
               "pyBigWig",
               "PyVCF"]

if platform.startswith("darwin"):
    bin_dir = "mac"
    libRGT = "librgt_mac.so"
else:
    bin_dir = "linux"
    libRGT = "librgt_linux.so"
    common_deps.append("ngslib")

tools_dictionary = {
    "core": (
        None,
        None,
        [],
        ["data/bin/" + bin_dir + "/bedToBigBed", "data/bin/" + bin_dir + "/bigBedToBed",
         "data/bin/" + bin_dir + "/wigToBigWig", "data/bin/" + bin_dir + "/bigWigMerge",
         "data/bin/" + bin_dir + "/bedGraphToBigWig"]
    ),
    "motifanalysis": (
        "rgt-motifanalysis",
        "rgt.motifanalysis.Main:main",
        ["Biopython>=1.64", "fisher>=0.1.4"],
        ["data/bin/" + bin_dir + "/bedToBigBed", "data/bin/" + bin_dir + "/bigBedToBed"]
    ),
    "hint": (
        "rgt-hint",
        "rgt.HINT.Main:main",
        ["scikit-learn>=0.14", "hmmlearn", "pyx==0.12.1"],
        []
    ),
    "THOR": (
        "rgt-THOR",
        "rgt.THOR.THOR:main",
        ["scikit-learn>=0.17.1", "hmmlearn", "matplotlib>=1.1.0", "mpmath", "HTSeq"],
        ["data/bin/" + bin_dir + "/wigToBigWig", "data/bin/" + bin_dir + "/bigWigMerge",
         "data/bin/" + bin_dir + "/bedGraphToBigWig"]
    ),
    "filterVCF": (
        "rgt-filterVCF",
        "rgt.filterVCF.filterVCF:main",
        [],
        []
    ),
    "viz": (
        "rgt-viz",
        "rgt.viz.Main:main",
        ["matplotlib>=1.1.0", "matplotlib_venn"],
        []
    ),
    "TDF": (
        "rgt-TDF",
        "rgt.tdf.Main:main",
        ["matplotlib>=1.1.0", "natsort"],
        []
    )
}


###################################################################################################
# Auxiliary Functions/Classes
###################################################################################################


# PassThroughOptionParser Class
class PassThroughOptionParser(OptionParser):
    """
    An unknown option pass-through implementation of OptionParser.
    When unknown arguments are encountered, bundle with largs and try again,
    until rargs is depleted.
    sys.exit(status) will still be called if a known argument is passed
    incorrectly (e.g. missing arguments or bad argument types, etc.)
    """

    def _process_args(self, largs, rargs, values):
        while rargs:
            try:
                OptionParser._process_args(self, largs, rargs, values)
            except (BadOptionError, AmbiguousOptionError), e:
                largs.append(e.opt_str)


# recursive_chown_chmod Function
def recursive_chown_chmod(path_to_walk, uid, gid, file_permission, path_permission):
    """
    Recursively applies chown from path.
    """
    for root_dir, directory_list, file_list in walk(path_to_walk):
        chown(root_dir, uid, gid)
        chmod(root_dir, path_permission)
        for f in file_list:
            current_complete_file = path.join(root_dir, f)
            chown(current_complete_file, uid, gid)
            chmod(current_complete_file, file_permission)


###################################################################################################
# Processing Input Arguments
###################################################################################################

# Parameters
rgt_data_base_name = "rgtdata"
usage_message = "python setup.py install [python options] [RGT options]"
version_message = "Regulatory Genomics Toolbox (RGT). Version: " + str(current_version)

# Initializing Option Parser
parser = PassThroughOptionParser(usage=usage_message, version=version_message)

# Parameter: RGT Data Location
param_rgt_data_location_name = "--rgt-data-path"
parser.add_option(param_rgt_data_location_name, type="string", metavar="STRING",
                  help="Path containing data used by RGT tool.",
                  dest="param_rgt_data_location", default=path.join(getenv('HOME'), rgt_data_base_name))

# Parameter: Tool
param_rgt_tool_name = "--rgt-tool"
parser.add_option(param_rgt_tool_name, type="string", metavar="STRING",
                  help=("The tool which will be installed. If this argument is not used, "
                        "then the complete package is installed. The current available options "
                        "are: " + ", ".join(tools_dictionary.keys()) + "; where 'core' means that "
                        "only the RGT python library will be installed with no further command-line "
                        "tools. You can also provide multiple tools in a list separated by comma."),
                  dest="param_rgt_tool", default=",".join(tools_dictionary.keys()))

# Processing Options
options, arguments = parser.parse_args()
if path.basename(options.param_rgt_data_location) != rgt_data_base_name:
    options.param_rgt_data_location = path.join(options.param_rgt_data_location, rgt_data_base_name)
if options.param_rgt_data_location[0] == "~":
    options.param_rgt_data_location = path.join(getenv('HOME'), options.param_rgt_data_location[2:])
options.param_rgt_tool = options.param_rgt_tool.split(",")

# Manually Removing Additional Options from sys.argv
new_sys_argv = []
for e in sys.argv:
    if param_rgt_data_location_name == e[:len(param_rgt_data_location_name)]:
        continue
    elif param_rgt_tool_name == e[:len(param_rgt_tool_name)]:
        continue
    new_sys_argv.append(e)

sys.argv = new_sys_argv

# Defining entry points
current_entry_points = {"console_scripts": []}
for tool_option in options.param_rgt_tool:
    if tool_option != "core":
        current_entry_points["console_scripts"].append(" = ".join(tools_dictionary[tool_option][:2]))

# Defining install requirements
current_install_requires = common_deps
for tool_option in options.param_rgt_tool:
    current_install_requires += tools_dictionary[tool_option][2]

###################################################################################################
# Creating Data Path
###################################################################################################

# Creating Data Path
if not path.exists(options.param_rgt_data_location):
    makedirs(options.param_rgt_data_location)

# Creating data.config
data_config_file_name = path.join(options.param_rgt_data_location, "data.config")
# if not os.path.isfile(data_config_file_name):
data_config_file = open(data_config_file_name, "w")
data_config_file.write("# Configuration file loaded at rgt startup. CAREFUL: any changes shall be overwritten\n"
                       "# whenever rgt is (re)installed. Use data.config.user for permanent changes.\n\n")

genome = "mm9"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_mm9.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.mm9\n"))
data_config_file.write("genes_Gencode: " + path.join(genome_dir, "genes_Gencode_mm9.bed\n"))
data_config_file.write("genes_RefSeq: " + path.join(genome_dir, "genes_RefSeq_mm9.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "gencode.vM1.annotation.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_mouse.txt\n\n"))
data_config_file.write("repeat_maskers: " + path.join(genome_dir, "repeat_maskers\n\n"))
genome = "mm10"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_mm10.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.mm10\n"))
data_config_file.write("genes_Gencode: " + path.join(genome_dir, "genes_Gencode_mm10.bed\n"))
data_config_file.write("genes_RefSeq: " + path.join(genome_dir, "genes_RefSeq_mm10.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "gencode.vM11.annotation.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_mouse.txt\n\n"))
genome = "hg19"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_hg19.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.hg19\n"))
data_config_file.write("genes_Gencode: " + path.join(genome_dir, "genes_Gencode_hg19.bed\n"))
data_config_file.write("genes_RefSeq: " + path.join(genome_dir, "genes_RefSeq_hg19.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "gencode.v19.annotation.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_human.txt\n\n"))
data_config_file.write("repeat_maskers: " + path.join(genome_dir, "repeat_maskers\n\n"))
genome = "hg38"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_hg38.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.hg38\n"))
data_config_file.write("genes_Gencode: " + path.join(genome_dir, "genes_Gencode_hg38.bed\n"))
data_config_file.write("genes_RefSeq: " + path.join(genome_dir, "genes_RefSeq_hg38.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "gencode.v24.annotation.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_human.txt\n\n"))
data_config_file.write("repeat_maskers: " + path.join(genome_dir, "repeat_maskers\n\n"))
genome = "zv9"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_zv9_ensembl_release_79.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.zv9\n"))
data_config_file.write("gene_regions: " + path.join(genome_dir, "genes_zv9.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "Danio_rerio.Zv9.79.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_zebrafish.txt\n\n"))
genome = "zv10"
genome_dir = path.join(options.param_rgt_data_location, genome)
data_config_file.write("[" + genome + "]\n")
data_config_file.write("genome: " + path.join(genome_dir, "genome_zv10_ensembl_release_84.fa\n"))
data_config_file.write("chromosome_sizes: " + path.join(genome_dir, "chrom.sizes.zv10\n"))
data_config_file.write("gene_regions: " + path.join(genome_dir, "genes_zv10.bed\n"))
data_config_file.write("annotation: " + path.join(genome_dir, "Danio_rerio.GRCz10.84.gtf\n"))
data_config_file.write("gene_alias: " + path.join(genome_dir, "alias_zebrafish.txt\n\n"))

data_config_file.write("[MotifData]\n")
data_config_file.write("pwm_dataset: motifs\n")
data_config_file.write("logo_dataset: logos\n")
data_config_file.write("repositories: jaspar_vertebrates,uniprobe_primary\n\n")
data_config_file.write("[HmmData]\n")
data_config_file.write("default_hmm_dnase: fp_hmms/dnase.hmm\n")
data_config_file.write("default_hmm_dnase_bc: fp_hmms/dnase_bc.hmm\n")
data_config_file.write("default_hmm_atac: fp_hmms/atac.hmm\n")
data_config_file.write("default_hmm_atac_bc: fp_hmms/atac.pkl\n")
data_config_file.write("default_hmm_histone: fp_hmms/histone.hmm\n")
data_config_file.write("default_hmm_dnase_histone: fp_hmms/dnase_histone.hmm\n")
data_config_file.write("default_hmm_dnase_histone_bc: fp_hmms/dnase_histone_bc.hmm\n")
data_config_file.write("default_hmm_atac_histone: fp_hmms/atac_histone.hmm\n")
data_config_file.write("default_hmm_atac_histone_bc: fp_hmms/atac_histone_bc.hmm\n")
data_config_file.write("default_bias_table_F_SH: fp_hmms/single_hit_bias_table_F.txt\n")
data_config_file.write("default_bias_table_R_SH: fp_hmms/single_hit_bias_table_R.txt\n")
data_config_file.write("default_bias_table_F_DH: fp_hmms/double_hit_bias_table_F.txt\n")
data_config_file.write("default_bias_table_R_DH: fp_hmms/double_hit_bias_table_R.txt\n")
data_config_file.write("default_bias_table_F_ATAC: fp_hmms/atac_bias_table_F.txt\n")
data_config_file.write("default_bias_table_R_ATAC: fp_hmms/atac_bias_table_R.txt\n\n")
data_config_file.write("[Library]\n")
data_config_file.write("path_triplexator: " + path.join(options.param_rgt_data_location, "lib/libtriplexator.so") + "\n")
data_config_file.write("path_c_rgt: " + path.join(options.param_rgt_data_location, "lib/" + libRGT) + "\n")

data_config_file.close()

# Creating data.config.user, but only if not already present
user_config_file_name = path.join(options.param_rgt_data_location, "data.config.user")
if not os.path.isfile(user_config_file_name):
    user_config_file = open(user_config_file_name, "w")

    user_config_file.write("# Here you can overwrite any property set in the data.config file. It shall not be\n"
                           "# be overwritten in any case, so if you are experiencing problems rename or remove this\n"
                           "# file. See data.config for how the file should be formatted.\n\n")
    genome = "self_defined"
    genome_dir = path.join(options.param_rgt_data_location, genome)
    user_config_file.write("# Template to add a genomic section.\n")
    user_config_file.write("#[" + genome + "]\n")
    user_config_file.write("#genome: undefined\n")
    user_config_file.write("#chromosome_sizes: undefined\n")
    user_config_file.write("#gene_regions: undefined\n")
    user_config_file.write("#annotation: undefined\n")
    user_config_file.write("#gene_alias: undefined\n\n")

# Creating data.config.path
script_dir = path.dirname(path.abspath(__file__))
data_config_path_file_name = path.join(script_dir, "rgt", "data.config.path")
data_config_path_file = open(data_config_path_file_name, "w")
data_config_path_file.write(data_config_file_name)
data_config_path_file.close()

# Copying data from package folder to installation folder
"""
Copy Files Dictionary:
  * Insert in the dictionary below a key + list in the format X:[Y1,Y2,...Yn], representing:
    - X: A path in the data folder structure.
    - Y1,Y2,...Yn: files/folders inside the path X to be copied.
"""
copy_files_dictionary = {
    ".": ["setupGenomicData.py", "setupLogoData.py"],
    "lib": ["libtriplexator.so", libRGT],
    "hg19": ["genes_Gencode_hg19.bed", "chrom.sizes.hg19", "alias_human.txt", "genes_RefSeq_hg19.bed"],
    "hg38": ["genes_Gencode_hg38.bed", "chrom.sizes.hg38", "alias_human.txt", "genes_RefSeq_hg38.bed"],
    "mm9": ["genes_Gencode_mm9.bed", "chrom.sizes.mm9", "alias_mouse.txt", "genes_RefSeq_mm9.bed"],
    "mm10": ["genes_Gencode_mm10.bed", "chrom.sizes.mm10", "alias_mouse.txt", "genes_RefSeq_mm10.bed"],
    "zv9": ["genes_zv9.bed", "chrom.sizes.zv9", "alias_zebrafish.txt"],
    "zv10": ["genes_zv10.bed", "chrom.sizes.zv10", "alias_zebrafish.txt"],
    "fp_hmms": ["dnase.hmm", "dnase_bc.hmm", "histone.hmm", "dnase_histone.hmm", "dnase_histone_bc.hmm",
                "single_hit_bias_table_F.txt", "single_hit_bias_table_R.txt", "atac.pkl", "atac.hmm", "atac_bc.hmm",
                "atac_bias_table_F.txt", "atac_bias_table_R.txt", "atac_histone.hmm", "atac_histone_bc.hmm",
                "double_hit_bias_table_F.txt", "double_hit_bias_table_R.txt", "H3K4me3_proximal.hmm"],
    "motifs": ["jaspar_vertebrates", "uniprobe_primary", "uniprobe_secondary", "hocomoco", "hocomoco.fpr",
               "jaspar_vertebrates.fpr", "uniprobe_primary.fpr", "uniprobe_secondary.fpr"],
    "fig": ["rgt_logo.gif", "style.css", "default_motif_logo.png", "jquery-1.11.1.js", "jquery.tablesorter.min.js",
            "tdf_logo.png", "viz_logo.png"],
}
for copy_folder in copy_files_dictionary.keys():
    copy_dest_path = path.join(options.param_rgt_data_location, copy_folder)
    if not path.exists(copy_dest_path): makedirs(copy_dest_path)
    for copy_file in copy_files_dictionary[copy_folder]:
        copy_source_file = path.join(script_dir, "data", copy_folder, copy_file)
        copy_dest_file = path.join(copy_dest_path, copy_file)
        if os.path.isfile(copy_source_file):
            copy(copy_source_file, copy_dest_file)
        else:
            dir_util.copy_tree(copy_source_file, copy_dest_file)

            # if not path.exists(copy_dest_file):
            # try: dir_util.copy_tree(copy_source_file, copy_dest_file)
            # except OSError as exc:
            #     if exc.errno == ENOTDIR:
            #         copy(copy_source_file, copy_dest_file)
            #     else:
            #         raise

###################################################################################################
# Setup Function
###################################################################################################

# Parameters
short_description = "Toolkit to perform regulatory genomics data analysis"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatics",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["ChIP-seq", "DNase-seq", "Peak Calling", "Motif Discovery", "Motif Enrichment", "HMM"]
author_list = ["Eduardo G. Gusmao", "Manuel Allhoff", "Joseph Chao-Chung Kuo", "Fabio Ticconi", "Ivan G. Costa"]
corresponding_mail = "software@costalab.org"
license_type = "GPL"
package_data_dictionary = {"rgt": [path.basename(data_config_path_file_name)]}

# External scripts
external_scripts = []
for tool_option in options.param_rgt_tool:
    for e in tools_dictionary[tool_option][3]:
        external_scripts.append(e)

# Fetching Additional Structural Files
readme_file_name = path.join(path.dirname(path.abspath(__file__)), "README.rst")

# Fetching Long Description
readme_file = open(readme_file_name, "r")
long_description = readme_file.read()
readme_file.close()

# Setup Function

setup(name="RGT",
      version=current_version,
      description=short_description,
      long_description=long_description,
      classifiers=classifiers_list,
      keywords=", ".join(keywords_list),
      author=", ".join(author_list),
      author_email=corresponding_mail,
      license=license_type,
      packages=find_packages(),
      package_data=package_data_dictionary,
      entry_points=current_entry_points,
      install_requires=current_install_requires,
      scripts=external_scripts,
      url="http://www.regulatory-genomics.org",
      download_url="https://github.com/CostaLab/reg-gen/archive/{0}.zip".format(current_version),
      platforms=supported_platforms)

###################################################################################################
# Termination
###################################################################################################

# Removing data.config.path
remove(data_config_path_file_name)

# Modifying Permissions when Running Superuser/Admin
# $SUDO_USER exists only if you are sudo, and returns the original user name
current_user = getenv("SUDO_USER")
default_file_permission = 0644
default_path_permission = 0755

if current_user:
    current_user_uid = getpwnam(current_user).pw_uid
    current_user_gid = getpwnam(current_user).pw_gid
    recursive_chown_chmod(options.param_rgt_data_location, current_user_uid, current_user_gid, default_file_permission,
                          default_path_permission)
