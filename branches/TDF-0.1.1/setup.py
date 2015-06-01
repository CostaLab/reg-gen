import sys # Cannot be changed
from os import walk, chown, chmod, path, getenv, makedirs, remove
from sys import platform, exit
from pwd import getpwnam
from shutil import copy, copytree
import distutils.dir_util
from errno import ENOTDIR
from optparse import OptionParser, BadOptionError, AmbiguousOptionError
from setuptools import setup, find_packages

"""
Installs the RGT tool with standard setuptools options and additional
options specific for RGT.

Authors: Eduardo G. Gusmao, Manuel Allhoff, Joseph Kuo, Ivan G. Costa.

Installs the RGT tool with standard setuptools options and additional
options specific for RGT.
"""

###################################################################################################
# Unsupported Platforms
###################################################################################################

supported_platforms = ["linux","linux2","darwin"]
if(platform not in supported_platforms):
  print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
  exit(0)

###################################################################################################
# Parameters
###################################################################################################

"""
Current Version Standard: X.Y.Z
X: Major RGT release.
Y: Major Specific-Tool release.
Z: Minor RGT/Specific-Tool bug corrections.
"""
current_version = "0.0.1"

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
tools_dictionary = {
"core": (
  None,
  None,
  [],
  []
),
#"motifanalysisold": (
#  "rgt-motifanalysisold",
#  "rgt.motifanalysisold.main:main",
#  ["numpy>=1.4.0","scipy>=0.7.0","Biopython>=1.60","pandas==0.7.1","fisher>=0.1.0","statsmodels>=0.4.0","HTML>=0.04","matplotlib>=1.1.0"],
#  []
#), 
"motifanalysis": (
  "rgt-motifanalysis",
  "rgt.motifanalysis.Main:main",
  ["numpy>=1.4.0","scipy>=0.7.0","Biopython>=1.64","pysam>=0.7.5","fisher>=0.1.4"],
  ["data/bin/bedToBigBed","data/bin/bigBedToBed"]
), 
"hint": (
  "rgt-hint",
  "rgt.HINT.Main:main",
  ["numpy>=1.4.0","scipy>=0.7.0","scikit-learn<=0.14","pysam>=0.7.5"],
  []
), 
"ODIN": (
  "rgt-ODIN",
  "rgt.ODIN.ODIN:main",
  ["hmmlearn", "scikit-learn", "numpy>=1.4.0", "scipy>=0.7.0", "pysam>=0.7.5", "HTSeq", "mpmath"],
  []
), 
"THOR": (
"rgt-THOR",
"rgt.THOR.THOR:main",
["hmmlearn", "scikit-learn", "numpy>=1.4.0", "scipy>=0.7.0", "pysam>=0.7.5", "HTSeq", "mpmath"],
[]
),                 
"filterVCF": (
"rgt-filterVCF",
"rgt.filterVCF.filterVCF:main",
["PyVCF", "numpy>=1.4.0", "scipy>=0.7.0"],
[]
),
"viz": (
  "rgt-viz",
  "rgt.viz.Main:main",
  ["numpy>=1.4.0","scipy>=0.7.0","matplotlib>=1.1.0"],
  []
),
"TDF": (
  "rgt-TDF",
  "rgt.triplex.Main:main",
  ["numpy>=1.4.0","scipy>=0.7.0","matplotlib>=1.1.0", "pysam>=0.7.5"],
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
                OptionParser._process_args(self,largs,rargs,values)
            except (BadOptionError,AmbiguousOptionError), e:
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
            current_complete_file = path.join(root_dir,f)
            chown(current_complete_file, uid, gid)
            chmod(current_complete_file, file_permission)

###################################################################################################
# Processing Input Arguments
###################################################################################################

# Parameters
rgt_data_base_name = "rgtdata"
usage_message = "python setup.py install [python options] [RGT options]"
version_message = "Regulatory Genomics Toolbox (RGT). Version: "+str(current_version)

# Initializing Option Parser
parser = PassThroughOptionParser(usage = usage_message, version = version_message)

# Parameter: RGT Data Location
param_rgt_data_location_name = "--rgt-data-path"
parser.add_option(param_rgt_data_location_name, type = "string", metavar="STRING",
                  help = "Path containing data used by RGT tool.",
                  dest = "param_rgt_data_location", default = path.join(getenv('HOME'),rgt_data_base_name))

# Parameter: Tool
param_rgt_tool_name = "--rgt-tool"
parser.add_option(param_rgt_tool_name, type = "string", metavar="STRING",
                  help = ("The tool which will be installed. If this argument is not used, "
                          "then the complete package is installed. The current available options"
                          "are: "+", ".join(tools_dictionary.keys())+"; where 'core' means that"
                          "only the RGT python library will be installed with no further command-line"
                          "tools. You can also provide multiple tools in a list separated by comma."), 
                  dest = "param_rgt_tool", default = ",".join(tools_dictionary.keys()))

# Processing Options
options, arguments = parser.parse_args()
if(path.basename(options.param_rgt_data_location) != rgt_data_base_name):
    options.param_rgt_data_location = path.join(options.param_rgt_data_location,rgt_data_base_name)
if(options.param_rgt_data_location[0] == "~"):
    options.param_rgt_data_location = path.join(getenv('HOME'),options.param_rgt_data_location[2:])
options.param_rgt_tool = options.param_rgt_tool.split(",")

# Manually Removing Additional Options from sys.argv
new_sys_argv = []
for e in sys.argv:
    if(param_rgt_data_location_name == e[:len(param_rgt_data_location_name)]): continue
    elif(param_rgt_tool_name == e[:len(param_rgt_tool_name)]): continue
    new_sys_argv.append(e)

sys.argv = new_sys_argv

# Defining entry points
current_entry_points = {"console_scripts" : []}
for tool_option in options.param_rgt_tool:
    if(tool_option != "core"):
        current_entry_points["console_scripts"].append(" = ".join(tools_dictionary[tool_option][:2]))

# Defining install requirements
current_install_requires = []
for tool_option in options.param_rgt_tool: current_install_requires += tools_dictionary[tool_option][2]

###################################################################################################
# Creating Data Path
###################################################################################################

# Creating Data Path
if not path.exists(options.param_rgt_data_location):
    makedirs(options.param_rgt_data_location)

# Creating data.config
data_config_file_name = path.join(options.param_rgt_data_location, "data.config")
data_config_file = open(data_config_file_name,"w")
data_config_file.write("[GenomeData]\n")
data_config_file.write("genome: genome.fa\n")
data_config_file.write("chromosome_sizes: chrom.sizes\n")
data_config_file.write("association_file: association_file.bed\n")
data_config_file.write("gencode_annotation: gencode_annotation.gtf\n")
data_config_file.write("gene_alias: alias.txt\n\n")
data_config_file.write("[MotifData]\n")
data_config_file.write("pwm_dataset: motifs\n")
data_config_file.write("logo_dataset: logos\n")
data_config_file.write("repositories: jaspar_vertebrates,uniprobe_primary\n\n")
data_config_file.write("[HmmData]\n")
data_config_file.write("default_hmm: fp_hmms/H3K4me3_proximal.hmm\n\n")
data_config_file.close()

# Creating data.config.path
script_dir = path.dirname(path.abspath(__file__))
data_config_path_file_name = path.join(script_dir,"rgt","data.config.path")
data_config_path_file = open(data_config_path_file_name,"w")
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
".": ["setupGenomicData.py","setupLogoData.py"],
"hg19": ["association_file.bed","chrom.sizes","alias.txt"],
"mm9": ["association_file.bed","chrom.sizes","alias.txt"],
"fp_hmms": ["H3K4me3_proximal.hmm"],
"motifs": ["jaspar_vertebrates", "uniprobe_primary", "uniprobe_secondary", "hocomoco", "hocomoco.fpr", "jaspar_vertebrates.fpr", "uniprobe_primary.fpr", "uniprobe_secondary.fpr"],
"fig": ["rgt_logo.gif","style.css","default_motif_logo.png","jquery-1.11.1.js","jquery.tablesorter.min.js","tdf_logo.png", "viz_logo.png"],
}
for copy_folder in copy_files_dictionary.keys():
    copy_dest_path = path.join(options.param_rgt_data_location,copy_folder)
    if not path.exists(copy_dest_path): makedirs(copy_dest_path)
    for copy_file in copy_files_dictionary[copy_folder]:
        copy_source_file = path.join(script_dir,"data",copy_folder,copy_file)
        copy_dest_file = path.join(copy_dest_path,copy_file)
        if not path.exists(copy_dest_file): 
            try: copytree(copy_source_file, copy_dest_file)
            except OSError as exc:
                if exc.errno == ENOTDIR: 
                    copy(copy_source_file, copy_dest_file)
                else: 
                    raise
    
###################################################################################################
# Setup Function
###################################################################################################

# Parameters
short_description = "Toolkit to perform regulatory genomics data analysis"
classifiers_list = ["Topic :: Scientific/Engineering :: Bio-Informatic",
                    "Topic :: Scientific/Engineering :: Artificial Intelligence"]
keywords_list = ["ChIP-seq","DNase-seq","Peak Calling","Motif Discovery","Motif Enrichment","HMM"]
author_list = ["Eduardo G. Gusmao","Manuel Allhoff","Joseph Kuo","Ivan G. Costa"]
corresponding_mail = "software@costalab.org"
license_type = "GPL"
package_data_dictionary = {"rgt": [path.basename(data_config_path_file_name)]}

# External scripts
external_scripts=[]
for tool_option in options.param_rgt_tool:
    for e in tools_dictionary[tool_option][3]: external_scripts.append(e)

# Fetching Additional Structural Files
readme_file_name = path.join(path.dirname(path.abspath(__file__)), "README.txt")

# Fetching Long Description
readme_file = open(readme_file_name,"r")
long_description =readme_file.read() + "nn"
readme_file.close()

# Setup Function
setup(name = "RGT",
      version = current_version,
      description = short_description,
      long_description = long_description,
      classifiers = classifiers_list,
      keywords = ", ".join(keywords_list),
      author = ", ".join(author_list),
      author_email = corresponding_mail,
      license = license_type,
      packages = find_packages(),
      package_data = package_data_dictionary,
      entry_points = current_entry_points,
      install_requires = current_install_requires,
      scripts = external_scripts
)

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

if(current_user):
    current_user_uid = getpwnam(current_user).pw_uid
    current_user_gid = getpwnam(current_user).pw_gid
    recursive_chown_chmod(options.param_rgt_data_location,current_user_uid,current_user_gid,default_file_permission,default_path_permission)


