
###################################################################################################
# Setup Binaries Script
#   This script will fetch the necessary binaries for all tools. Be aware that it will consider
#   only 64 bits system (linux and MAC OS X).
###################################################################################################

# Import
from os import system, path, walk, remove, mkdir
from sys import platform
from glob import glob
from optparse import OptionParser

# Optional Input
usage_message = "python setupData.py [options]"

# Initializing Option Parser
parser = OptionParser(usage = usage_message)
options, arguments = parser.parse_args()


###################################################################################################
# Parameters
###################################################################################################

# Current rgt data path
curr_dir = path.dirname(path.realpath(__file__))

# Platform
supported_platforms = ["linux","linux2","darwin"]
if(platform not in supported_platforms):
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Binary Data
###################################################################################################

# Fetching UCSC binaries
ucsc_util_list = ["bigBedToBed","bedToBigBed"]
if("linux" in platform): bin_root = "http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/"
else: bin_root = "http://hgdownload.cse.ucsc.edu/admin/exe/macOSX.x86_64/"
out_loc = path.join(curr_dir,"bin")
if(not path.exists(out_loc)): mkdir(out_loc)
for util_tool in ucsc_util_list:
    print "Downloading UCSC binary script "+util_tool
    if(not path.isfile(path.join(out_loc,util_tool))):
        system("wget "+bin_root+util_tool+" -P "+out_loc)
    print "OK"


