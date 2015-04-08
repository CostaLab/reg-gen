
###################################################################################################
# Setup Genomic Data Script
#   This script will fetch the necessary genomic data for all tools. Be aware that it will download
#   big data sets and it might be storage consuming. Also, this script will need an Internet
#   connection and may take several minutes to execute depending on the Internet speed.
###################################################################################################

# Import
from os import system, path, walk, remove, mkdir
from sys import platform
from glob import glob
from optparse import OptionParser
import gzip

# Optional Input
usage_message = "python setupGenomicData.py [options]"

# Initializing Option Parser
parser = OptionParser(usage = usage_message)

# Parameter: RGT Data Location
parser.add_option("--all", dest = "all", action = "store_true", default = False,
                  help = ("Fetch all data sets."))
parser.add_option("--hg19", dest = "hg19", action = "store_true", default = False,
                  help = ("Fetch human genome files."))
parser.add_option("--mm9", dest = "mm9", action = "store_true", default = False,
                  help = ("Fetch mouse files."))
parser.add_option("--hg19-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg19 genome (all chromosomes in the same file).",
                  dest = "hg19_genome_path", default = None)
parser.add_option("--mm9-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm9 genome (all chromosomes in the same file).",
                  dest = "mm9_genome_path", default = None)
parser.add_option("--hg19-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg19 GTF file.",
                  dest = "hg19_gtf_path", default = None)
parser.add_option("--mm9-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm9 GTF file.",
                  dest = "mm9_gtf_path", default = None)
options, arguments = parser.parse_args()
if(options.all):
    options.hg19 = True
    options.mm9 = True

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
# Genomic Data HG19
###################################################################################################

if(options.hg19):

    output_location = path.join(curr_dir,"hg19")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome.fa")
    if(options.hg19_genome_path):
        print "Creating symbolic link to HG19 genome"
        system("ln -s "+options.hg19_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
        chr_list = ["chr"+str(e) for e in range(1,23)+["X","Y","M"]]
        output_genome_file = open(output_genome_file_name,"w")
        for chr_name in chr_list:
            print "Downloading hg19 genome ("+chr_name+")"
            gz_file_name = path.join(output_location,chr_name+".fa.gz")
            if(path.isfile(gz_file_name)): remove(gz_file_name)
            system("wget "+gen_root_url+chr_name+".fa.gz -P "+output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write( gz_file.read() )
            gz_file.close()
            remove(gz_file_name)
            print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"gencode_annotation.gtf")
    if(options.hg19_gtf_path):
        print "Creating symbolic link to HG19 GTF"
        system("ln -s "+options.hg19_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        gtf_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"gencode.v19.annotation.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading hg19 GTF (gene annotation)"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print "OK"

###################################################################################################
# Genomic Data MM9
###################################################################################################

if(options.mm9):

    output_location = path.join(curr_dir,"mm9")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome.fa")
    if(options.mm9_genome_path):
        print "Creating symbolic link to MM9 genome"
        system("ln -s "+options.mm9_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/"
        chr_list = ["chr"+str(e) for e in range(1,20)+["X","Y","M"]]
        output_genome_file = open(output_genome_file_name,"w")
        for chr_name in chr_list:
            print "Downloading MM9 genome ("+chr_name+")"
            gz_file_name = path.join(output_location,chr_name+".fa.gz")
            if(path.isfile(gz_file_name)): remove(gz_file_name)
            system("wget "+gen_root_url+chr_name+".fa.gz -P "+output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write( gz_file.read() )
            gz_file.close()
            remove(gz_file_name)
            print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"gencode_annotation.gtf")
    if(options.mm9_gtf_path):
        print "Creating symbolic link to MM9 GTF"
        system("ln -s "+options.mm9_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        gtf_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"gencode.vM1.annotation.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading MM9 GTF (gene annotation)"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print "OK"


