
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
parser.add_option("--hg38", dest = "hg38", action = "store_true", default = False,
                  help = ("Fetch human genome files."))
parser.add_option("--mm9", dest = "mm9", action = "store_true", default = False,
                  help = ("Fetch mouse files."))
parser.add_option("--mm10", dest = "mm10", action = "store_true", default = False,
                  help = ("Fetch mouse files."))
parser.add_option("--zv9", dest = "zv9", action = "store_true", default = False,
                  help = ("Fetch zebrafish files."))
parser.add_option("--zv10", dest = "zv10", action = "store_true", default = False,
                  help = ("Fetch zebrafish files."))

parser.add_option("--hg19-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg19 genome (all chromosomes in the same file).",
                  dest = "hg19_genome_path", default = None)
parser.add_option("--hg38-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg38 genome (all chromosomes in the same file).",
                  dest = "hg38_genome_path", default = None)
parser.add_option("--mm9-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm9 genome (all chromosomes in the same file).",
                  dest = "mm9_genome_path", default = None)
parser.add_option("--mm10-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm9 genome (all chromosomes in the same file).",
                  dest = "mm10_genome_path", default = None)
parser.add_option("--zv9-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing zv9 genome (all chromosomes in the same file).",
                  dest = "zv9_genome_path", default = None)
parser.add_option("--zv10-genome-path", type = "string", metavar="STRING",
                  help = "Path to an already existing zv10 genome (all chromosomes in the same file).",
                  dest = "zv10_genome_path", default = None)

parser.add_option("--hg19-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg19 GTF file.",
                  dest = "hg19_gtf_path", default = None)
parser.add_option("--hg38-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing hg38 GTF file.",
                  dest = "hg38_gtf_path", default = None)
parser.add_option("--mm9-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm9 GTF file.",
                  dest = "mm9_gtf_path", default = None)
parser.add_option("--mm10-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing mm10 GTF file.",
                  dest = "mm10_gtf_path", default = None)
parser.add_option("--zv9-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing zv9 GTF file.",
                  dest = "zv9_gtf_path", default = None)
parser.add_option("--zv10-gtf-path", type = "string", metavar="STRING",
                  help = "Path to an already existing zv10 GTF file.",
                  dest = "zv10_gtf_path", default = None)

options, arguments = parser.parse_args()
if(options.all):
    options.hg19 = True
    options.hg38 = True
    options.mm9 = True
    options.mm10 = True
    options.zv9 = True
    options.zv10 = True

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
    output_genome_file_name = path.join(output_location,"genome_hg19.fa")
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
    gtf_output_file_name = path.join(output_location,"gencode.v19.annotation.gtf")
    if(options.hg19_gtf_path):
        print "Creating symbolic link to HG19 GTF"
        system("ln -s "+options.hg19_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        gtf_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"gencode.v19.annotation.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading hg19 GTF (gene annotation) from genode"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print "OK"

###################################################################################################
# Genomic Data HG38
###################################################################################################

if(options.hg38):

    output_location = path.join(curr_dir,"hg38")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome_hg38.fa")
    if(options.hg38_genome_path):
        print "Creating symbolic link to HG38 genome"
        system("ln -s "+options.hg38_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz"
        # chr_list = ["chr"+str(e) for e in range(1,23)+["X","Y","M"]]
        output_genome_file = open(output_genome_file_name,"w")
        # for chr_name in chr_list:
        print "Downloading hg38 genome"
        gz_file_name = path.join(output_location,"GRCh38.primary_assembly.genome.fa.gz")
        # if(path.isfile(gz_file_name)): remove(gz_file_name)
        system("wget "+gen_root_url+" -P "+output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_genome_file.write( gz_file.read() )
        gz_file.close()
        remove(gz_file_name)
        print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"gencode.v24.annotation.gtf")
    if(options.hg38_gtf_path):
        print "Creating symbolic link to HG38 GTF"
        system("ln -s "+options.hg38_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        gtf_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"gencode.v24.annotation.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading hg19 GTF (gene annotation) from genode"
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
    output_genome_file_name = path.join(output_location,"genome_mm9.fa")
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
    gtf_output_file_name = path.join(output_location,"gencode.vM1.annotation.gtf")
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

###################################################################################################
# Genomic Data MM10
###################################################################################################

if(options.mm10):

    output_location = path.join(curr_dir,"mm10")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome_mm10.fa")
    if(options.mm10_genome_path):
        print "Creating symbolic link to MM10 genome"
        system("ln -s "+options.mm10_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz"
        output_genome_file = open(output_genome_file_name, "w")
        print "Downloading mm10 genome"
        gz_file_name = path.join(output_location, "GRCm38.primary_assembly.genome.fa.gz")
        system("wget " + gen_root_url + " -P " + output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_genome_file.write(gz_file.read())
        gz_file.close()
        remove(gz_file_name)
        print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"gencode.vM11.annotation.gtf")
    if(options.mm10_gtf_path):
        print "Creating symbolic link to MM10 GTF"
        system("ln -s "+options.mm10_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        gtf_url = "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"gencode.vM11.annotation.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading MM10 GTF (gene annotation)"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print "OK"

###################################################################################################
# Genomic Data ZV9
###################################################################################################

if(options.zv9):

    output_location = path.join(curr_dir,"zv9")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome_zv9_ensembl_release_79.fa")
    if(options.zv9_genome_path):
        print "Creating symbolic link to ZV9 genome"
        system("ln -s "+options.zv9_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "ftp://ftp.ensembl.org/pub/release-79/fasta/danio_rerio/dna/"
        chr_list = [str(e) for e in range(1,25) ]+["MT"]
        output_genome_file = open(output_genome_file_name,"w")
        for chr_name in chr_list:
            print "Downloading ZV9 genome (chromosome "+chr_name+")"
            gz_file_name = path.join(output_location,"Danio_rerio.Zv9.dna.chromosome."+chr_name+".fa.gz")
            if(path.isfile(gz_file_name)): remove(gz_file_name)
            # Danio_rerio.Zv9.dna.chromosome.8.fa.gz
            system("wget "+gen_root_url+"Danio_rerio.Zv9.dna.chromosome."+chr_name+".fa.gz -P "+output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write( gz_file.read() )
            gz_file.close()
            remove(gz_file_name)
            print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"Danio_rerio.Zv9.79.gtf")
    if(options.zv9_gtf_path):
        print "Creating symbolic link to ZV9 GTF"
        system("ln -s "+options.zv9_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        # ftp://ftp.ensembl.org/pub/release-79/gtf/danio_rerio/Danio_rerio.Zv9.79.gtf.gz
        gtf_url = "ftp://ftp.ensembl.org/pub/release-79/gtf/danio_rerio/Danio_rerio.Zv9.79.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"Danio_rerio.Zv9.79.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading ZV9 GTF (gene annotation)"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        ## TODO add chr in each entry
        print "OK"


###################################################################################################
# Genomic Data ZV10
###################################################################################################

if(options.zv10):

    output_location = path.join(curr_dir,"zv10")
    if(not path.exists(output_location)): mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location,"genome_zv10_ensembl_release_84.fa")
    if(options.zv10_genome_path):
        print "Creating symbolic link to ZV10 genome"
        system("ln -s "+options.zv10_genome_path+" "+output_genome_file_name)
        print "OK"
    else:
        gen_root_url = "ftp://ftp.ensembl.org/pub/release-84/fasta/danio_rerio/dna/"
        chr_list = [str(e) for e in range(1,25) ]+["MT"]
        output_genome_file = open(output_genome_file_name,"w")
        for chr_name in chr_list:
            print "Downloading ZV10 genome (chromosome "+chr_name+")"
            gz_file_name = path.join(output_location,"Danio_rerio.GRCz10.dna.chromosome."+chr_name+".fa.gz")
            if(path.isfile(gz_file_name)): remove(gz_file_name)
            # Danio_rerio.GRCz10.dna.chromosome.8.fa.gz
            system("wget "+gen_root_url+"Danio_rerio.GRCz10.dna.chromosome."+chr_name+".fa.gz -P "+output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write( gz_file.read() )
            gz_file.close()
            remove(gz_file_name)
            print "OK"
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location,"Danio_rerio.GRCz10.84.gtf")
    if(options.zv10_gtf_path):
        print "Creating symbolic link to ZV10 GTF"
        system("ln -s "+options.zv10_gtf_path+" "+gtf_output_file_name)
        print "OK"
    else:
        # ftp://ftp.ensembl.org/pub/release-84/gtf/danio_rerio/Danio_rerio.GRCz10.84.gtf.gz
        gtf_url = "ftp://ftp.ensembl.org/pub/release-84/gtf/danio_rerio/Danio_rerio.GRCz10.84.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location,"Danio_rerio.GRCz10.84.gtf.gz")
        if(path.isfile(gtf_output_file_name_gz)): remove(gtf_output_file_name_gz)
        print "Downloading ZV10 GTF (gene annotation)"
        system("wget "+gtf_url+" -P "+output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name,"w")
        gtf_output_file.write( gz_file.read() )
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        # system("sed -i -e 's/^/chr/' "+gtf_output_file_name)
        ## TODO add chr in each entry
        print "OK"
